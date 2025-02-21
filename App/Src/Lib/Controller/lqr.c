#include <stddef.h>
#include "lqr.h"
#include "chassis.h"
#include "user_lib.h"
#include "math.h"
#include "remote.h"
#include "buzzer.h"
#include "filter.h"

#define BALANCE_WHEEL_R 0.06f //Æ½ºâ±øÂÖ×Ó°ë¾¶m

struct MovingAverageFilter theta_ddot_filter_L, theta_ddot_filter_R;

float init_wheel_K_L[6] = {-17.162366, -2.016433, -1.130678, -1.530085, 15.895931, 2.388100};
float init_wheel_K_R[6] = {-17.162366, -2.016433, -1.130678, -1.530085, 15.895931, 2.388100};
float init_joint_K_L[6] = {15.857528, 2.386820, 1.485220, 1.953273, 37.371067, 4.026670};
float init_joint_K_R[6] = {15.857528, 2.386820, 1.485220, 1.953273, 37.371067, 4.026670};

float wheel_K_L[6] = {0, 0, 0, 0, 0, 0};
float joint_K_L[6] = {0, 0, 0, 0, 0, 0};

float wheel_K_R[6] = {0, 0, 0, 0, 0, 0};
float joint_K_R[6] = {0, 0, 0, 0, 0, 0};

float wheel_fitting_factor[6][4] = {
        {-3964.789786,1391.433119,-304.422387,-9.911213},
        {188.923997,-48.441236,-10.961798,-3.519304},

        {2197.033570,-560.667080,42.546502,-9.557216},
        {5976.125503,-1591.508475,150.116594,-12.176106},

        {-2732.364398,845.352460,-160.751119,34.961334},
        {-1645.126724,461.185227,-56.408919,6.651699}
};
float joint_fitting_factor[6][4] = {
        {8784.509441,-2674.427092,311.122091,5.732456},
        {217.630936,-55.175163,1.322862,3.683477},

        {4475.867511,-1149.150365,87.535506,4.274444},
        {551.266460,-47.896455,-22.847299,6.847363},

        {-800.481139,-44.949756,90.264108,30.104826},
        {2.236401,-55.993238,21.385757,2.228422}
};

////ÎÞ×°¼×°å
//float wheel_fitting_factor[6][4] = {
//        {-3821.520996,1325.093326,-287.827024,-9.935187},
//        {242.908604,-63.719965,-8.846525,-3.608726},
//
//        {2176.974867,-556.618549,42.332514,-9.325974},
//        {5974.388056,-1592.564625,150.675500,-12.062120},
//
//        {-2709.457075,827.227880,-155.172135,34.669776},
//        {-1839.699731,513.018275,-61.478342,6.998230}
//};
//float joint_fitting_factor[6][4] = {
//        {8570.748116,-2600.218118,300.664063,5.869814},
//        {200.528456,-50.160108,0.910455,3.737808},
//
//        {4277.394108,-1099.902838,83.950809,4.416515},
//        {362.461098,-0.050670,-26.761287,6.978830},
//
//        {-897.748195,-4.637262,84.120334,29.714867},
//        {50.355116,-70.985863,23.455697,1.992841}
//};

/////////////////////////////////////////////////////////////////////////////////////////////////>>

void chassis_K_matrix_fitting(float L0, float K[6], const float KL[6][4]) {
  for (int i = 0; i < 6; i++) {
    K[i] = KL[i][0] * powf(L0, 3) + KL[i][1] * powf(L0, 2) + KL[i][2] * powf(L0, 1) + KL[i][3] * powf(L0, 0);
  }
}

static float cal_leg_theta(float phi0, float phi) {
  float theta = 0, alpha = 0;//alpha is the Angle at which the virtual joint motor is turned
  alpha = PI / 2 - phi0;

  if (alpha * phi < 0) {
    theta = ABS(alpha) - ABS(phi);
    if ((alpha > 0) && (phi < 0)) {
      theta *= -1;
    } else {

    }
  } else {
    theta = ABS(alpha) + ABS(phi);
    if ((alpha < 0) && (phi < 0)) {
    } else {
      theta *= -1;
    }
  }
  return theta;
}

void Matrix_multiply(int rows1, int cols1, float matrix1[rows1][cols1],
                     int rows2, int cols2, float matrix2[rows2][cols2],
                     float result[rows1][cols2]) {
  if (cols1 != rows2)
    return;

  // Perform matrix multiplication
  for (int i = 0; i < rows1; ++i) {
    for (int j = 0; j < cols2; ++j) {
      result[i][j] = 0;
      for (int k = 0; k < cols1; ++k) {
        result[i][j] += matrix1[i][k] * matrix2[k][j];
      }
    }
  }
}

static void state_variable_update(struct Leg *leg_L, struct Leg *leg_R, float phi, float phi_dot) {
  if (leg_L == NULL || leg_R == NULL) {
    return;
  }
  //theta_last
  leg_L->state_variable_feedback.theta_last = leg_L->state_variable_feedback.theta;
  leg_R->state_variable_feedback.theta_last = leg_R->state_variable_feedback.theta;

  //theta
  leg_L->state_variable_feedback.theta =
      cal_leg_theta(leg_L->vmc.forward_kinematics.fk_phi.phi0, phi);
  leg_R->state_variable_feedback.theta =
      cal_leg_theta(leg_R->vmc.forward_kinematics.fk_phi.phi0, phi);

  //theta_ddot
  leg_L->state_variable_feedback.theta_dot_last = leg_L->state_variable_feedback.theta_dot;
  leg_L->state_variable_feedback.theta_dot =
      (leg_L->state_variable_feedback.theta - leg_L->state_variable_feedback.theta_last)
          / (CHASSIS_PERIOD * 0.001f);
  float theta_ddot_raw_L =
      (leg_L->state_variable_feedback.theta_dot - leg_L->state_variable_feedback.theta_dot_last)
          / (CHASSIS_PERIOD * 0.001f);
  updateFilter(&theta_ddot_filter_L, theta_ddot_raw_L);
  leg_L->state_variable_feedback.theta_ddot = getFilteredValue(&theta_ddot_filter_L);

  leg_R->state_variable_feedback.theta_dot_last = leg_R->state_variable_feedback.theta_dot;
  leg_R->state_variable_feedback.theta_dot =
      (leg_R->state_variable_feedback.theta - leg_R->state_variable_feedback.theta_last)
          / (CHASSIS_PERIOD * 0.001f);
  float theta_ddot_raw_R =
      (leg_R->state_variable_feedback.theta_dot - leg_R->state_variable_feedback.theta_dot_last)
          / (CHASSIS_PERIOD * 0.001f);
  updateFilter(&theta_ddot_filter_R, theta_ddot_raw_R);
  leg_R->state_variable_feedback.theta_ddot = getFilteredValue(&theta_ddot_filter_R);

  //x
  if (get_chassis()->chassis_ctrl_info.v_m_per_s != 0
      || get_chassis()->is_chassis_offground == true) {
    leg_L->state_variable_feedback.x = 0;
    leg_R->state_variable_feedback.x = 0;
  } else {
    leg_L->state_variable_feedback.x =
        leg_L->state_variable_feedback.x
            + CHASSIS_PERIOD * 0.001f * ((float) get_wheel_motors()[0].angular_vel * BALANCE_WHEEL_R);
    leg_R->state_variable_feedback.x =
        leg_R->state_variable_feedback.x
            + CHASSIS_PERIOD * 0.001f * ((float) -get_wheel_motors()[1].angular_vel * BALANCE_WHEEL_R);

    VAL_LIMIT(leg_L->state_variable_feedback.x, -0.05, 0.05)
    VAL_LIMIT(leg_R->state_variable_feedback.x, -0.05, 0.05)
  }

  if (get_chassis()->chassis_ctrl_info.v_m_per_s != 0) {
    leg_L->state_variable_feedback.x_dot = leg_L->kalman_result[0];
    leg_R->state_variable_feedback.x_dot = leg_R->kalman_result[0];
  } else {
    leg_L->state_variable_feedback.x_dot = (float) get_wheel_motors()[0].angular_vel * BALANCE_WHEEL_R;
    leg_R->state_variable_feedback.x_dot = (float) -get_wheel_motors()[1].angular_vel * BALANCE_WHEEL_R;
  }
  leg_L->state_variable_feedback.x_dot_last = leg_L->state_variable_feedback.x_dot;
  leg_L->state_variable_feedback.x_ddot =
      (leg_L->state_variable_feedback.x_dot - leg_L->state_variable_feedback.x_dot_last)
          / (CHASSIS_PERIOD * 0.001f);

  leg_R->state_variable_feedback.x_dot_last = leg_R->state_variable_feedback.x_dot;
  leg_R->state_variable_feedback.x_ddot =
      (leg_R->state_variable_feedback.x_dot - leg_R->state_variable_feedback.x_dot_last)
          / (CHASSIS_PERIOD * 0.001f);

  //phi
  leg_L->state_variable_feedback.phi = phi;
  leg_L->state_variable_feedback.phi_dot = phi_dot;

  leg_R->state_variable_feedback.phi = phi;
  leg_R->state_variable_feedback.phi_dot = phi_dot;
}

static void state_variable_set(struct Chassis *chassis) {
  if (chassis == NULL) {
    return;
  }

  chassis->leg_L.state_variable_set_point.x = 0;
  chassis->leg_L.state_variable_set_point.x_dot = chassis->chassis_ctrl_info.v_m_per_s;
  chassis->leg_L.state_variable_set_point.theta = 0.05236f;
  chassis->leg_L.state_variable_set_point.theta_dot = 0;
  chassis->leg_L.state_variable_set_point.phi = 0;
  chassis->leg_L.state_variable_set_point.phi_dot = 0;

  chassis->leg_R.state_variable_set_point.x = 0;
  chassis->leg_R.state_variable_set_point.x_dot = chassis->chassis_ctrl_info.v_m_per_s;
  chassis->leg_R.state_variable_set_point.theta = 0.05236f;
  chassis->leg_R.state_variable_set_point.theta_dot = 0;
  chassis->leg_R.state_variable_set_point.phi = 0;
  chassis->leg_R.state_variable_set_point.phi_dot = 0;
}

static void state_variable_error(struct Leg *leg_L, struct Leg *leg_R) {
  if (leg_L == NULL || leg_R == NULL) {
    return;
  }

  leg_L->state_variable_error.x = leg_L->state_variable_feedback.x - leg_L->state_variable_set_point.x;
  leg_L->state_variable_error.x_dot = leg_L->state_variable_feedback.x_dot - leg_L->state_variable_set_point.x_dot;
  leg_L->state_variable_error.theta = leg_L->state_variable_feedback.theta - leg_L->state_variable_set_point.theta;
  leg_L->state_variable_error.theta_dot =
      leg_L->state_variable_feedback.theta_dot - leg_L->state_variable_set_point.theta_dot;
  leg_L->state_variable_error.phi = leg_L->state_variable_feedback.phi - leg_L->state_variable_set_point.phi;
  leg_L->state_variable_error.phi_dot =
      leg_L->state_variable_feedback.phi_dot - leg_L->state_variable_set_point.phi_dot;

  leg_R->state_variable_error.x = leg_R->state_variable_feedback.x - leg_R->state_variable_set_point.x;
  leg_R->state_variable_error.x_dot = leg_R->state_variable_feedback.x_dot - leg_R->state_variable_set_point.x_dot;
  leg_R->state_variable_error.theta = leg_R->state_variable_feedback.theta - leg_R->state_variable_set_point.theta;
  leg_R->state_variable_error.theta_dot =
      leg_R->state_variable_feedback.theta_dot - leg_R->state_variable_set_point.theta_dot;
  leg_R->state_variable_error.phi = leg_R->state_variable_feedback.phi - leg_R->state_variable_set_point.phi;
  leg_R->state_variable_error.phi_dot =
      leg_R->state_variable_feedback.phi_dot - leg_R->state_variable_set_point.phi_dot;
}

static void state_variable_out(struct Chassis *chassis) {
  if (chassis == NULL) {
    return;
  }
  if (chassis->is_chassis_offground) {
//    for (int i = 0; i < 6; i++) {
//      joint_K_L[i] *= 0.2f;
//      joint_K_R[i] *= 0.2f;
//    }
//    for (int i = 0; i < 6; i++) {
//      wheel_K_L[i] = 0;
//      wheel_K_R[i] = 0;
//    }
//    for (int i = 2; i < 6; i++) {
//      joint_K_L[i] = 0;
//      joint_K_R[i] = 0;
//    }
  }
//
//  if (chassis->chassis_ctrl_mode == CHASSIS_INIT) {
//
//    chassis->leg_L.state_variable_wheel_out.theta = chassis->leg_L.state_variable_error.theta * init_wheel_K_L[0];
//    chassis->leg_L.state_variable_wheel_out.theta_dot =
//        chassis->leg_L.state_variable_error.theta_dot * init_wheel_K_L[1];
//    chassis->leg_L.state_variable_wheel_out.x = chassis->leg_L.state_variable_error.x * init_wheel_K_L[2];
//    chassis->leg_L.state_variable_wheel_out.x_dot = chassis->leg_L.state_variable_error.x_dot * init_wheel_K_L[3];
//    chassis->leg_L.state_variable_wheel_out.phi = chassis->leg_L.state_variable_error.phi * init_wheel_K_L[4];
//    chassis->leg_L.state_variable_wheel_out.phi_dot = chassis->leg_L.state_variable_error.phi_dot * init_wheel_K_L[5];
//
//    chassis->leg_L.state_variable_joint_out.theta = chassis->leg_L.state_variable_error.theta * init_joint_K_L[0];
//    chassis->leg_L.state_variable_joint_out.theta_dot =
//        chassis->leg_L.state_variable_error.theta_dot * init_joint_K_L[1];
//    chassis->leg_L.state_variable_joint_out.x = chassis->leg_L.state_variable_error.x * init_joint_K_L[2];
//    chassis->leg_L.state_variable_joint_out.x_dot = chassis->leg_L.state_variable_error.x_dot * init_joint_K_L[3];
//    chassis->leg_L.state_variable_joint_out.phi = chassis->leg_L.state_variable_error.phi * init_joint_K_L[4];
//    chassis->leg_L.state_variable_joint_out.phi_dot = chassis->leg_L.state_variable_error.phi_dot * init_joint_K_L[5];
//
//    chassis->leg_R.state_variable_wheel_out.theta = chassis->leg_R.state_variable_error.theta * init_wheel_K_R[0];
//    chassis->leg_R.state_variable_wheel_out.theta_dot =
//        chassis->leg_R.state_variable_error.theta_dot * init_wheel_K_R[1];
//    chassis->leg_R.state_variable_wheel_out.x = chassis->leg_R.state_variable_error.x * init_wheel_K_R[2];
//    chassis->leg_R.state_variable_wheel_out.x_dot = chassis->leg_R.state_variable_error.x_dot * init_wheel_K_R[3];
//    chassis->leg_R.state_variable_wheel_out.phi = chassis->leg_R.state_variable_error.phi * init_wheel_K_R[4];
//    chassis->leg_R.state_variable_wheel_out.phi_dot = chassis->leg_R.state_variable_error.phi_dot * init_wheel_K_R[5];
//
//    chassis->leg_R.state_variable_joint_out.theta = chassis->leg_R.state_variable_error.theta * init_joint_K_R[0];
//    chassis->leg_R.state_variable_joint_out.theta_dot =
//        chassis->leg_R.state_variable_error.theta_dot * init_joint_K_R[1];
//    chassis->leg_R.state_variable_joint_out.x = chassis->leg_R.state_variable_error.x * init_joint_K_R[2];
//    chassis->leg_R.state_variable_joint_out.x_dot = chassis->leg_R.state_variable_error.x_dot * init_joint_K_R[3];
//    chassis->leg_R.state_variable_joint_out.phi = chassis->leg_R.state_variable_error.phi * init_joint_K_R[4];
//    chassis->leg_R.state_variable_joint_out.phi_dot = chassis->leg_R.state_variable_error.phi_dot * init_joint_K_R[5];


  chassis->leg_L.state_variable_wheel_out.theta = chassis->leg_L.state_variable_error.theta * wheel_K_L[0];
  chassis->leg_L.state_variable_wheel_out.theta_dot = chassis->leg_L.state_variable_error.theta_dot * wheel_K_L[1];
  chassis->leg_L.state_variable_wheel_out.x = chassis->leg_L.state_variable_error.x * wheel_K_L[2];
  chassis->leg_L.state_variable_wheel_out.x_dot = chassis->leg_L.state_variable_error.x_dot * wheel_K_L[3];
  chassis->leg_L.state_variable_wheel_out.phi = chassis->leg_L.state_variable_error.phi * wheel_K_L[4];
  chassis->leg_L.state_variable_wheel_out.phi_dot = chassis->leg_L.state_variable_error.phi_dot * wheel_K_L[5];

  chassis->leg_L.state_variable_joint_out.theta = chassis->leg_L.state_variable_error.theta * joint_K_L[0];
  chassis->leg_L.state_variable_joint_out.theta_dot = chassis->leg_L.state_variable_error.theta_dot * joint_K_L[1];
  chassis->leg_L.state_variable_joint_out.x = chassis->leg_L.state_variable_error.x * joint_K_L[2];
  chassis->leg_L.state_variable_joint_out.x_dot = chassis->leg_L.state_variable_error.x_dot * joint_K_L[3];
  chassis->leg_L.state_variable_joint_out.phi = chassis->leg_L.state_variable_error.phi * joint_K_L[4];
  chassis->leg_L.state_variable_joint_out.phi_dot = chassis->leg_L.state_variable_error.phi_dot * joint_K_L[5];

  chassis->leg_R.state_variable_wheel_out.theta = chassis->leg_R.state_variable_error.theta * wheel_K_R[0];
  chassis->leg_R.state_variable_wheel_out.theta_dot = chassis->leg_R.state_variable_error.theta_dot * wheel_K_R[1];
  chassis->leg_R.state_variable_wheel_out.x = chassis->leg_R.state_variable_error.x * wheel_K_R[2];
  chassis->leg_R.state_variable_wheel_out.x_dot = chassis->leg_R.state_variable_error.x_dot * wheel_K_R[3];
  chassis->leg_R.state_variable_wheel_out.phi = chassis->leg_R.state_variable_error.phi * wheel_K_R[4];
  chassis->leg_R.state_variable_wheel_out.phi_dot = chassis->leg_R.state_variable_error.phi_dot * wheel_K_R[5];

  chassis->leg_R.state_variable_joint_out.theta = chassis->leg_R.state_variable_error.theta * joint_K_R[0];
  chassis->leg_R.state_variable_joint_out.theta_dot = chassis->leg_R.state_variable_error.theta_dot * joint_K_R[1];
  chassis->leg_R.state_variable_joint_out.x = chassis->leg_R.state_variable_error.x * joint_K_R[2];
  chassis->leg_R.state_variable_joint_out.x_dot = chassis->leg_R.state_variable_error.x_dot * joint_K_R[3];
  chassis->leg_R.state_variable_joint_out.phi = chassis->leg_R.state_variable_error.phi * joint_K_R[4];
  chassis->leg_R.state_variable_joint_out.phi_dot = chassis->leg_R.state_variable_error.phi_dot * joint_K_R[5];

}

/*******************************************************************************
 *                                     LQR                                     *
 *******************************************************************************/

void lqr_ctrl(struct Chassis *chassis) {
  chassis_K_matrix_fitting(chassis->leg_L.vmc.forward_kinematics.fk_L0.L0 * 0.25f, wheel_K_L, wheel_fitting_factor);
  chassis_K_matrix_fitting(chassis->leg_L.vmc.forward_kinematics.fk_L0.L0 * 0.25f, joint_K_L, joint_fitting_factor);
  chassis_K_matrix_fitting(chassis->leg_R.vmc.forward_kinematics.fk_L0.L0 * 0.25f, wheel_K_R, wheel_fitting_factor);
  chassis_K_matrix_fitting(chassis->leg_R.vmc.forward_kinematics.fk_L0.L0 * 0.25f, joint_K_R, joint_fitting_factor);
  state_variable_update(&chassis->leg_L,
                        &chassis->leg_R,
                        chassis->imu_reference.pitch_angle,
                        chassis->imu_reference.pitch_gyro);
  state_variable_set(chassis);
  state_variable_error(&chassis->leg_L, &chassis->leg_R);
  state_variable_out(chassis);
}