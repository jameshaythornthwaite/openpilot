
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4363878931908636218) {
   out_4363878931908636218[0] = delta_x[0] + nom_x[0];
   out_4363878931908636218[1] = delta_x[1] + nom_x[1];
   out_4363878931908636218[2] = delta_x[2] + nom_x[2];
   out_4363878931908636218[3] = delta_x[3] + nom_x[3];
   out_4363878931908636218[4] = delta_x[4] + nom_x[4];
   out_4363878931908636218[5] = delta_x[5] + nom_x[5];
   out_4363878931908636218[6] = delta_x[6] + nom_x[6];
   out_4363878931908636218[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2367728689678461946) {
   out_2367728689678461946[0] = -nom_x[0] + true_x[0];
   out_2367728689678461946[1] = -nom_x[1] + true_x[1];
   out_2367728689678461946[2] = -nom_x[2] + true_x[2];
   out_2367728689678461946[3] = -nom_x[3] + true_x[3];
   out_2367728689678461946[4] = -nom_x[4] + true_x[4];
   out_2367728689678461946[5] = -nom_x[5] + true_x[5];
   out_2367728689678461946[6] = -nom_x[6] + true_x[6];
   out_2367728689678461946[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_748139986113288943) {
   out_748139986113288943[0] = 1.0;
   out_748139986113288943[1] = 0.0;
   out_748139986113288943[2] = 0.0;
   out_748139986113288943[3] = 0.0;
   out_748139986113288943[4] = 0.0;
   out_748139986113288943[5] = 0.0;
   out_748139986113288943[6] = 0.0;
   out_748139986113288943[7] = 0.0;
   out_748139986113288943[8] = 0.0;
   out_748139986113288943[9] = 1.0;
   out_748139986113288943[10] = 0.0;
   out_748139986113288943[11] = 0.0;
   out_748139986113288943[12] = 0.0;
   out_748139986113288943[13] = 0.0;
   out_748139986113288943[14] = 0.0;
   out_748139986113288943[15] = 0.0;
   out_748139986113288943[16] = 0.0;
   out_748139986113288943[17] = 0.0;
   out_748139986113288943[18] = 1.0;
   out_748139986113288943[19] = 0.0;
   out_748139986113288943[20] = 0.0;
   out_748139986113288943[21] = 0.0;
   out_748139986113288943[22] = 0.0;
   out_748139986113288943[23] = 0.0;
   out_748139986113288943[24] = 0.0;
   out_748139986113288943[25] = 0.0;
   out_748139986113288943[26] = 0.0;
   out_748139986113288943[27] = 1.0;
   out_748139986113288943[28] = 0.0;
   out_748139986113288943[29] = 0.0;
   out_748139986113288943[30] = 0.0;
   out_748139986113288943[31] = 0.0;
   out_748139986113288943[32] = 0.0;
   out_748139986113288943[33] = 0.0;
   out_748139986113288943[34] = 0.0;
   out_748139986113288943[35] = 0.0;
   out_748139986113288943[36] = 1.0;
   out_748139986113288943[37] = 0.0;
   out_748139986113288943[38] = 0.0;
   out_748139986113288943[39] = 0.0;
   out_748139986113288943[40] = 0.0;
   out_748139986113288943[41] = 0.0;
   out_748139986113288943[42] = 0.0;
   out_748139986113288943[43] = 0.0;
   out_748139986113288943[44] = 0.0;
   out_748139986113288943[45] = 1.0;
   out_748139986113288943[46] = 0.0;
   out_748139986113288943[47] = 0.0;
   out_748139986113288943[48] = 0.0;
   out_748139986113288943[49] = 0.0;
   out_748139986113288943[50] = 0.0;
   out_748139986113288943[51] = 0.0;
   out_748139986113288943[52] = 0.0;
   out_748139986113288943[53] = 0.0;
   out_748139986113288943[54] = 1.0;
   out_748139986113288943[55] = 0.0;
   out_748139986113288943[56] = 0.0;
   out_748139986113288943[57] = 0.0;
   out_748139986113288943[58] = 0.0;
   out_748139986113288943[59] = 0.0;
   out_748139986113288943[60] = 0.0;
   out_748139986113288943[61] = 0.0;
   out_748139986113288943[62] = 0.0;
   out_748139986113288943[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_1532464344372039942) {
   out_1532464344372039942[0] = state[0];
   out_1532464344372039942[1] = state[1];
   out_1532464344372039942[2] = state[2];
   out_1532464344372039942[3] = state[3];
   out_1532464344372039942[4] = state[4];
   out_1532464344372039942[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1532464344372039942[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1532464344372039942[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2578963947479186478) {
   out_2578963947479186478[0] = 1;
   out_2578963947479186478[1] = 0;
   out_2578963947479186478[2] = 0;
   out_2578963947479186478[3] = 0;
   out_2578963947479186478[4] = 0;
   out_2578963947479186478[5] = 0;
   out_2578963947479186478[6] = 0;
   out_2578963947479186478[7] = 0;
   out_2578963947479186478[8] = 0;
   out_2578963947479186478[9] = 1;
   out_2578963947479186478[10] = 0;
   out_2578963947479186478[11] = 0;
   out_2578963947479186478[12] = 0;
   out_2578963947479186478[13] = 0;
   out_2578963947479186478[14] = 0;
   out_2578963947479186478[15] = 0;
   out_2578963947479186478[16] = 0;
   out_2578963947479186478[17] = 0;
   out_2578963947479186478[18] = 1;
   out_2578963947479186478[19] = 0;
   out_2578963947479186478[20] = 0;
   out_2578963947479186478[21] = 0;
   out_2578963947479186478[22] = 0;
   out_2578963947479186478[23] = 0;
   out_2578963947479186478[24] = 0;
   out_2578963947479186478[25] = 0;
   out_2578963947479186478[26] = 0;
   out_2578963947479186478[27] = 1;
   out_2578963947479186478[28] = 0;
   out_2578963947479186478[29] = 0;
   out_2578963947479186478[30] = 0;
   out_2578963947479186478[31] = 0;
   out_2578963947479186478[32] = 0;
   out_2578963947479186478[33] = 0;
   out_2578963947479186478[34] = 0;
   out_2578963947479186478[35] = 0;
   out_2578963947479186478[36] = 1;
   out_2578963947479186478[37] = 0;
   out_2578963947479186478[38] = 0;
   out_2578963947479186478[39] = 0;
   out_2578963947479186478[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2578963947479186478[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2578963947479186478[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2578963947479186478[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2578963947479186478[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2578963947479186478[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2578963947479186478[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2578963947479186478[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2578963947479186478[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2578963947479186478[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2578963947479186478[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2578963947479186478[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2578963947479186478[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2578963947479186478[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2578963947479186478[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2578963947479186478[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2578963947479186478[56] = 0;
   out_2578963947479186478[57] = 0;
   out_2578963947479186478[58] = 0;
   out_2578963947479186478[59] = 0;
   out_2578963947479186478[60] = 0;
   out_2578963947479186478[61] = 0;
   out_2578963947479186478[62] = 0;
   out_2578963947479186478[63] = 1;
}
void h_25(double *state, double *unused, double *out_947661413387208961) {
   out_947661413387208961[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7546209857534697687) {
   out_7546209857534697687[0] = 0;
   out_7546209857534697687[1] = 0;
   out_7546209857534697687[2] = 0;
   out_7546209857534697687[3] = 0;
   out_7546209857534697687[4] = 0;
   out_7546209857534697687[5] = 0;
   out_7546209857534697687[6] = 1;
   out_7546209857534697687[7] = 0;
}
void h_24(double *state, double *unused, double *out_8368114139978590899) {
   out_8368114139978590899[0] = state[4];
   out_8368114139978590899[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6001824904607495755) {
   out_6001824904607495755[0] = 0;
   out_6001824904607495755[1] = 0;
   out_6001824904607495755[2] = 0;
   out_6001824904607495755[3] = 0;
   out_6001824904607495755[4] = 1;
   out_6001824904607495755[5] = 0;
   out_6001824904607495755[6] = 0;
   out_6001824904607495755[7] = 0;
   out_6001824904607495755[8] = 0;
   out_6001824904607495755[9] = 0;
   out_6001824904607495755[10] = 0;
   out_6001824904607495755[11] = 0;
   out_6001824904607495755[12] = 0;
   out_6001824904607495755[13] = 1;
   out_6001824904607495755[14] = 0;
   out_6001824904607495755[15] = 0;
}
void h_30(double *state, double *unused, double *out_3780962411196355005) {
   out_3780962411196355005[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6410255963779292765) {
   out_6410255963779292765[0] = 0;
   out_6410255963779292765[1] = 0;
   out_6410255963779292765[2] = 0;
   out_6410255963779292765[3] = 0;
   out_6410255963779292765[4] = 1;
   out_6410255963779292765[5] = 0;
   out_6410255963779292765[6] = 0;
   out_6410255963779292765[7] = 0;
}
void h_26(double *state, double *unused, double *out_8534016932641807622) {
   out_8534016932641807622[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6700224773109598113) {
   out_6700224773109598113[0] = 0;
   out_6700224773109598113[1] = 0;
   out_6700224773109598113[2] = 0;
   out_6700224773109598113[3] = 0;
   out_6700224773109598113[4] = 0;
   out_6700224773109598113[5] = 0;
   out_6700224773109598113[6] = 0;
   out_6700224773109598113[7] = 1;
}
void h_27(double *state, double *unused, double *out_6892963051267079214) {
   out_6892963051267079214[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3381221781696382211) {
   out_3381221781696382211[0] = 0;
   out_3381221781696382211[1] = 0;
   out_3381221781696382211[2] = 0;
   out_3381221781696382211[3] = 1;
   out_3381221781696382211[4] = 0;
   out_3381221781696382211[5] = 0;
   out_3381221781696382211[6] = 0;
   out_3381221781696382211[7] = 0;
}
void h_29(double *state, double *unused, double *out_6168003306232441360) {
   out_6168003306232441360[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1343810939479272313) {
   out_1343810939479272313[0] = 0;
   out_1343810939479272313[1] = 1;
   out_1343810939479272313[2] = 0;
   out_1343810939479272313[3] = 0;
   out_1343810939479272313[4] = 0;
   out_1343810939479272313[5] = 0;
   out_1343810939479272313[6] = 0;
   out_1343810939479272313[7] = 0;
}
void h_28(double *state, double *unused, double *out_2117460458383937735) {
   out_2117460458383937735[0] = state[5];
   out_2117460458383937735[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3807987766153760559) {
   out_3807987766153760559[0] = 0;
   out_3807987766153760559[1] = 0;
   out_3807987766153760559[2] = 0;
   out_3807987766153760559[3] = 0;
   out_3807987766153760559[4] = 0;
   out_3807987766153760559[5] = 1;
   out_3807987766153760559[6] = 0;
   out_3807987766153760559[7] = 0;
   out_3807987766153760559[8] = 0;
   out_3807987766153760559[9] = 0;
   out_3807987766153760559[10] = 0;
   out_3807987766153760559[11] = 0;
   out_3807987766153760559[12] = 0;
   out_3807987766153760559[13] = 0;
   out_3807987766153760559[14] = 1;
   out_3807987766153760559[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
