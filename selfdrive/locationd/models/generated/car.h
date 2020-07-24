/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4363878931908636218);
void inv_err_fun(double *nom_x, double *true_x, double *out_2367728689678461946);
void H_mod_fun(double *state, double *out_748139986113288943);
void f_fun(double *state, double dt, double *out_1532464344372039942);
void F_fun(double *state, double dt, double *out_2578963947479186478);
void h_25(double *state, double *unused, double *out_947661413387208961);
void H_25(double *state, double *unused, double *out_7546209857534697687);
void h_24(double *state, double *unused, double *out_8368114139978590899);
void H_24(double *state, double *unused, double *out_6001824904607495755);
void h_30(double *state, double *unused, double *out_3780962411196355005);
void H_30(double *state, double *unused, double *out_6410255963779292765);
void h_26(double *state, double *unused, double *out_8534016932641807622);
void H_26(double *state, double *unused, double *out_6700224773109598113);
void h_27(double *state, double *unused, double *out_6892963051267079214);
void H_27(double *state, double *unused, double *out_3381221781696382211);
void h_29(double *state, double *unused, double *out_6168003306232441360);
void H_29(double *state, double *unused, double *out_1343810939479272313);
void h_28(double *state, double *unused, double *out_2117460458383937735);
void H_28(double *state, double *unused, double *out_3807987766153760559);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
