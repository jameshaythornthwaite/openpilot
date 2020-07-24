/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_907173371217757873);
void inv_err_fun(double *nom_x, double *true_x, double *out_8081364877058298739);
void H_mod_fun(double *state, double *out_2045697321123555115);
void f_fun(double *state, double dt, double *out_5473203447822424204);
void F_fun(double *state, double dt, double *out_1156352550666527776);
void h_3(double *state, double *unused, double *out_4701478374013215574);
void H_3(double *state, double *unused, double *out_3890623566473753871);
void h_4(double *state, double *unused, double *out_6516659658683726097);
void H_4(double *state, double *unused, double *out_8117098049898543368);
void h_9(double *state, double *unused, double *out_3490194137804294281);
void H_9(double *state, double *unused, double *out_4981401289114065623);
void h_10(double *state, double *unused, double *out_3023605263989570279);
void H_10(double *state, double *unused, double *out_6948261743464078135);
void h_12(double *state, double *unused, double *out_6876824376853794091);
void H_12(double *state, double *unused, double *out_7435554915144550553);
void h_31(double *state, double *unused, double *out_348388349878683834);
void H_31(double *state, double *unused, double *out_906168875739015311);
void h_32(double *state, double *unused, double *out_529659947281651289);
void H_32(double *state, double *unused, double *out_4975449126032130713);
void h_13(double *state, double *unused, double *out_2263906346962012941);
void H_13(double *state, double *unused, double *out_403583952460359094);
void h_14(double *state, double *unused, double *out_3490194137804294281);
void H_14(double *state, double *unused, double *out_4981401289114065623);
void h_19(double *state, double *unused, double *out_3524830045566905425);
void H_19(double *state, double *unused, double *out_2385117127606726621);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);