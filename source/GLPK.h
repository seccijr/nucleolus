#include "common.h"
#include <glpk.h>

void GLPK(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &excess, double &prec, vector<bool> &unsettled, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, vector<vector<bool>> &A, double &t1, vector<double> &singleton_bounds, bool &nlsu);
void pivot(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<vector<bool>> &A, vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned int &piv, unsigned int &sr_count, unsigned short int &iter, vector<bool> &unsettled_p, vector<double> &singleton_bounds, double &epsi_old, bool &nlsu);
void step(vector<bool> &T, vector<bool> &T2, vector<bool> &unsettled, vector<bool> &unsettled_p, unsigned int &s, vector<vector<bool>> &A, double &epsi, vector<double> &excess, vector<double> &d, unsigned short int &n, vector<double> &x, vector<double> &singleton_bounds, bool &disp, double &prec);

void subroutine(vector<bool> &U, vector<bool> &U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2,
                vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n,
                unsigned int &tight_size, unsigned short int &tight2_size, unsigned short int &rank, bool &disp,
                vector<vector<bool>> &Asettled, unsigned int &sr_count, bool &u, unsigned int &s,
                vector<unsigned int> &T_coord, vector<unsigned int> &T2_coord, vector<bool> &unsettled,
                double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled, bool &nlsu);

void subr_upd(vector<vector<double>> &Arref, vector<bool> &J, unsigned int &i, unsigned short int &n, double &prec,
              vector<bool> &U, vector<bool> &U2, unsigned int &sumt, unsigned int &sumt2, vector<bool> &t,
              vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &tight_size,
              unsigned short int &tight2_size, unsigned short int &rank, vector<bool> &unsettled,
              vector<vector<bool>> &Asettled, bool &disp, unsigned int &s, vector<unsigned int> &T_coord,
              vector<unsigned int> &T2_coord, double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled,
              glp_prob *lp, int ia[], int ja[], double ar[], int &ia_index, int &q_constr_index,
              vector<int> &q_indices);

void BNF_mem(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &excess, double &prec,
             vector<bool> &unsettled, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t,
             vector<double> &x, vector<bool> &a, double &t1, vector<double> &singleton_bounds, bool &nlsu);

void pivot_mem(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<bool> &a,
               vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank,
               vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned int &piv,
               unsigned int &sr_count, unsigned short int &iter, vector<bool> &unsettled_p, vector<double> &singleton_bounds,
               double &epsi_old, bool &nlsu);

void step_mem(vector<bool> &T, vector<bool> &T2, vector<bool> &unsettled, vector<bool> &unsettled_p, unsigned int &s,
              vector<bool> &a, double &epsi, vector<double> &excess, vector<double> &d, unsigned short int &n,
              vector<double> &x, vector<double> &singleton_bounds, bool &disp, double &prec);

void tight_coal(vector<bool> &T, vector<double> &excess, double &epsi, double &prec, unsigned int &s,
                vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size);

void tight_coal2(vector<bool> &T2, vector<double> &x, vector<double> &singleton_bounds, double &prec, unsigned short int &n,
                 vector<unsigned int> &T2_coord, vector<bool> &unsettled_p, unsigned short int &t2_size);

void imprdir(vector<double> &d, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size,
             vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<bool> &U, vector<bool> &U2,
             unsigned short int &rank, vector<vector<bool>> &Asettled, bool &disp);
