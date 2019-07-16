#include "common.h"

void subr_upd_sg(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, IloExpr &q, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda, unsigned short int &rank, bool &disp, unsigned int &sett, vector<vector<bool>> &Asettled, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <bool> &v, vector<unsigned int> &T2_coord);
void subroutine_sg(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, unsigned int &sr, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <bool> &v, vector<unsigned int> &T2_coord);
void subr_upd(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, IloExpr &q, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda, unsigned short int &rank, bool &disp, unsigned int &sett, vector<vector<bool>> &Asettled, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <double> &v, vector<unsigned int> &T2_coord);
void subroutine(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, unsigned int &sr, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <double> &v, vector<unsigned int> &T2_coord);
void iteration_sg_mem(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds);
void SP_sg_mem(bool &disp, unsigned short int &n, vector<bool> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu);
void iteration_mem(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds, bool &nlsu);
void SP_mem(bool &disp, unsigned short int &n, vector<double> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu);
void iteration_sg(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<bool> &singleton_bounds, bool &nlsu);
void SP_sg(bool &disp, unsigned short int &n, vector<bool> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu);
void iteration(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds, bool &nlsu);
void SP(bool &disp, unsigned short int &n, vector<double> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu);
