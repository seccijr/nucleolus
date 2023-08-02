#include "gen_game.h"

int main()
{
	unsigned short int n = 0;
	unsigned short int type = 0;
	unsigned int seed = 0;
	bool disp = false;
	bool memo = false;
	bool nlsu = false;
	ifstream inp;
	cout << "Reading the input...";
	inp.open("output_game.txt");
	inp >> n >> type >> seed >> disp >> memo >> nlsu;
	inp.close();
	cout << "done!" << endl;
	if (seed == 0)
		seed = GetTickCount();
	srand(seed);
	unsigned int s = pow(2, n) - 2;
	vector<double> x(n, 0);
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	double prec = pow(10, -6);
	vector<double> excess(s, 0);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	unsigned short int iter = 0;
	unsigned int piv = 0;
	unsigned int sr = 0;
	double t = 0;
	vector<double> v(s + 1, 0);
	cout << "Generating game...";
	if (type == 1)
		type1(v, s, n);
	else if (type == 2)
		type2(v, s, n);
	else if (type == 4)
		type4(v, s, n);
	cout << "done!" << endl;
	ofstream res;
	res.open("v.txt", ofstream::out | ofstream::trunc);
	for (unsigned int i = 0; i < s + 1; i++)
		res << v[i] << endl;
	res.close();
	return 0;
}