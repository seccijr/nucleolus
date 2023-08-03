#include "gen_game.h"

int main(int argc, char *argv[])
{
	if (argc != 4)
	{
		return -1;
	}
	unsigned short int n = atoi(argv[1]);
	unsigned short int type = atoi(argv[2]);
	unsigned int seed = atoi(argv[3]);
	srand(seed);
	unsigned int s = pow(2, n) - 2;
	vector<double> x(n, 0);
	vector<double> v(s + 1, 0);
	if (type == 1)
		type1(v, s, n);
	else if (type == 2)
		type2(v, s, n);
	else if (type == 4)
		type4(v, s, n);
	for (unsigned int i = 0; i < s + 1; i++)
		cout << v[i] << endl;
	return 0;
}