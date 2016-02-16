#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

string generate_canonical_matrix(const int m, const int limit_number_generate) {

	srand (time(NULL));int number = rand() % 100 + 10;
	ostringstream convert;
	convert << number;

	string a = "../data/problem";
	string b = convert.str();
	string c = ".txt";
	string namefile = a + b + c;

	ofstream myfile(namefile.c_str());
	if (myfile.is_open()) {

		//cost coefficients
		for (int i = 0; i < m * 2; ++i)
		myfile << rand() % limit_number_generate << " ";

		myfile << endl << endl;

		//b coefficients
		for (int i = 0; i < m; ++i)
		myfile << rand() % 31 << " ";

		myfile << endl << endl;

		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < m * 2; ++j) {
				if (j < m) {
					(j == i) ? myfile << 1 << " " : myfile << 0 << " ";
				} else
				myfile << rand() % 41 +(-20) << " ";

			}
			myfile << endl;
		}

		myfile.close();
		return namefile;
	} else {
		cerr << "Unable to save file!" << endl;
		return 0;
	}

}
