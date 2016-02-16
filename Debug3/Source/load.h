/**
 * @file load.h
 * library for support
 *
 * @author  Riccardo Zanella, riccardozanella89@gmail.com
 * @version 2.0
 */

#ifndef SOURCE_LOAD_H_
#define SOURCE_LOAD_H_

// Includes:
#include <cstdio>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "cpxmacro.h"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

// error status and message buffer
int status;
char errmsg[BUF_SIZE];

const int NAME_SIZE = 512;
char name[NAME_SIZE];

const double epsilon_8_1 = 1.e-7L;
const double epsilon_8_2 = 1.e-7L;
const double epsilon_8_3 = 1.e-5L;
const double epsilon_8_4 = 1.e-6L;

//number of iteration of the problem
static int iter;

//flag for solve
static bool flag_find = true;

//number of variable
int N;

double integer = -CPX_INFBOUND;
//number of variable original
int Num_original_variables;

//number of variable original
int Num_original_constraints;

//index for e_k
int k;

//number of constraint
int num_constraint;

//slack_variable
int slack = 0;

//Coefficient cost
std::vector<double> c;

//price matrix
std::vector<std::vector<double> > A;

//known terms
std::vector<double> b;

//dual variables P_1 problem
std::vector<double> dual_varVals_P1;

//dual variables P_2 problem
std::vector<double> dual_varVals_P2;

//gamma
int gam;

//min_sol from P_1/P_2
double min_sol;


std::vector<double> int_var;

/**
 change sign matrix A
 @param  none
 @return void
 */
inline void change_sign_A() {
	for (unsigned int i = 0; i < A.size(); i++) {
		for (unsigned int j = 0; j < A[i].size(); j++) {
			if (A[i][j] != 0)
				A[i][j] = -A[i][j];
		}
	}
}

/**
 Method that load from file the problem (for example of format file see folder data)
 @param  (ifstream &) , object ifstream
 @return void
 */
void load_problem(ifstream &myfile) {

	bool flag_obb = true;
	bool flag_known = true;

	std::string line;
	std::string number;

	if (!myfile.is_open()) {
		cerr << "Unable to open file!!!!!" << endl;
		exit(1);
	}
	do {
		getline(myfile, line);
		//if aren't a comment or empty line
		if ((line[0] != '/' && line[1] != '/') && line.length() != 0) {

			std::istringstream is(line);

			int count = 0;
			while (is >> number) {
				count++;
			}

			N = count;
			Num_original_variables = N;

			c.resize(N);

			std::istringstream it(line);

			int i = 0;
			while (it >> number) {
				c[i] = atof(number.c_str());
				i++;
			}
			flag_obb = false;
		}

	} while (flag_obb);

	do {

		getline(myfile, line);
		//if aren't a comment or empty line
		if ((line[0] != '/' && line[1] != '/') && line.length() != 0) {

			std::istringstream is(line);
			int count = 0;
			while (is >> number) {
				count++;
			}

			num_constraint = count;
			Num_original_constraints = num_constraint;

			b.resize(num_constraint);

			std::istringstream it(line);

			int i = 0;
			while (it >> number) {

				b[i] = atof(number.c_str());
				i++;
			}

			flag_known = false;
		}
	} while (flag_known);

	//ALLOCATE MATRIX A
	A.resize(num_constraint);
	for (int i = 0; i < num_constraint; ++i) {
		A[i].resize(N);
	}

	int i = 0;
	int j = 0;
	while (getline(myfile, line)) {
		//if aren't a comment or empty line
		if ((line[0] != '/' && line[1] != '/') && line.length() != 0) {
			std::istringstream it(line);
			j = 0;
			if (i >= num_constraint) {
				cerr << "Matrix doesn't respect standard line";
				exit(1);
			}
			while (it >> number) {
				if (j >= N) {
					cerr << "Matrix doesn't respect standard columns";
					exit(1);
				}
				A[i][j] = atof(number.c_str());
				j++;

			}
			i++;

		}
	}

}

/**
 Method that print matrix A
 @param  none
 @return void
 */
void print_matrix() {

	cout << endl;
	cout << "Matrix A " << endl;

	for (unsigned int i = 0; i < A.size(); i++) {
		for (unsigned int j = 0; j < A[i].size(); j++) {

			cout << A[i][j] << " ";
		}
		cout << endl;
	}

}

/**
 Method that print vector_c
 @param  none
 @return void
 */
void print_vect_c() {

	cout << "Vector c" << endl;

	for (int i = 0; i < N; i++)
		cout << c[i] << " ";

	cout << endl;
}

/**
 Method that print vector_b
 @param  none
 @return void
 */
void print_vect_b() {

	cout << endl << "Vector b" << endl;
	for (int i = 0; i < num_constraint; i++)
		cout << b[i] << " ";

	cout << endl;
}

void print_u_variables() {
	cout << "u variables (the last is u_0): " << endl;
	for (std::vector<double>::const_iterator i = dual_varVals_P1.begin();
			i != dual_varVals_P1.end(); ++i)
		std::cout << *i << ' ';

	cout << endl;
}

void print_v_variables() {

	cout << "v variables (the last is v_0): " << endl;
	for (std::vector<double>::const_iterator i = dual_varVals_P2.begin();
			i != dual_varVals_P2.end(); ++i)
		std::cout << *i << ' ';

	cout << endl;
}




#endif /* SOURCE_LOAD_H_ */
