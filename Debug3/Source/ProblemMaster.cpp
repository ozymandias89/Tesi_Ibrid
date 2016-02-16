/*
 * ProblemMaster.cpp
 *
 *  Created on: 26 gen 2016
 *      Author: riccardo
 */

#include "ProblemMaster.h"

ProblemMaster::ProblemMaster(bool verbose) {
	this->verbose = verbose;

}

ProblemMaster::~ProblemMaster() {
	// TODO Auto-generated destructor stub
}

void ProblemMaster::setupLP(CEnv env, Prob lp) {

	{	// variables
		static const char* varType = NULL;
		double obj = 0.0;
		double lb = 0.0;
		double ub = CPX_INFBOUND;

		for (int i = 0; i < N; i++) {
			obj = c[i];
			snprintf(name, NAME_SIZE, "x_%i", i);
			char* varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
					&varName);

		}

	}

	// constraints

	{
		std::vector<int> idx;
		std::vector<double> coef;

		for (int i = 0; i < num_constraint; i++) {
			char sense = 'E';
			int matbeg = 0;
			double rhs = b[i];
			int nzcnt = 0;

			for (int iter = 0; iter < N; iter++) {

				if (A[i][iter] != 0) {
					idx.push_back(iter);
					coef.push_back(A[i][iter]);
					nzcnt++;
				}

			}

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}
	}

}

void ProblemMaster::remove_constraint(CEnv env, Prob lp, int constraint) {

	if (verbose) {
		print_matrix();
		print_vect_c();
		print_vect_b();
		cout << "Delete redundant constraint number: " << constraint << endl;
	}

	CHECKED_CPX_CALL(CPXdelrows, env, lp, constraint, constraint);

	if (constraint < Num_original_constraints) {
		b.erase(b.begin() + constraint);
		A[constraint].clear();
		A.erase(A.begin() + constraint);

		Num_original_constraints--;
		num_constraint--;

	} else {

		int column = (Num_original_variables + constraint
				- Num_original_constraints);
		CHECKED_CPX_CALL(CPXdelcols, env, lp, column, column);

		b.erase(b.begin() + constraint);

		for (unsigned int i = 0; i < A.size(); ++i) {
			A[i].erase(A[i].begin() + column);
		}

		A[constraint].clear();
		A.erase(A.begin() + constraint);

		c.erase(c.begin() + column);

		num_constraint--;
		N--;
	}

	if (verbose) {
		print_matrix();
		print_vect_c();
		print_vect_b();
	}

}

void ProblemMaster::step1(CEnv env, Prob lp) {

	if (verbose) {
		cout << endl;
		cout << "STEP 1:" << endl;
	}

	CHECKED_CPX_CALL(CPXlpopt, env, lp);
	int stat = CPXgetstat(env, lp);

	if (verbose)
		cout << endl << "Status problem " << stat << endl;

	//check feasible
	//----------------------------------------------------------------
	// Remove redundant constraint (if exist) else STOP CONDITION 1
	//----------------------------------------------------------------
	if (stat != CPX_STAT_INFEASIBLE) {
		bool flag_redundant;
		do {
			flag_redundant = false;
			int num_cols = CPXgetnumcols(env, lp);
			int num_rows = CPXgetnumrows(env, lp);
			double redlb[num_cols], redub[num_cols];
			for (int i = 0; i < num_cols; ++i) {
				redlb[i] = 0;
				redub[i] = CPX_INFBOUND;
			}
			int rstat[num_rows];
			CHECKED_CPX_CALL(CPXbasicpresolve, env, lp, redlb, redub,
					&rstat[0]);

			int i = 0;
			while (i < num_rows) {
				if (rstat[i] == -1) {
					flag_redundant = true;
					break;
				}
				i++;
			}

			//remove constraint
			if (flag_redundant) {
				if (verbose)
					cout << "Constraint index " << i << " is redundant "
							<< endl;
				remove_constraint(env, lp, i);
			} else if (!flag_redundant && verbose) {
				cout << "No detect redundant constraints. " << endl;
				cout << endl;
			}

		} while (flag_redundant);

	} else {
		cout << endl;
		cout << " Iteration number: " << iter << endl;
		throw std::runtime_error(" STOP CONDITION STEP 1 ");
	}

}

void ProblemMaster::print_objval(CEnv env, Prob lp) {

	double objval;
	CHECKED_CPX_CALL(CPXgetobjval, env, lp, &objval);
	std::cout << endl << "Obj val: " << objval << std::endl;
}

void ProblemMaster::set_var_P(CEnv env, Prob lp) {

	int cur_numcols = CPXgetnumcols(env, lp);

	varVals.clear();
	varVals.resize(cur_numcols);
	CHECKED_CPX_CALL(CPXgetx, env, lp, &varVals[0], 0, cur_numcols - 1);

}

void ProblemMaster::print_var_P(CEnv env, Prob lp) {

	cout << "PRIMAL VARIABLES: " << endl;
	int cur_numcols = CPXgetnumcols(env, lp);
	int surplus;
	status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, 0,
			cur_numcols - 1);
	int cur_colnamespace = -surplus; // the space needed to save the names

	// allocate memory
	char** cur_colname = (char **) malloc(sizeof(char *) * cur_numcols);
	char* cur_colnamestore = (char *) malloc(cur_colnamespace);

	// get the names
	CPXgetcolname(env, lp, cur_colname, cur_colnamestore, cur_colnamespace,
			&surplus, 0, cur_numcols - 1);

	//  print index, name and value of each column
	for (int i = 0; i < cur_numcols; i++)
		cout << cur_colname[i] << " = " << varVals[i] << endl;

	// free
	free(cur_colname);
	free(cur_colnamestore);

}

/**
 Method that set and print dual variables
 @param  (CEnv env, Prob lp, bool prob,  bool verbose), environment of the problem,
 problem , flag(true P1 problem, false P2 problem) and bool verbose
 @return void
 */
void ProblemMaster::set_and_print_var_D(CEnv env, Prob lp, bool prob) {

	if (verbose)
		cout << endl;
	int num_rows = CPXgetnumrows(env, lp);

	if (prob) {
		dual_varVals_P1.clear();
		dual_varVals_P1.resize(num_rows);
		CHECKED_CPX_CALL(CPXgetpi, env, lp, &dual_varVals_P1[0], 0,
				num_rows - 1);

		if (verbose) {
			cout << "DUAL VARIABLES (the last is u_0): " << endl;
			for (int i = 0; i < num_rows; i++) {
				cout << dual_varVals_P1[i] << " ";
			}
			cout << endl;
		}

	} else {
		dual_varVals_P2.clear();
		dual_varVals_P2.resize(num_rows);
		CHECKED_CPX_CALL(CPXgetpi, env, lp, &dual_varVals_P2[0], 0,
				num_rows - 1);
		if (verbose) {
			cout << "DUAL VARIABLES (the last is v_0): " << endl;
			for (int i = 0; i < num_rows; i++) {
				cout << dual_varVals_P2[i] << " ";
			}
			cout << endl;
		}
	}

}

int ProblemMaster::select_fractionar_var() {

	// Selects variable with maximal fractionary value
	int index = -1;
	double max_fractionary = 10e-6;

	for (int i = 0; i < Num_original_variables; i++) {
		const double value = varVals[i], fractionary = fabs(
				value - round(value));

		if (fractionary > max_fractionary) {
			max_fractionary = fractionary;
			index = i;
		}
	}

	return index;
}

void ProblemMaster::create_P1_prob(CEnv env, Prob lp, int index) {

	if (verbose)
		cout << endl << "SUB_PROBLEM P1" << endl;

	double rhs = floor(varVals[index]);

	char sense = 'L';
	int matbeg = 0;
	const int idx = index;
	const double coef = 1;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
			&coef, 0, 0);

	if (verbose)
		cout << "insert inequality x_" << index << " <= " << rhs << endl;

}

void ProblemMaster::create_P2_prob(CEnv env, Prob lp, int index) {

	if (verbose) {
		cout << endl;
		cout << "SUB_PROBLEM P2" << endl;
	}

	double rhs = floor(varVals[index]) + 1;
	char sense = 'G';
	int matbeg = 0;
	const int idx = index;
	const double coef = 1;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 1, &rhs, &sense, &matbeg, &idx,
			&coef, 0, 0);

	if (verbose)
		cout << "insert inequality x_" << index << " >= " << rhs << endl;

}

double* ProblemMaster::solve_P1_Problem(CEnv env, Prob lp, int index) {

	static double z[2];
	z[0] = CPX_INFBOUND;
	z[1] = CPX_INFBOUND;
	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	int stat = CPXgetstat(env, lp);
	if (verbose)
		cout << endl << "Status problem " << stat << endl;

	int cur_numrows = CPXgetnumrows(env, lp);

// print and set solution and create and resolve P_2 problem"
	if (stat != CPX_STAT_INFEASIBLE) {
		if (verbose)
			cout << "FEASIBLE " << endl;
		gam = floor(varVals[index]);
		print_objval(env, lp);
		set_var_P(env, lp);
		print_var_P(env, lp);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z[0]);
		set_and_print_var_D(env, lp, true);
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		if (verbose)
			cout << "delete last inequality " << endl;

		create_P2_prob(env, lp, index);

		z[1] = solve_P2_Problem(env, lp, index);

	} else {
		if (verbose)
			cout << "No solution for P1 problem exists.. " << endl;

		// add Slack variables
		static const char* varType = NULL;
		double obj = 0.0;
		double lb = 0.0;
		double ub = CPX_INFBOUND;
		snprintf(name, NAME_SIZE, "S_%i", slack);
		slack++;
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		//add cut
		std::vector<int> idx;
		std::vector<double> coef;

		double rhs = floor(varVals[index]) + 1;
		char sense = 'E';
		int matbeg = 0;

		// x_k
		idx.push_back(index);
		coef.push_back(1);

		// S
		idx.push_back(N);
		coef.push_back(-1);

		num_constraint++;
		N++;

		//add 0 to c
		c.push_back(0);

		//extend A matrix
		A.resize(num_constraint);
		for (int i = 0; i < num_constraint; i++)
			A[i].resize(N);

		for (int i = 0; i < N; i++) {
			if (i == index) {
				A[(num_constraint - 1)][i] = 1;
			} else if (i == N - 1) {
				A[(num_constraint - 1)][i] = -1;
			} else
				A[(num_constraint - 1)][i] = 0;
		}

		if (verbose)
			cout << "delete last inequality " << endl;
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 2, &rhs, &sense, &matbeg,
				&idx[0], &coef[0], 0, 0);

		b.push_back(rhs);

		if (verbose) {
			cout << "Resolve a new problem P1.. " << endl;
			cout << "add inequality x_" << index << " >= " << rhs << endl;
		}

		double sum = 0;
		sum += int_var[index] - rhs;

		cout << "SLACK: " << sum << endl;
		if (sum < -1.e-4) {
			cout << "Negative Slack!" << endl;
			exit(1);
		}


		solve(env, lp);

	}

	return z;

}

double ProblemMaster::solve_P2_Problem(CEnv env, Prob lp, int index) {

	double z = CPX_INFBOUND;

	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	int stat = CPXgetstat(env, lp);
	if (verbose)
		cout << endl << "Status problem " << stat << endl;

	int cur_numrows = CPXgetnumrows(env, lp);

	if (stat != CPX_STAT_INFEASIBLE) {
		if (verbose)
			cout << "FEASIBLE " << endl;
		print_objval(env, lp);
		set_var_P(env, lp);
		print_var_P(env, lp);
		CHECKED_CPX_CALL(CPXgetobjval, env, lp, &z);
		set_and_print_var_D(env, lp, false);
		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		if (verbose)
			cout << "delete last inequality " << endl;

	} else {
		if (verbose)
			cout << "No solution for P2 problem exists.. " << endl;

		// add Slack variables
		static const char* varType = NULL;
		double obj = 0.0;
		double lb = 0.0;
		double ub = CPX_INFBOUND;
		snprintf(name, NAME_SIZE, "S_%i", slack);
		slack++;
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		//add cut
		std::vector<int> idx;
		std::vector<double> coef;

		double rhs = floor(varVals[index]);
		char sense = 'E';
		int matbeg = 0;

		//x_k
		idx.push_back(index);
		coef.push_back(1);

		//S
		idx.push_back(N);
		coef.push_back(1);

		num_constraint++;
		N++;

		//add 0 to c
		c.push_back(0);

		//extend A matrix
		A.resize(num_constraint);
		for (int i = 0; i < num_constraint; i++)
			A[i].resize(N);

		for (int i = 0; i < N; i++) {
			if (i == index) {
				A[(num_constraint - 1)][i] = 1;
			} else if (i == N - 1) {
				A[(num_constraint - 1)][i] = 1;
			} else
				A[(num_constraint - 1)][i] = 0;
		}

		cout << "delete last inequality " << endl;

		CHECKED_CPX_CALL(CPXdelrows, env, lp, cur_numrows - 1, cur_numrows - 1);
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, 2, &rhs, &sense, &matbeg,
				&idx[0], &coef[0], 0, 0);

		if (verbose) {
			cout << "Resolve a new problem P2.. " << endl;
			cout << "add inequality x_" << index << " <= " << rhs << endl;
			cout << "Restart from step 1 with new problem master: " << endl;
		}

		b.push_back(rhs);

		double sum = 0;
		sum += rhs - int_var[index];

		if (sum < -1.e-4) {
			cout << "Negative Slack!" << endl;
			exit(1);
		}


		solve(env, lp);
	}

	return z;
}

void ProblemMaster::solve(CEnv env, Prob lp) {

	solve_integer_problem(env, lp, true);
	// --------------------------------------------------
	// 2. solve linear problem
	// --------------------------------------------------
	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	int stat = CPXgetstat(env, lp);

	// --------------------------------------------------
	// 3. STOP CONDITION
	// --------------------------------------------------
	if (stat == CPX_STAT_UNBOUNDED || stat == CPX_STAT_INFEASIBLE) {
		cout << endl << " STOP CONDITION STEP 3 " << endl;
		cout << " Iteration number: " << iter << endl;
		throw std::runtime_error("  STOP CONDITION STEP 3 ");
	}

	cout << endl << "PROBLEM MASTER:" << endl;

	// --------------------------------------------------
	// 4. print solution
	// --------------------------------------------------
	print_objval(env, lp);

	// --------------------------------------------------
	// 5. set number and value of variable
	//    (cur_numcols,varVals) and print these
	// --------------------------------------------------
	set_var_P(env, lp);
	print_var_P(env, lp);

	// --------------------------------------------------
	// 6. chose the best fractional variable
	// --------------------------------------------------
	int index = select_fractionar_var();

	// --------------------------------------------------------
	// 7. if x solution aren't integer create P1 and P2 problem
	// --------------------------------------------------------
	if (index != -1) {

		if (verbose) {
			cout << endl << "More fractional variable choose " << varVals[index]
					<< endl;

			cout << "Index of variable choose: " << index << endl;
		}

		//create problem P_1
		create_P1_prob(env, lp, index);

		// --------------------------------------------------------
		// 8. solve sub_problems (P_1 and P_2) return min solution
		// --------------------------------------------------------
		double* z = solve_P1_Problem(env, lp, index);

		/////////////////////////////////////////////////////

		// ------------------------------------------------
		// 9. only if both problems have solution else get
		//		the best solution and stop
		// ------------------------------------------------
		if (*z < CPX_INFBOUND && *(z + 1) < CPX_INFBOUND && flag_find) {
			flag_find = false;

			min_sol = std::min(*z, *(z + 1));
			k = index;
		}

	} else {
		cout << " Iteration number: " << iter << endl;
		throw std::runtime_error(
				" The last solution is the best integer solution. STOP CONDITION STEP 4 ");
	}

}


void ProblemMaster::add_constraint_R(CEnv env, Prob lp,
		std::set<std::vector<double> > R, long aggressivity) {

	//change sign matrix A
	change_sign_A();

	static const char* varType = NULL;
	double obj = 0.0;
	double lb = 0.0;
	double ub = CPX_INFBOUND;
	char sense = 'E';
	int matbeg = 0;
	int nzcnt = 0;

	std::vector<int> idx;
	std::vector<double> coef;

	std::vector<double> y_tilde;
	int iterator = N;

	double sum = 0;
	//for each element in a set insert a new constraint
	for (std::set<std::vector<double> >::iterator it = R.begin(); it != R.end();
			++it) {

		y_tilde = *it;

		double M = 0;
		double m = CPX_INFBOUND;

		for (int i = 0; i < iterator; i++) {

			if (fabs(y_tilde[i]) != 0) {
				if ((fabs(y_tilde[i])) > M) {
					M = (fabs(y_tilde[i]));
				}
				if ((fabs(y_tilde[i])) < m) {
					m = (fabs(y_tilde[i]));
				}

			}
		}

		double aggress = (M/m);

		if (verbose)
		cout<< endl << "AGGRESSIVITY: " << M / m << endl;

		if ( aggress < aggressivity ){
			//cout<< endl << "AGGRESSIVITY: " << M / m << "AGGIUNTO" << endl;

		sum = 0;
		nzcnt = 0;

		// add Slack variables
		snprintf(name, NAME_SIZE, "S_%i", slack);
		slack++;
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		//-S
		coef.push_back(-1);
		idx.push_back(N);
		nzcnt++;

		//a_T * x
		for (int i = 0; i < iterator; i++) {

			if (y_tilde[i] != 0) {
				sum += (y_tilde[i] * int_var[i]);
				//	cout << "###### " << y_tilde[i] << "   " << int_var[i]  << endl;

				coef.push_back(y_tilde[i]);
				idx.push_back(i);
				nzcnt++;
			}
		}

		//beta
		sum = -sum;
		double rhs = y_tilde[iterator];
		sum += y_tilde[iterator];

		sum = -sum;

		cout << "SLACK: " << sum << endl;
		if (sum < -1.e-4){
			cout << "Negative Slack!" << endl;
			exit(1);
		}

		num_constraint++;
		N++;

		//add 0 to c
		c.push_back(0);

		//extend A matrix
		A.resize(num_constraint);
		for (int i = 0; i < num_constraint; i++)
			A[i].resize(N);

		for (int i = 0; i < iterator; i++)
			A[(num_constraint - 1)][i] = y_tilde[i];

		A[(num_constraint - 1)][N - 1] = -1;

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
				&matbeg, &idx[0], &coef[0], 0, 0);

		b.push_back(rhs);

		idx.clear();
		coef.clear();
	}
	}

}

void ProblemMaster::set_integer_variable(CEnv env, Prob lp, bool verbose) {

	if (verbose)
		cout << "PRIMAL VARIABLES INTEGER: " << endl;
	int cur_numcols = CPXgetnumcols(env, lp);

	int_var.clear();
	int_var.resize(cur_numcols);
	CHECKED_CPX_CALL(CPXgetx, env, lp, &int_var[0], 0, cur_numcols - 1);

	int surplus;
	status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, 0,
			cur_numcols - 1);
	int cur_colnamespace = -surplus; // the space needed to save the names

	// allocate memory
	char** cur_colname = (char **) malloc(sizeof(char *) * cur_numcols);
	char* cur_colnamestore = (char *) malloc(cur_colnamespace);

	// get the names
	CPXgetcolname(env, lp, cur_colname, cur_colnamestore, cur_colnamespace,
			&surplus, 0, cur_numcols - 1);

	if (verbose) {
		//  print index, name and value of each column
		for (int i = 0; i < cur_numcols; i++) {
			cout << cur_colname[i] << " = " << int_var[i] << endl;
		}
	}
	// free
	free(cur_colname);
	free(cur_colnamestore);
}

void ProblemMaster::solve_integer_problem(CEnv env, Prob lp, bool verbose) {

	//change problem to integer
	CHECKED_CPX_CALL(CPXchgprobtype, env, lp, CPXPROB_MILP);

	//set varibles to integer
	char ctype[N];
	for (int i = 0; i < N; i++) {
		if (i < Num_original_variables) {
			ctype[i] = CPX_INTEGER;
		} else
			ctype[i] = CPX_CONTINUOUS;
	}
	//const_cast<char *>(ctype);
	CHECKED_CPX_CALL(CPXcopyctype, env, lp, ctype);

	CHECKED_CPX_CALL(CPXmipopt, env, lp);
	int stat = CPXgetstat(env, lp);

	if (stat != CPXMIP_INFEASIBLE) {
		if (verbose) {
			cout << endl << "Problem solved to integer: " << endl;
			print_objval(env, lp);
			double risult;
			CHECKED_CPX_CALL(CPXgetobjval, env, lp, &risult);
			if (integer == -CPX_INFBOUND)
				integer = risult;
			else if (fabs(integer - risult) > 1.e-3L) {
				cerr << "We have lose one integer solution: " << endl;
				cout << "Number of constraints added: "
						<< CPXgetnumrows(env, lp) - Num_original_constraints
						<< endl;
				exit(1);
			}
			set_integer_variable(env, lp, verbose);
		}

	} else {
		cout << endl;
		cout << " Integer problem not resolvable! " << endl;
		cout << " Iteration number: " << iter << endl;
		CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);
		// free allocate memory
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(const_cast<cpxenv **>(&env));
		exit(0);
	}

	CHECKED_CPX_CALL(CPXchgprobtype, env, lp, CPXPROB_LP);
}

void ProblemMaster::save(CEnv env, Prob lp) {

	CHECKED_CPX_CALL(CPXwriteprob, env, lp, "../data/problem.lp", 0);

}
