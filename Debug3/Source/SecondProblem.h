/*
 * SecondProblem.h
 *
 *  Created on: 04 nov 2015
 *      Author: riccardo
 */

#ifndef SOURCE_SECONDPROBLEM_H_
#define SOURCE_SECONDPROBLEM_H_

// Includes:
#include "load.h"

using namespace std;

class SecondProblem {
public:
	SecondProblem(bool verbose = false);
	virtual ~SecondProblem();

	// PREDICATES:
	void set_up_variable(CEnv env, Prob lp);

	/**
	 Method that set the second problem
	 @param  (CEnv env, Prob lp)
	 @return void
	 */
	void setupSP(CEnv env, Prob lp);

	/** Method that evaluate the vector r as sum of rows of C (matrix of dual problem (13-19) + 21) and set r vector
	 @param  bool verbose
	 @return none
	 */
	void evaluate_rT();

	/**
	 Auxiliar method, set variables result from second problem and insert y tilde in set R
	 @param  (CEnv env, Prob lp)
	 @return none
	 */
	void set_solution(CEnv env, Prob lp);

	/**
	 Method that solve the second problem (13-19) + 21)
	 @param  CEnv env, Prob lp, bool verbose
	 @return bool infeasibility of problem
	 */
	bool solve(CEnv env, Prob lp);

	void add_constraint(CEnv env, Prob lp, int constraint_to_add);

	/**
	 Method that creates tight C'y >= d subsystem and adds to dual problem
	 @param  CEnv env, Prob lp, bool verbose
	 @return void
	 */
	void step8_1(CEnv env, Prob lp);

	/**
	 Method that add constraint r_T * y = r_T * y
	 @param  CEnv env, Prob lp, bool verbose
	 @return void
	 */
	void step8_2(CEnv env, Prob lp);

	/**
	 Test if y_tilde == y_barr
	 @param  none
	 @return bool
	 */
	bool y_tilde_EQ_y_bar();

	void save (CEnv env, Prob lp);

	/**
	 Method print u
	 @param  none
	 @return none
	 */
	void print_u();
	/**
	 Method print v
	 @param  none
	 @return none
	 */
	void print_v();
	/**
	 Method print r
	 @param  none
	 @return none
	 */
	void print_r();
	void print_c();
	/**
	 Method print a
	 @param  none
	 @return none
	 */
	void print_a();
	/**
	 Method print u0
	 @param  none
	 @return none
	 */
	void print_u0();
	/**
	 Method print v0
	 @param  none
	 @return none
	 */
	void print_v0();
	/**
	 Method print beta
	 @param  none
	 @return none
	 */
	void print_beta();
	/**
	 Method print y_tilde
	 @param  none
	 @return none
	 */
	void print_y_tilde();

	void print_y_bar();

	// ATTRIBUTES:
	set<vector<double> > R;
	set<int> satisfy_constraint_list;
	vector<double> cost;
	vector<double> u;
	vector<double> v;
	vector<double> a;
	double u0;
	double v0;
	double beta;
	vector<double> rt;
	vector<double> y_tilde;
	bool verbose;

};

#endif /* SOURCE_SECONDPROBLEM_H_ */
