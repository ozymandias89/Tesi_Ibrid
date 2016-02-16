/*
 * ProblemMaster.h
 *
 *  Created on: 26 gen 2016
 *      Author: riccardo
 */

#ifndef SOURCE_PROBLEMMASTER_H_
#define SOURCE_PROBLEMMASTER_H_

#include "load.h"

class ProblemMaster {
public:
	ProblemMaster(bool verbose = false);
	virtual ~ProblemMaster();


// PREDICATES:

	/**
	 Method that set the primary problem
	 @param  (CEnv env, Prob lp)
	 @return void
	 */
	void setupLP(CEnv env, Prob lp);
	void remove_constraint(CEnv env, Prob lp, int constraint);
	void step1(CEnv env, Prob lp);
	/**
	 Method that print object function
	 @param  (CEnv env, Prob lp, bool verbose), environment of the problem and problem and verbose
	 @return void
	 */
	void print_objval(CEnv env, Prob lp);
	/**
	 Method that set and print primal variable
	 @param  (CEnv env, Prob lp, bool verbose), environmant of the problem, problem and verbose
	 @return void
	 */
	void set_var_P(CEnv env, Prob lp);
	void print_var_P(CEnv env, Prob lp);
	void set_and_print_var_D(CEnv env, Prob lp, bool prob);
	/**
	 Method that chooses the best x fractional variable
	 @param  (vector<double>)
	 @return int, return index of higher variable functional (-1 if no variable is fractional)
	 */
	int select_fractionar_var();

	void create_P1_prob(CEnv env, Prob lp, int index);
	void create_P2_prob(CEnv env, Prob lp, int index);

	double* solve_P1_Problem(CEnv env, Prob lp, int index);
	double solve_P2_Problem(CEnv env, Prob lp, int index);
	void solve(CEnv env, Prob lp);
	/**
	 Method that add constraints in set R to the first problem.
	 @param  (CEnv env, Prob lp, std::set<std::vector<double> > R, long aggressivity) environment of the problem,
	 problem ,set R and aggressivity filter
	 @return void
	 */
	void add_constraint_R(CEnv env, Prob lp, std::set<std::vector<double> > R, long aggressivity);
	void set_integer_variable(CEnv env, Prob lp, bool verbose);
	void solve_integer_problem(CEnv env, Prob lp, bool verbose);
	void save (CEnv env, Prob lp);

// ATTRIBUTES:

bool verbose;
//primary variables
std::vector<double> varVals;

};

#endif /* SOURCE_PROBLEMMASTER_H_ */
