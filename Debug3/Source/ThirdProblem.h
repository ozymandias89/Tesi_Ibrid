/*
 * Thirdproblem.h
 *
 *  Created on: 14 nov 2015
 *      Author: riccardo
 */

#ifndef SOURCE_THIRDPROBLEM_H_
#define SOURCE_THIRDPROBLEM_H_

// Includes:
#include "load.h"

class ThirdProblem {
public:
	ThirdProblem(vector<double> y_til, vector<double> cost,
			bool verbose = false);
	virtual ~ThirdProblem();

	// PREDICATES:

	/**
	 Print vector
	 @param  (vector<template T> vector)
	 @return void
	 */
	template<typename T>
	void print_vector(vector<T> vector);

	/**
	 Calculate y bar minimum y tilde
	 @param  (vector<double> vector)
	 @return void
	 */
	void y_bar_MIN_y_tilde(vector<double> c);

	/**
	 Method that set the third problem
	 @param set<int> constraints
	 @return void
	 */
	void setup(set<int> constraints);

	/**
	 solve third problem and calculate lambda, return bool infeasible problem
	 @param  set<int> constraints
	 @return bool
	 */
	bool solve(set<int> constraints);

	/**
	 update y_bar step 8.4
	 @param  (CEnv env, Prob lp, vector<double>& c)
	 @return void
	 */
	void update_y_bar(vector<double>& c);

	// ATTRIBUTES:
	vector<double> y_tilde;
	vector<double> t;
	vector<double> result;
	vector<char> sense;
	double ub, lb;
	int constraint_to_add;
	double lambda;

	bool verbose;
};

#endif /* SOURCE_THIRDPROBLEM_H_ */
