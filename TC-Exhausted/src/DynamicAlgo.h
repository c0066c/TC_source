/*
 * DynamicAlgo.h
 *
 *  Created on: Sep 12, 2012
 *      Author: anas
 */
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

#ifndef DYNAMICALGO_H_
#define DYNAMICALGO_H_

class DynamicAlgo {
public:
	// Constructor
	DynamicAlgo(float alpha,double alphaScale,double RScale,int D,int sigma,int delta);

	// Destructor
	~DynamicAlgo();

	/*
	// Returns the minimum value of Ri
	double minR(int i);

	// Returns the maximum value of Ri
	double maxR(int i);
	*/

	// Calculate Rmin(i) and Rmax(i) for all i values and store them into minR and maxR arrays
	// Note: Rmax is stored in maxR[0]
	void minMaxR();

	// Generates (or reads) the values of R(i,j), C(i,j,e), P(i,j,e)
	void getData();

	// Retruns the integration
	double integ(int j);

	// Retruns the sum of reliability levels for the functions from u to n with versions predetermined by the user
	double rho(int u);

	// Calculates Hj (for j=1,2,...,K) for the current values of i,t and r, returns the maximum value
	// and stores the index of the minimum value into minInd variable
	double minH();

	/*
	// Stores the index of the version that achieves the minimum expected RT penalty into minInd
	// for the current values of i(n and 1),t and r, and returns the maximum value
	double minJ();
	*/

	// Fills the dynamic programming table
	double fillTable();
	int special_fun0;
	int special_fun2;
	int special_fun;
	int special_fun3;
	int special_fun4;

private:

	float alpha;
	int n,i,j,t,sigma,delta,K,D,minInd,Rmax,Dmax; //K:max no. of versions, minInd:the index of the
	double ***G,***C,***P,**R, *H,*minR,*maxR,alphaScale,RScale,tableScale,r;
	int *theta,***jStar,*size; // size: number of versions for each function

	int DEBUG;



};

#endif /* DYNAMICALGO_H_ */
