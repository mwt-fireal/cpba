#pragma once

#ifndef SCATTER_H
#define SCATTER_H

#include <string>
#include <map>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

using namespace std;

class Scatter
{
public:
	static gsl_matrix *PSTAR_W;
	static gsl_matrix *PSTAR_B;
	static gsl_matrix *PSTAR_A;
	static gsl_matrix *PSTAR_V;

	static map<double, double*> spcache;
	static map<string, double> slcache;

	//static double getScatPwr(double CTavg, double Ez);
	static double getSPratio(double CTavg, double Ez);
	static double scatLength(string mat);
	static double calcTdM(string mat, double Ez, double Eo);
};

#endif
