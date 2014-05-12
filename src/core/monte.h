#pragma once

#ifndef MONTE_H
#define MONTE_H

using namespace std;

#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <boost/format.hpp>

class MonteCarlo
{
public:
	int xDim;
	int zDim;
	gsl_vector *x;
	gsl_vector *z;

	gsl_matrix *D;

	MonteCarlo(int Eo, int FS, int ThetaDeg, int Hstep);

private:
	void init(int Eo, int FS, int ThetaDeg, int Hstep);
};

#endif
