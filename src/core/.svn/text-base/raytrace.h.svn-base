#pragma once

using namespace std;

#include "cpba.h"
#include "util.h"
#include "data.h"
#include "interp.h"

#include <gsl/gsl_vector.h>
#include <map>
#include <string>

#ifndef RAYTRACE_H
#define RAYTRACE_H

class zROI
{
public:
	int M;
	double x1, x2, z1, z2;
	gsl_vector *L;
	gsl_vector *ZL;
	
	zROI(int dim)
	{
		this->M = 0;
		this->L = gsl_vector_calloc(dim);
		this->ZL = gsl_vector_calloc(dim);
	}
};

class RayTrace
{
public:
	static void trace();

	static double range(double Energy, string material);
	static bool inROI(gsl_vector *roiX, gsl_vector *roiZ, double xo, double zo);
	static void zoLineTest(zROI *rz, double xo, double zo, int zLineSign);
};

#endif
