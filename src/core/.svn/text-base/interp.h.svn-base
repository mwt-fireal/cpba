#pragma once

// INTERP_H might actually be defined somewhere.
#ifndef CPBA_INTERP_H
#define CPBA_INTERP_H

using namespace std;

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_bspline.h>

#include "util.h"

class BilinearInterpolation
{
public:
	static vector<gsl_matrix *> meshgrid(gsl_vector *x, gsl_vector *z);
	static void get_neigbour_indices(gsl_vector* inArr,double x ,int& lowerX ,int& upperX);
	static double evalulate_2d_constant_bilinear(gsl_vector* xAxis,gsl_vector* yAxis, gsl_matrix* zSurface, double xcoord ,double ycoord);
	static gsl_matrix *evalulate_2d_bilinear(gsl_vector* xAxis,gsl_vector* yAxis, gsl_matrix* zSurface, gsl_matrix *xi, gsl_matrix *yi);
};

class CubicSplineInterpolation
{
public:
	
	static double evaluate_constant_spline(int vec_size, double *x, double *y, double xi);
	static double evaluate_constant_spline(int vec_size, gsl_vector *x, gsl_vector *y, double xi);

	static void evaluate_spline(gsl_vector *yi, double *x, double *y, double *xi);
	static void evaluate_spline(gsl_vector *yi, gsl_vector *x, gsl_vector *y, gsl_vector *xi);
};

class LinearInterpolation
{
public:
	
	static double evaluate_constant_linear(int vec_size, double *x, double *y, double xi);
	static double evaluate_constant_linear(int vec_size, gsl_vector *x, gsl_vector *y, double xi);

	static void evaluate_linear(gsl_vector *yi, double *x, double *y, double *xi);
	static void evaluate_linear(gsl_vector *yi, gsl_vector *x, gsl_vector *y, gsl_vector *xi);
};

#endif
