#include "interp.h"
#include "util.h"
#include "out.h"

using namespace std;

gsl_matrix *BilinearInterpolation::evalulate_2d_bilinear(gsl_vector* xAxis,gsl_vector* yAxis, gsl_matrix* zSurface, gsl_matrix *xi, gsl_matrix *yi)
{
	gsl_matrix *xyim = gsl_matrix_calloc(xi->size1, xi->size2);

	unsigned int i, j;
	for(i = 0; i < xi->size1; i++)
	{
		for(j = 0; j < xi->size2; j++)
		{
			double xc = gsl_matrix_get(xi, i, j);
			double yc = gsl_matrix_get(yi, i, j);
			double xyi = BilinearInterpolation::evalulate_2d_constant_bilinear(xAxis, yAxis, zSurface, xc, yc);

			//don't allow values less than zero.
			xyi = xyi < 0 ? 0 : xyi;

			gsl_matrix_set(xyim, i, j, xyi);
		}
	}

	return xyim;
}

void BilinearInterpolation::get_neigbour_indices(gsl_vector* inArr, double x, int& lowerX, int& upperX)
{
	int N = inArr->size;
	if (x <= gsl_vector_get(inArr,0))
	{
		lowerX = 1;
		upperX = 1;
	}
	else if (x >= gsl_vector_get(inArr,N-1)) 
	{
		lowerX = N;
		upperX = N;
	}
	else
	{
		for (int i = 2; i<=N; i++)
		{
			if (x < gsl_vector_get(inArr,i-1)) 
			{
				lowerX = i - 1;
				upperX = i;
				break;
			}
			else if (x == gsl_vector_get(inArr,i-1))
			{
				lowerX = i;
				upperX = i;
				break;
			}
		}
	}
}

double BilinearInterpolation::evalulate_2d_constant_bilinear(gsl_vector* xAxis, gsl_vector* yAxis, gsl_matrix* zSurface, double xcoord ,double ycoord)
{
	double fxy;

	int nx = xAxis->size;
	int ny = yAxis->size;
	int lx;
	int ux;

	BilinearInterpolation::get_neigbour_indices(xAxis, xcoord, lx, ux);

	int ly;
	int uy;

	BilinearInterpolation::get_neigbour_indices( yAxis, ycoord, ly, uy);

	double fQ11 = gsl_matrix_get(zSurface,lx-1, ly-1);
	double fQ21 = gsl_matrix_get(zSurface,ux-1, ly-1);
	double fQ12 = gsl_matrix_get(zSurface,lx-1, uy-1);
	double fQ22 = gsl_matrix_get(zSurface,ux-1, uy-1);

	//if point exactly found on a node do not interpolate
	if ((lx == ux) && (ly == uy))
	{
		return fQ11;
	}

	double x = xcoord;
	double y = ycoord;

	double x1 = gsl_vector_get(xAxis,lx-1);
	double x2 = gsl_vector_get(xAxis,ux-1);
	double y1 = gsl_vector_get(yAxis,ly-1);
	double y2 = gsl_vector_get(yAxis,uy-1);

	//if xcoord lies exactly on an xAxis node do linear interpolation
	if (lx == ux)
	{
		return fQ11 + (fQ12 - fQ11) * (y - y1) / (y2 - y1);
	}

	//if ycoord lies exactly on an xAxis node do linear interpolation
	
	if (ly == uy)
	{
		return fQ11 + (fQ22 - fQ11) * (x - x1) / (x2 - x1);
	}

	//perform interpolation.
	fxy = fQ11 * (x2 - x) * (y2 - y);
	fxy = fxy + fQ21 * (x - x1) * (y2 - y);
	fxy = fxy + fQ12 * (x2 - x) * (y - y1);
	fxy = fxy + fQ22 * (x - x1) * (y - y1);
	fxy = fxy / ((x2 - x1) * (y2 - y1));

	return fxy;
}

vector<gsl_matrix *> BilinearInterpolation::meshgrid(gsl_vector *x, gsl_vector *z)
{
	//we take the rows from the length of z, and the columns from the length of x.

	gsl_matrix *zz = gsl_matrix_calloc(z->size, x->size);
	gsl_matrix *xx = gsl_matrix_calloc(z->size, x->size);

	unsigned int i, j;
	for(i = 0; i < z->size; i++)
	{
		for(j = 0; j < x->size; j++)
		{
			gsl_matrix_set(xx, i, j, gsl_vector_get(x, j));
			gsl_matrix_set(zz, i, j, gsl_vector_get(z, i));
		}
	}

	vector<gsl_matrix *> c;

	c.push_back(xx);
	c.push_back(zz);

	return c;
}

double CubicSplineInterpolation::evaluate_constant_spline(int vec_size, gsl_vector *x, gsl_vector *y, double xi)
{
	double *xA = cutil::vec_to_array(x);
	double *yA = cutil::vec_to_array(y);

	double val = CubicSplineInterpolation::evaluate_constant_spline(vec_size, xA, yA, xi);

	/* clamp to zero if the result is erroneous */
	if(!cutil::is_finite(val))
	{
		val = 0;
	}

	free(xA);
	free(yA);

	return val;
}


double CubicSplineInterpolation::evaluate_constant_spline(int vec_size, double *x, double *y, double xi)
{
	gsl_interp_accel *acc;
    	gsl_spline *spline;

	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_cspline, vec_size);
	gsl_spline_init(spline, x, y, vec_size);

	double tvalue;
	
	gsl_spline_eval_e(spline, xi, acc, &tvalue);

	//don't allow negative values, NAN, or non-finite numbers.
	if(tvalue == GSL_NAN || tvalue < 0 || !cutil::is_finite(tvalue))
	{
		tvalue = 0;
	}

	gsl_interp_accel_free(acc);
	gsl_spline_free(spline);

	return tvalue;
}

void CubicSplineInterpolation::evaluate_spline(gsl_vector *yi, gsl_vector *x, gsl_vector *y, gsl_vector *xi)
{
	double *xA = cutil::vec_to_array(x);
	double *yA = cutil::vec_to_array(y);
	double *xiA = cutil::vec_to_array(xi);

	CubicSplineInterpolation::evaluate_spline(yi, xA, yA, xiA);

	free(xA);
	free(yA);
	free(xiA);
}

void CubicSplineInterpolation::evaluate_spline(gsl_vector *yi, double *x, double *y, double *xi)
{
	gsl_interp_accel *acc;
    gsl_spline *spline;

	int vec_size = sizeof(x) / sizeof(double);

	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_cspline, vec_size);
	gsl_spline_init(spline, x, y, vec_size);

	for(unsigned int j = 0; j < yi->size; j++)
	{
		double tvalue;
		gsl_spline_eval_e(spline, xi[j], acc, &tvalue);

		//don't allow negative values, NAN, or non-finite numbers.
		if(tvalue == GSL_NAN || tvalue < 0 || !cutil::is_finite(tvalue))
		{
			tvalue = 0;
		}

		gsl_vector_set(yi, j, tvalue);
	}

	gsl_interp_accel_free(acc);
	gsl_spline_free(spline);
}


double LinearInterpolation::evaluate_constant_linear(int vec_size, gsl_vector *x, gsl_vector *y, double xi)
{
	double *xA = cutil::vec_to_array(x);
	double *yA = cutil::vec_to_array(y);

	double val = LinearInterpolation::evaluate_constant_linear(vec_size, xA, yA, xi);

	free(xA);
	free(yA);

	return val;
}


double LinearInterpolation::evaluate_constant_linear(int vec_size, double *x, double *y, double xi)
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, vec_size);
	gsl_interp_init(interp, x, y, vec_size);

	double tvalue;
	
	gsl_interp_eval_e(interp, x, y, xi, acc, &tvalue);

	//don't allow negative values, NAN, or non-finite numbers.
	if(tvalue == GSL_NAN || tvalue < 0 || !cutil::is_finite(tvalue))
	{
		//Output::print(boost::format("invalid value found for evaluate_constant_linear, value is %1%") % tvalue);
		tvalue = 0;
	}

	gsl_interp_accel_free(acc);
	gsl_interp_free(interp);

	return tvalue;
}

void LinearInterpolation::evaluate_linear(gsl_vector *yi, gsl_vector *x, gsl_vector *y, gsl_vector *xi)
{
	double *xA = cutil::vec_to_array(x);
	double *yA = cutil::vec_to_array(y);
	double *xiA = cutil::vec_to_array(xi);

	LinearInterpolation::evaluate_linear(yi, xA, yA, xiA);

	free(xA);
	free(yA);
	free(xiA);
}

void LinearInterpolation::evaluate_linear(gsl_vector *yi, double *x, double *y, double *xi)
{
	int vec_size = sizeof(x) / sizeof(double);

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, vec_size);
	gsl_interp_init(interp, x, y, vec_size);

	for(unsigned int j = 0; j < yi->size; j++)
	{
		double tvalue;
		
		gsl_interp_eval_e(interp, x, y, xi[j], acc, &tvalue);

		//don't allow negative values, NAN, or non-finite numbers.
		if(tvalue == GSL_NAN || tvalue < 0 || !cutil::is_finite(tvalue))
		{
			tvalue = 0;
		}

		gsl_vector_set(yi, j, tvalue);
	}

	gsl_interp_accel_free(acc);
	gsl_interp_free(interp);
}
