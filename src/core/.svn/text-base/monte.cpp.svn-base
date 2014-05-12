#include "monte.h"
#include "util.h"
#include "out.h"
#include "cpba.h"

#include <boost/format.hpp>

using namespace std;

MonteCarlo::MonteCarlo(int Eo, int FS, int ThetaDeg, int Hstep)
{
	this->init(Eo, FS, ThetaDeg, Hstep);
}

void MonteCarlo::init(int Eo, int FS, int ThetaDeg, int Hstep)
{
	int rows;
	FILE *fptr;
	Output::print("MonteCarlo::init() started");

	boost::format df("%1%/MC/Dose_E%2%_FS%3%_TILT%4%_STEP%5%.txt");
	df % CPBA::INPUT_PATH % Eo % FS % ThetaDeg % Hstep;

	boost::format xf("%5%/MC/x_E%1%_FS%2%_TILT%3%_STEP%4%.txt");
	xf % Eo % FS % ThetaDeg % Hstep % CPBA::INPUT_PATH;

	boost::format zf("%5%/MC/z_E%1%_FS%2%_TILT%3%_STEP%4%.txt");
	zf % Eo % FS % ThetaDeg % Hstep % CPBA::INPUT_PATH;

	int i;

	this->D = cutil::read_matrix_from_file(df);

	double a;

	fptr = fopen(xf.str().c_str(), "r");
	rows = cutil::fast_count_lines(fptr);
	rewind(fptr);

	this->x = gsl_vector_calloc(rows);

	for(i = 0; i < rows; i++)
	{
		fscanf(fptr, "%lf", &a);
		gsl_vector_set(this->x, i, a);
	}

	fclose(fptr);

	fptr = fopen(zf.str().c_str(), "r");
	rows = cutil::fast_count_lines(fptr);
	rewind(fptr);

	this->z = gsl_vector_calloc(rows);

	for(i = 0; i < rows; i++)
	{
		fscanf(fptr, "%lf", &a);
		gsl_vector_set(this->z, i, a);
	}

	fclose(fptr);

	this->xDim = this->x->size;
	this->zDim = this->z->size;

	Output::print("MonteCarlo::init() finished");

}
