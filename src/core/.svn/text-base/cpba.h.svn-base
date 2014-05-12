// libcpba.h

#pragma once

using namespace std;

// master include file for libcpba

#include "data.h"
#include "init.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <boost/program_options.hpp>


class CPBAParameters
{
public:
	double Eo;
	double WX;
	double FMCS;
	double ThetaDeg;
	double Hstep;
	double DX;
	double DZ;
};

class CPBA
{
public:
	static PBAData *data;
	
	static string INPUT_PATH;
	static boost::program_options::variables_map PROG_ARGS;

	gsl_matrix *PDD;
	gsl_matrix *Zeff;
	gsl_matrix *W1;
	gsl_matrix *W2;
	gsl_vector *BP;
	gsl_vector *BN1;
	gsl_vector *BN2;
	gsl_matrix *SIGhNUC1;
	gsl_matrix *SIGhNUC2;
	gsl_matrix *Da;
	gsl_matrix *Db;
	gsl_matrix *pNorm;
	gsl_matrix *nNorm;
	gsl_matrix *Zp;
	gsl_matrix *Zn;

	bool materialCreated;

	CPBA(CPBAParameters *cp);
	
	void init(CPBAMaterialParameters *mp);

	gsl_vector *extractPDD(gsl_matrix *D);
	void execute();

	void calculateDose();

	void finalizeData();

	void save();

	
};
