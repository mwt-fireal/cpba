#include "scatter.h"
#include "data.h"
#include "cpba.h"
#include "util.h"
#include "interp.h"
#include "mat.h"
#include "out.h"
#include <map>
#include <string>

using namespace std;

gsl_matrix *Scatter::PSTAR_W;
gsl_matrix *Scatter::PSTAR_B;
gsl_matrix *Scatter::PSTAR_A;
gsl_matrix *Scatter::PSTAR_V;

map<double, double*> Scatter::spcache;
map<string, double> Scatter::slcache;

double Scatter::scatLength(string mat)
{
	map<string, double>::iterator it;

	it = Scatter::slcache.find(mat);

	if(it != Scatter::slcache.end())
	{
		//Output::print(boost::format("Scatter::scatLength cache hit for material %1%, returning %2%") % mat % it->second);
		return it->second;
	}

	double alpha = 1/137.03604;
	double Na = 6.022045E23;
	double re = 2.8179402894E-13;

	Material *m = (Material *)CPBA::data->MaterialDB->mdb->find(mat)->second;

	double invScatLen = 0.0;

	map<string, string> mc;

	for(int i = 0; i < m->numElements; i++)
	{
		mc = m->elements.at(i);
		double f = atof(mc["f"].c_str());
		double Z = atof(mc["Z"].c_str());
		double A = atof(mc["A"].c_str());

		invScatLen = invScatLen + f * (alpha * Na * pow(re, 2) * pow(Z, 2) / A * (2 * log(33219 * pow((A*Z), -1/3)) - 1));
	}

	if(invScatLen == 0)
	{
		return 0;
	}

	double Xs = 1 / (m->rho * invScatLen);

	//Output::print(boost::format("Scatter::scatLength cache MISS for mat %1%, Xs is %2%") % mat % Xs);
	Scatter::slcache.insert(pair<string, double>(mat, Xs));
	//Scatter::slcache[mat] = Xs;

	return Xs;
}

double Scatter::calcTdM(string mat, double Ez, double Eo)
{
	double mpc2 = 938.2796;
	double Es = 15;
	double Ecutoff = 0.5435;

	double SPWR = 0.0;

	if(mat.compare("VACUUM") != 0)
	{
		double Xs = Scatter::scatLength(mat);

		double pvo = (pow((Eo+mpc2), 2) - pow(mpc2, 2)) / (Eo+mpc2);
		double pv = (pow((Ez+mpc2), 2) - pow(mpc2, 2)) / (Ez+mpc2);

		double fdM = 0.0;

		if(Ez > Ecutoff && Ez != Eo)
		{
			fdM = .5244 + .1975 * log10(1.0 - pow((pv/pvo), 2)) + .232 * log10(pv) - .0098 * log10(pv) * log10(1 - pow((pv/pvo), 2));
		}
		else if(Ez <= Ecutoff)
		{
			Ez = Ecutoff;
			pv = (pow((Ez+mpc2), 2) - pow(mpc2, 2)) / (Ez+mpc2);
			fdM = .5244 + .1975 * log10(1.0 - pow((pv/pvo), 2)) + .232 * log10(pv) - .0098 * log10(pv) * log10(1 - pow((pv/pvo), 2));
		}
		else if(Ez == Eo)
		{
			Ez = 0.95 * Eo;
			pv = (pow((Ez+mpc2), 2) - pow(mpc2, 2)) / (Ez+mpc2);
			fdM = .5244 + .1975 * log10(1.0 - pow((pv/pvo), 2)) + .232 * log10(pv) - .0098 * log10(pv) * log10(1 - pow((pv/pvo), 2));
		}

		SPWR = fdM * pow((Es/pv), 2) * 1 / Xs;
	}

	return SPWR;
}

/*

double Scatter::getScatPwr(double CTavg, double Ez)
{
	double CTn[4] = {0.0, 100.0, 500.0, 900.0};

	gsl_matrix *W = cutil::read_matrix_from_file("Input Data/scatWATER.txt");
	gsl_matrix *B = cutil::read_matrix_from_file("Input Data/scatCOMPACTBONE.txt");
	gsl_matrix *A = cutil::read_matrix_from_file("Input Data/scatAIR.txt");

	double bone_rho = ((Material *)CPBA::data->MaterialDB->mdb->find("COMPACTBONE")->second)->rho;
	double air_rho = ((Material *)CPBA::data->MaterialDB->mdb->find("AIR")->second)->rho;

	double Wval, Bval, Aval = 0.0;

	if(Ez > gsl_matrix_get(W, 0, 0))
	{
		Wval = LinearInterpolation::evaluate_constant_linear(W->size1, cutil::get_matrix_col(W, 0),cutil::get_matrix_col(W, 1), Ez);
	}
	else
	{
		Wval = gsl_matrix_get(W, 0, 1);
	}

	if(Ez > gsl_matrix_get(B, 0, 0))
	{
		Bval = LinearInterpolation::evaluate_constant_linear(B->size1, cutil::get_matrix_col(B, 0), cutil::get_matrix_col(B, 1), Ez) * bone_rho;
	}
	else
	{
		Bval = gsl_matrix_get(B, 0, 1) * bone_rho;
	}

	if(Ez > gsl_matrix_get(A, 0, 0))
	{
		Aval = LinearInterpolation::evaluate_constant_linear(A->size1, cutil::get_matrix_col(A, 0),cutil::get_matrix_col(A, 1), Ez) * air_rho;
	}
	else
	{
		Aval = gsl_matrix_get(A, 0, 1) * air_rho;
	}

	double LASP[4] = {0.0, Aval, Wval, Bval};

	double T = CubicSplineInterpolation::evaluate_constant_spline(4, CTn, LASP, CTavg);

	gsl_matrix_free(W);
	gsl_matrix_free(B);
	gsl_matrix_free(A);

	return T;
}

*/

double Scatter::getSPratio(double CTavg, double Ez)
{
	map<double, double*>::iterator it;

	if(!cutil::is_finite(Ez))
	{
		Ez = 0;
	}

	it = Scatter::spcache.find(Ez);

	double *LCSP;
	double CTn[4] = {0.0, 100.0, 500.0, 900.0};

	if(it == Scatter::spcache.end())
	{
		double bone_rho = ((Material *)CPBA::data->MaterialDB->mdb->find("COMPACTBONE")->second)->rho;
		double air_rho = ((Material *)CPBA::data->MaterialDB->mdb->find("AIR")->second)->rho;

		if(Scatter::PSTAR_W == (gsl_matrix *)NULL)
		{
			Scatter::PSTAR_W = cutil::read_matrix_from_file(boost::format("%1%/pstarWATER.txt") % CPBA::INPUT_PATH);
		}

		if(Scatter::PSTAR_B == (gsl_matrix *)NULL)
		{
			Scatter::PSTAR_B = cutil::read_matrix_from_file(boost::format("%1%/pstarCOMPACTBONE.txt") % CPBA::INPUT_PATH);
		}

		if(Scatter::PSTAR_A == (gsl_matrix *)NULL)
		{
			Scatter::PSTAR_A = cutil::read_matrix_from_file(boost::format("%1%/pstarAIR.txt") % CPBA::INPUT_PATH);
		}

		//Output::print(boost::format("CTavg: %1%, Ez: %2%") % CTavg % Ez);
		//Output::print_matrix("W", W);
		//Output::print_matrix("B", B);
		//Output::print_matrix("A", A);

		double Wval, Bval, Aval = 0.0;

		gsl_vector *c1, *c2;

		if(Ez > gsl_matrix_get(Scatter::PSTAR_W, 0, 0))
		{
			c1 = cutil::get_matrix_col(Scatter::PSTAR_W, 0);
			c2 = cutil::get_matrix_col(Scatter::PSTAR_W, 1);
			Wval = LinearInterpolation::evaluate_constant_linear(Scatter::PSTAR_W->size1, c1, c2, Ez);
			gsl_vector_free(c1);
			gsl_vector_free(c2);

		}
		else
		{
			Wval = gsl_matrix_get(Scatter::PSTAR_W, 0, 1);
		}

		if(Ez > gsl_matrix_get(Scatter::PSTAR_B, 0, 0))
		{
			c1 = cutil::get_matrix_col(Scatter::PSTAR_B, 0);
			c2 = cutil::get_matrix_col(Scatter::PSTAR_B, 1);
			Bval = LinearInterpolation::evaluate_constant_linear(Scatter::PSTAR_B->size1, c1, c2, Ez) * bone_rho;
			gsl_vector_free(c1);
			gsl_vector_free(c2);

		}
		else
		{
			Bval = gsl_matrix_get(Scatter::PSTAR_B, 0, 1) * bone_rho;
		}

		if(Ez > gsl_matrix_get(Scatter::PSTAR_A, 0, 0))
		{
			c1 = cutil::get_matrix_col(Scatter::PSTAR_A, 0);
			c2 = cutil::get_matrix_col(Scatter::PSTAR_A, 1);
			Aval = LinearInterpolation::evaluate_constant_linear(Scatter::PSTAR_A->size1, c1, c2, Ez) * air_rho;
			gsl_vector_free(c1);
			gsl_vector_free(c2);
		}
		else
		{
			Aval = gsl_matrix_get(Scatter::PSTAR_A, 0, 1) * air_rho;
		}

		//Output::print(boost::format("Wval: %1%, Bval: %2%, Aval: %3%") % Wval % Bval % Aval);

		
		LCSP = (double *)malloc(sizeof(double) * 4);

		LCSP[0] = 0.0;
		LCSP[1] = Aval/Wval;
		LCSP[2] = 1.0;
		LCSP[3] = Bval/Wval;

		//Output::print(boost::format("Scatter::getSPratio cache MISS for Ez %1%, LCSP is [%2%, %3%, %4%, %5%]") % Ez % LCSP[0] % LCSP[1] % LCSP[2] % LCSP[3]);
		Scatter::spcache.insert(pair<double, double*>(Ez, LCSP));
	}
	else
	{
		LCSP = it->second;
		//Output::print(boost::format("Scatter::getSPratio cache HIT for Ez %1%, LCSP is [%2%, %3%, %4%, %5%]") % Ez % LCSP[0] % LCSP[1] % LCSP[2] % LCSP[3]);
	}

	
	double rho = CubicSplineInterpolation::evaluate_constant_spline(4, CTn, LCSP, CTavg);
	
	/*
	if(rho > 0 && rho < 1)
	{
		Output::print(boost::format("Scatter::getSPratio::Ez %1%, CTavg %6%, LCSP is [%2%, %3%, %4%, %5%] == rho %7%") % Ez % LCSP[0] % LCSP[1] % LCSP[2] % LCSP[3] % CTavg % rho);
	}
	*/
	
	//Output::print(boost::format("returning spratio of: %1%") % rho);
	return rho;
}
