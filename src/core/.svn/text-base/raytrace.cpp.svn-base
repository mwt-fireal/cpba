#include <math.h>

#include "raytrace.h"
#include "scatter.h"
#include "interp.h"
#include "out.h"

#include <algorithm>
#include <sstream>
#include <cctype>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <time.h>
#include <map>
#include <string>

void RayTrace::trace()
{
	int i;
	CPBA::data->RAY = new Ray(CPBA::data->CTMAT->NZ, CPBA::data->CTMAT->NX);

	//shortcut!
	Ray *RAY = CPBA::data->RAY;

	double Ro = RayTrace::range(CPBA::data->BEAM->Eo, "WATER");

	//Output::print(boost::format("Ro is: %1%") % Ro);

	double Ez, E2, M0, M1, M2, SM0, SM1, SM2;
	
	double current_ctn;
	
	for(i = 0; i < CPBA::data->CTMAT->NX; i++)
	{
		//fprintf('        --> %2.1f%% complete ...\n',i/s.CTMAT.NX*100);
		Output::print(boost::format("--> Raytrace %1% out of %2%") % (i+1) % CPBA::data->CTMAT->NX);

		gsl_matrix_set(RAY->Energy, 0, i, CPBA::data->BEAM->Eo);
		gsl_matrix_set(RAY->Energy2, 0, i, CPBA::data->BEAM->Eo);

		Ez = E2 = CPBA::data->BEAM->Eo;
		M0 = M1 = M2 = 0.0;
		SM0 = SM1 = SM2 = 0.0;
		
		//Output::print(boost::format("RIO SIZE: %1%") % CPBA::data->ROI->size());

		for(int j = 0; j < CPBA::data->CTMAT->NZ; j++)
		{
			for(unsigned int k = 0; k < CPBA::data->ROI->size(); k++)
			{
				ROIData *r = CPBA::data->ROI->at(k);
				//Output::print(boost::format("r->XCon: %1%") % r->XCon);
				//Output::print(boost::format("r->ZCon: %1%") % r->ZCon);
				//Output::print(boost::format("CTMAT->X[%1%]: %2%") % i % gsl_vector_get(CPBA::data->CTMAT->X, i));
				//Output::print(boost::format("CTMAT->Z[%1%]: %2%") % i % gsl_vector_get(CPBA::data->CTMAT->Z, j));
				if(RayTrace::inROI(r->XCon, r->ZCon, gsl_vector_get(CPBA::data->CTMAT->X, i), gsl_vector_get(CPBA::data->CTMAT->Z, j)))
				{					
					gsl_matrix_set(RAY->CTN, j, i, r->CTN);
					current_ctn = gsl_matrix_get(RAY->CTN, j, i);
					//Output::print(boost::format("current ctn: %1%") % current_ctn);

					RAY->MAT[j][i] = r->MAT;
					
					double ratio = Scatter::getSPratio(current_ctn, Ez);
					
					//Output::print(boost::format("spRatio for ctn %1%, Ez %2%, is %3%") % current_ctn % Ez % ratio);

					double zff = CPBA::data->CTMAT->DZ/2 * ratio;
					gsl_matrix_set(RAY->Zeff, j, i, zff);
				}
				else
				{
					//Output::print(boost::format("NOT inROI for [%1%][%2%]") % j % i);
				}
			}
			
			//Output::print_matrix("CTN", RAY->CTN);
			//Output::print_matrix("Zeff", RAY->Zeff);

			// Find material density
			double rho = ((Material *)CPBA::data->MaterialDB->mdb->find(RAY->MAT[j][i])->second)->rho;
			gsl_matrix_set(RAY->rho, j, i, rho);

			if(j > 0)
			{
				double CTavg = .5 * (current_ctn + gsl_matrix_get(RAY->CTN, j-1, i));
				double ratio = Scatter::getSPratio(CTavg, Ez);
				
				//Output::print(boost::format("******spRatio for ctn %1%, Ez %2%, is %3%") % CTavg % Ez % ratio);
				
				gsl_matrix_set(RAY->Zeff, j, i, gsl_matrix_get(RAY->Zeff, j-1, i) + CPBA::data->CTMAT->DZ * ratio);
			}

			gsl_matrix *pstar;

			if(RAY->MAT[j][i] == "WATER")
			{
				if(Scatter::PSTAR_W == (gsl_matrix *)NULL)
				{
					Scatter::PSTAR_W = cutil::read_matrix_from_file(boost::format("%1%/pstarWATER.txt") % CPBA::INPUT_PATH);
				}

				pstar = Scatter::PSTAR_W;
			}
			else if(RAY->MAT[j][i] == "COMPACTBONE")
			{
				if(Scatter::PSTAR_B == (gsl_matrix *)NULL)
				{
					Scatter::PSTAR_B = cutil::read_matrix_from_file(boost::format("%1%/pstarCOMPACTBONE.txt") % CPBA::INPUT_PATH);
				}

				pstar = Scatter::PSTAR_B;
			}
			else if(RAY->MAT[j][i] == "AIR")
			{
				if(Scatter::PSTAR_A == (gsl_matrix *)NULL)
				{
					Scatter::PSTAR_A = cutil::read_matrix_from_file(boost::format("%1%/pstarAIR.txt") % CPBA::INPUT_PATH);
				}

				pstar = Scatter::PSTAR_A;
			}
			else if(RAY->MAT[j][i] == "VACUUM")
			{
				if(Scatter::PSTAR_V == (gsl_matrix *)NULL)
				{
					Scatter::PSTAR_V = cutil::read_matrix_from_file(boost::format("%1%/pstarVACUUM.txt") % CPBA::INPUT_PATH);
				}

				pstar = Scatter::PSTAR_V;
			}
			else
			{
				Output::print(boost::format("RayTrace::trace() -- ERROR: unknown material type %1%") % RAY->MAT[j][i]);
				exit(1);
			}
			
			if(Scatter::PSTAR_W == (gsl_matrix *)NULL)
			{
				Scatter::PSTAR_W = cutil::read_matrix_from_file(boost::format("%1%/pstarWATER.txt") % CPBA::INPUT_PATH);
			}

			gsl_matrix *pstar2 = Scatter::PSTAR_W;

			double Sz = LinearInterpolation::evaluate_constant_linear(pstar->size1, cutil::get_matrix_col(pstar, 0), cutil::get_matrix_col(pstar, 1), Ez) * rho;
			double S2 = LinearInterpolation::evaluate_constant_linear(pstar2->size1, cutil::get_matrix_col(pstar2, 0), cutil::get_matrix_col(pstar2, 1), E2);

			gsl_matrix_set(RAY->Sz, j, i, Sz);
			gsl_matrix_set(RAY->Sz2, j, i, S2);

			if(RAY->MAT[j][i].compare("VACUUM") != 0)
			{
				Ez = Ez - Sz * CPBA::data->CTMAT->DZ;
				E2 = E2 - S2 * CPBA::data->CTMAT->DZ;
			}

			Ez = Ez < 0 ? 0 : Ez;
			E2 = E2 < 0 ? 0 : E2;

			gsl_matrix_set(RAY->Energy, j, i, Ez);
			gsl_matrix_set(RAY->Energy2, j, i, E2);

			if(RAY->MAT[j][i].compare("VACUUM") != 0)
			{
				gsl_matrix_set(RAY->SPratio, j, i, Scatter::getSPratio(current_ctn, Ez) / rho);
			}
			else
			{
				gsl_matrix_set(RAY->SPratio, j, i, 0);
			}

			double SIGTHZ = Scatter::calcTdM(RAY->MAT[j][i], Ez, CPBA::data->BEAM->Eo);
			double SIGTH2;

			if(RAY->MAT[j][i].compare("VACUUM") != 0)
			{
				SIGTH2 = Scatter::calcTdM("WATER", Ez, CPBA::data->BEAM->Eo);
			}
			else
			{
				SIGTH2 = 0.0;
			}

			gsl_matrix_set(RAY->SIGTH, j, i, SIGTHZ);
			gsl_matrix_set(RAY->SIGTH2, j, i, SIGTH2);

			M2 = M2 - 2 * M1 + M0 + SIGTHZ * CPBA::data->CTMAT->DZ;
			M1 = M1 - M0 - SIGTHZ * CPBA::data->CTMAT->DZ;
			M0 = M0 + SIGTHZ * CPBA::data->CTMAT->DZ;

			SM2 = SM2 - 2 * SM1 + SM0 + SIGTH2 * CPBA::data->CTMAT->DZ;
			SM1 = SM1 - SM0 - SIGTH2 * CPBA::data->CTMAT->DZ;
			SM0 = SM0 + SIGTH2 * CPBA::data->CTMAT->DZ;

			double A0, A1, A2, B;

			A0 = (0.5 * M0);
			A1 = (-0.5 * (M1 + M0/2) * CPBA::data->CTMAT->DZ);
			A2 = (0.5 * (M2 + M1 + M0/3) * pow(CPBA::data->CTMAT->DZ, 2));
			B = A0 * A2 - pow(A1, 2);

			gsl_matrix_set(RAY->A0, j, i, A0);
			gsl_matrix_set(RAY->A1, j, i, A1);
			gsl_matrix_set(RAY->A2, j, i, A2);
			gsl_matrix_set(RAY->B, j, i, B);

			if(gsl_matrix_get(RAY->Zeff, j, i) <= Ro)
			{
				gsl_matrix_set(RAY->SIGX1, j, i, sqrt((M0/3 + M1 + M2) * pow(CPBA::data->CTMAT->DZ, 2)));
				gsl_matrix_set(RAY->SIGX2, j, i, sqrt((SM0/3 + SM1 + SM2) * pow(CPBA::data->CTMAT->DZ, 2)));
			}
			else
			{
				gsl_matrix_set(RAY->SIGX1, j, i, gsl_matrix_get(RAY->SIGX1, j-1, i));
				gsl_matrix_set(RAY->SIGX2, j, i, gsl_matrix_get(RAY->SIGX2, j-1, i));
			}
		}
	}

	/*
	string file;
	
	boost::format ff("DebugData");
	
	boost::filesystem::path out_dir = ff.str();
	
	if(boost::filesystem::exists(out_dir))
	{
		boost::filesystem::remove_all(out_dir);
	}
	
	boost::filesystem::create_directories(out_dir);
	
	file = boost::str(boost::format("%1%/RayTrace_Zeff.txt") % ff.str());
	Output::print(file);
	cutil::write_data_to_file(file, RAY->Zeff);
	*/
	
	//Output::print_matrix("RAY->rho", RAY->rho);
	//Output::print_matrix("RAY->CTN", RAY->CTN);
	//Output::print_matrix("RAY->Zeff", RAY->Zeff);
	//Output::print_matrix("RAY->Sz", RAY->Sz);
	//Output::print_matrix("RAY->Sz2", RAY->Sz2);
	//Output::print_matrix("RAY->Energy", RAY->Energy);
	//Output::print_matrix("RAY->Energy2", RAY->Energy2);
	//Output::print_matrix("RAY->SPratio", RAY->SPratio);
	//Output::print_matrix("RAY->SIGTH", RAY->SIGTH);
	//Output::print_matrix("RAY->SIGTH2", RAY->SIGTH2);
	//Output::print_matrix("RAY->A0", RAY->A0);
	//Output::print_matrix("RAY->A1", RAY->A1);
	//Output::print_matrix("RAY->A2", RAY->A2);
	//Output::print_matrix("RAY->B", RAY->B);
	//Output::print_matrix("RAY->SIGX1", RAY->SIGX1);
	//Output::print_matrix("RAY->SIGX2", RAY->SIGX2);
}

bool RayTrace::inROI(gsl_vector *roiX, gsl_vector *roiZ, double xo, double zo)
{
	int i;
	bool inside = false;
	int zLineSign = 0;

	zROI *zm = new zROI(20);

	for(unsigned int i = 0; i < roiX->size; i++)
	{
		int i2 = (i % roiX->size) + 1;

		if(i2 == roiX->size)
		{
			i2 = 0;
		}

		zm->x1 = gsl_vector_get(roiX, i);
		zm->x2 = gsl_vector_get(roiX, i2);

		zm->z1 = gsl_vector_get(roiZ, i);
		zm->z2 = gsl_vector_get(roiZ, i2);

		if(zm->x2 < zm->x1)
		{
			if(xo <= zm->x1 && xo >= zm->x2)
			{
				zLineSign = -1;
				//Output::print("zLineSign set to -1"); 
				RayTrace::zoLineTest(zm, xo, zo, zLineSign);
			}
		}
		else
		{
			if(xo <= zm->x2 && xo >= zm->x1)
			{
				//Output::print("zLineSign set to 1");
				zLineSign = 1;
				RayTrace::zoLineTest(zm, xo, zo, zLineSign);
			}

		}
	}

	//Output::print(boost::format("zLineSign: %1%") % zLineSign);

	if(zLineSign != 0)
	{
		int iset = 0;

		//Output::print_vector("zm->ZL",zm->ZL);
		//Output::print_vector("zm->L",zm->L);
		//Output::print(boost::format("zm->M: %1%") % zm->M);

		for(i = 1; i < zm->M; i++)
		{
			double ZL_1 = gsl_vector_get(zm->ZL, i);
			double ZL_1m = gsl_vector_get(zm->ZL, i-1);
			if(ZL_1m < ZL_1)
			{
				gsl_vector_set(zm->ZL, i, ZL_1m);
				gsl_vector_set(zm->ZL, i-1, ZL_1);
				double Ltemp = gsl_vector_get(zm->L, i);
				gsl_vector_set(zm->L, i, gsl_vector_get(zm->L, i-1));
				gsl_vector_set(zm->L, i-1, Ltemp);
				iset = 1;
			}
		}

		//Output::print_vector("zm->ZL",zm->ZL);
		//Output::print_vector("zm->L",zm->L);
		//Output::print(boost::format("zm->M: %1%") % zm->M);

		for(i = 0; i < zm->M; i+=2)
		{
			if(gsl_vector_get(zm->L, i) && gsl_vector_get(zm->L, i+1))
			{
				//Output::print("inside is now set to true");
				//Output::print(boost::format("two variables are: %1%, %2%") % gsl_vector_get(zm->L, i) % gsl_vector_get(zm->L, i+1));
				return true;
			}
			else
			{
				inside = false;
			}
		}
	}
	else
	{
		inside = false;
	}

	return inside;
}

void RayTrace::zoLineTest(zROI *zm, double xo, double zo, int zLineSign)
{
	double x21 = zm->x2 - zm->x1;

	if(x21 == 0)
	{
		x21 = 0.0001;
	}

	double zLine = zm->z1 + (zm->z2 - zm->z1) * (xo - zm->x1) / x21;

	int zoBelow = zo >= zLine ? 1 : -1;

	zm->M++;

	int index = zm->M - 1;

	gsl_vector_set(zm->ZL, index, zLine);
	gsl_vector_set(zm->L, index, 1);

	if((zLineSign * zoBelow) == -1)
	{
		gsl_vector_set(zm->L, index, 0);
	}
}

double RayTrace::range(double Energy, string material)
{
	string rangeFile;
	gsl_matrix *m;

	transform(material.begin(), material.end(), material.begin(), (int(*)(int))toupper);

	Output::print(material);

	if(material == "WATER")
	{
		if(Scatter::PSTAR_W == (gsl_matrix *)NULL)
		{
			Scatter::PSTAR_W = cutil::read_matrix_from_file(boost::format("%1%/pstarWATER.txt") % CPBA::INPUT_PATH);
		}

		m = Scatter::PSTAR_W;
	}

	if(material == "VACUUM")
	{
		if(Scatter::PSTAR_V == (gsl_matrix *)NULL)
		{
			Scatter::PSTAR_V = cutil::read_matrix_from_file(boost::format("%1%/pstarVACUUM.txt") % CPBA::INPUT_PATH);
		}

		m = Scatter::PSTAR_V;
	}


	if(material == "COMPACTBONE")
	{
		if(Scatter::PSTAR_B == (gsl_matrix *)NULL)
		{
			Scatter::PSTAR_B = cutil::read_matrix_from_file(boost::format("%1%/pstarCOMPACTBONE.txt") % CPBA::INPUT_PATH);
		}

		m = Scatter::PSTAR_B;
	}

	if(material == "AIR")
	{
		if(Scatter::PSTAR_A == (gsl_matrix *)NULL)
		{
			Scatter::PSTAR_A = cutil::read_matrix_from_file(boost::format("%1%/pstarAIR.txt") % CPBA::INPUT_PATH);
		}

		m = Scatter::PSTAR_A;
	}

	gsl_vector *v1 = cutil::get_matrix_col(m, 0);
	gsl_vector *v4 = cutil::get_matrix_col(m, 2);
	
	/*
    % Return range in g/cm^2. The interpolation method is described in
    % detail in the energy.m function.
	*/
	//Output::print_vector("v1", v1);
	//Output::print_vector("v4", v4);

	return exp(CubicSplineInterpolation::evaluate_constant_spline(v1->size, cutil::vec_operation_log_to_array(v1), cutil::vec_operation_log_to_array(v4), log(Energy)));
}
