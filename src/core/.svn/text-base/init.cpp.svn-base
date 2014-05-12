#include "init.h"
#include "cpba.h"
#include "util.h"
#include "out.h"
#include <limits>
#include <cmath>
#include <boost/format.hpp>
#include <gsl/gsl_blas.h>

using namespace std;

void PBAInitialization::createSlab(double xslab1, double xslab2, double thickness, double depth, int CTNUM_SLAB, string material_SLAB, int CTNUM_BK, string material_BK)
{
	Output::print("PBAInitialization::createSlab()");
	double xs[5];

	if(xslab1 == numeric_limits<double>::min() && xslab2 != numeric_limits<double>::max())
	{
		//Left single-sided laterally infinite slab
		xs[0] = CPBA::data->AXIS->Xmin;
		xs[1] = CPBA::data->AXIS->Xmin;
		xs[2] = xslab2;
		xs[3] = xslab2;
		xs[4] = CPBA::data->AXIS->Xmin;
	}
	else if(xslab1 != numeric_limits<double>::min() && xslab2 == numeric_limits<double>::max())
	{
		//Right single-sided laterally infinite slab
		xs[0] = xslab1;
		xs[1] = xslab1;
		xs[2] = CPBA::data->AXIS->Xmax;
		xs[3] = CPBA::data->AXIS->Xmax;
		xs[4] = xslab1;
	}
	else if(xslab1 == numeric_limits<double>::min() && xslab2 == numeric_limits<double>::max())
	{
		//Laterally infinite slab
		//double xs[5] = {CPBA::data->AXIS->Xmin, CPBA::data->AXIS->Xmin, CPBA::data->AXIS->Xmax, CPBA::data->AXIS->Xmax,CPBA::data->AXIS->Xmin};
		xs[0] = CPBA::data->AXIS->Xmin;
		xs[1] = CPBA::data->AXIS->Xmin;
		xs[2] = CPBA::data->AXIS->Xmax;
		xs[3] = CPBA::data->AXIS->Xmax;
		xs[4] = CPBA::data->AXIS->Xmin;
	}
	else
	{
		//Right parallelpiped
		//double xs[5] = {sx, sx, ex, ex, sx};
		xs[0] = xslab1;
		xs[1] = xslab1;
		xs[2] = xslab2;
		xs[3] = xslab2;
		xs[4] = xslab1;
	}

	double zs[5] = {depth+thickness, depth, depth, depth+thickness, depth+thickness};

	double xw[5] = {CPBA::data->AXIS->Xmin, CPBA::data->AXIS->Xmin, CPBA::data->AXIS->Xmax, CPBA::data->AXIS->Xmax, CPBA::data->AXIS->Xmin};
	double zw[5] = {CPBA::data->AXIS->Zmax, CPBA::data->AXIS->Zmin, CPBA::data->AXIS->Zmin, CPBA::data->AXIS->Zmax, CPBA::data->AXIS->Zmax};

	ROIData *rr;
	
	rr = new ROIData;
	rr->NPts = 5;
	rr->XCon = cutil::array_to_vec(xw, 5);
	rr->ZCon = cutil::array_to_vec(zw, 5);
	rr->CTN = CTNUM_BK;
	rr->MAT = material_BK;

	//Output::print_vector("XCon", rr->XCon);
	//Output::print_vector("ZCon", rr->ZCon);

	CPBA::data->ROI->push_back(rr);

	rr = new ROIData;

	rr->NPts = 5;
	rr->XCon = cutil::array_to_vec(xs, 5);
	rr->ZCon = cutil::array_to_vec(zs, 5);
	
	//Output::print_vector("rr->XCon", rr->XCon);
	//Output::print_vector("rr->ZCon", rr->ZCon);
	
	rr->CTN = CTNUM_SLAB;
	rr->MAT = material_SLAB;

	CPBA::data->ROI->push_back(rr);

	CPBA *cpba = (CPBA *)CPBA::data->cpba;
	cpba->materialCreated = true;

	Output::print(boost::format("SLAB::MATERIAL: %1%") % material_SLAB);
	Output::print(boost::format("BACKGROUND::MATERIAL: %1%") % material_BK);
	Output::print("PBAInitialization::createSlab() finished");
}

void PBAInitialization::createPhantom(int CTNUM, string material)
{
	vector<double> xs;
	vector<double> zs;

	if(CPBA::data->PARAM->Hstep != 0)
	{
		/*
		 % If there is a step gradient designed for the surface contour,
            % design a step phantom
            fprintf('        -- Step surface designed for patient surface contour\n');
            xs = [s.AXIS.Xmin,s.AXIS.Xmin,s.AXIS.X(s.AXIS.NX/2),s.AXIS.X(s.AXIS.NX/2),s.AXIS.Xmax,s.AXIS.Xmax,s.AXIS.Xmin];
            zs = [s.AXIS.Zmax,s.AXIS.Zmin,s.AXIS.Zmin,s.PARAM.Hstep,s.PARAM.Hstep,s.AXIS.Zmax,s.AXIS.Zmax];
		*/

		
		xs.push_back(CPBA::data->AXIS->Xmin);
		xs.push_back(CPBA::data->AXIS->Xmin);
		//we have to shift one index to the left to match with the matlab version.
		xs.push_back(gsl_vector_get(CPBA::data->AXIS->X, ((int)CPBA::data->AXIS->NX/2)-1));
		xs.push_back(gsl_vector_get(CPBA::data->AXIS->X, ((int)CPBA::data->AXIS->NX/2)-1));
		xs.push_back(CPBA::data->AXIS->Xmax);
		xs.push_back(CPBA::data->AXIS->Xmax);
		xs.push_back(CPBA::data->AXIS->Xmin);

		zs.push_back(CPBA::data->AXIS->Zmax);
		zs.push_back(CPBA::data->AXIS->Zmin);
		zs.push_back(CPBA::data->AXIS->Zmin);
		zs.push_back(CPBA::data->PARAM->Hstep);
		zs.push_back(CPBA::data->PARAM->Hstep);
		zs.push_back(CPBA::data->AXIS->Zmax);
		zs.push_back(CPBA::data->AXIS->Zmax);
	}
	else
	{
		xs.push_back(CPBA::data->AXIS->Xmin);
		xs.push_back(CPBA::data->AXIS->Xmin);
		xs.push_back(CPBA::data->AXIS->Xmax);
		xs.push_back(CPBA::data->AXIS->Xmax);
		xs.push_back(CPBA::data->AXIS->Xmin);

		zs.push_back(CPBA::data->AXIS->Zmax);
		zs.push_back(CPBA::data->AXIS->Zmin);
		zs.push_back(CPBA::data->AXIS->Zmin);
		zs.push_back(CPBA::data->AXIS->Zmax);
		zs.push_back(CPBA::data->AXIS->Zmax);
	}

	ROIData *rr = new ROIData;

	rr->NPts = xs.size();
	rr->XCon = cutil::stl_vector_to_gsl_vector(&xs);
	rr->ZCon = cutil::stl_vector_to_gsl_vector(&zs);

	//Output::print_vector("rr->XCon", rr->XCon);
	//Output::print_vector("rr->ZCon", rr->ZCon);

	rr->CTN = CTNUM;
	rr->MAT = material;

	CPBA::data->ROI->push_back(rr);

	CPBA *cpba = (CPBA *)CPBA::data->cpba;
	cpba->materialCreated = true;
}

void PBAInitialization::processROIs()
{
	gsl_matrix *m;
	gsl_matrix *m2;

	gsl_vector *sm;

	for(unsigned int i = 0; i < CPBA::data->ROI->size(); i++)
	{
		ROIData *rr = CPBA::data->ROI->at(i);
		double xs1[4] = {gsl_vector_min(rr->XCon), gsl_vector_min(rr->XCon), gsl_vector_max(rr->ZCon), gsl_vector_max(rr->ZCon)};
		rr->Xcorner = cutil::array_to_vec(xs1, 4);

		double xs2[4] = {gsl_vector_max(rr->ZCon), gsl_vector_min(rr->ZCon), gsl_vector_min(rr->ZCon), gsl_vector_max(rr->ZCon)};
		rr->Zcorner = cutil::array_to_vec(xs2, 4);

		//Output::print_vector("rr->Xcorner", rr->Xcorner);
		//Output::print_vector("rr->Zcorner", rr->Zcorner);

		m = gsl_matrix_alloc(2, 2);
		sm = gsl_vector_calloc(2);

		gsl_vector_set(sm, 0, cos(CPBA::data->PARAM->Theta));
		gsl_vector_set(sm, 1, sin(CPBA::data->PARAM->Theta));
		
		gsl_matrix_set_row(m, 0, sm);

		gsl_vector_free(sm);

		sm = gsl_vector_calloc(2);

		gsl_vector_set(sm, 0, -sin(CPBA::data->PARAM->Theta));
		gsl_vector_set(sm, 1, cos(CPBA::data->PARAM->Theta));

		gsl_matrix_set_row(m, 1, sm);

		gsl_vector_free(sm);

		gsl_matrix *result = gsl_matrix_alloc(2, rr->ZCon->size);

		m2 = gsl_matrix_alloc(2, rr->ZCon->size);

		gsl_matrix_set_row(m2, 0, rr->XCon);
		gsl_matrix_set_row(m2, 1, rr->ZCon);

		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m, m2, 0.0, result);

		rr->XCon = cutil::get_matrix_row(result, 0);
		rr->ZCon = cutil::get_matrix_row(result, 1);

		//Output::print_vector("rr->XCon", rr->XCon);
		//Output::print_vector("rr->ZCon", rr->ZCon);

		gsl_matrix_free(m2);
	}
}

void PBAInitialization::createAxis()
{
	Output::print("PBAInitialization::createAxis()");
	Grid *axis = new Grid();

	PBAData *local = CPBA::data;

	int xo = 10;
	int zo = 4;

	axis->DX = CPBA::data->PARAM->DX;
	axis->DZ = CPBA::data->PARAM->DZ;

	double Xmin = floor((CPBA::data->COLL->Xmin) /(cos(CPBA::data->PARAM->Theta)) - xo - (CPBA::data->BEAM->Rp + CPBA::data->PARAM->Hstep + zo) * sin(CPBA::data->PARAM->Theta));
	double Xmax = ceil(CPBA::data->COLL->Xmax + xo);
	double Zmin = 0.0;
	double Zmax = ceil((CPBA::data->BEAM->Rp + CPBA::data->PARAM->Hstep + zo) * cos(CPBA::data->PARAM->Theta));

	PBAInitialization::cartesianGrid(axis, Xmin, Xmax, Zmin, Zmax, axis->DX, axis->DZ);

	vector<int> ii = cutil::vec_operation_greater_than_index(axis->Z, axis->DZ/2);

	vector<double> iz;

	double dz_factor = axis->DZ / 2.0;
	iz.push_back(axis->DZ/2);

	for(unsigned int i = 0; i < ii.size(); i++)
	{
		iz.push_back(gsl_vector_get(axis->Z, ii.at(i)) + dz_factor);
	}

	gsl_vector_free(axis->Z);
	axis->Z = cutil::stl_vector_to_gsl_vector(&iz);
	
	axis->Zmin = gsl_vector_min(axis->Z);
	axis->Zmax = gsl_vector_max(axis->Z);
	axis->Xmin = gsl_vector_min(axis->X);
	axis->Xmax = gsl_vector_max(axis->X);
	axis->NZ = axis->Z->size;
	axis->NX = axis->X->size;

	//Output::print(boost::format("NX: %1%, NZ: %2%") % axis->NX % axis->NZ);


	double xx[5] = {axis->Xmin, axis->Xmin, axis->Xmax, axis->Xmax, axis->Xmin};
	double zz[5] = {axis->Zmax, axis->Zmin, axis->Zmin, axis->Zmax, axis->Zmax};

	axis->Xcorner = cutil::array_to_vec(xx, 5);
	axis->Zcorner = cutil::array_to_vec(zz, 5);

	//Output::print_vector("Xcorner", axis->Xcorner);
	//Output::print_vector("Zcorner", axis->Zcorner);

	CPBA::data->AXIS = axis;

	Output::print("PBAInitialization::createAxis() finished");
}

void PBAInitialization::cartesianGrid(Grid *axis, double Xmin, double Xmax, double Zmin, double Zmax, double DX, double DZ)
{
	double Xst;
	double Xend;
	double Zst;
	double Zend;

	gsl_vector *Xleft;
	gsl_vector *Xright;
	gsl_vector *Zleft;
	gsl_vector *Zright;
	gsl_vector *Xr;
	gsl_vector *Zl;
	gsl_vector *Xl;
	gsl_vector *Zr;

	if(Xmin < 0)
	{
		Xst = ceil(Xmin/DX) * DX + DX / 2.0;
	}
	else
	{
		Xst = floor(Xmax/DX) * DX + DX / 2.0;
	}

	if(Xmax > 0)
	{
		Xend = floor(Xmax/DX) * DX - DX / 2.0;
	}
	else
	{
		Xend = ceil(Xmax/DX) * DX - DX / 2.0;
	}


	if(Zmin < 0)
	{
		Zst = ceil(Zmin/DZ) * DZ + DZ / 2.0;
	}
	else
	{
		Zst = floor(Zmin/DZ) * DZ + DZ / 2.0;
	}

	if(Zmax > 0)
	{
		Zend = floor(Zmax/DZ) * DZ - DZ / 2.0;
	}
	else
	{
		Zend = ceil(Zmax/DZ) * DZ - DZ / 2.0;
	}

	if(abs(Xmin) >= abs(Xmax))
	{
		//Output::print("1");
		Xleft = cutil::MakeVectorSequence<double>(Xst, DX, 0);
		Xr = cutil::vec_operation_copy(Xleft);
		gsl_vector_reverse(Xr);
		cutil::vec_operation_inline_negate(Xr);

		Xright = cutil::vec_operation_less_than_equal(Xr, abs(Xend), true);
	}
	else
	{
		//Output::print("2");
		Xright = cutil::MakeVectorSequence<double>(0, DX, Xend);
		Xl = cutil::vec_operation_copy(Xright);
		gsl_vector_reverse(Xl);
		cutil::vec_operation_inline_negate(Xl);
		Xleft = cutil::vec_operation_less_than_equal(Xl, abs(Xst), true);
	}

	//Output::print_vector("Xleft", Xleft);
	//Output::print_vector("Xright", Xright);

	if(abs(Zmin) >= abs(Zmax))
	{
		Zleft = cutil::MakeVectorSequence<double>(Zst, DZ, 0);
		Zr = cutil::vec_operation_copy(Zleft);
		gsl_vector_reverse(Zr);
		cutil::vec_operation_inline_negate(Zr);

		Zright = cutil::vec_operation_less_than_equal(Zr, abs(Zend), true);
	}
	else
	{
		Zright = cutil::MakeVectorSequence<double>(0, DZ, Zend);
		Zl = cutil::vec_operation_copy(Zright);
		gsl_vector_reverse(Zl);
		cutil::vec_operation_inline_negate(Zl);
		Zleft = cutil::vec_operation_less_than_equal(Zl, abs(Zst), true);
	}

	//Output::print_vector("Zleft", Zleft);
	//Output::print_vector("Zright", Zright);

	vector<double> XT;
	vector<double> ZT;

	unsigned int i;

	//Output::print(boost::format("left_spoke: %1%, right_spoke: %2%") % ((-DX)/2.0) % (DX/2.0));

	gsl_vector *xl_temp = cutil::vec_operation_less_than_not_zero(Xleft, (-DX)/2.0);
	gsl_vector *xr_temp = cutil::vec_operation_greater_than_not_zero(Xright, DX/2.0);

	double val;
	if(xl_temp != (gsl_vector *)NULL)
	{
		for(i = 0; i < xl_temp->size; i++)
		{
			val = gsl_vector_get(xl_temp, i);
			if(val < -(DX/2))
			{
				//Output::print(boost::format("pushing back %1%, its less than %2%") % val % (-DX/2));
				XT.push_back(val);
			}
		}
	}

	XT.push_back(-DX/2);
	XT.push_back(DX/2);

	if(xr_temp != (gsl_vector *)NULL)
	{
		for(i = 0; i < xr_temp->size; i++)
		{
			val = gsl_vector_get(xr_temp, i);
			if(val > (DX/2))
			{
				XT.push_back(val);
			}
		}
	}

	gsl_vector *X = cutil::stl_vector_to_gsl_vector(&XT);

	//Output::print_vector("X", X);

	// Z 
	gsl_vector *zl_temp = cutil::vec_operation_less_than_not_zero(Zleft, (-DZ)/2.0);
	gsl_vector *zr_temp = cutil::vec_operation_greater_than_not_zero(Zright, DZ/2.0);
	
	if(zl_temp != (gsl_vector *)NULL)
	{
		for(i = 0; i < zl_temp->size; i++)
		{
			val = gsl_vector_get(zl_temp, i);

			if(val < -(DZ/2))
			{
				ZT.push_back(val);
			}
		}
	}

	ZT.push_back(-DX/2);
	ZT.push_back(DX/2);

	if(zr_temp != (gsl_vector *)NULL)
	{
		for(i = 0; i < zr_temp->size; i++)
		{
			val = gsl_vector_get(zr_temp, i);

			if(val > (DZ/2))
			{
				ZT.push_back(val);
			}
		}
	}

	gsl_vector *Z = cutil::stl_vector_to_gsl_vector(&ZT);

	axis->X = X;
	axis->Z = Z;
	
}

void PBAInitialization::createGrids()
{
	PBAData *dd = CPBA::data;

	int i;
	Grid *CTMAT = new Grid;
	DoseMatrix *DOSEMAT = new DoseMatrix;

	CTMAT->Zmin = 200;
	CTMAT->Zmax = -200;
	CTMAT->Xmin = 100;
	CTMAT->Xmax = -100;

	ROIData *rr = CPBA::data->ROI->at(0);

	for(i = 0; i < rr->NPts; i++)
	{
		double xc = gsl_vector_get(rr->XCon, i);
		double zc = gsl_vector_get(rr->ZCon, i);

		if(xc < CTMAT->Xmin)
		{
			CTMAT->Xmin = xc;
		}

		if(xc > CTMAT->Xmax)
		{
			CTMAT->Xmax = xc;
		}

		if(zc < CTMAT->Zmin)
		{
			CTMAT->Zmin = zc;
		}

		if(zc > CTMAT->Zmax)
		{
			CTMAT->Zmax = zc;
		}
	}
	
	//Output::print(boost::format("CTMAT-Zmin: %1%") % CTMAT->Zmin);

	double XWID = CPBA::data->COLL->WX + 2.0;

	DOSEMAT->DZ = CPBA::data->AXIS->DZ;
	DOSEMAT->DX = CPBA::data->AXIS->DX;

	DOSEMAT->NX = (int)(2 * floor(XWID/DOSEMAT->DX));

	DOSEMAT->X = gsl_vector_calloc(DOSEMAT->NX);

	for(i = 0; i < DOSEMAT->NX; i++)
	{
		double part = ((i+1) - ( (DOSEMAT->NX - 1.0) / 2.0 ) - 1.0);
		gsl_vector_set(DOSEMAT->X, i, part * DOSEMAT->DX);
	}

	DOSEMAT->NZ = (int)floor((CTMAT->Zmax - CTMAT->Zmin) / DOSEMAT->DZ);
	
	double a = 0.0;

	while(a >= CTMAT->Zmin)
	{
		a = a - DOSEMAT->DZ / 2.0;
	}
	
	//Output::print(boost::format("a is: %1%") % a);

	DOSEMAT->Z = gsl_vector_calloc(DOSEMAT->NZ);

	for(i = 0; i < DOSEMAT->NZ; i++)
	{
		gsl_vector_set(DOSEMAT->Z, i, i * DOSEMAT->DZ + a);
	}

	if(CPBA::data->PARAM->ThetaDeg == 0)
	{
		gsl_vector_add_constant(DOSEMAT->Z, 0 - gsl_vector_get(DOSEMAT->Z, 0));
		gsl_vector_add_constant(DOSEMAT->Z, DOSEMAT->DZ / 2);
	}
	
	//Output::print_vector("DOSEMAT-Z", DOSEMAT->Z);

	DOSEMAT->Xmin = CTMAT->Xmin = gsl_vector_min(DOSEMAT->X);
	DOSEMAT->Xmax = CTMAT->Xmax = gsl_vector_max(DOSEMAT->X);
	DOSEMAT->Zmin = CTMAT->Zmin = gsl_vector_min(DOSEMAT->Z);
	DOSEMAT->Zmax = CTMAT->Zmax = gsl_vector_max(DOSEMAT->Z);

	CTMAT->DX = DOSEMAT->DX;
	CTMAT->DZ = DOSEMAT->DZ;
	CTMAT->NX = DOSEMAT->NX;
	CTMAT->NZ = DOSEMAT->NZ;

	CTMAT->X = cutil::vec_operation_copy(DOSEMAT->X);
	CTMAT->Z = cutil::vec_operation_copy(DOSEMAT->Z);

	CPBA::data->CTMAT = CTMAT;
	CPBA::data->DOSEMAT = DOSEMAT;
}
