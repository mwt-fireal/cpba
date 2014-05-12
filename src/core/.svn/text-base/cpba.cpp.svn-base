// This is the main DLL file.

#include <time.h>

#include "cpba.h"
#include "util.h"
#include "monte.h"
#include "interp.h"
#include "raytrace.h"
#include "init.h"
#include "out.h"
#include "scatter.h"
#include "pinnacle.h"

#include <sstream>
#include <cmath>

#include <gsl/gsl_sf.h>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

PBAData *CPBA::data;
string CPBA::INPUT_PATH;
boost::program_options::variables_map CPBA::PROG_ARGS;

CPBA::CPBA(CPBAParameters *cp)
{
	//double Eo, double WX, double FMCS, double ThetaDeg, double Hstep, double DX, double DZ)
	this->materialCreated = false;

	CPBA::data = new PBAData;
	CPBA::data->cpba = (void *)this;
	CPBA::data->MaterialDB = new PBAMaterials((boost::format("%1%/mat.txt") % CPBA::INPUT_PATH).str());
	CPBA::data->BEAM = new Beam;
	CPBA::data->COLL = new Collimator;
	CPBA::data->PARAM = new PBAParameters;
	CPBA::data->ROI = new vector<ROIData*>;
	CPBA::data->GRID = new Grid;
	CPBA::data->AXIS = new Grid;
	
	CPBA::data->DOSEMAT = new DoseMatrix;

	CPBA::data->BEAM->Eo = cp->Eo;
	CPBA::data->BEAM->Rp = 2.2e-3 * pow(cp->Eo, 1.77);

	CPBA::data->COLL->WX = cp->WX;

	/*
	COLL1.Xmin = -COLL1.WX/2;
    COLL1.Xmax = COLL1.WX/2;
	*/

	CPBA::data->COLL->Xmin = 0 - (CPBA::data->COLL->WX / 2.0);
	CPBA::data->COLL->Xmax = CPBA::data->COLL->WX / 2.0;

	CPBA::data->PARAM->FMCS = cp->FMCS;
	CPBA::data->PARAM->ThetaDeg = cp->ThetaDeg;
	CPBA::data->PARAM->Theta = cp->ThetaDeg * M_PI/180.0;
	CPBA::data->PARAM->Hstep = cp->Hstep;
	CPBA::data->PARAM->DX = cp->DX;
	CPBA::data->PARAM->DZ = cp->DZ;
}

void CPBA::init(CPBAMaterialParameters *mp)
{
	Scatter::PSTAR_W = (gsl_matrix *)NULL;
	Scatter::PSTAR_B = (gsl_matrix *)NULL;
	Scatter::PSTAR_A = (gsl_matrix *)NULL;

	Output::print("CPBA::init()");
	
	PBAInitialization::createAxis();

	switch(mp->type)
	{
		//void PBAInitialization::createSlab(double sx, double ex, double thickness, double depth, int CTNUM_SLAB, string material_SLAB, int CTNUM_BK, string material_BK)
		case CPBAMaterialParameters::CPBA_SLAB:
			Output::print("CPBA::CPBA_SLAB Created");
			PBAInitialization::createSlab(mp->sx, mp->ex, mp->thickness, mp->depth, mp->CTNUM_SLAB, mp->MATERIAL_SLAB, mp->CTNUM_BK, mp->MATERIAL_BK);
			break;
		//void PBAInitialization::createPhantom(int CTNUM, string material)
		case CPBAMaterialParameters::CPBA_PHANTOM:
			Output::print("CPBA::CPBA_FLAT_PHANTOM Created");
			PBAInitialization::createPhantom(mp->CTNUM_PHANTOM, mp->MATERIAL_PHANTOM);
			break;
	}

	;

	PBAInitialization::processROIs();

	PBAInitialization::createGrids();

	Output::print("CPBA::init() finished");
}

void CPBA::execute()
{
	clock_t start;
	double elapsed;
	Output::print("CPBA::execute()");

	start = clock(); 
	RayTrace::trace();

	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Raytracing: %1% seconds") % elapsed);

	start = clock(); 
	this->calculateDose();

	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Dose Caclulations: %1% seconds") % elapsed);

	//this->finalizeData();

	Output::print("CPBA::execute() finished");
}
     
void CPBA::calculateDose()
{
	DoseMatrix *DM = CPBA::data->DOSEMAT;

	int NX = DM->NX;
	int NZ = DM->NZ;
	
	Output::print(boost::format("DoseGrid[%1%][%2%]") % NZ % NX);

	DM->Dose = gsl_matrix_calloc(NZ, NX);
	DM->PDose = gsl_matrix_calloc(NZ, NX);
	DM->NDose = gsl_matrix_calloc(NZ, NX);
	DM->Fluence = gsl_matrix_calloc(NZ, NX);
	DM->PFluence = gsl_matrix_calloc(NZ, NX);
	DM->NFluence = gsl_matrix_calloc(NZ, NX);

	this->PDD = gsl_matrix_calloc(NZ, NX);
	this->Zeff = gsl_matrix_calloc(NZ, NX);
	this->W1 = gsl_matrix_calloc(NZ, NX);
	this->W2 = gsl_matrix_calloc(NZ, NX);
	this->SIGhNUC1 = gsl_matrix_calloc(NZ, NX);
	this->SIGhNUC2 = gsl_matrix_calloc(NZ, NX);
	this->pNorm = gsl_matrix_calloc(NZ, NX);
	this->nNorm = gsl_matrix_calloc(NZ, NX);
	this->Zp = gsl_matrix_calloc(NZ, NX);
	this->Zn = gsl_matrix_calloc(NZ, NX);
	this->BP = gsl_vector_calloc(NZ);
	this->BN1 = gsl_vector_calloc(NZ);
	this->BN2 = gsl_vector_calloc(NZ);

	gsl_matrix_scale(CPBA::data->RAY->SIGX1, CPBA::data->PARAM->FMCS);
	gsl_matrix_scale(CPBA::data->RAY->SIGX2, CPBA::data->PARAM->FMCS);

	MonteCarlo *mc = new MonteCarlo((int)floor(CPBA::data->BEAM->Eo + 0.5), (int)floor(CPBA::data->COLL->WX + 0.5), 0, 0);

	gsl_vector *PDDmc = this->extractPDD(mc->D);

	CPBA::data->PARAM->pbIndicies = new list<double>;

	Grid *CTMAT = CPBA::data->CTMAT;

	int midX;

	// be off by one on the negative side as compared to the matlab version since C indexes from 0:N-1 not 1:N
	if(DM->NX % 2 == 1)
	{
		midX = (DM->NX / 2);
	}
	else
	{
		midX = (DM->NX / 2) - 1;
	}
	
	//Output::print(boost::format("midX: %1%") % midX);

	for(unsigned int i = 0; i < CTMAT->X->size; i++)
	{
		double val = gsl_vector_get(CTMAT->X, i);

		if(val > CPBA::data->COLL->Xmin && val < CPBA::data->COLL->Xmax)
		{
			CPBA::data->PARAM->pbIndicies->push_back(i);
		}
	}

	CPBA::data->PARAM->NPB = CPBA::data->PARAM->pbIndicies->size();
	
	//Output::print_vector("pbIndicies", cutil::stl_list_to_gsl_vector(CPBA::data->PARAM->pbIndicies));
	//Output::print(boost::format("NPB: %1%") % CPBA::data->PARAM->NPB);
	
	//previous version. this is now Nels's verison, which has more columns in it i think.
	//gsl_matrix *halo = cutil::read_matrix_from_file(boost::format("%1%/Halo Files/%2%_haloOriginal.txt") % CPBA::INPUT_PATH % (int)CPBA::data->BEAM->Eo);
	gsl_matrix *halo = 
	cutil::read_matrix_from_file(boost::format("%1%/Halo Files/%2%_haloOriginal%3%cm.txt") % CPBA::INPUT_PATH % (int)CPBA::data->BEAM->Eo % (int)CPBA::data->COLL->WX);
	
	//Output::print_matrix("halo", halo);

	int PBCNT = 0;
	bool repeatFlag = true;
	double sigVal = 4;

	list<double>::iterator xip;

	//begin.
	for(int ip = 0; ip < CPBA::data->CTMAT->NX; ip++)
	{
		double ctX = gsl_vector_get(CTMAT->X, ip);

		if(ctX < CPBA::data->COLL->Xmax && ctX > CPBA::data->COLL->Xmin)
		{
			PBCNT++;

			Output::print(boost::format("--> Simulating pencil beam %1% of %2%") % PBCNT % CPBA::data->PARAM->NPB);

			for(int j = 0; j < NZ; j++)
			{
				//Output::print(boost::str(boost::format("**J**: %1%") % j);
				double dmZ = gsl_vector_get(DM->Z, j);
				
				gsl_vector *yi_vec = cutil::get_matrix_col(CPBA::data->RAY->Zeff, ip);
				
				
				if(dmZ <= CTMAT->Zmax && dmZ >= CTMAT->Zmin)
				{
					double output = LinearInterpolation::evaluate_constant_linear(CTMAT->Z->size, CTMAT->Z, yi_vec, dmZ);
					gsl_matrix_set(this->Zeff, j, ip, output);
				}
				else if(dmZ < CTMAT->Zmin)
				{
					//Output::print(boost::format("Zeff 2 -> [%1%][%2%]") % j % ip);
					gsl_matrix_set(this->Zeff, j, ip, gsl_matrix_min(CPBA::data->RAY->Zeff));
				}
				else if(dmZ > CTMAT->Zmax)
				{
					//Output::print(boost::format("Zeff 3 -> [%1%][%2%]") % j % ip);
					gsl_matrix_set(this->Zeff, j, ip, gsl_matrix_max(CPBA::data->RAY->Zeff));
				}

				if(j < DM->NZ && repeatFlag == true)
				{
					ROIData *cr = CPBA::data->ROI->at(0);

					if(RayTrace::inROI(cr->XCon, cr->ZCon, ctX, gsl_vector_get(DM->Z, j+1)))
					{
						gsl_matrix_set(this->PDD, j, ip, gsl_vector_get(PDDmc, 0));
						repeatFlag = false;
					}
				}
				
				gsl_vector_free(yi_vec);

				//find PDD and other parameters.
				double zf = gsl_matrix_get(this->Zeff, j, ip);
				double mc_z_begin = gsl_vector_get(mc->z, 0);
				double mc_z_end = gsl_vector_get(mc->z, mc->z->size-1);

				

				if(zf != 0)
				{
					if(zf <= mc_z_end && zf >= mc_z_begin)
					{
						//Output::print(boost::str(boost::format("PDD 1 -> [%1%][%2%]") % j % ip);
						gsl_matrix_set(this->PDD, j, ip, CubicSplineInterpolation::evaluate_constant_spline(mc->z->size, mc->z, PDDmc, zf));
					}
					else if(zf < mc_z_begin)
					{
						//Output::print(boost::str(boost::format("PDD 2 -> [%1%][%2%]") % j % ip);
						gsl_matrix_set(this->PDD, j, ip, gsl_vector_get(PDDmc, 0));
					}
					else if(zf > mc_z_end)
					{
						//Output::print(boost::str(boost::format("PDD 3 -> [%1%][%2%]") % j % ip);
						gsl_matrix_set(this->PDD, j, ip, gsl_vector_get(PDDmc, PDDmc->size - 1));
					}
					
					/* NEW NELS CODE */

					if(zf <= gsl_matrix_get(halo, halo->size1 - 1, 0) && zf >= gsl_matrix_get(halo, 0, 0))
					{
						//Output::print(boost::str(boost::format("W,SIG 1 -> [%1%][%2%]") % j % ip);
						gsl_matrix_set(this->W1, j, ip, CubicSplineInterpolation::evaluate_constant_spline(halo->size1, cutil::get_matrix_col(halo, 0), cutil::get_matrix_col(halo, 1), zf));
						gsl_matrix_set(this->W2, j, ip, CubicSplineInterpolation::evaluate_constant_spline(halo->size1, cutil::get_matrix_col(halo, 0), cutil::get_matrix_col(halo, 3), zf));
						
						gsl_matrix_set(this->SIGhNUC1, j, ip, CubicSplineInterpolation::evaluate_constant_spline(halo->size1, cutil::get_matrix_col(halo, 0), cutil::get_matrix_col(halo, 2), zf));
						gsl_matrix_set(this->SIGhNUC2, j, ip, CubicSplineInterpolation::evaluate_constant_spline(halo->size1, cutil::get_matrix_col(halo, 0), cutil::get_matrix_col(halo, 4), zf));
					}
					else if(zf < gsl_matrix_get(halo, 0, 0))
					{
						//Output::print(boost::str(boost::format("W,SIG 2 -> [%1%][%2%]") % j % ip);
						gsl_matrix_set(this->W1, j, ip, gsl_matrix_get(halo, 0, 1));
						gsl_matrix_set(this->W2, j, ip, gsl_matrix_get(halo, 0, 3));
						
						gsl_matrix_set(this->SIGhNUC1, j, ip, gsl_matrix_get(halo, 0, 2));
						gsl_matrix_set(this->SIGhNUC2, j, ip, gsl_matrix_get(halo, 0, 4));
					}
					else if(zf > gsl_matrix_get(halo, halo->size1 - 1, 0))
					{
						//Output::print(boost::str(boost::format("W,SIG 3 -> [%1%][%2%]") % j % ip);
						gsl_matrix_set(this->W1, j, ip, gsl_matrix_get(halo, halo->size1 - 1, 1));
						gsl_matrix_set(this->W2, j, ip, gsl_matrix_get(halo, halo->size1 - 1, 3));
						
						gsl_matrix_set(this->SIGhNUC1, j, ip, gsl_matrix_get(halo, halo->size1 - 1, 2));
						gsl_matrix_set(this->SIGhNUC2, j, ip, gsl_matrix_get(halo, halo->size1 - 1, 4));
					}
					
					cutil::make_finite_non_neg(this->SIGhNUC1);
					cutil::make_finite_non_neg(this->SIGhNUC2);
					cutil::make_finite_non_neg(this->W1);
					cutil::make_finite_non_neg(this->W2);


					//ignore these tests for now.
				}

				// this shows up three or four times.
				//double RR_SS_1 = sqrt(pow(gsl_matrix_get(CPBA::data->RAY->SIGX1, j, ip), 2) + pow(gsl_matrix_get(this->SIGhNUC, j, ip), 2));
				double RR_SS_2 = sqrt(pow(gsl_matrix_get(CPBA::data->RAY->SIGX2, j, ip), 2) + pow(gsl_matrix_get(this->SIGhNUC1, j, ip), 2));

				//double SIGTEST = sigVal * RR_SS_1 + CTMAT->DX;
				
				/* NELS says to make this 100 */
				double SIGTEST = 100.0;
				
				for(int i = 0; i < NX; i++)
				{
					double x_i = gsl_vector_get(DM->X, i);
					double x_ip = gsl_vector_get(CTMAT->X, ip);

					double Xa = gsl_vector_get(CTMAT->X, ip) + (CTMAT->DX / 2.0) - gsl_vector_get(DM->X, i);
					double Xb = gsl_vector_get(CTMAT->X, ip) - (CTMAT->DX / 2.0) - gsl_vector_get(DM->X, i);

					double SIGP = gsl_matrix_get(CPBA::data->RAY->SIGX1, j, ip);
					double SIGPw = gsl_matrix_get(CPBA::data->RAY->SIGX2, j, ip);
					double SIGhN1 = gsl_matrix_get(this->SIGhNUC1, j, ip);
					double SIGhN2 = gsl_matrix_get(this->SIGhNUC2, j, ip);
					double SIGN1 = sqrt(pow(SIGP, 2) + pow(SIGhN1, 2));
					double SIGN2 = sqrt(pow(SIGP, 2) + pow(SIGhN2, 2));
					
					if(abs((double)(x_i - x_ip)) <= SIGTEST)
					{
						
						double FSCATp = .5 * (gsl_sf_erf(0.7071 * (Xa / SIGP)) - gsl_sf_erf(0.7071 * (Xb / SIGP)));
						double FSCATn1 = .5 * (gsl_sf_erf(0.7071 * (Xa / SIGN1)) - gsl_sf_erf(0.7071 * (Xb / SIGN1)));
						//double FSCATn2 = .5 * (gsl_sf_erf(0.7071 * (Xa / SIGN2)) - gsl_sf_erf(0.7071 * (Xb / SIGN2)));
						//double FSCATn2 = 1 / M_PI * (atan((x_ip + (CPBA::data->CTMAT->DX / 2.0) - gsl_matrix_get(DM-X, i)
						
						// jesus, i think i have this right? D:
						double FSn2_a = 1 / M_PI;
						double FSn2_b = (x_ip + (CPBA::data->CTMAT->DX / 2.0) - x_i) / gsl_matrix_get(this->SIGhNUC2, j, ip);
						double FSn2_c = (x_ip - (CPBA::data->CTMAT->DX / 2.0) - x_i) / gsl_matrix_get(this->SIGhNUC2, j, ip);
						
						double FSCATn2 = FSn2_a * (atan(FSn2_b) - atan(FSn2_c));
						
						if(!cutil::is_finite(FSCATp))
						{
							FSCATp = 0;
						}
						
						if(!cutil::is_finite(FSCATn1))
						{
							FSCATn1 = 0;
						}
						
						if(!cutil::is_finite(FSCATn2))
						{
							FSCATn2 = 0;
						}
						
						gsl_matrix_set(DM->PFluence, j, i, gsl_matrix_get(DM->PFluence, j, i) + FSCATp);
						gsl_matrix_set(DM->NFluence, j, i, gsl_matrix_get(DM->NFluence, j, i) + (FSCATn1 + FSCATn2));
						
						double bp_a = gsl_sf_erf(CPBA::data->COLL->WX / (2 * sqrt(2.0) * gsl_matrix_get(CPBA::data->RAY->SIGX2, j, ip)));
						double bn1_a = gsl_sf_erf(CPBA::data->COLL->WX / (2 * sqrt(2.0) * RR_SS_2));
						
						if(!cutil::is_finite(bp_a))
						{
							bp_a = 1;
						}
						
						if(!cutil::is_finite(bn1_a))
						{
							bn1_a = 1;
						}
						
						
						gsl_vector_set(this->BP, j, 1 / bp_a);
						gsl_vector_set(this->BN1, j, 1 / bn1_a);
						double BN2_T = 1 / (( 2 / M_PI ) * atan(CPBA::data->COLL->WX / ( 2 * gsl_matrix_get(this->SIGhNUC2, j, ip))));
						gsl_vector_set(this->BN2, j, BN2_T);

						//Output::print(boost::str(boost::format("BP: %1%, BN: %2%") % BP % BN);

						double pdose = gsl_matrix_get(DM->PDose, j, i) + (1 - gsl_matrix_get(this->W1, j, ip) - gsl_matrix_get(this->W2, j, ip)) 
						* FSCATp * gsl_matrix_get(this->PDD, j, ip) * gsl_vector_get(this->BP, j) * gsl_matrix_get(CPBA::data->RAY->SPratio, j, ip);
							
						double ndose = gsl_matrix_get(DM->NDose, j, i) + gsl_matrix_get(this->W1, j, ip)
							* FSCATn1 * gsl_matrix_get(this->PDD, j, ip) * gsl_vector_get(this->BN1, j) 
							* gsl_matrix_get(CPBA::data->RAY->SPratio, j, ip) + gsl_matrix_get(this->W2, j, ip) 
						* FSCATn2 * gsl_matrix_get(this->PDD, j, ip) * gsl_vector_get(this->BN2, j) * gsl_matrix_get(CPBA::data->RAY->SPratio, j, ip);
							
						double dose = pdose + ndose;

						if(!cutil::is_finite(pdose))
						{
							pdose = 0;
						}

						if(!cutil::is_finite(ndose))
						{
							ndose = 0;
						}

						if(!cutil::is_finite(dose))
						{
							dose = 0;
						}

						gsl_matrix_set(DM->PDose, j, i, pdose);
						gsl_matrix_set(DM->NDose, j, i, ndose);
						gsl_matrix_set(DM->Dose, j, i, dose);
						
					}
				}
			}
		}
	}

	//Output::print_matrix("DM->PDose", DM->PDose);
	//Output::print_matrix("DM->NDose", DM->NDose);
	//Output::print_matrix("DM->Dose", DM->Dose);
	//Output::print_matrix("DM->PFluence", DM->PFluence);
	//Output::print_matrix("DM->NFluence", DM->NFluence);
	//Output::print_matrix("this->PDD", this->PDD);
	//Output::print_matrix("this->W", this->W);
	//Output::print_matrix("this->SIGhNUC", this->SIGhNUC);
	//Output::print_matrix("this->Zeff", this->Zeff);
	
	Output::print("CPBA::calculateDose() finished");
}

void CPBA::finalizeData()
{
	double scale_factor;
	vector<gsl_matrix *> xx_zz_c;
	Output::print("CPBA::finalizeData()");
	xx_zz_c = BilinearInterpolation::meshgrid(CPBA::data->DOSEMAT->X, CPBA::data->DOSEMAT->Z);

	gsl_matrix *xx = xx_zz_c.at(0);
	gsl_matrix *zz = xx_zz_c.at(1);

	//Output::print_matrix("xx", xx);
	//Output::print_matrix("zz", zz);

	this->Da = gsl_matrix_calloc(CPBA::data->DOSEMAT->Dose->size1, CPBA::data->DOSEMAT->Dose->size2);
	gsl_matrix_memcpy(this->Da, CPBA::data->DOSEMAT->Dose);

	scale_factor = (1 / gsl_matrix_max(CPBA::data->DOSEMAT->Dose)) * 100.0;

	gsl_matrix_scale(this->Da, scale_factor);

	//xm, zm (from montecarlo data)
	MonteCarlo *mc = new MonteCarlo((int)CPBA::data->BEAM->Eo, (int)CPBA::data->COLL->WX, (int)CPBA::data->PARAM->ThetaDeg, (int)CPBA::data->PARAM->Hstep);

	if(CPBA::data->PARAM->ThetaDeg == 0)
	{
		//2DInterpolation::evalulate_2d_bilinear(gsl_vector* xAxis,gsl_vector* yAxis, gsl_matrix* zSurface, gsl_matrix *xi, gsl_matrix *yi)
		this->Db = BilinearInterpolation::evalulate_2d_bilinear(mc->x, mc->z, mc->D, xx, zz);
		
		scale_factor = (1 / gsl_matrix_max(this->Db)) * 100.0;

		gsl_matrix_scale(this->Db, scale_factor);

	}
	else
	{
		gsl_matrix *xr = gsl_matrix_calloc(CPBA::data->CTMAT->NZ, CPBA::data->CTMAT->NX);
		gsl_matrix *zr = gsl_matrix_calloc(CPBA::data->CTMAT->NZ, CPBA::data->CTMAT->NX);

		for(unsigned int j = 0; j < xr->size1; j++)
		{
			for(unsigned int i = 0; i < xr->size2; i++)
			{
				double xxr = gsl_vector_get(CPBA::data->CTMAT->X, i) * cos(0.0 - CPBA::data->PARAM->Theta) + gsl_vector_get(CPBA::data->CTMAT->Z, j) * sin(0.0 - CPBA::data->PARAM->Theta);
				double zzr = (0 - gsl_vector_get(CPBA::data->CTMAT->X, i)) * sin(0.0 - CPBA::data->PARAM->Theta) + gsl_vector_get(CPBA::data->CTMAT->Z, j) * cos(0.0 - CPBA::data->PARAM->Theta);

				gsl_matrix_set(xr, j, i, xxr);
				gsl_matrix_set(zr, j, i, zzr);
			}
		}

		this->Db = BilinearInterpolation::evalulate_2d_bilinear(mc->x, mc->z, mc->D, xr, zr);
		
		scale_factor = gsl_matrix_max(this->Db) * 100.0;

		gsl_matrix_scale(this->Db, scale_factor);

		gsl_matrix_free(xr);
		gsl_matrix_free(zr);
	}

	//Output::print_matrix("this->Da", this->Da);
	//Output::print_matrix("this->Db", this->Db);

	Output::print("CPBA::finalizeData() finished");
}

gsl_vector *CPBA::extractPDD(gsl_matrix *D)
{
	int xDim = D->size2;
	int midX;
	gsl_vector *PDD;
	gsl_vector *PDD2;

	if(xDim % 2 == 1)
	{
		// we need to be off by one on the negative side as compared to the matlab version since matlab indexes from 1:N and 
		// C indexes to 0:N-1
		midX = (xDim - 1) / 2;
		PDD = cutil::get_matrix_col(D, midX);
	}
	else
	{
		//be off by one to match matlab.
		midX = (xDim / 2) - 1;

		PDD = cutil::get_matrix_col(D, midX);
		PDD2 = cutil::get_matrix_col(D, midX+1);
		//Output::print_vector("PDD", PDD);
		//Output::print_vector("PDD2", PDD2);

		for(unsigned int i = 0; i < PDD->size; i++)
		{
			double val1 = gsl_vector_get(PDD, i);
			double val2 = gsl_vector_get(PDD2, i);

			gsl_vector_set(PDD, i, 0.5 * (val1 + val2));
		}

		//Output::print_vector("PDD", PDD);

	}

	return PDD;
}

void CPBA::save()
{
	string file;
	
	int pinnacle_mode = CPBA::PROG_ARGS["pinnacle"].as<int>();
	
	if(pinnacle_mode == 0)
	{

		/* we need to write out our files */
		boost::format ending("E%1%_FS%2%_TILT%3%_STEP%4%");
		ending % (int)CPBA::data->BEAM->Eo % (int)CPBA::data->COLL->WX % (int)CPBA::data->PARAM->ThetaDeg % (int)CPBA::data->PARAM->Hstep;

		boost::format ff("Output/PBA_%1%");
		ff % ending.str();

		boost::filesystem::path out_dir = ff.str();

		if(boost::filesystem::exists(out_dir))
		{
			boost::filesystem::remove_all(out_dir);
		}

		boost::filesystem::create_directories(out_dir);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_Dose" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->Dose);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_x" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->X);
		
		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_z" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->Z);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_A0" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->A0);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_A1" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->A1);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_B" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->B);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_CTN" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->CTN);
		
		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_energy" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->Energy);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_energy2" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->Sz);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_NDOSE" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->NDose);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_PDOSE" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->PDose);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_PBindices" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, cutil::stl_list_to_gsl_vector(CPBA::data->PARAM->pbIndicies));

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_rho" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->rho);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_scatPwr" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->SIGTH);
		
		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_scatPwrH2O" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->SIGTH2);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_SIG" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->SIGX1);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_SIGH2O" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->SIGX2);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_SPratio" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->SPratio);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_Zeff" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->RAY->Zeff);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_W1nuc" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->W1);
		
		
		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_W2nuc" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->W2);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_SIGnuc1" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->SIGhNUC1);
		
		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_SIGnuc2" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->SIGhNUC2);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_TotFluence" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->Fluence);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_pFluence" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->PFluence);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_nFluence" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, CPBA::data->DOSEMAT->NFluence);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_PDD" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->PDD);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_Zp" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->Zp);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_Zn" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->Zn);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_pNorm" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->pNorm);

		file = boost::str(boost::format("%1%/%2%_%3%.txt") % ff.str() % "PBA_nNorm" % ending);
		Output::print(file);
		cutil::write_data_to_file(file, this->nNorm);
	}
	
	if(pinnacle_mode == 1)
	{
	
		// pinnacle output files.
		
		//CPBA::Z_DEPTH = vm["zdepth"].as<double>();
		
		string root_path = CPBA::PROG_ARGS["pinnaclepath"].as<string>();
		
		boost::filesystem::path pinn_dir = root_path;
		
		if(!boost::filesystem::exists(pinn_dir))
		{
			boost::filesystem::create_directories(pinn_dir);
		}
		
		file = boost::str(boost::format("%1%/CPBA.Script") % root_path);
		Output::print(file);
		Pinnacle::write_store(file);
		
		
		file = boost::str(boost::format("%1%/PBA_PinnacleDoseGrid.bin") % root_path);
		Output::print(file);
		Pinnacle::write_dosegrid(file);
	}

	/*
	'PBA_Dose',s.DOSEMAT.Dose
	'PBA_x',s.DOSEMAT.X,
	'PBA_z's.DOSEMAT.Z,
	'MC_Dose',Dbsave,
	'PBA_A0',s.RAY.A0
	'PBA_A1',s.RAY.A1,
	'PBA_A2',s.RAY.A2
	'PBA_B',s.RAY.B,
	'PBA_CTN',s.RAY.CTN,
	'PBA_energy',s.RAY.Energy
	'PBA_energy2',s.RAY.Sz
	'PBA_NDOSE',s.DOSEMAT.NDose
	'PBA_PDOSE',s.DOSEMAT.PDose
	'PBA_PBindices',s.PARAM.pbIndices
	'PBA_rho',s.RAY.rho
	'PBA_scatPwr',s.RAY.SIGTH
	'PBA_scatPwrH2O',s.RAY.SIGTH2
	'PBA_SIG',s.RAY.SIGX1
	'PBA_SIGH2O',s.RAY.SIGX2
	'PBA_SPratio',s.RAY.SPratio
	'PBA_Zeff',s.RAY.Zeff
	'PBA_Wnuc',W
	'PBA_SIGnuc',SIGhNUC
	'PBA_TotFluence',s.RAY.Zeff
	'PBA_pFluence',s.DOSEMAT.Fluence
	'PBA_nFluence',s.DOSEMAT.NFluence
	'PBA_PDD',PDD
	'PBA_Zp',Zp
	'PBA_Zn',Zn
	'PBA_pNorm',pNorm
	'PBA_nNorm',nNorm
	*/
}
