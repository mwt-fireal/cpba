#pragma once

#ifndef DATA_H
#define DATA_H

using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <vector>
#include <list>
#include <string>

#include "mat.h"
#include "util.h"

class Beam
{
public:
	double Eo;
	double Rp;
};

class Collimator
{
public:
	double WX;
	double Xmin;
	double Xmax;

	Collimator()
	{
		this->WX = 0.0;
		this->Xmin = 0.0;
		this->Xmax = 0.0;
	}
};

class PBAParameters
{
public:
	double Theta;
	double ThetaDeg;
	double FMCS;
	double Hstep;
	double DX;
	double DZ;
	double NPB;
	
	// may change.
	list<double> *pbIndicies;
};

class ROIData
{
public:
	double NCont;
	gsl_vector *XCon;
	gsl_vector *ZCon;
	double NPts;
	double CTN;
	string MAT;
	gsl_vector *Xcorner;
	gsl_vector *Zcorner;

	ROIData()
	{
		this->NCont = 0;
		this->NPts = 0;
		this->CTN = 0;
	}
};

class Grid
{
public:
	int NX;
	double Xmin;
	double Xmax;

	//probably wrong
	gsl_vector *X;

	int NZ;
	double Zmin;
	double Zmax;

	gsl_vector *Z;

	gsl_vector *Xcorner;
	gsl_vector *Zcorner;
	double DX;
	double DZ;

	Grid()
	{
		this->NX = 0;
		this->Xmin = 0.0;
		this->Xmax = 0.0;
		
		this->NZ = 0;
		this->Zmin = 0.0;
		this->Zmax = 0.0;

		this->DX = 0.0;
		this->DZ = 0.0;

		this->X = NULL;
		this->Z = NULL;
		this->Xcorner = NULL;
		this->Zcorner = NULL;
	}

	static Grid* duplicate(Grid *oo)
	{
		Grid *ng = new Grid;
		ng->NX = oo->NX;
		ng->Xmin = oo->Xmin;
		ng->Xmax = oo->Xmax;

		ng->X = cutil::vec_operation_copy(oo->X);
		ng->NZ = oo->NZ;
		ng->Zmin = oo->Zmin;
		ng->Zmax = oo->Zmax;

		ng->Z = cutil::vec_operation_copy(oo->Z);
		ng->Xcorner = cutil::vec_operation_copy(oo->Xcorner);
		ng->Zcorner = cutil::vec_operation_copy(oo->Zcorner);
		ng->DX = oo->DX;
		ng->DZ = oo->DZ;

		return ng;
	}
};

class Ray
{
public:
	gsl_matrix *CTN;
	vector<vector<string> > MAT;
	gsl_matrix *Zeff;
	gsl_matrix *rho;
	gsl_matrix *Energy;
	gsl_matrix *Energy2;

	gsl_matrix *Sz;
	gsl_matrix *Sz2;
	gsl_matrix *SPratio;
	gsl_matrix *A0;
	gsl_matrix *A1;
	gsl_matrix *A2;
	gsl_matrix *B;
	gsl_matrix *SIGTH;
	gsl_matrix *SIGTH2;
	gsl_matrix *SIGX1;
	gsl_matrix *SIGX2;

	Ray(int NZ, int NX)
	{
		this->CTN = gsl_matrix_calloc(NZ, NX);
		this->Zeff = gsl_matrix_calloc(NZ, NX);
		this->rho = gsl_matrix_calloc(NZ, NX);
		this->Energy = gsl_matrix_calloc(NZ, NX);
		this->Energy2 = gsl_matrix_calloc(NZ, NX);
		this->SPratio = gsl_matrix_calloc(NZ, NX);
		this->Sz = gsl_matrix_calloc(NZ, NX);
		this->Sz2 = gsl_matrix_calloc(NZ, NX);
		this->A0 = gsl_matrix_calloc(NZ, NX);
		this->A1 = gsl_matrix_calloc(NZ, NX);
		this->A2 = gsl_matrix_calloc(NZ, NX);
		this->B = gsl_matrix_calloc(NZ, NX);
		this->SIGTH = gsl_matrix_calloc(NZ, NX);
		this->SIGTH2 = gsl_matrix_calloc(NZ, NX);
		this->SIGX1 = gsl_matrix_calloc(NZ, NX);
		this->SIGX2 = gsl_matrix_calloc(NZ, NX);

		this->MAT.resize(NZ);
		for(int i = 0; i < NZ; i++)
		{
			this->MAT.at(i).resize(NX);

			//this->MAT[i] = new vector<string>;
			//this->MAT[i]->resize(NX);

			for(int j = 0; j < NX; j++)
			{
				this->MAT[i][j] = "VACUUM";
			}
		}
	}
};

class DoseMatrix
{
public:
	gsl_matrix *Dose;
	gsl_matrix *PDose;
	gsl_matrix *NDose;
	gsl_matrix *Fluence;
	gsl_matrix *PFluence;
	gsl_matrix *NFluence;

	int NX;
	double Xmin;
	double Xmax;

	gsl_vector *X;

	int NZ;
	double Zmin;
	double Zmax;

	gsl_vector *Z;

	gsl_vector *Xcorner;
	gsl_vector *Zcorner;
	double DX;
	double DZ;

	DoseMatrix()
	{
		this->NX = 0;
		this->Xmin = 0.0;
		this->Xmax = 0.0;
		
		this->NZ = 0;
		this->Zmin = 0.0;
		this->Zmax = 0.0;

		this->DX = 0.0;
		this->DZ = 0.0;

		this->X = NULL;
		this->Z = NULL;
		this->Xcorner = NULL;
		this->Zcorner = NULL;

		this->Dose = NULL;
		this->PDose = NULL;
		this->NDose = NULL;
		this->Fluence = NULL;
		this->PFluence = NULL;
		this->NFluence = NULL;
	}

	static DoseMatrix* duplicate(Grid *oo)
	{
		DoseMatrix *ng = new DoseMatrix;
		ng->NX = oo->NX;
		ng->Xmin = oo->Xmin;
		ng->Xmax = oo->Xmax;

		ng->X = cutil::vec_operation_copy(oo->X);
		ng->NZ = oo->NZ;
		ng->Zmin = oo->Zmin;
		ng->Zmax = oo->Zmax;

		ng->Z = cutil::vec_operation_copy(oo->Z);
		ng->Xcorner = cutil::vec_operation_copy(oo->Xcorner);
		ng->Zcorner = cutil::vec_operation_copy(oo->Zcorner);
		ng->DX = oo->DX;
		ng->DZ = oo->DZ;

		return ng;
	}
};

class PBAData
{
public:
	clock_t start_time;
	void *cpba;
	PBAMaterials *MaterialDB;
	Beam *BEAM;
	Collimator *COLL;
	PBAParameters *PARAM;
	vector<ROIData *> *ROI;
	Grid *GRID;
	Grid *AXIS;
	Ray *RAY;
	DoseMatrix *DOSEMAT;
	Grid *CTMAT;
};

#endif
