#pragma once

#include <string>

#ifndef INIT_H
#define INIT_H

using namespace std;

#include <gsl/gsl_matrix.h>

#include "data.h"

class CPBAMaterialParameters
{
public:
	static const int CPBA_SLAB = 0;
	static const int CPBA_PHANTOM = 1;

	int type;

	//SLAB
	double sx;
	double ex;
	double thickness;
	double depth;
	
	int CTNUM_SLAB;
	string MATERIAL_SLAB;
	
	int CTNUM_BK;
	string MATERIAL_BK;

	//PHANTOM
	int CTNUM_PHANTOM;
	string MATERIAL_PHANTOM;
};

class PBAInitialization
{
public:
	static void createSlab(double sx, double ex, double thickness, double depth, int CTNUM_SLAB, string material_SLAB, int CTNUM_BK, string material_BK);
	static void createPhantom(int CTNUM, string material);
	static void processROIs();
	static void createAxis();
	static void cartesianGrid(Grid *axis, double Xmin, double Xmax, double Zmin, double Zmax, double DX, double DZ);
	static void createGrids();
};


#endif
