//
//  pinnacle.cpp
//  cpba_toolset
//
//  Created by Michael Thomas on 1/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "pinnacle.h"
#include "cpba.h"
#include "out.h"

#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

using namespace std;

void Pinnacle::write_store(string filename)
{
	/*
	 At .TestDimX = Float {
	 Value = 521;
	 };
	 At .TestDimY = Float {
	 Value = 120;
	 };
	*/
	
	boost::format ss;
	
	ofstream file;
	
	if(boost::filesystem::exists(filename))
	{
		boost::filesystem::remove(filename);
	}
	
	file.open(filename.c_str());
		
	file << boost::str(boost::format("TrialList.Current.DoseGrid.Dimension.X=%1%;\n") % (CPBA::data->DOSEMAT->NX));
	file << boost::str(boost::format("TrialList.Current.DoseGrid.Dimension.Y=%1%;\n") % (CPBA::data->DOSEMAT->NZ));
	file << boost::str(boost::format("TrialList.Current.DoseGrid.Origin.X=%1%;\n") % gsl_vector_get(CPBA::data->DOSEMAT->X, 0));
	file << boost::str(boost::format("TrialList.Current.DoseGrid.Origin.Y=%1%;\n\n\n") % gsl_vector_get(CPBA::data->DOSEMAT->Z, 0));	
	
	file.close();
}


void Pinnacle::write_dosegrid(string filename)
{	
	if(boost::filesystem::exists(filename))
	{
		boost::filesystem::remove(filename);
	}

	ofstream binary_file(filename.c_str(), ios::out | ios::binary);
	
	int z_depth = CPBA::PROG_ARGS["zdepth"].as<int>();
	
	for(int y = 0; y < z_depth; y++)
	{
		for(int z = 0; z < CPBA::data->DOSEMAT->NZ; z++)
		{
			for(int x = 0; x < CPBA::data->DOSEMAT->NX; x++)
			{
				float vv = (float)gsl_matrix_get(CPBA::data->DOSEMAT->Dose, z, x);
				binary_file.write((char *)&vv, sizeof(float));
			}
		}
	}
	
	binary_file.close();
}


/* used by dose_scale */
dose_volume *Pinnacle::read_and_scale_dosegrid(string filename, int z, int y, int x, double scale_factor)
{
	if(!boost::filesystem::exists(filename))
	{
		Output::print(boost::format("Pinnacle::read_and_scale_dosegrid()::error, %1% does not exist or is not readable") % filename);
		exit(1);
	}
	
	ifstream binary_file(filename.c_str(), ios::in | ios::binary);
	
	int z_depth = CPBA::PROG_ARGS["z"].as<int>();
	int y_depth = CPBA::PROG_ARGS["y"].as<int>();
	int x_depth = CPBA::PROG_ARGS["x"].as<int>();
	
	//array_type A(boost::extents[3][4][2]); 
	
	dose_volume *DV = new dose_volume();
	DV->resize(z_depth);
	
	float vv;
	
	for(int z = 0; z < z_depth; z++)
	{
		gsl_matrix *m = gsl_matrix_calloc(y_depth, x_depth);
		
		DV->at(z) = m;
		
		for(int y = 0; y < y_depth; y++)
		{
			for(int x = 0; x < x_depth; x++)
			{
				binary_file.read((char *)&vv, sizeof(float));
				vv = vv * scale_factor;
				gsl_matrix_set(m, y, x, vv);
			}
		}
	}
	
	binary_file.close();
	
	return DV;
}

void Pinnacle::write_computed_dosegrid(dose_volume *DV, string filename)
{	
	if(boost::filesystem::exists(filename))
	{
		boost::filesystem::remove(filename);
	}
	
	ofstream binary_file(filename.c_str(), ios::out | ios::binary);
	
	int z_depth = DV->size();
	
	for(int z = 0; z < z_depth; z++)
	{
		gsl_matrix *m = (gsl_matrix *)DV->at(z); 
		
		for(int y = 0; y < m->size1; y++)
		{
			for(int x = 0; x < m->size2; x++)
			{
				float vv = (float)gsl_matrix_get(m, y, x);
				binary_file.write((char *)&vv, sizeof(float));
			}
		}
	}
	
	binary_file.close();
	
}
