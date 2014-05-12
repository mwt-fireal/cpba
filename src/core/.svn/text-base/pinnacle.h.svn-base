//
//  pinnacle.h
//  cpba_toolset
//
//  Created by Michael Thomas on 1/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef cpba_toolset_pinnacle_h
#define cpba_toolset_pinnacle_h

#include<boost/multi_array.hpp>
#include<string>
#include<vector>
#include "cpba.h"
#include "data.h"

using namespace std;

typedef vector<gsl_matrix *> dose_volume;

class Pinnacle
{
public:
	static void write_store(string filename);
	static void write_dosegrid(string filename);
	
	/* used by dose_scale */
	static dose_volume *read_and_scale_dosegrid(string filename, int z, int y, int x, double scale_factor);
	static void write_computed_dosegrid(dose_volume *DV, string filename);
};

#endif
