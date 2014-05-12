#include "mat.h"
#include "util.h"
#include "out.h"

#include <map>
#include <list>
#include <vector>
#include <sstream>
#include <boost/filesystem.hpp>

using namespace std;

PBAMaterials::PBAMaterials(string sourceFile)
{
	this->mdb = new map<string, Material*>;
	
	boost::filesystem::path p = sourceFile;
	
	if(!boost::filesystem::exists(p))
	{
		Output::print(boost::format("cutil::read_matrix_from_file::error, matrix file %1% does not exist or is not readable") % sourceFile);
		exit(1);
	}

	//Output::print(sourceFile);
	FILE *fptr = fopen(sourceFile.c_str(), "r");

	char *line = cutil::fast_get_line(fptr);

	//bool in_material = false;

	/* laziness */
	char buf1[100];
	char buf2[100];
	int ibuf1 = 0;
	int ibuf2 = 0;
	double dbuf1 = 0.0;
	double dbuf2 = 0.0;

	char obuf[1000];

	while(line != (char *)NULL)
	{
		if(line[0] != '-' && line[0] > 32)
		{
			// new symbol

			sscanf(line, "%s %s %lf %d", buf1, buf2, &dbuf1, &ibuf1);

			int num_elements = ibuf1;

			Material *m = new Material(buf1, buf2, dbuf1, ibuf1);

			for(int i = 0; i < num_elements; i++)
			{
				line = cutil::fast_get_line(fptr);
				//Output::print(line);

				/* H 2 1 1.007940 0.111894 */
				sscanf(line, "%s %d %d %lf %lf", buf1, &ibuf1, &ibuf2, &dbuf1, &dbuf2);

				map<string, string> mc;

				mc.insert(pair<string, string>("symbol", string(buf1)));
				mc.insert(pair<string, string>("elNum", cutil::tostring(ibuf1)));
				mc.insert(pair<string, string>("Z", cutil::tostring(ibuf2)));
				mc.insert(pair<string, string>("A", cutil::tostring(dbuf1)));
				mc.insert(pair<string, string>("f", cutil::tostring(dbuf2)));

				m->elements.push_back(mc);
			}

			sprintf(obuf, "inserted symbol: %s\n", m->symbol.c_str());
			//Output::print(obuf);
			this->mdb->insert(pair<string, Material*>(m->symbol, m));
		}

		free(line);
		line = cutil::fast_get_line(fptr);
	}

	if(line != (char *)NULL)
	{
		free(line);
	}
}
