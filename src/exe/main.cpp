#include <math.h>
#include <time.h>

#include "cpba.h"
#include "util.h"
#include "monte.h"
#include "interp.h"
#include "raytrace.h"
#include "init.h"
#include "out.h"

#include <sstream>

#include <gsl/gsl_sf.h>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

using namespace std;

int main(int argc, char **argv)
{
    int opt;
    
	#ifdef _WINDOWS
		_CrtSetDbgFlag(_CRTDBG_CHECK_DEFAULT_DF);
	#endif
	
	clock_t start = clock();
	double elapsed;
    
	namespace po = boost::program_options;

	po::options_description cmdline_options;
	
	po::options_description generic("Generic Options");
	po::options_description settings("CPBA Settings");
	po::options_description slab("SLAB PHANTOM Settings");
	po::options_description flat("FLAT PHANTOM Settings");
	
	generic.add_options()
		("help", "produce help message")
	("zdepth", po::value<int>()->default_value(5), "pinnacle z depth")
	("pinnaclepath", po::value<string>()->default_value("/tmp/cpba"), "pinnacle store path")
	("pinnacle", po::value<int>()->default_value(0), "operate in pinnacle mode (save only pinnacle files)")
	("inputpath", po::value<string>()->default_value("./InputData"), "path to input directory");
	
	settings.add_options()
		("eo", po::value<double>()->default_value(100.0), "Energy (MeV)")
		("wx", po::value<double>()->default_value(4.0), "Field with at final collimator (cm)")
		("fmcs", po::value<double>()->default_value(1.0), "MCS correction factor")
		("thetadeg", po::value<double>()->default_value(0.0), "Angle between axis and phantom surface normal (degrees)")
		("hstep", po::value<double>()->default_value(4.0), "Height of surface irregularity (cm)")
		("dx", po::value<double>()->default_value(0.5), "Pencil beam width (cm)")
		("dz", po::value<double>()->default_value(0.5), "Pencil beam step size (cm)")
	("material", po::value<string>()->default_value("SLAB"), "Material phantom type (SLAB/FLAT)");
	
	slab.add_options()
		("smallest_x", po::value<string>()->default_value("-Inf"), "Smallest x-value of slab (cm) [+/- Inf allowed]")
		("largest_x", po::value<string>()->default_value("+Inf"), "Largest x-value of slab (cm) [+/-Inf allowed]")
		("thickness", po::value<double>()->default_value(2.0), "Thickness of SLAB PHANTOM (cm)")
		("depth", po::value<double>()->default_value(2.0), "Depth of SLAB PHANTOM (cm)")
		("ctnum_slab", po::value<int>()->default_value(900), "CT # of SLAB PHANTOM (100=AIR, 500=WATER, 900=BONE)")
		("ctnum_bk", po::value<int>()->default_value(500), "CT # of SLAB PHANTOM background (100=AIR, 500=WATER, 900=COMPACTBONE)")
		("material_slab", po::value<string>()->default_value("COMPACTBONE"), "SLAB PHANTOM material (AIR, WATER, COMPACTBONE)")
	("material_bk", po::value<string>()->default_value("WATER"), "SLAB PHANTOM background material (AIR, WATER, COMPACTBONE)");

	flat.add_options()
		("ctnum_flat", po::value<int>()->default_value(500), "CT # of FLAT PHANTOM (100=AIR, 500=WATER, 900=COMPACTBONE)")
	("material_flat", po::value<string>()->default_value("WATER"), "FLAT PHANTOM material (AIR, WATER, COMPACTBONE)");
    
	cmdline_options.add(generic).add(settings).add(slab).add(flat);

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
	po::notify(vm);
    
	if (vm.count("help"))
	{
		cout << cmdline_options << endl;
		return 1;
	}

	char *input_env = getenv("CPBA_INPUT_PATH");

	boost::filesystem::path p;

	if(input_env != NULL)
	{
		Output::print(boost::format("CPBA_INPUT_PATH: %1%") % input_env);
		p = string(input_env);
	}
	else
	{
		p = vm["inputpath"].as<string>();
	}

	
	CPBA::PROG_ARGS = vm;
	
	CPBA::INPUT_PATH = boost::filesystem::absolute(p).generic_string();
	
	Output::print(boost::format("CPBA::main(): Input Path --> %1%") % CPBA::INPUT_PATH);
	
	/*
	 150 
	 10 
	 1 
	 0 
	 0 
	 5.0
	 5.0
	 SLAB:-Inf:+Inf:2:2:500:WATER:500:WATER
	 */
	
	/*
	 //SIMPLE SIMULATION
	 CPBAParameters *cp = new CPBAParameters();
	 cp->Eo = 150;
	 cp->WX = 10;
	 cp->FMCS = 1;
	 cp->ThetaDeg = 0; //45
	 cp->Hstep = 0;
	 cp->DX = 5.0;
	 cp->DZ = 5.0;
	 */
	
	
	/*
	 //STANDARD SIMULATION
	 CPBAParameters *cp = new CPBAParameters();
	 cp->Eo = 150;
	 cp->WX = 4; //FS
	 cp->FMCS = 1;
	 cp->ThetaDeg = 0;
	 cp->Hstep = 0;
	 cp->DX = 0.1;
	 cp->DZ = 0.1;
	 */
	
	//EXTENDED SIMULATION
	CPBAParameters *cp = new CPBAParameters();
	cp->Eo = vm["eo"].as<double>();
	cp->WX = vm["wx"].as<double>();
	cp->FMCS = vm["fmcs"].as<double>();
	cp->ThetaDeg = vm["thetadeg"].as<double>();
	cp->Hstep = vm["hstep"].as<double>();
	cp->DX = vm["dx"].as<double>();
	cp->DZ = vm["dz"].as<double>();
	
	
	CPBA *cpba = new CPBA(cp);
	
	CPBA::data->start_time = start;
	
	CPBAMaterialParameters *mp = new CPBAMaterialParameters();
	
	/*
	 mp->type = CPBAMaterialParameters::CPBA_SLAB;
	 mp->sx = numeric_limits<double>::min();
	 mp->ex = numeric_limits<double>::max();
	 mp->thickness = 2.0;
	 mp->depth = 2.0;
	 mp->CTNUM_SLAB = 500;
	 mp->MATERIAL_SLAB = "WATER";
	 mp->CTNUM_BK = 500;
	 mp->MATERIAL_BK = "WATER";
	 */
	
	//PHANTOM
	
	//Output::print(boost::format("smallest_x: %1%") % vm["smallest_x"].as<string>());
	//Output::print(boost::format("largest_x: %1%") % vm["largest_x"].as<string>());
	
	//Output::print(boost::format("smallest_x, doing compare against -Inf, %1%") % vm["smallest_x"].as<string>().compare("-Inf"));
	//Output::print(boost::format("largest_x, doing compare against +Inf, %1%") % vm["largest_x"].as<string>().compare("+Inf"));
	
	if(vm["material"].as<string>().compare("FLAT") == 0)
	{
		mp->type = CPBAMaterialParameters::CPBA_PHANTOM;
		mp->CTNUM_PHANTOM = vm["ctnum_flat"].as<int>();
		mp->MATERIAL_PHANTOM = vm["material_flat"].as<string>();
	}
	else if(vm["material"].as<string>().compare("SLAB") == 0)
	{
		mp->type = CPBAMaterialParameters::CPBA_SLAB;
		
		if(vm["smallest_x"].as<string>().compare("-Inf") == 0)
		{
			Output::print("smallest_x set to numeric_limits::min");
			mp->sx = numeric_limits<double>::min();
		}
		else if(vm["smallest_x"].as<string>().compare("+Inf") == 0)
		{
			Output::print("smallest_x set to numeric_limits::max");
			mp->sx = numeric_limits<double>::max();
		}
		else
		{
			mp->sx = vm["smallest_x"].as<double>();
		}
		
		if(vm["largest_x"].as<string>().compare("-Inf") == 0)
		{
			Output::print("largest_x set to numeric_limits::min");
			mp->ex = numeric_limits<double>::min();
		}
		else if(vm["largest_x"].as<string>().compare("+Inf") == 0)
		{
			Output::print("largest_x set to numeric_limits::max");
			mp->ex = numeric_limits<double>::max();
		}
		else
		{
			mp->ex = vm["largest_x"].as<double>();
		}
		
		mp->thickness = vm["thickness"].as<double>();
		mp->depth = vm["depth"].as<double>();
		mp->CTNUM_SLAB = vm["ctnum_slab"].as<int>();
		mp->MATERIAL_SLAB = vm["material_slab"].as<string>();
		mp->CTNUM_BK = vm["ctnum_bk"].as<int>();
		mp->MATERIAL_BK = vm["material_bk"].as<string>();
	}
	else
	{
		Output::print(boost::format("Error: unknown material type \"%1%\"") % vm["material"].as<string>());
		exit(1);
	}
	
	
	cpba->init(mp);
	
	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Initialization: %1% seconds") % elapsed);
	
	cpba->execute();
	
	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Execution: %1% seconds") % elapsed);
	
	cpba->save();
	
	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Total Execution Time: %1% seconds") % elapsed);
}
