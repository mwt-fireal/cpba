#include "out.h"
#include "cpba.h"
#include "pinnacle.h"

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
	
	po::options_description generic("Dose Volume Scaling Options");
	
	generic.add_options()
		("help", "produce help message")
		("inputfile", po::value<string>()->default_value("/tmp/DoseVolume"), "pinnacle dose volume file")
		("outputdir", po::value<string>()->default_value("/tmp/cpba"), "pinnacle store path")
		("index", po::value<string>()->default_value("0001"), "unique write index")
		("factor", po::value<double>()->default_value(2.0), "scale factor")
		("x", po::value<int>()->default_value(10), "X dimension")
		("y", po::value<int>()->default_value(10), "Y dimension")
	("z", po::value<int>()->default_value(10), "Z dimension");
	
	cmdline_options.add(generic);
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
	po::notify(vm);
    
	if (vm.count("help"))
	{
		cout << cmdline_options << endl;
		return 1;
	}
	
	boost::filesystem::path p = vm["inputfile"].as<string>();
	
	if(!boost::filesystem::exists(p))
	{
		Output::print(boost::format("error, DoseVolume file %1% does not exist or is not readable") % p.c_str());
		exit(1);
	}
	
	CPBA::PROG_ARGS = vm;
	
	Output::print(boost::format("CPBA::dose_scale(): Input Path --> %1%") % vm["inputfile"].as<string>());
	Output::print(boost::format("CPBA::dose_scale(): Scale Factor --> %1%") % vm["factor"].as<double>());
	
	int x = vm["x"].as<int>();
	int y = vm["y"].as<int>();
	int z = vm["z"].as<int>();
	
	double scale_factor = vm["factor"].as<double>();
	
	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Initialization: %1% seconds") % elapsed);

	/* scale */
	dose_volume *DV = Pinnacle::read_and_scale_dosegrid(vm["inputfile"].as<string>(), z, y, x, scale_factor);
	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Scaling: %1% seconds") % elapsed);

	/* write */
	string outfile = boost::str(boost::format("%1%/dv.%2%") % vm["outputdir"].as<string>() % vm["index"].as<string>());
	Pinnacle::write_computed_dosegrid(DV, outfile);
	elapsed = ((double)clock() - start) / CLOCKS_PER_SEC;
	Output::print(boost::format("--->TIMING: Total Execution Time: %1% seconds") % elapsed);
	
}
