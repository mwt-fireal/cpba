#ifndef __PENCILDEBUGH__
#define __PENCILDEBUGH__

using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <boost/format.hpp>

#include <string>

class Output
{
public:
	static void print(string message);
	static void print(const char *message);
	static void print(boost::format message);

	static void print_only(string message);
	static void print_only(boost::format message);

	static void print_vector(const char *prefix, gsl_vector *v);
	static void print_vector(string prefix, gsl_vector *v);
	static void print_vector(boost::format prefix, gsl_vector *v);

	static void print_matrix(const char *prefix, gsl_matrix *m);
	static void print_matrix(string prefix, gsl_matrix *m);
	static void print_matrix(boost::format prefix, gsl_matrix *m);
};

#endif
