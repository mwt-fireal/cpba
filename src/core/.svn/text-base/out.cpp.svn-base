#include <iostream>
#include <string>
#include <sstream>

#include "out.h"

using namespace std;

#ifdef _WINDOWS 
	#include <Windows.h>
#endif



void Output::print_vector(const char *prefix, gsl_vector *vector)
{
	Output::print_vector(string(prefix), vector);
}

void Output::print_vector(boost::format prefix, gsl_vector *vector)
{
	Output::print_vector(boost::str(prefix), vector);
}

void Output::print_vector(string prefix, gsl_vector *vector)
{
	stringstream ss;
	stringstream f_ss;
	for(unsigned int i = 0; i < vector->size; i++)
	{
		if(i > 0)
		{
			ss << " ";
		}

		ss << gsl_vector_get(vector, i);
	}

	f_ss << prefix << ": " << ss.str();

	Output::print(f_ss.str());
}

void Output::print_matrix(const char *prefix, gsl_matrix *m)
{
	Output::print_matrix(string(prefix), m);
}

void Output::print_matrix(boost::format prefix, gsl_matrix *m)
{
	Output::print_matrix(boost::str(prefix), m);
}

void Output::print_matrix(string prefix, gsl_matrix *m)
{
	Output::print(boost::format("%1%[%2%][%3%]") % prefix % m->size1 % m->size2);

	for(unsigned int i = 0; i < m->size1; i++)
	{
		for(unsigned int j = 0; j < m->size2; j++)
		{
			Output::print_only(boost::format("%1%") % gsl_matrix_get(m, i, j));

			if(j+1 < m->size2)
			{
				Output::print_only(" ");
			}
		}

		Output::print_only("\r\n");
	}
}

void Output::print(const char *message)
{
	Output::print(string(message));
}

void Output::print(boost::format message)
{
	Output::print(boost::str(message));
}

void Output::print(string message)
{
	stringstream ss;
	ss << "DEBUG::" << message;

	#ifdef _WINDOWS 
		OutputDebugString(ss.str().c_str());
		OutputDebugString("\r\n");
	#else
		cout << ss.str() << endl;
	#endif
}

void Output::print_only(string message)
{
	#ifdef _WINDOWS 
		OutputDebugString(message.c_str());
	#else
		cout << message;
	#endif
}

void Output::print_only(boost::format message)
{
	Output::print_only(boost::str(message));
}


