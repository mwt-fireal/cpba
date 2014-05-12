#pragma once

using namespace std;

#include<stdio.h>

#include<vector>
#include<list>
#include<string>
#include<sstream>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>

#include<boost/format.hpp>

class cutil
{
public:

	//templated.
	//all static template functions must be declared inline to allow for correct preprocessor expansion 
	//to satisfy any missing symbol errors.

	template<typename ST, typename T>
	static ST MakeSequence(T begin, T step, T end)
	{
		ST vals;
		T i;
	
		i = begin;
		while(i <= end)
		{
			vals.push_back(i);
			i += step;
		}

		return vals;
	}

	template<typename T>
	static T* MakeArraySequence(T begin, T step, T end)
	{
		//first we need to figure out how many elements.

		int count = 0;
		T i;
		int index;

		for(i = begin; i <= end; i+=step)
		{
			count++;
		}

		T *ary = (T *)malloc(sizeof(T) * ++count);

		index = 0;
		T current = begin;

		while(current < end)
		{
			ary[index] = current;
			current += step;
			index++;
		}

		return ary;
	}


	template<typename T>
	static gsl_vector *MakeVectorSequence(T begin, T step, T end)
	{
		//first we need to figure out how many elements.
		T i;
		int count = 0;
		T current = begin;
		int index = 0;

		for(i = begin; i <= end; i+=step)
		{
			count++;
		}

		//char obuf[1000];
		//sprintf(obuf, "from %lf to %lf, with step of %lf, found count of %d", begin, end, step, count);
		//Output::print(obuf);

		gsl_vector *v = gsl_vector_alloc(count);

		while(current < end)
		{
			gsl_vector_set(v, index, current);
			current += step;
			index++;
		}

		return v;
	}

	template<typename T>
	static int numel(vector<T> ary)
	{
		return ary.size();
	}


	template<typename T>
	static int numel(vector< vector<T> > ary)
	{
		int sz = 0;

		int dim = ary.size();

		for(int i = 0; i < dim; i++)
		{
			sz += ((vector<T>)ary[i]).size();
		}

		return sz;
	}

	template<typename T>
	static int numel(vector<T> ary, int begin, int end)
	{
		if(end == -1)
		{
			return ary.size() - begin;
		}
		else
		{
			if(end > ary.size())
			{
				return (ary.size() - begin);
			}

			return (end - begin) + 1;
		}
	}

	template<typename T>
	static string tostring(T val)
	{
		stringstream ss;
		ss << val;

		return ss.str();
	}

	template<typename T>
	static int numel(vector< vector<T> > ary, int begin, int end)
	{
		int sz = 0;

		int dim = ary.size();

		if(end == -1)
		{
			end = dim;
		}

		for(int i = begin; i < end; i++)
		{
			sz += ((vector<T>)ary[i]).size();
		}

		return sz;
	}

	static gsl_vector *stl_vector_to_gsl_vector(vector<double>* vec);
	static gsl_vector *stl_list_to_gsl_vector(list<double> *vec);
	
	//non templated.
	static double *vec_to_array(gsl_vector *vec);
	static gsl_vector *array_to_vec(double *ary, int size);
	static void vec_operation_log(gsl_vector *v);
	static void vec_operation_abs(gsl_vector *v);
	static gsl_vector *vec_operation_abs_copy(gsl_vector *v);
	static gsl_vector *vec_operation_copy(gsl_vector *v);
	static gsl_vector *vec_operation_chop(gsl_vector *v, int start, int end);
	static gsl_vector *vec_operation_negate(gsl_vector *v);
	static void vec_operation_inline_negate(gsl_vector *v);
	static double *vec_operation_log_to_array(gsl_vector *v);
	static gsl_vector *vec_operation_not_zero(gsl_vector *v);
	static gsl_vector *vec_operation_add_constant(gsl_vector *v, double scale_factor);
	static int vec_operation_greater_than(gsl_vector *v, double value);
	static vector<int> vec_operation_greater_than_index(gsl_vector *v, double value);
	static vector<int> vec_operation_less_than_index(gsl_vector *v, double value);
	static int vec_operation_less_than(gsl_vector *v, double value);
	static gsl_vector *vec_operation_equal(gsl_vector *v, double value);
	static gsl_vector *vec_operation_less_than_equal(gsl_vector *v, double value, bool abs);
	static gsl_vector *vec_operation_less_than_not_zero(gsl_vector *v, double value);
	static gsl_vector *vec_operation_greater_than_not_zero(gsl_vector *v, double value);
	static vector<gsl_vector *> *vec_operation_unique(gsl_vector *v);
	static int fast_count_lines(FILE *fptr);
	static char *fast_get_line(FILE *fptr);
	static int scan_columns(FILE *fptr);
	
	static void make_finite_non_neg(gsl_matrix *m);
	static void make_finite_non_neg(gsl_vector *v);

	static gsl_matrix *read_matrix_from_file(string filename);
	static gsl_matrix *read_matrix_from_file(boost::format filename);

	static void write_matrix_to_file(boost::format filename, gsl_matrix *m);
	static void write_matrix_to_file(string filename, gsl_matrix *m);

	static void write_vector_to_file(boost::format filename, gsl_vector *v);
	static void write_vector_to_file(string filename, gsl_vector *v);

	static void write_data_to_file(boost::format filename, gsl_vector *v);
	static void write_data_to_file(string filename, gsl_vector *v);
	static void write_data_to_file(boost::format filename, gsl_matrix *m);
	static void write_data_to_file(string filename, gsl_matrix *m);

	static gsl_vector *get_matrix_row(gsl_matrix *M, int row);
	static gsl_vector *get_matrix_col(gsl_matrix *M, int col);

	static bool is_finite(double v);
};
