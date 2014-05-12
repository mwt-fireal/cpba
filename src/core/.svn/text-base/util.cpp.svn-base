#include<stdio.h>

#include<vector>

#include "util.h"

#include "out.h"

#include<boost/format.hpp>
#include<boost/filesystem.hpp>
#include<fstream>
#include<cmath>

using namespace std;

int cutil::fast_count_lines(FILE *fptr)
{
	//rewind file, make sure we're at the beginning.
	rewind(fptr);

	int count = 0;
	char *line = cutil::fast_get_line(fptr);

	while(line != (char *)NULL)
	{
		count++;
		line = cutil::fast_get_line(fptr);
		//Output::print(line);
	}
	
	//Output::print(boost::format("fast count lines returns: %1%") % count);

	return count;
}

char *cutil::fast_get_line(FILE *fptr)
{
	
	//char db[1000];
	char a, b, c;

	int current_size = 50;
	//we'll start with enough memory for 50 characters.
	char *da = (char *)malloc(sizeof(char) * current_size);

	int size  = 0;

	if(feof(fptr))
	{
		return (char *)NULL;
	}

	while(!feof(fptr))
	{
		b = fgetc(fptr);

		if(b < 0)
		{
			if(size > 0)
			{
				return da;
			}
			else
			{
				return (char *)NULL;
			}
		}

		if(size+1 == current_size)
		{
			current_size += 50;
			da = (char *)realloc((void *)da, sizeof(char) * current_size);
		}

		// if its \n, and the previous character was not \r (windows is \r\n), count it.
		// if its \r, count it regardless of whether a \n follows it, or its just a \r by itself.
		if((int) b == 13)
		{
			//lookahead
			c = fgetc(fptr);
			
			if((int)c != 10)
			{
				fputc(c, fptr);
			}
			
			da[size++] = '\0';
			return da;
		}
		else if((int)b == 10)
		{
			da[size++] = '\0';
			return da;
		}
		else
		{
			da[size++] = b;
		}

		// store our previous character 
		a = b;
	}

	return da;
}

int cutil::scan_columns(FILE *fptr)
{
	rewind(fptr);

	char *line = cutil::fast_get_line(fptr);

	//scan for spaces.

	int pos = 0;

	char c = line[pos];

	int columns = 0;

	bool in_seperator = false;

	//\0 indicates end of line)
	while(c != '\0')
	{
		//Output::print(boost::format("character is: %1%, as ord: %2%") % c % (int)c);
		if((int)c == 32 || (int)c == 9)
		{
			//Output::print("determined that it is a whitespace character, setting in separator");
			in_seperator = true;
		}
		else
		{
			//begin our first column. this happens if there is leading whitespace before the first real character.
			if(columns == 0)
			{
				in_seperator = false;
				columns++;
			}

			if(in_seperator)
			{
				in_seperator = !in_seperator;
				//Output::print("determined that is is NOT a whitespace character, incrementing columns, resetting in_seperator");
				columns++;
			}
		}

		c = line[++pos];
	}

	free((void *)line);

	return columns;
}

gsl_matrix *cutil::read_matrix_from_file(boost::format filename)
{
	return cutil::read_matrix_from_file(boost::str(filename));
}

gsl_matrix *cutil::read_matrix_from_file(string filename)
{
	boost::filesystem::path p = filename;
	
	if(!boost::filesystem::exists(p))
	{
		Output::print(boost::format("cutil::read_matrix_from_file::error, matrix file %1% does not exist or is not readable") % filename);
		exit(1);
	}
	
	FILE *fptr = fopen(filename.c_str(), "r");

	rewind(fptr);

	int rows = cutil::fast_count_lines(fptr);
	rewind(fptr);
	int columns = cutil::scan_columns(fptr);
	rewind(fptr);

	gsl_matrix *M = gsl_matrix_alloc(rows, columns);

	//gsl_matrix_fread(fptr, M);

	double a;

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++)
		{
			fscanf(fptr, "%lf", &a);
			gsl_matrix_set(M, i, j, a);
		}
	}

	fclose(fptr);

	return M;
}

void cutil::write_matrix_to_file(boost::format filename, gsl_matrix *m)
{
	cutil::write_matrix_to_file(boost::str(filename), m);
}

void cutil::write_matrix_to_file(string filename, gsl_matrix *m)
{
	ofstream fptr(filename.c_str());
	boost::format ff("%1%");

	if(!fptr.is_open())
	{
		Output::print(boost::format("error: could not open %1% for writing") % filename);
		exit(1);
	}

	for(unsigned int i = 0; i < m->size1; i++)
	{
		for(unsigned int j = 0; j < m->size2; j++)
		{
			ff % gsl_matrix_get(m, i, j);

			fptr << ff;

			if(j+1 < m->size2)
			{
				fptr << " ";
			}
		}

		fptr << endl;
	}

	fptr.close();
}

void cutil::write_vector_to_file(boost::format filename, gsl_vector *v)
{
	cutil::write_vector_to_file(boost::str(filename), v);
}

void cutil::write_vector_to_file(string filename, gsl_vector *v)
{
	ofstream fptr(filename.c_str());
	boost::format ff("%1%");

	if(!fptr.is_open())
	{
		Output::print(boost::format("error: could not open %1% for writing") % filename);
		exit(1);
	}

	for(unsigned int i = 0; i < v->size; i++)
	{
		ff % gsl_vector_get(v, i);

		fptr << ff << endl;
	}

	fptr.close();
}

void cutil::write_data_to_file(boost::format filename, gsl_vector *v)
{
	cutil::write_vector_to_file(filename, v);
}
void cutil::write_data_to_file(string filename, gsl_vector *v)
{
	cutil::write_vector_to_file(filename, v);
}

void cutil::write_data_to_file(boost::format filename, gsl_matrix *m)
{
	cutil::write_matrix_to_file(filename, m);
}

void cutil::write_data_to_file(string filename, gsl_matrix *m)
{
	cutil::write_matrix_to_file(filename, m);
}

double *cutil::vec_to_array(gsl_vector *vec)
{
	string msg;

	double *ary = (double *)malloc(sizeof(double) * vec->size);

	for(unsigned int i = 0; i < vec->size; i++)
	{
		ary[i] = gsl_vector_get(vec, i);
	}

	return ary;
}

gsl_vector *cutil::array_to_vec(double *ary, int size)
{
	gsl_vector *vv = gsl_vector_calloc(size);

	for(int i = 0; i < size; i++)
	{
		gsl_vector_set(vv, i, ary[i]);
	}

	return vv;
}

gsl_vector *cutil::vec_operation_add_constant(gsl_vector *v, double scale_factor)
{
	gsl_vector *v_scale = gsl_vector_calloc(v->size);
	gsl_vector_memcpy(v_scale, v);

	gsl_vector_add_constant(v_scale, scale_factor);

	return v_scale;
}

void cutil::vec_operation_log(gsl_vector *v)
{
	for(unsigned int i = 0; i < v->size; i++)
	{
		gsl_vector_set(v, i, log(gsl_vector_get(v, i)));
	}
}

gsl_vector *cutil::vec_operation_copy(gsl_vector *v)
{
	gsl_vector *vv = gsl_vector_calloc(v->size);

	gsl_vector_memcpy(vv, v);

	return vv;
}

gsl_vector *cutil::vec_operation_chop(gsl_vector *v, int start, int end)
{
	gsl_vector *vv = gsl_vector_calloc(end-start);

	for(int i = start; i <= end; i++)
	{
		gsl_vector_set(vv, i, gsl_vector_get(v, i));
	}

	return vv;
}

void cutil::vec_operation_abs(gsl_vector *v)
{
	for(unsigned int i = 0; i < v->size; i++)
	{
		gsl_vector_set(v, i, fabs(gsl_vector_get(v, i)));
	}
}


gsl_vector *cutil::stl_vector_to_gsl_vector(vector<double>* vec)
{
	gsl_vector *v;

	if(vec->size() == 0)
	{
		return (gsl_vector *)NULL;
	}

	v = gsl_vector_calloc(vec->size());

	vector<double>::iterator it;

	int i = 0;
	for(it = vec->begin(); it != vec->end(); it++)
	{
		gsl_vector_set(v, i++, (double)*(it));
	}

	return v;
}


gsl_vector *cutil::stl_list_to_gsl_vector(list<double>* vec)
{
	gsl_vector *v = gsl_vector_calloc(vec->size());

	list<double>::iterator it;

	int i = 0;
	for(it = vec->begin(); it != vec->end(); it++)
	{
		gsl_vector_set(v, i++, (double)*(it));
	}

	return v;
}

gsl_vector *cutil::vec_operation_abs_copy(gsl_vector *v)
{
	gsl_vector *vv = gsl_vector_calloc(v->size);
	for(unsigned int i = 0; i < v->size; i++)
	{
		gsl_vector_set(vv, i, fabs(gsl_vector_get(v, i)));
	}

	return vv;
}

double *cutil::vec_operation_log_to_array(gsl_vector *v)
{
	double *vary = (double *)malloc(sizeof(double) * v->size);

	for(unsigned int i = 0; i < v->size; i++)
	{
		vary[i] = log(gsl_vector_get(v, i));
	}

	return vary;
}

int cutil::vec_operation_greater_than(gsl_vector *v, double value)
{
	//O(2n) ... 
	int count = 0;
	for(unsigned int i = 0; i < v->size; i++)
	{
		if(gsl_vector_get(v, i) > value)
		{
			count++;
		}
	}

	return count;
}

vector<int> cutil::vec_operation_greater_than_index(gsl_vector *v, double value)
{
	vector<int> ii;
	//O(2n) ... 
	for(unsigned int i = 0; i < v->size; i++)
	{
		if(gsl_vector_get(v, i) > value)
		{
			ii.push_back(i);
		}
	}

	return ii;
}

vector<int> cutil::vec_operation_less_than_index(gsl_vector *v, double value)
{
	vector<int> ii;
	//O(2n) ... 
	for(unsigned int i = 0; i < v->size; i++)
	{
		if(gsl_vector_get(v, i) < value)
		{
			ii.push_back(i);
		}
	}

	return ii;
}

int cutil::vec_operation_less_than(gsl_vector *v, double value)
{
	//O(2n) ... 
	int count = 0;
	for(unsigned int i = 0; i < v->size; i++)
	{
		if(gsl_vector_get(v, i) < value)
		{
			count++;
		}
	}

	return count;
}

gsl_vector *cutil::vec_operation_less_than_equal(gsl_vector *v, double value, bool use_abs)
{
	vector<double> *nv = new vector<double>;

	for(unsigned int i = 0; i < v->size; i++)
	{
		double va = use_abs == true ? abs(gsl_vector_get(v, i)) : gsl_vector_get(v, i);

		if(va <= value)
		{
			nv->push_back(va);
		}
	}


	return cutil::stl_vector_to_gsl_vector(nv);

}

gsl_vector *cutil::vec_operation_less_than_not_zero(gsl_vector *v, double value)
{
	vector<double> *nv = new vector<double>;

	//char obuf[100];

	for(unsigned int i = 0; i < v->size; i++)
	{
		double va = gsl_vector_get(v, i);

		if(va != 0 && va < value)
		{
			if(abs(va - value) < 0.0001)
			{
				continue;
			}

			nv->push_back(va);
			//Output::print(boost::format("cutil::vec_operation_less_than_not_zero -- pushing back %1%, its less than %2%") % va % value);
		}
	}

	return cutil::stl_vector_to_gsl_vector(nv);
}

gsl_vector *cutil::vec_operation_greater_than_not_zero(gsl_vector *v, double value)
{
	vector<double> *nv = new vector<double>;

	//char obuf[100];

	for(unsigned int i = 0; i < v->size; i++)
	{
		double va = gsl_vector_get(v, i);

		//sprintf(obuf, "va = %lf", va);
		//Output::print(obuf);

		if(va != 0 && va > value)
		{
			if(abs(va - value) < 0.0001)
			{
				continue;
			}

			nv->push_back(va);
			//sprintf(obuf, "PUSHED BACK: va = %lf, value of %lf", va, value);
			//Output::print(obuf);
		}
	}

	return cutil::stl_vector_to_gsl_vector(nv);
}

gsl_vector *cutil::vec_operation_negate(gsl_vector *v)
{
	gsl_vector *a = gsl_vector_calloc(v->size);

	for(unsigned int i = 0; i < a->size; i++)
	{
		gsl_vector_set(a, i, 0 - (gsl_vector_get(v, i)));
	}

	return a;
}

void cutil::vec_operation_inline_negate(gsl_vector *v)
{
	for(unsigned int i = 0; i < v->size; i++)
	{
		gsl_vector_set(v, i, 0 - (gsl_vector_get(v, i)));
	}
}

//returns the indexes of those elements that are equal.
gsl_vector *cutil::vec_operation_equal(gsl_vector *v, double value)
{
	//O(2n) ... 
	int count = 0;
	unsigned int i;
	for(i = 0; i < v->size; i++)
	{
		if(gsl_vector_get(v, i) == value)
		{
			count++;
		}
	}

	string rr;

	gsl_vector *vv = gsl_vector_calloc(count);

	count = 0;
	for(i = 0; i < v->size; i++)
	{
		if(gsl_vector_get(v, i) == value)
		{
			gsl_vector_set(vv, count, i);
			count++;
		}
	}

	return vv;
}

vector<gsl_vector *> *cutil::vec_operation_unique(gsl_vector *v)
{
	unsigned int i, j;
	list<double> original;

	for(i = 0; i < v->size; i++)
	{
		original.push_back(gsl_vector_get(v, i));
	}

	//[b1]
	list<double> b1 = original;
	b1.sort();
	b1.unique();

	//[m1] (O(2n))

	/*

	A = 
   1   1   5   6   2   3   3   9   8   6   2   4

	[b1, m1, n1] = unique(A, 'first') 
b1 = 
   1   2   3   4   5   6   8   9 
m1 =  
   1   5   6  12   3   4   9   8 
n1 = 
   1   1   5   6   2   3   3   8   7   6   2   4

   */

	Output::print("b1 is created..");

	list<double> m1;
	list<double>::iterator b1_it;
	list<double>::iterator original_it;

	int index;
	for(b1_it = b1.begin(); b1_it != b1.end(); b1_it++)
	{
		index = 0;
		original_it = original.begin();

		while(true)
		{
			if(*b1_it == *original_it)
			{
				m1.push_back(index);
				break;
			}

			original_it++;
			index++;
		}
	}
	
	Output::print("m1 is created..");

	list<double> n1;

	for(original_it = original.begin(); original_it != original.end(); original_it++)
	{
		j = 0;

		for(b1_it = b1.begin(); b1_it != b1.end(); b1_it++)
		{
			if(*original_it == *b1_it)
			{
				//Output::print(boost::format("pushing back: %1%") % j);
				n1.push_back(j);
			}

			j++;
		}
	}
	
	Output::print("n1 is created..");

	vector<gsl_vector *> *ar = new vector<gsl_vector *>;

	ar->push_back(cutil::stl_list_to_gsl_vector(&b1));
	ar->push_back(cutil::stl_list_to_gsl_vector(&m1));
	ar->push_back(cutil::stl_list_to_gsl_vector(&n1));

	return ar;
}


gsl_vector *cutil::vec_operation_not_zero(gsl_vector *v)
{
	//O(2n) ... 
	int count = 0;
	unsigned int i;
	for(i = 0; i < v->size; i++)
	{
		if(gsl_vector_get(v, i) != 0)
		{
			count++;
		}
	}

	gsl_vector *nv = gsl_vector_calloc(count);

	double elem;

	for(i = 0; i < v->size; i++)
	{
		elem = gsl_vector_get(v, i);

		if(elem != 0)
		{
			gsl_vector_set(nv, count, elem);
			count++;
		}
	}

	return nv;
}

void cutil::make_finite_non_neg(gsl_matrix *m)
{
	for(int i = 0; i < m->size1; i++)
	{
		for(int j = 0; j < m->size2; j++)
		{
			double ans = gsl_matrix_get(m, i, j);
			if(!cutil::is_finite(ans) || ans < 0)
			{
				gsl_matrix_set(m, i, j, 0);
			}
		}
	}
}

void cutil::make_finite_non_neg(gsl_vector *v)
{
	for(int j = 0; j < v->size; j++)
	{
		double ans = gsl_vector_get(v, j);
		if(!cutil::is_finite(ans) || ans < 0)
		{
			gsl_vector_set(v, j, 0);
		}
	}
}

gsl_vector *cutil::get_matrix_row(gsl_matrix *M, int i)
{
	gsl_vector *v = gsl_vector_calloc(M->size2);

	gsl_matrix_get_row(v, M, i);

	return v;
}

gsl_vector *cutil::get_matrix_col(gsl_matrix *M, int i)
{
	gsl_vector *v = gsl_vector_calloc(M->size1);

	gsl_matrix_get_col(v, M, i);

	return v;
}

bool cutil::is_finite(double v)
{
	//this looks like it will always return true, but if v is NaN or INF, this will always 
	//evaluate to false.
	if(v == v)
	{
		return true;
	}

	//Output::print(boost::format("FAILED IS_FINITE on value %1%") % v);
	return false;
}
