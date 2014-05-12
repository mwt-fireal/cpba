#pragma once

/* use PBAMATH because MATH is probably defined by something, somewhere */
#ifndef __PBAMATH__
#define __PBAMATH__

using namespace std;

#include<map>
#include<list>
#include<vector>
#include<string>

class MaterialComponent
{
public:
	string symbol;
	double elNum;
	double Z;
	double A;
	double f;
};


class Material
{
public:
	string symbol;
	string formula;
	double rho;
	int numElements;

	vector<map<string, string> > elements;
	vector<map<string, string> >::iterator elem_it;

	Material(string symbol, string formula, double rho, int numElements)
	{
		this->symbol = symbol;
		this->formula = formula;
		this->rho = rho;
		this->numElements = numElements;
	}
};

class PBAMaterials
{
public:
	map<string, Material*> *mdb;
	PBAMaterials(string source_file);
};

#endif
