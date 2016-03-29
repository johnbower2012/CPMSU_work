#ifndef VECTORCLASS_H
#define VECTORCLASS_H

#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>

using namespace std;

class vector{
	public:
		int length;
		double* array;
		vector();
		vector(int);

		void alloc(int);
		void del();

		vector operator+(const vector&);
		vector operator-(const vector&);
		double operator*(const vector&);
		vector operator/(const vector&);
};

#endif

