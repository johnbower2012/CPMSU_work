#ifndef PROJECT4_LIBRARY_H
#define PROJECT4_LIBRARY_H

#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<random>
#include<chrono>
#include "time.h"

using namespace std;

ofstream spinfile;

void array_alloc(double*& a, int length);
void array_delete(double*& a);
void matrix_alloc(int**& A, int rows, int columns);
void matrix_delete(int**& A, int rows);
void matrix_print(int**& A, int rows, int columns);

void array_alloc(double*& a, int length){
	a = new double[length];
}
void array_delete(double*& a){
	delete[] a;
}
void matrix_alloc(int**& A, int rows, int columns){
	int i, j;
	A = new int*[rows];
	for(i=0;i<rows;i++){
		A[i] = new int[columns];
	}
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			A[i][j] = 0.0;
		}
	}
}
void matrix_delete(int**& A, int rows){
	int i;
	for(i=0;i<rows;i++){
		delete[] A[i];
	}
	delete[] A;
}
void matrix_print(int**& A, int rows, int columns){
	int i, j;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			cout << setw(10) << A[i][j];
		}
		cout << endl;
	}
}

				energy += delta_e;
				magnetization += spin[x][y]*2;
				counter++;
			}
		}
		else{
			spin[x][y] *= -1;
			energy += delta_e;
			magnetization += spin[x][y]*2;
			counter++;
		}

		energy_mean += energy;
		energy_sq_mean += energy*energy;
		magnetization_mean += magnetization;
		magnetization_sq_mean += magnetization*magnetization;
		magnetization_mean_absv += abs(magnetization);

		cycles++;
	}

	energy_mean /= ((double) cycles);
	magnetization_mean /= ((double) cycles);
	magnetization_mean_absv /= ((double) cycles);
	energy_sq_mean /= ((double) cycles);
	magnetization_sq_mean /= ((double) cycles);

	heat_capacity = beta*beta*(energy_sq_mean - energy_mean*energy_mean);
	susceptibility = (magnetization_sq_mean - magnetization_mean_absv*magnetization_mean_absv);

	energy_mean /= n_sq;
	magnetization_mean /= n_sq;
	magnetization_mean_absv /= n_sq;
	heat_capacity /= n_sq;
	susceptibility /= n_sq;

	array_delete(exponential); 
}

#endif
