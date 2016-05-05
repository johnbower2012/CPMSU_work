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


/*
void apriori_energy_mean(double& energy, double T, double J);
void apriori_magnetization_mean_absv(double& magnetization, double T, double J);
void spin_energy_change(int**& spin, double& energy, double J, int x, int y, int rows, int columns);
void spin_evolve(int**& spin, double T, double J, int n, int m, double& energy_mean, double& heat_capacity, double& magnetization_mean, double& magnetization_mean_absv, double& susceptibility, int& counter);
void spin_evolve_write(int**& spin, double T, double J, int n, char* outfilename, int m, int printcount, double& energy_mean, double& heat_capacity, double& magnetization_mean, double& magnetization_mean_absv, double& susceptibility);


void apriori_energy_mean(double& energy, double T, double J){
	double argument = 8.0*J/T;
	energy = -2.0*J*sinh(argument)/(3.0+cosh(argument));
}
void apriori_magnetization_mean_absv(double& magnetization, double T, double J){
	double argument = 8.0*J/T;
	magnetization = (2.0*exp(argument) + 4.0)/(3.0 + cosh(argument))/4.0;
}

void spin_energy_change(int**& spin, double& energy, double J, int x, int y, int rows, int columns){
	energy = 2*J*spin[x][y]*(spin[(x+1)%rows][y] + spin[(x+rows-1)%rows][y] + spin[x][(y+1)%columns] + spin[x][(y+columns-1)%columns]);
}
void spin_evolve(int**& spin, double T, double J, int n, int m, double& energy_mean, double& heat_capacity, double& magnetization_mean, double& magnetization_mean_absv, double& susceptibility, int& counter){
	double v1, beta, delta_e, old_energy, n_sq;
	double energy, energy_sq_mean;
	double magnetization, magnetization_sq_mean;
	double* exponential;
	energy_mean = 0.0; magnetization_mean = 0.0; magnetization_mean_absv = 0.0; 
	energy_sq_mean = 0.0; magnetization_sq_mean = 0.0;
	heat_capacity = 0.0; susceptibility = 0.0;
	n_sq = n*n;
	beta = 1.0/T;

	int cycles, i, expo, x, y;
	counter = 0; cycles = 0;

	unsigned seed;

	//Allocate array for expon values
	array_alloc(exponential, 2);
	for(i=0;i<2;i++){
		exponential[i] = exp(-J*beta*(i+1)*4);
	}

	//Random number generators
	seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	uniform_real_distribution<double> double_dist(0.0,1.0);
	uniform_int_distribution<int> int_dist(0,n-1);

	//calculate initial energy and magn.
	spin_energy_total(spin, energy, J, n, n);
	spin_magnetization(spin, magnetization, n, n);
	
	for(i=0;i<m;i++){
		old_energy = energy;
		x = int_dist(generator); y = int_dist(generator);
		spin_energy_change(spin, delta_e, J, x, y, n, n);
		v1 = double_dist(generator);
		if(delta_e>0){
			expo = (delta_e/4.0/J-1);
			if(v1>exponential[expo]){
				//Faster to check for > than to check <=;
			}
			else{
				spin[x][y] *= -1;
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
void spin_evolve_write(int**& spin, double T, double J, int n, char* outfilename, int m, int printcount, double& energy_mean, double& heat_capacity, double& magnetization_mean, double& magnetization_mean_absv, double& susceptibility){
	double v1, beta, delta_e, old_energy, n_sq;
	double energy_sq_mean, magnetization_sq_mean;
	double* exponential;
	energy_mean = 0.0; magnetization_mean = 0.0; magnetization_mean_absv = 0.0; 
	energy_sq_mean = 0.0; magnetization_sq_mean = 0.0;
	heat_capacity = 0.0; susceptibility = 0.0;
	n_sq = n*n;
	beta = 1.0/T;

	int cycles, counter, i, expo, x, y;
	counter = 0; cycles = 0;

	unsigned seed;

	//Allocate array for expon values
	array_alloc(exponential, 2);
	for(i=0;i<2;i++){
		exponential[i] = exp(-J*beta*(i+1)*4);
	}

	//Random number generators
	seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	uniform_real_distribution<double> double_dist(0.0,1.0);
	uniform_int_distribution<int> int_dist(0,n-1);
	
	spinfile.open(outfilename);
	spinfile << "#c" << setw(15) << "AC" << setw(15) << "E" << setw(15) <<  "M" << setw(15) << "abs(M)";
	spinfile << setw(15) << "Cv" << setw(15) << "Khi" << endl;
	for(i=0;i<m;i++){
		old_energy = energy;
		x = int_dist(generator); y = int_dist(generator);
		spin_energy_change(spin, delta_e, J, x, y, n, n);
		v1 = double_dist(generator);
		if(delta_e>0){
			expo = (delta_e/4.0/J-1);
			if(v1>exponential[expo]){
				//Faster to check for > than to check <=;
			}
			else{
				spin[x][y] *= -1;
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

		if(i%printcount == 0){
			spinfile << i << setw(15) << counter << setw(15) << energy_mean/n_sq/((double) cycles);
			spinfile << setw(15) << magnetization_mean/n_sq/((double) cycles);
			spinfile << setw(15) << magnetization_mean_absv/n_sq/((double) cycles);
			heat_capacity = beta*beta*(energy_sq_mean/((double) cycles) - energy_mean*energy_mean/((double) cycles)/((double) cycles))/n_sq;
			susceptibility = beta*(magnetization_sq_mean/((double) cycles) - magnetization_mean_absv*magnetization_mean_absv/((double) cycles)/((double) cycles))/n_sq;
			spinfile << setw(15) << heat_capacity << setw(15) << susceptibility << endl;
		}
	}
	spinfile.close();

	energy_mean /= ((double) cycles);
	magnetization_mean /= ((double) cycles);
	magnetization_mean_absv /= ((double) cycles);
	energy_sq_mean /= ((double) cycles);
	magnetization_sq_mean /= ((double) cycles);

	heat_capacity = beta*beta*(energy_sq_mean - energy_mean*energy_mean);
	susceptibility = beta*(magnetization_sq_mean - magnetization_mean_absv*magnetization_mean_absv);

	energy_mean /= n_sq;
	magnetization_mean /= n_sq;
	magnetization_mean_absv /= n_sq;
	heat_capacity /= n_sq;
	susceptibility /= n_sq;

	array_delete(exponential); 
}
*/

#endif
