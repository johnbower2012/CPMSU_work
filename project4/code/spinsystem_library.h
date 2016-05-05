#ifndef SPINSYSTEM_LIBRARY_H
#define SPINSYSTEM_LIBRARY_H

#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<random>
#include<chrono>
#include "time.h"

using namespace std;

class spin_system{
	public:
		int** spin;
		int rows, columns;
		double J, T;
		double energy, magnetization;

		spin_system();
		spin_system(int,int,double,double);
		spin_system(int,int,double,double,int);

		void print();
		void energy_total();
		void magnetization_total();
		void randomize();
		void align_up();
		void align_down();
};
class spin_evolution{
	public:
		spin_system spintoddler;
		int* probability;
		double exponential[2];
		double energy_mean, energy_sq_mean;
		double heat_capacity, e_variance;
		double magn_mean, magn_mean_abs, magn_sq_mean;
		double susceptibility, m_variance;
		
		spin_evolution(spin_system&);

		void energy_change_check(double&,int&,int,int);
		void evolve_test(bool&,double,int);
		void evolve(double,int,int);
		void update_statistics();
		void reset_statistics();
		void final_form(int);
		void evolution(int&,int);
		void print();
};

spin_evolution::spin_evolution(spin_system& spinbaby){
	int i, j;	
	spintoddler.rows = spinbaby.rows;
	spintoddler.columns = spinbaby.columns;
	spintoddler.J = spinbaby.J;
	spintoddler.T = spinbaby.T;
	spintoddler.energy = spinbaby.energy;
	spintoddler.magnetization = spinbaby.magnetization;
	
	spintoddler.spin = new int*[spintoddler.rows];
	for(i=0;i<spintoddler.rows;i++){
		spintoddler.spin[i] = new int[spintoddler.columns];
	}

	for(i=0;i<spintoddler.rows;i++){
		for(j=0;j<spintoddler.columns;j++){
			spintoddler.spin[i][j] = spinbaby.spin[i][j];
		}
	}
	probability = new int[4*spintoddler.rows*spintoddler.columns+1];
	for(i=0;i<2*spintoddler.rows*spintoddler.columns+1;i++){
		probability[i] = 0;
	}

	exponential[0] = exp(-4.0*spintoddler.J/spintoddler.T);
	exponential[1] = exp(-8.0*spintoddler.J/spintoddler.T);

	energy_mean = 0.0; energy_sq_mean = 0.0;
	heat_capacity = 0.0; e_variance = 0.0;
	magn_mean = 0.0; magn_mean_abs = 0.0; magn_sq_mean = 0.0;
	susceptibility = 0.0; m_variance = 0.0;
}
void spin_evolution::energy_change_check(double& delta_e, int& index, int x, int y){
	int n, m;
	double sum;
	n = spintoddler.rows;
	m = spintoddler.columns;
	sum = (spintoddler.spin[(x+1)%n][y] + spintoddler.spin[(x+n-1)%n][y] + spintoddler.spin[x][(y+1)%m] + spintoddler.spin[x][(y+m-1)%m]);
	delta_e = 0;
	if(sum==0){}
	else{
		delta_e = 2.0*((double) spintoddler.spin[x][y])*sum;
		index = abs(sum)/2 - 1;
	}
}	
void spin_evolution::evolve_test(bool& test, double delta_e, int index){
	double random;
	unsigned seed;

	test = false;

	seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	uniform_real_distribution<double> double_dist(0.0,1.0);

	if(delta_e>0){
		random = double_dist(generator);
		if(random > exponential[index]){
			//faster to check > than both < and =
		}
		else{
			test = true;
		}
	}
	else{
		test = true;
	}
}		
void spin_evolution::evolve(double delta_e, int x, int y){
	spintoddler.spin[x][y] *= -1;
	spintoddler.magnetization += 2*spintoddler.spin[x][y];
	spintoddler.energy += delta_e;
}
void spin_evolution::update_statistics(){
	int index;
	energy_mean += spintoddler.energy;
	energy_sq_mean += spintoddler.energy*spintoddler.energy;
	magn_mean += spintoddler.magnetization;
	magn_mean_abs += abs(spintoddler.magnetization);
	magn_sq_mean += spintoddler.magnetization*spintoddler.magnetization;
	index = spintoddler.energy/2 + spintoddler.rows*spintoddler.columns;
	probability[index] += 1;
}
void spin_evolution::reset_statistics(){
	int i;	
	energy_mean = 0.0;
	magn_mean = 0.0;
	magn_mean_abs = 0.0;
	energy_sq_mean = 0.0;
	magn_sq_mean = 0.0;

	e_variance = 0.0;
	heat_capacity = 0.0;
	m_variance = 0.0;
	susceptibility = 0.0;

	for(i=0;i<2*spintoddler.rows*spintoddler.columns+1;i++){
		probability[i] = 0;
	}

	exponential[0] = exp(-4.0*spintoddler.J/spintoddler.T);
	exponential[1] = exp(-8.0*spintoddler.J/spintoddler.T);
}
void spin_evolution::final_form(int cycles){
	double T, spins, cycles_d;
	T = spintoddler.T;
	spins = spintoddler.rows*spintoddler.columns;
	cycles_d = cycles;

	energy_mean /= cycles_d;
	magn_mean /= cycles_d;
	magn_mean_abs /= cycles_d;
	energy_sq_mean /= cycles_d;
	magn_sq_mean /= cycles_d;

	e_variance = (energy_sq_mean - energy_mean*energy_mean)/spins;
	heat_capacity = e_variance/T/T;
	m_variance = (magn_sq_mean - magn_mean*magn_mean)/spins;
	susceptibility = (magn_sq_mean - magn_mean_abs*magn_mean_abs)/T/spins;

	energy_mean /= spins;
	magn_mean /= spins;
	magn_mean_abs /= spins;
	energy_sq_mean /= spins;
	magn_sq_mean /= spins;
}
void spin_evolution::evolution(int& accstates, int cycles){
	bool test;

	int i, x, y, index;
	accstates = 0;

	double delta_e;

	unsigned seed;

	//Random number generators
	seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	uniform_int_distribution<int> int_x(0,spintoddler.rows-1);
	uniform_int_distribution<int> int_y(0,spintoddler.columns-1);

	for(i=0;i<cycles;i++){
		x = int_x(generator); y = int_y(generator);
		energy_change_check(delta_e, index, x, y);
		evolve_test(test, delta_e, index);
		if(test==true){
			evolve(delta_e,x,y);
			accstates++;
		}
		update_statistics();
	}
	final_form(cycles);
}
void spin_evolution::print(){
	int i;
	spintoddler.print();
	cout << "<E>/N" << endl;
	cout << setw(15) << energy_mean << " +/- " << pow(e_variance,0.5) << endl;
	cout << "<M>/N" << endl;
	cout <<	setw(15) << magn_mean << " +/- " << pow(m_variance,0.5) << endl;
	cout << endl;
	cout << "<|M|>/N" << endl;
	cout << setw(15) << magn_mean_abs << endl;
	cout << "Cv/N" << endl;
	cout << setw(15) << heat_capacity << endl;
	cout << "Khi/N" << endl;
	cout << setw(15) << susceptibility << endl;
	cout << endl;
}

/********************
	SPIN_SYSTEM
********************/

spin_system::spin_system(){
	rows = 0; columns = 0;
	energy = 0;
	magnetization = 0;
	J = 0;
	T = 0;
}
spin_system::spin_system(int x, int y, double J0, double T0){
	int i, j;

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	uniform_int_distribution<int> distribution(0,1);

	J = J0; T = T0;
	rows = x; columns = y;

	spin = new int*[rows];
	for(i=0;i<rows;i++){
		spin[i] = new int[columns];
	}

	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			spin[i][j] = distribution(generator);
			if(spin[i][j]==0){
				spin[i][j] = -1;
			}
		}
	}
	magnetization_total();
	energy_total();
}
spin_system::spin_system(int x, int y, double J0, double T0, int align){
	int i, j;

	J = J0; T = T0;
	rows = x; columns = y;

	spin = new int*[rows];
	for(i=0;i<rows;i++){
		spin[i] = new int[columns];
	}

	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			spin[i][j] = align;
		}
	}
	magnetization = rows*columns*align;
	energy = -rows*columns*2.0*J;
}
void spin_system::print(){
	int i, j;
	cout << endl;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			cout << setw(5) << spin[i][j];
		}
		cout << endl;
	}
	cout << "energy" << endl;
	cout << setw(15) << energy << endl;
	cout << "magnetization" << endl;
	cout << setw(15) << magnetization << endl;
	cout << endl;
}
void spin_system::energy_total(){
	int i, j;
	energy = 0;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			energy += -spin[i][j]*(spin[(i+1)%rows][j] + spin[(i+rows-1)%rows][j] + spin[i][(j+1)%columns] + spin[i][(j+columns-1)%columns]);
		}
	}
	energy /= 2.0;
}
void spin_system::magnetization_total(){
	int i, j;
	magnetization = 0.0;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			magnetization += spin[i][j];
		}
	}
}
void spin_system::randomize(){
	int i, j;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator (seed);
	uniform_int_distribution<int> distribution(0,1);
	magnetization = 0;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			spin[i][j] = distribution(generator);
			if(spin[i][j]==0){
				spin[i][j] = -1;
			}
			magnetization += spin[i][j];
		}
	}
	energy_total();
}
void spin_system::align_up(){
	int i, j;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			spin[i][j] = 1;
		}
	}
	magnetization = rows*columns;
	energy = -rows*columns*2;
}
void spin_system::align_down(){
	int i, j;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			spin[i][j] = -1;
		}
	}
	magnetization = -rows*columns;
	energy = -rows*columns*2;
}


#endif
