#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<random>
#include<chrono>
#include "time.h"
#include "project4_library.h"
#include "spinsystem_library.h"

using namespace std;

ofstream ofile;

int main(int argc, char* argv[]){
	int i, n, accstates, cycles, align, printcount, cycpow, choice;
	double T, J, time, current_T;
	char* outfilename, *probfilename;
	clock_t start, finish;

	/**************************************************
	Enter required information via command line
		n = grid size
		T = temperature
			choice = 0; choice = 1
				T = exact temperature
				T = number of temp steps between
					Tmin and Tmax
		J = energy
		align
			= 0		random orientation
			= 1		fully aligned (up)
			= 2		keep alignment from previous
						computation
		outfilename = output file name
		cycles = power of 2 for number of MC cycles
					So 	5 ---> 2^5
						6 ---> 2^6
						etc
		choice (prob_write_temp) 
			= 0 	probability distribution file
			= 1		loop through cycles count
						2^0 ---> 2^entry
			= 2		loop through temperature
						Enter Tmin, Tmax
						T = steps from Tmin to Tmax
						Tmin + i*(Tmax-Tmin)/T
						Tmin ... Tmax in T+1 steps
	**************************************************/
	if(argc<8){
		cout << "Bad usage. Enter also 'n T J aligned_init_0_1_2 outfilename cycles prob_write_temp_0_1_2' on same line." << endl;
		exit(1);
	}
	else{
		n = atoi(argv[1]);
		T = atof(argv[2]);
		J = atof(argv[3]);
		align = atoi(argv[4]);
		outfilename = argv[5];
		cycles = atoi(argv[6]);
		choice = atoi(argv[7]);
	}


	/**************************************************
	Create and initialize nxn spin grid; calculate
		energy, magnetization, etc.
	**************************************************/
	spin_system spinbaby(n,n,J,T);
	if(align==1){
		spinbaby.align_up();
	}

	/**************************************************
	Create spin_evolution from the spin_system
	**************************************************/
	spin_evolution spinadult(spinbaby);


	/**************************************************
	Evolve according to choice value
		0 ---> probability
		1 ---> loop through cycle counts
		2 ---> loop through temperature values
				in T+1 steps, including Tmin, Tmax
	**************************************************/
	if(choice==0){
		start = clock();

		cycpow = pow(2,cycles);

		/**************************************************
			evolve through 2^cycles
			track the number of accepted states
			Note: we do this so the equilibrium is reached
					before computation
		**************************************************/
		spinadult.evolution(accstates, cycpow);

		cout << endl;
		cout << cycles << setw(15) << accstates;

		spinadult.print();

		finish = clock();
		time = (finish - start)/((double) CLOCKS_PER_SEC);

		cout << "time" << endl;
		cout << setw(15) << time << endl;

		/**************************************************
			reset statistics in spin_evolution now
				that equilibrium has been reached
		**************************************************/
		spinadult.reset_statistics();

		cycpow = pow(2,cycles);
		spinadult.evolution(accstates, cycpow);

		ofile.open(outfilename);
		ofile << setw(15) << "#E" << setw(15) << "E/N" << setw(15) << "Prob" <<  setw(15) << "count" << endl;
		for(i=0;i<2*n*n + 1; i++){
			ofile << setw(15) << (i-n*n)*2 << setw(15) << (double) ((i-n*n)*2)/(double) (n*n);
			ofile << setw(15) << (double) spinadult.probability[i]/(double) cycpow << setw(15) << spinadult.probability[i] << endl;
		}
		ofile.close();
	}
	else if(choice==1){
		ofile.open(outfilename);
		ofile << setw(15) << "#i" << setw(15) << "cyc" << setw(15) << "accs" << setw(15) << "time";
		ofile << setw(15) << "<E>" << setw(15) << "e_v" << setw(15) << "<M>" << setw(15) << "m_v";
		ofile << setw(15) << "<|M|>" << setw(15) << "m_abs_v" << setw(15) << "Cv" << setw(15) << "Khi" << endl;
		for(i=0;i<cycles+1;i++){
			/**************************************************
				Align according to choice for each loop
				reset stastics before calculations begin
			**************************************************/
			if(align==0){
				spinadult.spintoddler.randomize();
			}
			else if(align==1){
				spinadult.spintoddler.align_up();
			}
			spinadult.reset_statistics();

			start = clock();

			/**************************************************
				Evolve over 2^cycles MC
			**************************************************/
			cycpow = pow(2,i);
			spinadult.evolution(accstates, cycpow);

			finish = clock();
			time = (finish - start)/((double) CLOCKS_PER_SEC);

			ofile << setw(15) << i << setw(15) << cycpow << setw(15) << accstates << setw(15) << time;
			ofile << setw(15) << spinadult.energy_mean << setw(15) << spinadult.e_variance;
			ofile << setw(15) << spinadult.magn_mean << setw(15) << spinadult.m_variance;
			ofile << setw(15) << spinadult.magn_mean_abs << setw(15) << spinadult.m_abs_variance;
			ofile << setw(15) << spinadult.heat_capacity << setw(15) << spinadult.susceptibility << endl;
			cout << setw(15) << i << setw(15) << cycpow << setw(15) << time << " seconds" << endl;
		}
		ofile.close();
	}
	else if(choice==2){
		ofile.open(outfilename);
		ofile << "#T" << setw(15) << "time" << setw(15) << "accs";
		ofile << setw(15) << "<E>" << setw(15) << "e_v" << setw(15) << "<M>" << setw(15) << "m_v";
		ofile << setw(15) << "<|M|>" << setw(15) << "m_abs_v" << setw(15) << "Cv" << setw(15) << "Khi" << endl;

		/**************************************************
			input Tmin and Tmax
		**************************************************/
		double Tm, TM;
		cout << "Enter Tmin: ";
		cin >> Tm;
		cout << "Enter Tmax: ";
		cin >> TM;

		/**************************************************
			align as chosen
		**************************************************/
		cycpow = pow(2,cycles);
		if(align==0){
			spinadult.spintoddler.randomize();
		}
		else if(align==1){
			spinadult.spintoddler.align_up();
		}
		else if(align==2){
			spinadult.evolution(accstates, cycpow);
		}

		for(i=0;i<T+1;i++){
			/**************************************************
				Calculate T for this iteration
				reset statistics
				evolve
			**************************************************/
			start = clock();
			current_T = Tm + ((TM-Tm)/T)*(double)i;
			spinadult.spintoddler.T = current_T;
			spinadult.reset_statistics();
			spinadult.evolution(accstates, cycpow);
			finish = clock();
			time = (finish - start)/((double) CLOCKS_PER_SEC);
			ofile << current_T << setw(15) << time << setw(15) << accstates;
			ofile << setw(15) << spinadult.energy_mean << setw(15) << spinadult.e_variance;
			ofile << setw(15) << spinadult.magn_mean << setw(15) << spinadult.m_variance;
			ofile << setw(15) << spinadult.magn_mean_abs << setw(15) << spinadult.m_abs_variance;
			ofile << setw(15) << spinadult.heat_capacity << setw(15) << spinadult.susceptibility << endl;
			cout << setw(15) << current_T << setw(15) << accstates << setw(15) << time << " seconds" << endl;
		}
		ofile.close();
	}

	return 0;
}
