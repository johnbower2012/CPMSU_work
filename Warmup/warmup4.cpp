#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>

using namespace std;

double fcompder1(double x, double h);
double fcompder2(double x, const double& h);
double ferror(double exact, double calc);

ofstream ofile;

int main (int argc, char* argv[]){
// ns = number of steps  (need at least 2); isl = initial step length; fsl = final step length
// x = value to be computed; exact = analystic solution
	int ns;
	double isl, fsl, x, exact; 
	double* sind;
	char* ofilename;
	char* efilename;

	x = sqrt(2);
	exact = ((double) 1/(double) 3);

	if (argc<=5){
		cout << "Bad usage. Enter ofilename valueofx numberofsteps initialstep finalstep\n";	
		exit(1);
	}	
	else{
		ofilename = (argv[1]);
		efilename = (argv[2]);
		ns = atoi(argv[3]);
		isl = atof(argv[4]);
		fsl = atof(argv[5]);
		sind = new (nothrow) double [ns];
	}
//Check for proper dynamic memory allocation
	if (sind==nullptr) {
		cout << "Error allocating step memory. Terminating...\n";
		exit(1);
	}
	else {
//array of step values
		for(int i=0; i<ns; i++){
			sind[i]= pow(10.0,(log10(isl) + (log10(fsl)-log10(isl))*((double) i)/((double) (ns-1))));
		}
//Array of derivative values
		double* compder1 = new (nothrow) double [ns];
		double* compder2 = new (nothrow) double [ns];
		if(compder1==nullptr||compder2==nullptr){
			cout << "Error allocating derivative memory. Terminating...\n";
			exit(1);
		}
		else{
//opening outfile
			ofile.open(ofilename);
			ofile << setprecision(10) << scientific;
			ofile << "# h" << setw(20) << "d1" << setw(20) << "d2\n";
			for(int i=0; i<ns; i++){
				compder1[i]=fcompder1(x, sind[i]);
				compder2[i]=fcompder2(x, sind[i]);
				ofile << sind[i] << setw(20);
				ofile << compder1[i] << setw(20);
				ofile << compder2[i] << '\n';
			}
			ofile.close();

			double* error1 = new (nothrow) double [ns];
			double* error2 = new (nothrow) double [ns];
			if(error1==nullptr||error2==nullptr){
				cout << "Error allocating error memory. Terminating...\n";
				exit(1);
			}
			else{
				ofile.open(efilename);
				ofile << setprecision(10) << scientific;
				ofile << "# h" << setw(20) << "d1" << setw(20) << "d2\n";
				for(int i=0; i<ns; i++){
					error1[i]=ferror(exact,compder1[i]);
					error2[i]=ferror(exact,compder2[i]);
					ofile << log10(sind[i]) << setw(20);
					ofile << log10(error1[i]) << setw(20);
					ofile << log10(error2[i]) << '\n';
				}
			}
			delete[] sind;
			delete[] compder1;
			delete[] compder2;
			delete[] error1;
			delete[] error2;
		}		
	}
	return 0;
}
		
double fcompder1(double x, double h){
	double value;
	value = ((atan(x+h)-atan(x))/h);
	return value;
}

double fcompder2(double x, const double& h){
	double value;
	value = ((atan(x+h)-atan(x-h))/(2*h));
	return value;
}		

double ferror(double exact, double calc){
	double value;
	value = (fabs(exact-calc)/exact);
	return value;
}

