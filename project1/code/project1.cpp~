/* 
Here we have a second order differential equation, solved numerically by use of grid points. 
We have -u"(x) = f(x), x from (0,1), u(0)=u(1)=0, f(x) = 100exp(-10x). 
We solve this by transforming it into a linear algebra problem, wherein -(v_i+1 + v_i-1 - 2*v_i)/(h^2) = f_i. 
We first solve it as a general tridiagonal matrix with arrays a, b, c, then for our specific case where a_ii=2,
a_(|i-j|=1)=1, and a_(|i-j|>1)=0.

A = ( b_1 c_1  0   0   0  ...  0   0  )
    ( a_2 b_2 c_2  0   0  ...  0   0  )
    (  0  a_3 b_3 c_3  0  ...  0   0  )
    ( ... ... ... ... ... ... ... ... )
    (  0   0   0   0  ...  0  a_n b_n )
*/


#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include "time.h"

using namespace std;

ofstream ofile;

int main(int argc, char* argv[]){
	char* outfilename;
	int n, i, yesno, j;
	double h, time1, time2;
	double* a, *b, *c, *v1, *v2, *func;
	clock_t start, finish;

	if(argc<4){
		cout << "Bad usage. Include also in command line gridpoint_count outfilename '1 or 0 (for error)'";
		exit (1);
	}
	else{
		n = atoi(argv[1]);
		outfilename = argv[2];
		yesno = atoi(argv[3]);
		ofile.open(outfilename);
		h=(1.0/(n+1));
	}

//Dynamically allocate memory, with check.	
	a = new (nothrow) double[n];
	b = new (nothrow) double[n];
	c = new (nothrow) double[n];
	v1 = new (nothrow) double[n];
	v2 = new (nothrow) double[n];
	func = new (nothrow) double[n];

	if(a==nullptr||b==nullptr||c==nullptr||v1==nullptr||v2==nullptr||func==nullptr){
		cout << "Error allocating dynamic memory. Terminating...\n";
		exit(1);
	}

//Initialize array values.
	for(i=0;i<n;i++){
		a[i]=(-1.0);
		b[i]=(2.0);
		c[i]=(-1.0);
		func[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
	}

/*
  Start clock for general tridiagonal matrix
  Create upper diagonal matrix
  Back substitutions
  End clock
*/
	start = clock();

	for(i=1;i<n;i++){
		double temp;
		temp = a[i]/b[i-1];
		b[i] -= temp*c[i-1];
		func[i] -= temp*func[i-1];
	}
	v1[n-1]=func[n-1]/b[n-1];
	for(i=n-2;i>=0;i--){
		v1[i] = (func[i]-c[i]*v1[i+1])/b[i];
	}

	finish = clock();
	time1 = ( (finish - start)/(double) CLOCKS_PER_SEC);

//redfine b values for our specific tridiagonal, reinitialize func[i] values.
	for(i=0;i<n;i++){
		b[i]=(double) (i+2)/(double) (i+1);
		func[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
	}

/*
  Start clock for our specific tridiagonal
  Form A into upper diagonal
  Back substitutions
  End clock
*/
	start = clock();

	for(i=1;i<n;i++){
		func[i] += func[i-1]/b[i-1];
	}
	v2[n-1]=func[n-1]/b[n-1];
	for(i=n-2;i>=0;i--){
		v2[i] = (func[i]+v2[i+1])/b[i];
	}

	finish = clock();
	time2 = ((finish-start)/(double) CLOCKS_PER_SEC);

//Write data file
	ofile << "#_x_i" << setw(15) << "v1(x_i)" << setw(15) << "v2(x_i)" << setw(15) << "u(x_i)" << setw(15) << "log(error1)" << setw(15) << "log(error2)" << endl;
	ofile.precision(8);
	ofile << "1.0000000" << setw(15) << "0.0000000" << setw(15) << "0.0000000" << endl;
	for(i=n-1;i>=0;i--){
		double exact, x, error1, error2;
		x = (i+1)*h;
		exact = 1.0 - (1.0 - exp(-10.0))*x-exp(-10.0*x);
		error1 = log10(abs(v1[i]-exact)/exact);
		error2 = log10(abs(v2[i]-exact)/exact);
		ofile << x << setw(15) << v1[i] << setw(15) << v2[i] << setw(15) << exact << setw(15) << error1 << setw(15) << error2 << endl;
	}
	ofile << "0.0000000" << setw(15) << "0.0000000" << setw(15) << "0.000000" << endl;

	ofile.close();
	
	cout.precision(5);
	cout << scientific;
	cout << "v1: " << time1 << "s" << endl;
	cout << "v2: " << time2 << "s" << endl;

	delete[] a, b, c, v1, v2, func;

//Compute error for varying step counts, if desired
	if(yesno==1){
		ofile.open("error.dat");
		ofile.precision(8);
		ofile << "#_n" << setw(15) << "log(h)" << setw(20) << "max_log(error)" << endl;
		double exact, x, temp, error;
		int k;
		for(j=0;j<5;j++){
			for(k=2;k<=10;k++){
				n = k*pow(10,j);
				double* b = new (nothrow) double[n];
				double* func = new (nothrow) double[n];
				double* v = new (nothrow) double[n];
				if(b==nullptr||func==nullptr||v==nullptr){
					cout << "Unable to allocate error computation memory for " << n << " grid points. Terminating..." << endl;
					exit(1);
				}
				h = 1.0/((double) (n+1));
				for(i=0;i<n;i++){
					b[i]=(double) (i+2)/(double) (i+1);
					func[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
				}
				for(i=1;i<n;i++){
					func[i] += func[i-1]/b[i-1];
				}
				v[n-1]=func[n-1]/b[n-1];
				for(i=n-2;i>=0;i--){
					v[i] = (func[i]+v[i+1])/b[i];
				}
				for(i=0;i<n;i++){
					error=-20.0;
					temp = 0.0;
					x = (i+1)*h;
					exact = 1.0 - (1.0 - exp(-10.0))*x-exp(-10.0*x);
					temp = log10(abs(v[i]-exact)/exact);
					if(temp>error) error = temp;
				}
				ofile << n << setw(20-j) << log10(h) << setw(15) << error << endl;

				delete[] b;
				delete[] func;
				delete[] v;
			}
		}
	}

	

	return 0;
}	
