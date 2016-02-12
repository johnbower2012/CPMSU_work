/***************************************************************
This file is designed to solve a second order differential
equation by discretizing the domain into x_i's and forming a 
linear algebra problem. Given -u"(x) = f(x), u(1)=u(0)=0, we 
approximate -u"(x_i) = -(v_i+1 + v_i-1 - 2*v_i)/h^2 = f(x_i), and 
solve A*v=f, where A is a tridiagonal matrix, v is an array of 
our numerical solution, and f is an array of our f(x_i). Each has 
dimensionality n, where n is our number of grid points. We 
define h, our step size, as 1/(n+1) so as not to include the end 
points, which we fix by definition.

Here we simplify the problem due to A's tridiagonal form, first
solving for a general tridiagonal matrix, then recognizing the 
inherent symmetry to reduce the number of floating points 
operatings required by half. We speed the process by allocating 
space for only three arrays to represent the diagonal and 
off-diagonal elements, as shown below. 

A = ( b_1 c_1  0   0   0  ...  0   0  )
    ( a_2 b_2 c_2  0   0  ...  0   0  )
    (  0  a_3 b_3 c_3  0  ...  0   0  )
    ( ... ... ... ... ... ... ... ... )
    (  0   0   0   0  ...  0  a_n b_n )

****************************************************************/

#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include "time.h"

using namespace std;

ofstream ofile;

void Array_alloc(double*& a, int length);
void Matrix_alloc(double**& a, int rows, int columns);
void Array_dealloc(double*& a);
void Matrix_dealloc(double**& a, int rows);
void gen_tridmat_equ_solve(int n, double*& a, double*& b, double*& c, double*& v, double*& func);
void secder_tridmat_equ_solve(int n, double*& b, double*& v, double*& func);

int main(int argc, char* argv[]){
	char* outfilename;
	int n, i, j, k, error_yesno, timer_data_nofile, size;
	int sizefile = 0, sizeerror=0;
	double h, time1, time2;
	double* a, *b, *c, *v1, *v2, *func;
	clock_t start, finish;

	if(argc<5){
		cout << "Bad usage. Enter also grid_point_count outfilename error_yesno ('1' or '0') timer_data_nofile on same ('2','1', or '0') line." << endl;
		exit (1);
	}
	else{
		size = atoi(argv[1]);
		outfilename = argv[2];
		error_yesno = atoi(argv[3]);
		timer_data_nofile = atoi(argv[4]);
		if(timer_data_nofile==2||timer_data_nofile==1){
			ofile.open(outfilename);
		}
	}

/***************************************************************
Here the process of filtering the input data occurs.
If timer_data_nofile = 2, we created the appropriate limits for
the time file, more information below. If timer_data_nofile = 1, 
we set the forloop for one pass and merely set n to the desired 
size. If timer_data_nofile equals 0, we forgo file creation 
altogether for the loop. If error_yesno = 1, we set the 
appropriate size for sizeerror for the error loop. If error_yesno
equals 0, we forgo sizing sizeerror altogether since the loop 
will not be executed.
****************************************************************/
	
	if(timer_data_nofile==2){
		do{
			sizefile++;
			j = size/pow(10.0,sizefile);
		} while(j>1);
		if(j==1) sizefile*=2;
		if(j==0) sizefile = 2*sizefile - 1;
	}
	else if(timer_data_nofile==1){
		n=size;
		sizefile=1;
	}
	else if(timer_data_nofile==0){
		sizefile=0;
	}
	if(error_yesno==1){
		do{
			sizeerror++;
			j=size/pow(10.0,sizeerror);
		} while (j>1);
	}

/***************************************************************
Note the structure presented here. If a time file is called for, 
we perform the loops at n = 5, 10, 50, 100, 500, 100, etc until 
the first number in this series greater than or equal to the 
desired matrix size, writing each respective n and time to the
entered outfilename. If a numerical data file is called for, we
simply perform one loop and output the data itself.
****************************************************************/

	for(k=1;k<=sizefile;k++){
		if(timer_data_nofile==2){
			if(k%2==0) n=pow(10.0,k/2);
			else n=pow(10.0,(k+1)/2)/2;
		}
	
		h=(1.0/(double) (n+1));

//Dynamically allocate memory
	
		Array_alloc(a,n); Array_alloc(b,n); Array_alloc(c,n); Array_alloc(v1,n); Array_alloc(v2,n); Array_alloc(func, n);

//Initialize array values

		for(i=0;i<n;i++){
			a[i]=(-1.0);
			b[i]=(2.0);
			c[i]=(-1.0);
			func[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
		}

/********************************************
  Start clock for general tridiagonal matrix
  Create upper diagonal matrix
  Back substitutions
  End clock
*********************************************/

		start = clock();
	
		gen_tridmat_equ_solve(n, a, b, c, v1, func);

		finish = clock();
		time1 = ( (finish - start)/(double) CLOCKS_PER_SEC);

//redfine b values for our specific tridiagonal, reinitialize func[i] values

		for(i=0;i<n;i++){
			b[i]=(double) (i+2)/(double) (i+1);
			func[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
		}

/********************************************
  Start clock for our specific tridiagonal
  Form A into upper diagonal
  Back substitutions
  End clock
*********************************************/

		start = clock();

		secder_tridmat_equ_solve(n, b, v2, func);
	
		finish = clock();
		time2 = ((finish-start)/(double) CLOCKS_PER_SEC);

//Write data file

		if(timer_data_nofile==2){
			ofile.precision(6);
			if(k==1){
				ofile << "#_n" << setw(15) << "log10(n)" << setw(15) << "time_gen" << setw(15);
				ofile << "log10(time_gen)" << setw(15) << "time_our" << setw(15) << "log10(time_our)" << endl;	
			}						
			ofile << n << setw(15) << log10(n) << setw(15) << time1 << setw(15)<< log10(time1) << setw(15) << time2 << setw(15) << log10(time2) << endl;
		}
		else if(timer_data_nofile==1){
			ofile << "#_x_i" << setw(15) << "v1(x_i)" << setw(15) << "v2(x_i)" << setw(15) << "u(x_i)" << endl;
			ofile.precision(8);		
			ofile << "1" << setw(15) << "0" << setw(15) << "0" << endl;
			for(i=n-1;i>=0;i--){
				double exact, x;
				x = (i+1)*h;
				exact = 1.0 - (1.0 - exp(-10.0))*x-exp(-10.0*x);
				ofile << x << setw(15) << v1[i] << setw(15) << v2[i] << setw(15) << exact << endl;
			}
			ofile << "0" << setw(15) << "0" << setw(15) << "0" << endl;
		}
		
		Array_dealloc(a); Array_dealloc(b); Array_dealloc(c);
		Array_dealloc(v1); Array_dealloc(v2); Array_dealloc(func);
	}
	
	if(timer_data_nofile==2||timer_data_nofile==1){
		ofile.close();
	}

//Write time to screen if creating a data file
	
	if(timer_data_nofile==1){
		cout.precision(5);
		cout << scientific;
		cout << "\tv1: " << time1 << " seconds" << endl;
		cout << "\tv2: " << time2 << " seconds" << endl;
	}

//Compute error for varying step counts, if desired

	if(error_yesno==1){
		ofile.open("trid_error.dat");
		ofile.precision(8);
		ofile << "#_n" << setw(15) << "log(h)" << setw(20) << "max_log(error1)" << setw(15) << "max_log(error2)" << endl;
		double exact, x, temp1, temp2, error1, error2;
		int m;
		for(m=0;m<sizeerror;m++){
			for(k=2;k<=10;k++){
				n = k*pow(10,m);
				double* a, *b1, *b2, *c, *func1, *func2, *v1, *v2;
				Array_alloc(a,n); Array_alloc(b1,n); Array_alloc(b2,n); Array_alloc(c,n);
				Array_alloc(func1,n); Array_alloc(func2,n); Array_alloc(v1,n); Array_alloc(v2,n);
				
				h = 1.0/((double) (n+1));
				for(i=0;i<n;i++){
					a[i]=(-1.0);
					b1[i]=(2.0);
					c[i]=(-1.0);
					func1[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
					b2[i]=(double) (i+2)/(double) (i+1);
					func2[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
				}	

				gen_tridmat_equ_solve(n, a, b1, c, v1, func1);
				secder_tridmat_equ_solve(n, b2, v2, func2);
				
				error1=-20.0;
				error2=-20.0;
				for(i=0;i<n;i++){
					temp1 = 0.0;
					temp2 = 0.0;
					x = (double) (i+1)*h;
					exact = 1.0 - (1.0 - exp(-10.0))*x-exp(-10.0*x);
					temp1 = log10(abs(v1[i]-exact)/exact);
					temp2 = log10(abs(v2[i]-exact)/exact);
					if(temp1>error1) error1 = temp1;
					if(temp2>error2) error2 = temp2;
				}
				ofile << n << setw(20-m) << log10(h) << setw(15) << error1 << setw(15) << error2 << endl;

				Array_dealloc(a); Array_dealloc(b1); Array_dealloc(b2); Array_dealloc(c);
				Array_dealloc(func1); Array_dealloc(func2); Array_dealloc(v1); Array_dealloc(v2);
			}
		}
	}


	return 0;
}	

void Array_alloc(double*& a, int length){
	a = new double[length];
}

void Matrix_alloc(double**& a, int rows, int columns){
	a = new double* [rows];
	for(int count=0;count<rows;count++){
		a[count] = new double [columns];
	}
}

void Array_dealloc(double*& a){
	delete[] a;
}

void Matrix_dealloc(double**& a, int rows){
	for(int i=0;i<rows;i++){
		delete[] a[i];
	}
	delete[] a;
}

void gen_tridmat_equ_solve(int n, double*& a, double*& b, double*& c, double*& v, double*& func){
	int i;
	for(i=1;i<n;i++){
		double temp;
		temp = a[i]/b[i-1];
		b[i] -= temp*c[i-1];
		func[i] -= temp*func[i-1];
	}
	v[n-1]=func[n-1]/b[n-1];
	for(i=n-2;i>=0;i--){
		v[i] = (func[i]-c[i]*v[i+1])/b[i];
	}
}

void secder_tridmat_equ_solve(int n, double*& b, double*& v, double*& func){
	int i;	
	for(i=1;i<n;i++){
		func[i] += func[i-1]/b[i-1];
	}
	v[n-1]=func[n-1]/b[n-1];
	for(i=n-2;i>=0;i--){
		v[i] = (func[i]+v[i+1])/b[i];
	}
}

