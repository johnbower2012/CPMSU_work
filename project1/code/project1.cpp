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


/*
Functions to condense the main() function. 
The first is "general tridiagonal matrix equation solver."
The second is "second derivative ----," made for our specific form of the second derivative.
*/
void Array_alloc(double*& a, int length);
void Matrix_alloc(double**& a, int rows, int columns);
void Array_dealloc(double*& a);
void Matrix_dealloc(double**& a, int rows);
void gentridmat_equ_solve(int n, double*& a, double*& b, double*& c, double*& v, double*& func);
void secder_tridmat_equ_solve(int n, double*& b, double*& v, double*& func);

int main(int argc, char* argv[]){
	char* outfilename;
	int n, i, j, k, error_yesno, timerfile_yesno, size;
	int sizetimefile = 0;
	double h, time1, time2;
	double* a, *b, *c, *v1, *v2, *func;
	clock_t start, finish;

	if(argc<5){
		cout << "Bad usage. Include also in command line gridpoint_count outfilename '1 or 0 (for error)' '1 or 0' for time file" << endl;
		exit (1);
	}
	else{
		size = atoi(argv[1]);
		outfilename = argv[2];
		error_yesno = atoi(argv[3]);
		timerfile_yesno = atoi(argv[4]);
		ofile.open(outfilename);
	}

//Here we adjust our planned executions for the timefile forloop	
	do{
		sizetimefile++;
		j = size/pow(10.0,sizetimefile);
	} while(j>1);
	if(j==1) sizetimefile*=2;
	if(j==0) sizetimefile = 2*sizetimefile - 1;
/***************************************************************
Note the structure presented here. If a time file is called for, 
we perform the loops at n = 5, 10, 50, 100, 500, 100, etc until 
the first number in this series greater than or equal to the 
desired matrix size, writing each respective n and time to the
entered outfilename. If a numerical data file is called for, we
simply perform one loop and output the data itself.
****************************************************************/
	for(k=1;k<=sizetimefile;k++){
		if(timerfile_yesno==0){
			k=sizetimefile;
			n = size;
		}
		else if(timerfile_yesno==1){
			if(k%2==0) n=pow(10.0,k/2);
			else n=pow(10.0,(k+1)/2)/2;
		}
	
		h=(1.0/(double) (n+1));

//Dynamically allocate memory, with check.	
		Array_alloc(a, n); Array_alloc(b,n); Array_alloc(c,n); Array_alloc(v1,n); Array_alloc(v2,n); Array_alloc(func, n);

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
	
		gentridmat_equ_solve(n, a, b, c, v1, func);

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

		secder_tridmat_equ_solve(n, b, v2, func);
	
		finish = clock();
		time2 = ((finish-start)/(double) CLOCKS_PER_SEC);

//Write data file
		if(timerfile_yesno==0){
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
		else if(timerfile_yesno==1){
			ofile.precision(6);			
			ofile << n << setw(15) << time1 << setw(15) << time2 << endl;
		}
		
	Array_dealloc(a); Array_dealloc(b); Array_dealloc(c); Array_dealloc(v1); Array_dealloc(v2); Array_dealloc(func);
	}

	ofile.close();
	
	if(timerfile_yesno==0){
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
		int k;
		for(j=0;j<7;j++){
			for(k=2;k<=10;k++){
				n = k*pow(10,j);
				double* a, *b1, *b2, *c, *func1, *func2, *v1, *v2;
				Array_alloc(a,n); Array_alloc(b1,n); Array_alloc(b2,n); Array_alloc(c,n);
				Array_alloc(func1,n); Array_alloc(func2,n); Array_alloc(v1,n); Array_alloc(v2,n);
				if(a==nullptr||b1==nullptr||b2==nullptr||c==nullptr||func1==nullptr||func2==nullptr||v1==nullptr||v2==nullptr){
					cout << "Unable to allocate error computation memory for " << n << " grid points. Terminating..." << endl;
					exit(1);
				}
				h = 1.0/((double) (n+1));
				for(i=0;i<n;i++){
					a[i]=(-1.0);
					b1[i]=(2.0);
					c[i]=(-1.0);
					func1[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
					b2[i]=(double) (i+2)/(double) (i+1);
					func2[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
				}	

				gentridmat_equ_solve(n, a, b1, c, v1, func1);
				secder_tridmat_equ_solve(n, b2, v2, func2);
				
		
				error1=-20.0;
				error2=-20.0;
				for(i=0;i<n;i++){
					temp1 = 0.0;
					temp2 = 0.0;
					x = (i+1)*h;
					exact = 1.0 - (1.0 - exp(-10.0))*x-exp(-10.0*x);
					temp1 = log10(abs(v1[i]-exact)/exact);
					temp2 = log10(abs(v2[i]-exact)/exact);
					if(temp1>error1) error1 = temp1;
					if(temp2>error2) error2 = temp2;
				}
				ofile << n << setw(20-j) << log10(h) << setw(15) << error1 << setw(15) << error2 << endl;

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

void gentridmat_equ_solve(int n, double*& a, double*& b, double*& c, double*& v, double*& func){
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

