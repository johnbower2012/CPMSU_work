/*Note that LU decomp is not possible for all matrices.*/

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
void Matrix_LUdecomp(int size, double**& a, double**& l, double**& u);
void Matrix_Lyf_solve(int size, double**& l, double*& y, double*& func);
void Matrix_Uvy_solve(int size, double**& u, double*& v, double*& y);
void Matrix_LUdecomp_solve_complete(int size, double**& a, double*& v, double*& func);

int main(int argc, char* argv[]){
	int i, j, k, n, size, timerfile_yesno, error_yesno;
	int sizetimefile = 0;
	double h, time;
	double* v_lu, *y, *func;
	double** a, **u, **l;
	char* outfilename;
	clock_t start, finish;
//Desired matrix size--output file name--'1' or '0' for timefile or numerical data file, respectively
	if(argc<5){
		cout << "Bad usage. Enter also matrix_size outfilename error_yesno timerfile_yesno ('1' or '0' for last two) on same line." << endl;
		exit(1);
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
	
//Allocate the arrays and matrices
		Array_alloc(v_lu, n); Array_alloc(y, n); Array_alloc(func, n);
		Matrix_alloc(a, n, n); Matrix_alloc(u, n, n); Matrix_alloc(l, n, n);

//Initialize func[i] values && a[i][j]
		for(i=0;i<n;i++){
			func[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
		}
		for(i=0;i<n;i++){
			for(j=0;j<n;j++){
				if(i==j) a[i][j]=2.0;
				else if((i+1)==j||(i-1)==j) a[i][j]=-1.0;
				else a[i][j]=0.0;
			}
		}


/*********************
Begin computation time
**********************/
		start = clock();

//Here we conduct the decomp. Calculate the first (second) row (column) of u (l), second (third) row (column) of u (l), and so forth to the n-1th row, then nth column, then nth row.
		Matrix_LUdecomp(n, a, l, u);
//Calculate y from L*y=func	
		Matrix_Lyf_solve(n, l, y, func);
//Calculcate v_lu from U*v_lu=y
		Matrix_Uvy_solve(n, u, v_lu, y);

		finish = clock();
/*******************
End computation time
********************/

		time = ((finish-start)/(double) CLOCKS_PER_SEC);

		if(timerfile_yesno==0)	cout << '\t' << time << " seconds" << endl;

//Write to file
		if(timerfile_yesno==0){
			for(i=0;i<n;i++){
				ofile << h*(double) (i+1) << setw(20) << v_lu[i] << endl;
			}
		}
		else if(timerfile_yesno==1){
			ofile.precision(6);
			ofile << n << setw(15) << time << endl;
		}

		Array_dealloc(y); Array_dealloc(v_lu); Array_dealloc(func);
		Matrix_dealloc(a,n); Matrix_dealloc(u,n); Matrix_dealloc(l,n);
	}
	
	ofile.close();

/**************************
Calculate error if desired
***************************/


	if(error_yesno==1){
		ofile.open("LUerror.dat");
		ofile.precision(8);
		ofile << "#_n" << setw(15) << "log(h)" << setw(20) << "max_log(error)" << endl;
		double exact, x, temp, error;
		int m;
		for(m=0;m<3;m++){
			for(k=2;k<=10;k++){
				n=k*pow(10,m);
				h=(1.0/(double) (n+1));
				Array_alloc(v_lu, n); Array_alloc(func, n);
				Matrix_alloc(a, n, n);
				
				for(i=0;i<n;i++){
					func[i]=pow(h, 2.0)*100.0*exp(-10.0*h*((double) (i+1)));
				}
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						if(i==j) a[i][j]=2.0;
						else if((i+1)==j||(i-1)==j) a[i][j]=-1.0;
						else a[i][j]=0.0;
					}
				}

				Matrix_LUdecomp_solve_complete(n, a, v_lu, func);

				error=-20.0;
				for(i=0;i<n;i++){
					temp = 0.0;
					x = (i+1)*h;
					exact = 1.0 - (1.0 - exp(-10.0))*x-exp(-10.0*x);
					temp = log10(abs(v_lu[i]-exact)/exact);
					if(temp>error) error = temp;
				}

				ofile << n << setw(20-m) << log10(h) << setw(15) << error << endl;

				Array_dealloc(v_lu); Array_dealloc(func);
				Matrix_dealloc(a, n);
			}
		}
	}


	return 0;
}




/**************************
Begin function definitions
***************************/

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

void Matrix_LUdecomp(int size, double**& a, double**& l, double**& u){
	int i, j, k;	
	for(i=0; i<size;i++){
		for(j=0;j<size;j++){
			l[i][j]=0.0;
			u[i][j]=0.0;
		}
		l[i][i]=1.0;
	}
	for(i=0; i<size;i++){
		double temp;
		for(j=i;j<size;j++){
			temp = 0.0;
			for(k=0;k<i;k++){
				temp +=l[i][k]*u[k][j];
			}
			u[i][j] = a[i][j] - temp;
		}
		for(j=i+1;j<size;j++){
			temp = 0.0;
			for(k=0;k<i;k++){
				temp +=l[j][k]*u[k][i];
			}
			 l[j][i] = (a[j][i] - temp)/u[i][i];
		}
	}
}

void Matrix_Lyf_solve(int size, double**& l, double*& y, double*& func){
	int i, j;	
	for(i=0;i<size;i++){
		double temp = 0.0;
		for(j=0;j<i;j++){
			temp += l[i][j]*y[j];
		}
		y[i]=func[i]-temp;
	}
}

void Matrix_Uvy_solve(int size, double**& u, double*& v, double*& y){
	int i, j;	
	for(i=size-1;i>=0;i--){
		double temp = 0.0;
		for(j=i+1;j<size;j++){
			temp += u[i][j]*v[j];
		}
		v[i] = (y[i]-temp)/u[i][i];
	}
}

void Matrix_LUdecomp_solve_complete(int size, double**& a, double*& v, double*& func){
	int i, j, k;
	double** l, **u;
	double* y;
	l = new double* [size];
	for(int i=0;i<size;i++){
		l[i] = new double [size];
	}	
	u = new double* [size];
	for(int i=0;i<size;i++){
		u[i] = new double [size];
	}
	y = new double[size];

	for(i=0; i<size;i++){
		for(j=0;j<size;j++){
			l[i][j]=0.0;
			u[i][j]=0.0;
		}
		l[i][i]=1.0;
	}

	for(i=0; i<size;i++){
		double temp;
		for(j=i;j<size;j++){
			temp = 0.0;
			for(k=0;k<i;k++){
				temp +=l[i][k]*u[k][j];
			}
			u[i][j] = a[i][j] - temp;
		}
		for(j=i+1;j<size;j++){
			temp = 0.0;
			for(k=0;k<i;k++){
				temp +=l[j][k]*u[k][i];
			}
			 l[j][i] = (a[j][i] - temp)/u[i][i];
		}
	}

	for(i=0;i<size;i++){
		double temp = 0.0;
		for(j=0;j<i;j++){
			temp += l[i][j]*y[j];
		}
		y[i]=func[i]-temp;
	}
	for(i=size-1;i>=0;i--){
		double temp = 0.0;
		for(j=i+1;j<size;j++){
			temp += u[i][j]*v[j];
		}
		v[i] = (y[i]-temp)/u[i][i];
	}

	for(int i=0;i<size;i++){
		delete[] l[i];
		delete[] u[i];
	}
	delete[] l;
	delete[] u;
	delete[] y;
}
