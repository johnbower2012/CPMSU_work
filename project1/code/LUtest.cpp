/*A test to implement my LU decomp methodology. Note that LU decomp is not possible for all matrices.*/

#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include "time.h"

using namespace std;

ofstream ofile;

int main(int argc, char* argv[]){
	int i, j, k, n;
	double h, time;
	double* v_lu, *y, *func;
	char* outfilename;
	clock_t start, finish;

	if(argc<3){
		cout << "Bad usage. Enter also matrix_size out_file_name on same line." << endl;
		exit(1);
	}
	else{
		n = atoi(argv[1]);
		outfilename = argv[2];
		h=(1.0/(n+1));
	}

	ofile.open(outfilename);

	v_lu = new double[n]; y = new double[n]; func = new double[n];
	double** a = new double* [n], **u = new double*[n], **l = new double*[n], **a2 = new double*[n];
	for(i=0; i<n;i++){
		a[i] = new double[n];
		u[i] = new double[n];
		l[i] = new double[n];
		a2[i] = new double[n];
	}

//Initialize func[i] values && a[i][j], u[i][j], l[i][j], a2[i][j]
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
	for(j=0;j<n;j++){
		for(i=0;i<n;i++){
			l[i][j] = 0.0;
			u[i][j] = 0.0;
			a2[i][j] = 0.0;
		}
		l[j][j] = 1.0;
	}
//Start time;
	start = clock();

//Here we conduct the decomp. Calculate the first (second) row (column) of u (l), second (third) row (column) of u (l), and so forth to the n-1th row, then nth column, then nth row.
	for(i=0; i<n;i++){
		double temp;
		for(j=i;j<n;j++){
			temp = 0.0;
			for(k=0;k<i;k++){
				temp +=l[i][k]*u[k][j];
			}
			u[i][j] = a[i][j] - temp;
		}
		for(j=i+1;j<n;j++){
			temp = 0.0;
			for(k=0;k<i;k++){
				temp +=l[j][k]*u[k][i];
			}
			 l[j][i] = (a[j][i] - temp)/u[i][i];
		}
	}
//Calculate y from L*y=func	
	for(i=0;i<n;i++){
		double temp = 0.0;
		for(j=0;j<i;j++){
			temp += l[i][j]*y[j];
		}
		y[i]=func[i]-temp;
	}
//Calculcate v_lu from U*v_lu=y
	for(i=n-1;i>=0;i--){
		double temp = 0.0;
		for(j=i+1;j<n;j++){
			temp += u[i][j]*v_lu[j];
		}
		v_lu[i] = (y[i]-temp)/u[i][i];
	}
//End clock
	finish = clock();
	time = ((finish-start)/(double) CLOCKS_PER_SEC);
	cout << time << endl;
//Write to file
	for(i=0;i<n;i++){
		ofile << h*(double) (i+1) << setw(20) << v_lu[i] << endl;
	}
	
	ofile.close();

	for(i=0;i<n;i++){
		delete[] a[i];
		delete[] u[i];
		delete[] l[i];
		delete[] a2[i];
	}
	delete[] a;
	delete[] u;
	delete[] l;
	delete[] a2;
	delete[] y;
	delete[] func;
	delete[] v_lu;

	return 0;
}
				
