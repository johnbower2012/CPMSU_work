#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>

using namespace std;

void array_alloc(double*& a, int length);
void matrix_alloc(double**& A, int rows, int columns);
void array_delete(double*& a);
void matrix_delete(double**& A, int rows);
void symmat_offdiag_max(double**& A, int size, double& a_max, int& row, int& column);
void jacobi_simrot_eigen_solver(double**& A, double**& R, int size, double tolerance);
void matrix_print(double**& A, int rows, int columns);
void matrix_print_mult(double**& A, double**& B, int rowsA, int columnsA, int columnsB);

int main(int argc, char* argv[]){
	double** A, **R, **A2;
	double a_max, a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	double tau, cos, sin, tan, tolerance, var; 
	int i, j, k, l, rows, columns;
	a_max = 1;
	tolerance = 1e-10;
//Check proper entry arguments	
	if(argc<2){
		cout << "Bad usage. Enter also matrix_size on same line." << endl;
		exit(1);
	}
	else{
		rows = atoi(argv[1]);
		columns = rows;
	}
//allocate matrices
	matrix_alloc(A, rows, columns);
	matrix_alloc(A2, rows, columns);
	matrix_alloc(R, rows, columns);
//Initialize matrices
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			if(i==j){
				A[i][j] = 2.0; 
				A2[i][j] = 2.0; 
				R[i][j] = 1.0;			
			}
			else{
				A[i][j] = 1.0;
				A2[i][j] = 1.0;
				R[i][j] = 0.0;
			}
		}
	}
//conduct similarity rotations until 
	jacobi_simrot_eigen_solver(A,R,rows, tolerance);

	cout.precision(6);
	cout << "\nMatrix A:\n";
	matrix_print(A2,rows,columns);
	cout << "\nMatrix D:\n";
	matrix_print(A,rows,columns);
	cout << scientific;
	cout << "\nMatrix R:\n";
	matrix_print(R,rows,columns);
	cout << "\nMatrix lambdaR:\n";
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			cout << A[j][j]*R[i][j] << '\t';
		}
		cout << '\n';
	}
	cout << "\nMatrix A*R:\n";
	matrix_print_mult(A2,R,rows,columns,columns);
	cout << '\n';
	for(i=0;i<rows;i++){
		var = 0.0;
		for(j=0;j<rows;j++){
			var += R[i][j]*R[i][j];
		}
		var = sqrt(var);
		cout << "Unit test for eigenvector" << i+1 << ": " << var << '\n';
	}
	cout << '\n';
	matrix_delete(A2,rows);
	matrix_delete(A,rows);
	matrix_delete(R,rows);	
		

	return 0;
}





void array_alloc(double*& a, int length){
	a = new double[length];
}
void matrix_alloc(double**& A, int rows, int columns){
	A = new double*[rows];
	for(int i=0;i<rows;i++){
		A[i] = new double[columns];
	}
}
void array_delete(double*& a){
	delete[] a;
}
void matrix_delete(double**& A, int rows){
	for(int i=0;i<rows;i++){
		delete[] A[i];
	}
	delete[] A;
}
void symmat_offdiag_max(double**& A, int size, double& a_max, int& row, int& column){
	double temp = 0.0;
	for(int i=0;i<size;i++){
		for(int j=i+1;j<size;j++){
			temp = fabs(A[i][j]);
			if(temp>a_max){
				a_max = temp;
				row = i;
				column = j;
			}
		}
	}
}
void jacobi_simrot_eigen_solver(double**& A, double**& R, int size, double tolerance){
	double a_ik, a_il, a_kk, a_ll, r_ik, r_il;
	double tau, tan, cos, sin;	
	int i, j, k, l;	
	double a_max = 0.0;
	double temp = 0.0;
	for(i=0;i<size;i++){
		for(j=i+1;j<size;j++){
			temp = fabs(A[i][j]);
			if(temp>a_max){
				a_max = temp;
				k = i;
				l = j;
			}
		}
	}
	while(a_max>tolerance){
		tau = (A[l][l] - A[k][k])/(2*A[k][l]);
		if(tau>=0) tan = -tau + sqrt(tau*tau+1);
		else if(tau<0) tan = -tau - sqrt(tau*tau+1);
		cos = 1.0/sqrt(tan*tan + 1.0);
		sin = tan*cos;
	
		a_kk = A[k][k];
		a_ll = A[l][l];
		A[k][k] = a_kk*cos*cos + a_ll*sin*sin - 2.0*sin*cos*A[k][l];
		A[l][l] = a_kk*sin*sin + a_ll*cos*cos + 2.0*sin*cos*A[k][l];
		A[k][l] = 0.0;
		A[l][k] = 0.0;
		for(i=0;i<size;i++){
			if(i!=k && i!=l){
				a_ik = A[i][k];
				a_il = A[i][l];
				A[i][k] = a_ik*cos - a_il*sin;
				A[k][i] = A[i][k];
				A[i][l] = a_il*cos + a_ik*sin;
				A[l][i] = A[i][l];
			}
			r_ik = R[i][k];
			r_il = R[i][l];
			R[i][k] = r_ik*cos - r_il*sin;
			R[i][l] = r_ik*sin + r_il*cos;
		}
		a_max = 0.0;			
		for(i=0;i<size;i++){
			for(j=i+1;j<size;j++){
				temp = fabs(A[i][j]);
				if(temp>a_max){
					a_max = temp;
					k = i;
					l = j;
				}
			}
		}
	}
}
void matrix_print(double**& A, int rows, int columns){
	for(int i=0;i<rows;i++){
		for(int j=0;j<columns;j++){
			cout << A[i][j] << '\t';
		}
		cout << '\n';
	}
}
void matrix_print_mult(double**& A, double**& B, int rowsA, int columnsA, int columnsB){
	double var;	
	for(int i=0;i<rowsA;i++){
		for(int j=0;j<columnsB;j++){
			var = 0.0;
			for(int p=0;p<columnsA;p++){
				var += A[i][p]*B[p][j];
			}
			cout << var << '\t';
		}
		cout << '\n';
	}
}



