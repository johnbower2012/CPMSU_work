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
void array_print(double*& a, int length);
void matrix_print(double**& A, int rows, int columns);
void matrix_print_mult(double**& A, double**& B, int rowsA, int columnsA, int columnsB);
void matrix_diag_sort(double**& A, double**&B, int size);

ofstream ofile;

int main(int argc, char* argv[]){
	double** A, **R, **A2, **eigenvalues;
	char* outfilename;
	double tolerance, var, var1, var2;
	double h, h_sq, rho, rho_sq, rho_max, rho_min;
	double temp = 0.0;
	int i, j, k, rows, columns, n, eigvec_yesno;
	tolerance = 1e-10;
//Check proper entry arguments
	if(argc<6){
		cout << "Bad usage. Enter also rho_max rho_min matrix_size outfilename of_1_0 on same line." << endl;
		exit(1);
	}
	else{
		rho_max = atof(argv[1]);
		rho_min = atof(argv[2]);
		n = atoi(argv[3]);
		outfilename = argv[4];
		eigvec_yesno = atoi(argv[5]);
		rows = n-1;
		columns = rows;
		h = (rho_max - rho_min)/((double) n);
		h_sq = (h*h);
		var1 = 2.0/h_sq;
		var2 = -1.0/h_sq;
	}
//allocate matrices
	matrix_alloc(eigenvalues,rows, 2);
	matrix_alloc(A, rows, columns);
	matrix_alloc(A2, rows, columns);
	matrix_alloc(R, rows, columns);
//Initialize matrices
	for(i=0;i<rows;i++){
		rho = rho_min + (i+1)*h;
		rho_sq = rho*rho;
		eigenvalues[i][0] = 0.0;
		eigenvalues[i][1] = 0;
		for(j=0;j<columns;j++){
			if(i==j){
				A[i][j] = var1 + rho_sq;
				A2[i][j] = A[i][j];
				R[i][j] = 1.0;
			}
			else if(i==j+1||i==j-1){
				A[i][j] = var2;
				A2[i][j] = A[i][j];
			}
			else{
				A[i][j] = 0.0;
				A2[i][j] = 0.0;
				R[i][j] = 0.0;
			}
		}
	}	

//conduct similarity rotations until all offdiag elements < tolerance.
	jacobi_simrot_eigen_solver(A,R,rows,tolerance);
//Sort the eigenvalues from A from smallest to largest and note original positions
	matrix_diag_sort(A,eigenvalues,rows);

//Open output file
	ofile.open(outfilename);
	ofile.precision(6);
	if(eigvec_yesno==1){
		ofile << "#_i" << setw(15) << "eigenvalues" << setw(12) << "rho_i";
		for(i=0;i<rows;i++){
			ofile << setw(14) << "v_" << i;
		}
		ofile << endl;
		for(i=0;i<rows;i++){
			ofile << fixed << i << setw(15) << eigenvalues[i][0] << setw(15) << scientific << rho_min + h*(i+1);
			for(j=0;j<columns;j++){
				k=eigenvalues[j][1];
				ofile << setw(15) << R[i][k];
			}
			ofile << endl;
		}
	}
	else if(eigvec_yesno==0){
		ofile << "#_i" << setw(15) << "eigenvalues" << setw(15) << "j" << endl;
		for(i=0;i<rows;i++){
			ofile << i << setw(15) << eigenvalues[i][0] << setw(15) << eigenvalues[i][1] << endl;
		}
	}
	ofile.close();
	
	matrix_delete(eigenvalues,rows);
	matrix_delete(A2,rows);
	matrix_delete(A,rows);
	matrix_delete(R,rows);	
		

	return 0;
}


/******************************
Begin function definitions
*******************************/


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
void array_print(double*& a, int length){
	for(int i=0;i<length;i++){
		cout << a[i] << endl;
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
void matrix_diag_sort(double**& A, double**&B, int size){
	int i, j, q;	
	for(i=0;i<size;i++){
		if(i==0){
			B[0][0]=A[0][0];
			B[0][1] = 0;
		}
		else{
			if(A[i][i]>=B[i-1][0]){
				B[i][0]=A[i][i];
				B[i][1]=i;
			}
			else{
				for(j=0;j<i;j++){
					if(A[i][i]<B[j][0]){
						for(q=i;q>j;q--){
							B[q][0] = B[q-1][0];
							B[q][1] = B[q-1][1];
						}
						B[j][0] = A[i][i];
						B[j][1] = i;
						break;
					}
				}
			}
		}
	}	
}


/**********************************************************
Print matrices to screen to test for reliability visually
**********************************************************/

/***************************************************************************************
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
	cout << fixed;
	cout << "\nEigenvalues:\n";
	for(i=0;i<rows;i++){
		cout << lambda[i] << endl;
	}
	cout << endl;
****************************************************************************************/
