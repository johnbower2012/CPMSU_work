#include<iostream>
#include<iomanip>

using namespace std;

inline void threeDarray_alloc(double***& a, int d1, int d2, int d3);
inline void threeDarray_delete(double***& a, int d1, int d2);

int main(int argc, char* argv[]){
	double*** a;
	int d1, d2, d3;
	int i, j, k;
	d1 = atoi(argv[1]);
	d2 = atoi(argv[2]);
	d3 = atoi(argv[3]);

	threeDarray_alloc(a, d1, d2, d3);

	cout << endl;

	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				a[i][j][k] = i + j + k;
				cout << '\t' << a[i][j][k];
			}
			cout << endl;
		}
		cout << endl;
	}

	threeDarray_delete(a, d1, d2);

	return 0;
}

inline void threeDarray_alloc(double***& a, int d1, int d2, int d3){
	int i, j, k;	
	a = new double**[d1];
	for(i=0;i<d1;i++){
		a[i] = new double*[d2];
	}
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			a[i][j] = new double[d3];
		}
	}
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				a[i][j][k] = 0.0;
			}
		}
	}
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				a[i][j][k] = 0.0;
			}
		}
	}
}
inline void threeDarray_delete(double***& a, int d1, int d2){
	int i, j;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			delete[] a[i][j];
		}
	}
	for(i=0;i<d1;i++){
		delete[] a[i];
	}
	delete[] a;
}
