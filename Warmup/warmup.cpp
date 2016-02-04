#include<iostream>
#include<cmath>

using namespace std;

int main (int argc, char* argv[]){
// ns = number of steps  (need at least 2); isl = initial step length; fsl = final step length
// x = value to be computed
	int ns;
	double isl, fsl, x; 
	double* sind;
	x = atof(argv[1]);
	ns = atoi(argv[2]);
	isl = atof(argv[3]);
	fsl = atof(argv[4]);
	sind = new (nothrow) double [ns];
//Check for proper dynamic memory allocation
	if (sind==nullptr) cout << "Error allocating memory. Terminating...\n";
	else {
		for(int i=0; i<ns; i++){
			sind[i]= isl + (fsl-isl)*((double) i)/((double) (ns-1));
		}

/* Test for functionality of dynamnic step length array

		for(int i=0; i<ns; i++){
			cout << "Step Length " << i << ": " << sind[i] << '\n';
		}
*/

//Array of derivative values
		double* compder1;
		double* compder2;
		compder1 = new (nothrow) double [ns];
		compder2 = new (nothrow) double [ns];
		if(compder1==nullptr||compder2==nullptr) cout << "Error allocating memory. Terminating...\n";
		else{
			for(int i=0; i<ns; i++){
				compder1[i]=((atan(x+sind[i])-atan(x))/sind[i]);
				compder2[i]=((atan(x+sind[i])-atan(x-sind[i]))/(2*sind[i]));
			}
			for(int i=0; i<ns; i++){
				cout.precision(8);
				cout << scientific;
				cout << "h" << i+1 << ": " << sind[i] << '\t';
				cout << "value1: " << compder1[i] << '\t';
				cout << "value2: " << compder2[i] << '\n';
			}
			delete[] sind;
			delete[] compder1;
			delete[] compder2;
		}		
	}
	return 0;
}
		
		

