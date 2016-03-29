#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include "vectorclass.h"

using namespace std;

int main(){
	int n = 3;
	vector a(3);
	vector b(3);
	vector c(3);
	a.array[0] = 1;
	b.array[1] = 1;

	c = a/b;

	for(int i = 0; i<n; i++){
		cout << a.array[i] << '\t';
	}
	cout << endl;
	for(int i = 0; i<n; i++){
		cout << b.array[i] << '\t';
	}
	cout << endl;
	for(int i = 0; i<n; i++){
		cout << c.array[i] << '\t';
	}
	cout << endl;
	a.del();
	b.del();
	c.del();
}
