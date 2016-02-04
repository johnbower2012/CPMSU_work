#include<iostream>
#include<cmath>

using namespace std;

int main(){
	double x, s, d;
	cin>>x;
	cin>>s;
	double a[10];
	double b[10];
	double c[10];
	d = log10(s)-log10(x);
	for(int i=0; i<10; i++){
		a[i] = x + (s-x)*i/9;
		b[i] = log10(x) + d*i/9;
		c[i] = pow(10.0,b[i]);
		cout << a[i] << '\t' << b[i] << '\t' << c[i] << '\n';
	}
}
