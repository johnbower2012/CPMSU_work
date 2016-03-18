#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>

using namespace std;

class vector{
		int length;
		double* info;
	public:
		vector();
		vector(int length);
		void alloc(int length);
		void del();
		~vector();

		vector operator+(const vector&);
};

int main(int argc, char* argv[]){}


vector::vector(){
	length = 0;
	info = NULL;
}
vector::vector(int length){
	info = new double[length];
	for(int i=0;i<length;i++){
		info[i] = 0.0;
	}
}
void vector::alloc(int length){
	info = new double[length];
	for(int i=0;i<length;i++){
		info[i]=0.0;
	}
}
void vector::del(){
	delete[] info;
}
vector vector::operator+(const vector& a){
	vector temp;
