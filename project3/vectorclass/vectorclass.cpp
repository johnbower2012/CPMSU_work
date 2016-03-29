#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include "vectorclass.h"

using namespace std;

vector::vector(){
	length = 0;
	array = nullptr;
}
vector::vector(int size){
	length = size;
	array = new double[length];
	for(int i=0;i<length;i++){
		array[i] = 0.0;
	}
}
void vector::alloc(int size){
	length = size;
	array = new double[length];
	for(int i=0;i<length;i++){
		array[i]=0.0;
	}
}
void vector::del(){
	delete[] array;
	length = 0;
}
vector vector::operator+(const vector& a){
	int i;	
	vector temp(length);	
	for(i=0;i<length;i++){
		temp.array[i] = array[i] + a.array[i];
	}
	return temp; 
	temp.del();
}
vector vector::operator-(const vector& a){
	int i;
	vector temp(length);
	for(i=0;i<length;i++){
		temp.array[i] = array[i] - a.array[i];
	}
	return temp;
	temp.del();
}
double vector::operator*(const vector& a){
	int i;
	double temp = 0.0;
	for(i=0;i<length;i++){
		temp += array[i]*a.array[i];
	}
	return temp;
}
vector vector::operator/(const vector& a){
	int i, mod13, mod23;
	vector temp(3);
	for(i=0;i<3;i++){
		mod13 = (i+1)%3;
		mod23 = (i+2)%3;
		temp.array[i] = array[mod13]*a.array[mod23] - array[mod23]*a.array[mod13];
	}
	return temp;
	temp.del();
}

/*****************************
Problem functions/operations
******************************

vector operator=(const vector&);

vector vector::operator=(const vector& a){
	int i;
	vector temp(length);
	for(i=0;i<length;i++){
		temp.array[i] = a.array[i];
	}
	return temp;
	temp.del();
}

*****************************/


