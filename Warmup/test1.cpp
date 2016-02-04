#include<stdafx.h>
#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<string>
#include<windows.h>

using namespace std;

HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
COORD CursorPosition;

void gotoXY (int x, int y){
CursorPosition.X=x;
CursorPosition.Y=y;
SetConsoleCursorPosition(console,CursorPosition);
}

void gotoXY (int x, int y, string text){
CursorPosition.X=x;
CursorPosition.Y=y;
SetConsoleCursorPosition(console,CursorPosition);
cout << text;
}

void WaitKey(){
while(_kbhit()) _getch();
_getch();
while(_knhit()) _getch();
}

int main(int){
	int x[]={0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90};
	int y[]={0,7,11,14,11,19,9,9,9,10,10,11,11,10,9,9};
	string line(80,'_');
	for(int a=0; a< 14;a++){
		gotoXY(x[a],y[a],"$");
		gotoXY(x[a],16);
		cout << a+1;
	}
	gotoXY(0,15,line);
	gotoXY(0,25);
	WaitKey();
	return 0;
}
