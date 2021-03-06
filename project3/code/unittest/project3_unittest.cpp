#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include "time.h"
#include "project3_library_unittest.h"

using namespace std;

ofstream ofile, ofile2;

int main(int argc, char* argv[]){
	string periname, aphname;
	stringstream number;
	char* outfilename;
	int i, j, k, steps, method, outputtype, bodies;
	int stepsperyear;
	int dec_or_inc, aphguess, periguess;
	double tmax, time, tolerance;
	clock_t start, finish;

	if(argc<7){
		cout << "Bad usage. Enter also 'outfilename method_choice output_type steps tmax tolerance' on same line." << endl;
		exit(1);
	}
	else{
		outfilename = argv[1];
		method = atoi(argv[2]);
		outputtype = atoi(argv[3]);
		steps = atoi(argv[4]);
		tmax = atof(argv[5]);
		if(outputtype==2||outputtype==3||outputtype==4){
			tmax = atoi(argv[5]);
		}
		tolerance = atof(argv[6]);
		
	}
	double** output, **aphelion, **perihelion;

	massivebody sun(1,3.771551748320805E-03, 1.938413234187417E-03, -1.625928558000791E-04,365*3.428095941073785E-08, 365*6.978886434512168E-06, 365*-9.372671992938156E-09);
	massivebody mercury(1.6601e-7,3.566221110752382E-01,-1.449153604767920E-01,-4.453344798939488E-02, 365*5.324168987021533E-03,365*2.725352157519689E-02, 365*1.737877336062238E-03);
//	massivebody sun(1,0,0,0,0,0,0);
//	massivebody mercury(1.6601e-7,0.3075, 0, 0, 0, 12.44, 0);
	massivebody venus(2.4478383e-6,4.456189829066335E-01,-5.759198980424926E-01,-3.358283335789598E-02,365*1.593371764525070E-02, 365*1.222110108236829E-02,365*-7.520174121955455E-04);
	massivebody earth(3.003489e-6,-9.906650404586314E-01,4.353612431574581E-02,-1.569714899841466E-04,365*-9.933032846586867E-04,365*-1.724423592380582E-02,365*2.880574493748607E-07);
	massivebody moon(1.230004e-8, -9.917828800119005E-01, 4.587463373520608E-02, -3.563694480647205E-04, 365*-1.530892843498182E-03, 365*-1.746686609878291E-02, 365*2.803200805250705E-05);
	massivebody mars(3.227151e-7, -1.394909885799664E+00, -7.759974369033253E-01, 1.786251125355763E-02, 365*7.324673346484570E-03, 365*-1.102624283521118E-02, 365*-4.109846883854566E-04);
	massivebody jupiter(9.5479194e-4,-5.321136962878863E+00, 1.055810040982731E+00, 1.146123452783326E-01, 365*-1.556597232697437E-03, 365*-7.046207863842619E-03, 365*6.409351264039102E-05);
	massivebody saturn(2.858860e-4,-3.332098484988519E+00, -9.442038663142483E+00, 2.967846224795282E-01, 365*4.954645896896637E-03, 365*-1.873255191158998E-03, 365*-1.647257228743735E-04);
	massivebody uranus(4.366244e-5, 1.876841090026478E+01, 6.826065082612249E+00, -2.177966843356428E-01, 365*-1.372937790621792E-03, 365*3.512867479374193E-03, 365*3.086802162915850E-05);
	massivebody neptune(5.151389e-5, 2.803900494548452E+01, -1.054870089186826E+01, -4.289565554838171E-01, 365*1.084650993604757E-03, 365*2.957157649376530E-03, 365*-8.562727609311126E-05);
	massivebody pluto(6.5812e-9, 8.773368933896196E+00, -3.186331328356860E+01, 8.718065633574812E-01, 365*3.100891963853092E-03, 365*1.939401372093854E-04, 365*-9.194995916567601E-04);
	massivebody roguestar(1,-500,0,0,50,0,0);

	massivesystem solarsystem;
	solarsystem.add(sun);
	solarsystem.add(mercury);
	solarsystem.add(venus);
	solarsystem.add(earth);
	solarsystem.add(moon);
	solarsystem.add(mars);
	solarsystem.add(jupiter);
	solarsystem.add(saturn);
	solarsystem.add(uranus);
	solarsystem.add(neptune);
	solarsystem.add(pluto);
//	solarsystem.add(roguestar);

	bodies = solarsystem.massivebody_count;

	matrix_alloc(output, steps+1, bodies*8+1);

	//method choice
	if(method==1){
		start = clock();
		solarsystem.RK4(output,steps,tmax, tolerance);
		finish = clock();
	}
	else if(method==2){
		start = clock();
		solarsystem.verlet(output,steps,tmax, tolerance);
		finish = clock();
	}
	else if(method==3){
		double* mass = new double[bodies];
		start = clock();
		solarsystem.initialize(output, mass, steps);
		verlet(output,mass,bodies,steps,tmax, tolerance);
		finish = clock();
		delete[] mass;
	}
	else if(method==4){
		double* mass = new double[bodies];
		start = clock();
		solarsystem.initialize(output, mass, steps);
		verlet_relcor(output,mass,bodies,steps,tmax);
		finish = clock();
		delete[] mass;
	}

	//outputfile choice
	ofile.open(outfilename);
	ofile.precision(8);
	if(outputtype==1){
		for(i=0;i<steps+1;i++){
			for(j=0;j<bodies*8+1;j++){
				ofile << setw(15) << output[i][j];
			}
			ofile << endl;
		}
	}
	else if(outputtype==2){
		stepsperyear=steps/tmax;
		for(i=0;i<tmax+1;i++){
			for(j=0;j<bodies;j++){
				ofile << setw(15) << output[i*stepsperyear][j*8+1];
				ofile << setw(15) << output[i*stepsperyear][j*8+2];
				ofile << setw(15) << output[i*stepsperyear][j*8+3];
			}
			ofile << endl;
		}
	}
	else if(outputtype==3){
		matrix_alloc(aphelion,bodies,2);
		matrix_alloc(perihelion,bodies,2);
		stepsperyear=steps/tmax;
		for(i=0;i<tmax+1;i++){
			for(j=0;j<bodies*8+1;j++){
				ofile << setw(15) << output[i*stepsperyear][j];
			}
			ofile << endl;
		}
		helionstates(output, aphelion, perihelion, bodies, steps);
		for(i=0;i<bodies;i++){
				cout << setw(15) << aphelion[i][0] << setw(15) << aphelion[i][1];
				cout << setw(30) << perihelion[i][0] << setw(15) << perihelion[i][1];
				cout << endl;
		}
		matrix_delete(aphelion,bodies);
		matrix_delete(perihelion,bodies);
	}
	else if(outputtype==4){
		stepsperyear=steps/tmax;
		for(i=0;i<tmax+1;i++){
			for(j=0;j<bodies*8+1;j++){
				ofile << setw(15) << output[i*stepsperyear][j];
			}
			ofile << endl;
		}
		for(k=0;k<bodies;k++){
			periname = "perihelion_body_";
			aphname = "aphelion_body_";
			aphguess = 1;
			periguess = aphguess;
			number << k;
			periname += number.str() + ".dat";
			aphname += number.str() + ".dat";
			number.str("");
			matrix_alloc(aphelion,aphguess,5);
			matrix_alloc(perihelion,periguess,5);
			dec_or_inc=1;
			helionstates_dynamic(aphelion, aphguess, perihelion, periguess, output,steps,k,20,dec_or_inc);
			
			ofile2.open(periname);
			ofile2.precision(8);
			for(i=0;i<periguess;i++){
				for(j=0;j<5;j++){
					ofile2 << setw(15) << perihelion[i][j];
				}
				ofile2 << endl;
			}
			ofile2.close();
			ofile2.open(aphname);
			ofile2.precision(8);
			for(i=0;i<aphguess;i++){
				for(j=0;j<5;j++){
					ofile2 << setw(15) << aphelion[i][j];
				}
				ofile2 << endl;
			}
			ofile2.close();

			matrix_delete(aphelion,aphguess);
			matrix_delete(perihelion,periguess);
		}
	}

	time = (finish - start)/((double) CLOCKS_PER_SEC);
	cout << "time " << time << " seconds" << endl;
	ofile.close();

	matrix_delete(output,steps+1);
	return 0;
}
