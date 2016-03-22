#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>

using namespace std;

void array_alloc(double*& a, int length);
void array_delete(double*& a);
void matrix_alloc(double**&, int, int);
void matrix_delete(double*&, int);
void rel_coord(double**&, double*, int);


class massivebody{
	public:
		double mass;
		double position[3];
		double velocity[3];
		massivebody();
		massivebody(double, double, double, double, double, double, double);
		massivebody(massivebody, double, double, double, double, double, double, double);
};

class massivesystem{
	public:
		int massivebody_count;
		massivebody* system;
		double mass_total;
		double composition[3];
		double comvelocity[3];
		massivesystem();
		~massivesystem();
		void add(massivebody guy);
		void relpositions(massivesystem&, double**&);
		void forces(massivesystem&, double*&, const double**&);
		void initialize(massivesystem&, double**&, double*&, int);
		void print();
};

void verlet(double**&, double*&, int, int, double);
void verlet(massivesystem&, double**&, int, double);
void RK4(massivesystem&, double**&, int, double);

ofstream ofile;

int main(int argc, char* argv[]){
	char* outfilename;
	int steps, i, j, method;
	double tmax, bodies;
	if(argc<5){
		cout << "Bad usage. Enter also 'outfilename method_choice steps tmax' on same line." << endl;
		exit(1);
	}
	else{
		outfilename = argv[1];
		method = atoi(argv[2]);
		steps = atoi(argv[3]);
		tmax = atof(argv[4]);
	}
	double** output;

	massivebody sun(1,0,0,0,0,0,0);
	massivebody mercury(1.6596e-7,0,0.3871,0,1.590*2.0*M_PI,0,0);
	massivebody venus(2.4471e-6,0.723,0,0,0,-1.176*2.0*M_PI,0);
	massivebody earth(3.00251e-6,1,0,0,0,2.0*M_PI,0);
	massivebody moon(earth, 3.69415e-8, 2.57e-3, 0, 0, 0, 0, 0.034*2.0*M_PI);
	massivebody mars(3.2262e-7,0,0,1.5236,0.808*2.0*M_PI,0,0);
	massivebody jupiter(9.95424e-4,5.4,0,0,0,0,-0.438*2.0*M_PI);
	massivebody saturn(2.8573e-4,0,0,9.537,0,-0.324*2.0*M_PI,0);
	massivebody uranus(4.3638e-5,19.2,0,0,0,0.228*2.0*M_PI,0);
	massivebody neptune(5.1488e-5,0,30.070,0,0,0,0.182*2.0*M_PI);
	massivebody pluto(6.5812e-9,0,0,39.52,0,0.157*2.0*M_PI,0);
	massivebody roguestar(1,-500,0,0,50,0,0);

	massivesystem solarsystem;
	solarsystem.add(sun);
/*	solarsystem.add(mercury);
	solarsystem.add(venus); */
	solarsystem.add(earth);
/*	solarsystem.add(moon);
	solarsystem.add(mars);
	solarsystem.add(jupiter);
	solarsystem.add(saturn);
	solarsystem.add(neptune);
	solarsystem.add(uranus);
	solarsystem.add(pluto);*/
	//solarsystem.add(roguestar);

	if(method==1){
		verlet(solarsystem,output,steps,tmax);
	}
	else if(method==2){
		RK4(solarsystem,output,steps,tmax);
	}
	else if(method==3){
		double* mass = new double[solarsystem.massivebody_count];
		solarsystem.initialize(solarsystem, output, mass, steps);
		verlet(output,mass,solarsystem.massivebody_count,steps,tmax);
		delete[] mass;
	}

	ofile.open(outfilename);
	bodies = solarsystem.massivebody_count;
	for(i=0;i<steps+1;i++){
		for(j=0;j<bodies*6+1;j++){
			ofile << setw(15) << output[i][j];
		}
		ofile << endl;
	}
	ofile.close();

	for(i=0;i<steps+1;i++){
		delete[] output[i];
	}
	delete[] output;
	return 0;
}


void array_alloc(double*& a, int length){
	a = new double[length];
}
void array_delete(double*& a){
	delete[] a;
}
void matrix_alloc(double**& a, int rows, int columns){
	a = new double*[rows];
	for(int i=0;i<rows;i++){
		a[i] = new double[columns];
	}
}
void matrix_delete(double**& a, int rows){
	for(int i=0;i<rows;i++){
		delete a[i];
	}
	delete[] a;
}

massivebody::massivebody(){
	mass = 0.0;
	position[0] = 0.0;
	position[1] = 0.0;
	position[2] = 0.0;
	velocity[0] = 0.0;
	velocity[1] = 0.0;
	velocity[2] = 0.0;
}
massivebody::massivebody(double Mass, double x, double y, double z, double vx, double vy, double vz){
	mass = Mass;
	position[0] = x;
	position[1] = y;
	position[2] = z;
	velocity[0] = vx;
	velocity[1] = vy;
	velocity[2] = vz;
}
massivebody::massivebody(massivebody primary, double Mass, double x, double y, double z, double vx, double vy, double vz){
	mass = Mass;
	position[0] = primary.position[0] + x;
	position[1] = primary.position[1] + y;
	position[2] = primary.position[2] + z;
	velocity[0] = primary.velocity[0] + vx;
	velocity[1] = primary.velocity[1] + vy;
	velocity[2] = primary.velocity[2] + vz;
}
massivesystem::massivesystem(){
	massivebody_count = 0;
	system = nullptr;
	mass_total = 0;
	for(int i=0;i<3;i++){
		composition[i]=0.0;
		comvelocity[i]=0.0;
	}
}
massivesystem::~massivesystem(){
	delete[] system;
}
void massivesystem::add(massivebody guy){
	double mass_old = mass_total;
	massivebody* temp = new massivebody[massivebody_count+1];
	for(int i=0;i<massivebody_count;i++){
		temp[i].mass=system[i].mass;
		temp[i].position[0]=system[i].position[0];
		temp[i].position[1]=system[i].position[1];
		temp[i].position[2]=system[i].position[2];
		temp[i].velocity[0]=system[i].velocity[0];
		temp[i].velocity[1]=system[i].velocity[1];
		temp[i].velocity[2]=system[i].velocity[2];
	}
	temp[massivebody_count].mass=guy.mass;
	temp[massivebody_count].position[0]=guy.position[0];
	temp[massivebody_count].position[1]=guy.position[1];
	temp[massivebody_count].position[2]=guy.position[2];
	temp[massivebody_count].velocity[0]=guy.velocity[0];
	temp[massivebody_count].velocity[1]=guy.velocity[1];
	temp[massivebody_count].velocity[2]=guy.velocity[2];

	massivebody_count++;

	delete[] system;
	mass_total += guy.mass;
	for(int i=0;i<3;i++){
		composition[i] = (composition[i]*mass_old + guy.position[i]*guy.mass)/mass_total;
		comvelocity[i] = (comvelocity[i]*mass_old + guy.velocity[i]*guy.mass)/mass_total;
	}
	system = temp;
}
void massivesystem::print(){
	cout << "mass=" << setw(15) << mass_total << endl;
	for(int i=0;i<3;i++){
		cout << "pos[" << i<< "]=" << setw(15) << composition[i] << endl;
	}
	for(int i=0;i<3;i++){
		cout << "vel[" << i<< "]=" << setw(15) << comvelocity[i] << endl;
	}
}

void rel_coord(double**& relcoord, double** output, int i, int bodies){
	int m, n;
	for(m=0;m<bodies;m++){
		for(n=0;n<bodies;n++){
			if(m!=n){
				relcoord[m][n*4] = output[i][m*6+1] - output[i][n*6+1];
				relcoord[m][n*4+1] = output[i][m*6+2] - output[i][n*6+2];
				relcoord[m][n*4+2] = output[i][m*6+3] - output[i][n*6+3];
				relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
			}
		}
	}	
}
void massivesystem::relpositions(massivesystem& sol, double**& relcoord){
	int i, j;	
	int bodies = sol.massivebody_count;
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies;j++){
			if(i!=j){
				relcoord[i][j*4] = sol.system[i].position[0] - sol.system[j].position[0];
				relcoord[i][j*4+1] = sol.system[i].position[1] - sol.system[j].position[1];
				relcoord[i][j*4+2] = sol.system[i].position[2] - sol.system[j].position[2];
				relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
			}
		}
	}
}
void massivesystem::forces(massivesystem& sol, double*& force, const double**& relcoord){
	int j, k;
	int bodies = sol.massivebody_count;
	double rcubed, fourpisq;
	fourpisq = 4.0*M_PI*M_PI;
	for(j=0;j<bodies;j++){
		force[j*3] = 0.0;
		force[j*3+1] = 0.0;
		force[j*3+2] = 0.0;
	}
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
				force[j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
				force[j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[j*3] *= fourpisq;
		force[j*3+1] *= fourpisq;
		force[j*3+2] *= fourpisq;
	}
}
void massivesystem::initialize(massivesystem& sol, double**& output, double*& mass, int steps){
	int i, j;
	int bodies = sol.massivebody_count;
	output = new double*[steps+1];
	mass = new double[bodies];
	for(int i=0;i<steps+1;i++){
		output[i] = new double[bodies*6+1];
	}
	for(j=0;j<bodies;j++){
		mass[j] = sol.system[j].mass;
		output[0][j*6+1] = sol.system[j].position[0];
		output[0][j*6+2] = sol.system[j].position[1];
		output[0][j*6+3] = sol.system[j].position[2];
		output[0][j*6+4] = sol.system[j].velocity[0];
		output[0][j*6+5] = sol.system[j].velocity[1];
		output[0][j*6+6] = sol.system[j].velocity[2];
	}
	for(i=1;i<steps+1;i++){
		for(j=0;j<bodies;j++){
			output[i][j*6+1] = 0;
			output[i][j*6+2] = 0;
			output[i][j*6+3] = 0;
			output[i][j*6+4] = 0;
			output[i][j*6+5] = 0;
			output[i][j*6+6] = 0;
		}
	}
}
void verlet(double**& output, double*& mass, int bodies, int steps, double tmax){
	int i, j, k, m, n;
	double h, halfh, halfhsq, rcubed, fourpisq;
	h = (tmax-0.0)/((double) (steps));
	halfh = 0.5*h;
	halfhsq = halfh*h;	
	fourpisq = 4.0*M_PI*M_PI;
	//relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
	//technically double counting--consider revising
	double** relcoord = new double*[bodies];
	for(i=0;i<bodies;i++){
		relcoord[i] = new double[bodies*4];
	}
	//force provides fx,fy,fz for each body at step i and i+1
	double** force = new double*[2];
	for(i=0;i<2;i++){
		force[i] = new double[bodies*3];
	}
	//initialize all matrices to zero
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*4;j++){
			relcoord[i][j]=0.0;
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<bodies*3;j++){
			force[i][j]=0.0;
		}
	}
	//reinitialize relcoord to initial conditions
	for(m=0;m<bodies;m++){
		for(n=0;n<bodies;n++){
			if(m!=n){
				relcoord[m][n*4] = output[0][m*6+1] - output[0][n*6+1];
				relcoord[m][n*4+1] = output[0][m*6+2] - output[0][n*6+2];
				relcoord[m][n*4+2] = output[0][m*6+3] - output[0][n*6+3];
				relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
			}
		}
	}
	//initialize forces
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[0][j*3] += -mass[k]*relcoord[j][k*4]/rcubed;
				force[0][j*3+1] += -mass[k]*relcoord[j][k*4+1]/rcubed;
				force[0][j*3+2] += -mass[k]*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[0][j*3] *= fourpisq;
		force[0][j*3+1] *= fourpisq;
		force[0][j*3+2] *= fourpisq;
	}
	//begin solution for-loop
	for(i=0;i<steps;i++){
		//Set next time value
		output[i+1][0] = output[i][0] + h;
		//calculate x,y,z		
		for(j=0;j<bodies;j++){
			output[i+1][j*6+1] = output[i][j*6+1] + h*output[i][j*6+4] + halfhsq*force[0][j*3];
			output[i+1][j*6+2] = output[i][j*6+2] + h*output[i][j*6+5] + halfhsq*force[0][j*3+1];
			output[i+1][j*6+3] = output[i][j*6+3] + h*output[i][j*6+6] + halfhsq*force[0][j*3+2];
		}
		//new relcoord
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*6+1] - output[i+1][n*6+1];
					relcoord[m][n*4+1] = output[i+1][m*6+2] - output[i+1][n*6+2];
					relcoord[m][n*4+2] = output[i+1][m*6+3] - output[i+1][n*6+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//calculate force_i+1
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -mass[k]*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -mass[k]*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -mass[k]*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
		}
		//calculate vx, vy, vz
		for(j=0;j<bodies;j++){
			output[i+1][j*6+4] = output[i][j*6+4] + halfh*(force[1][j*3] + force[0][j*3]);
			output[i+1][j*6+5] = output[i][j*6+5] + halfh*(force[1][j*3+1] + force[0][j*3+1]);
			output[i+1][j*6+6] = output[i][j*6+6] + halfh*(force[1][j*3+2] + force[0][j*3+2]);
		}
		//force_i for next loop is force_i+1 for this loop. No need for double computation.
		for(j=0;j<bodies;j++){
			force[0][j*3] = force[1][j*3];
			force[0][j*3+1] = force[1][j*3+1];
			force[0][j*3+2] = force[1][j*3+2];
		}
	}
	//delete dross--keep output
	for(i=0;i<bodies;i++){
		delete[] relcoord[i];
	}
	delete[] relcoord;
	for(i=0;i<2;i++){
		delete[] force[i];
	}
	delete[] force;

}
void verlet(massivesystem& sol, double**& output, int steps, double tmax){
	int i, j, k, m, n, bodies;
	bodies = sol.massivebody_count;
	double h, halfh, halfhsq, rcubed, fourpisq;
	h = (tmax-0.0)/((double) (steps));
	halfh = 0.5*h;
	halfhsq = halfh*h;	
	fourpisq = 4.0*M_PI*M_PI;
	//output to provide t and x, y, z, vx, vy, vz for each body at each step i
	output = new double*[steps+1];	
	for(int i=0;i<steps+1;i++){
		output[i] = new double[bodies*6+1];
	}
	//relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
	//technically double counting--consider revising
	double** relcoord = new double*[bodies];
	for(i=0;i<bodies;i++){
		relcoord[i] = new double[bodies*4];
	}
	//force provides fx,fy,fz for each body at step i and i+1
	double** force = new double*[2];
	for(i=0;i<2;i++){
		force[i] = new double[bodies*3];
	}
	//initialize all matrices to zero
	for(i=0;i<steps+1;i++){
		for(j=0;j<bodies*6+1;j++){
			output[i][j]=0.0;
		}
	}
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*4;j++){
			relcoord[i][j]=0.0;
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<bodies*3;j++){
			force[i][j]=0.0;
		}
	}

	//reinitialize relcoord to initial conditions
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies;j++){
			if(i!=j){
				relcoord[i][j*4] = sol.system[i].position[0] - sol.system[j].position[0];
				relcoord[i][j*4+1] = sol.system[i].position[1] - sol.system[j].position[1];
				relcoord[i][j*4+2] = sol.system[i].position[2] - sol.system[j].position[2];
				relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
			}
		}
	}
	//initialize forces
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[0][j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
				force[0][j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
				force[0][j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[0][j*3] *= fourpisq;
		force[0][j*3+1] *= fourpisq;
		force[0][j*3+2] *= fourpisq;
	}
	//initialize initial conditions for output
	for(j=0;j<bodies;j++){
		output[0][j*6+1] = sol.system[j].position[0];
		output[0][j*6+2] = sol.system[j].position[1];
		output[0][j*6+3] = sol.system[j].position[2];
		output[0][j*6+4] = sol.system[j].velocity[0];
		output[0][j*6+5] = sol.system[j].velocity[1];
		output[0][j*6+6] = sol.system[j].velocity[2];
	}

	//begin solution for-loop
	for(i=0;i<steps;i++){
		//Set next time value
		output[i+1][0] = output[i][0] + h;
		//calculate x,y,z		
		for(j=0;j<bodies;j++){
			output[i+1][j*6+1] = output[i][j*6+1] + h*output[i][j*6+4] + halfhsq*force[0][j*3];
			output[i+1][j*6+2] = output[i][j*6+2] + h*output[i][j*6+5] + halfhsq*force[0][j*3+1];
			output[i+1][j*6+3] = output[i][j*6+3] + h*output[i][j*6+6] + halfhsq*force[0][j*3+2];
		}
		//new relcoord
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*6+1] - output[i+1][n*6+1];
					relcoord[m][n*4+1] = output[i+1][m*6+2] - output[i+1][n*6+2];
					relcoord[m][n*4+2] = output[i+1][m*6+3] - output[i+1][n*6+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//calculate force_i+1
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
		}
		//calculate vx, vy, vz
		for(j=0;j<bodies;j++){
			output[i+1][j*6+4] = output[i][j*6+4] + halfh*(force[1][j*3] + force[0][j*3]);
			output[i+1][j*6+5] = output[i][j*6+5] + halfh*(force[1][j*3+1] + force[0][j*3+1]);
			output[i+1][j*6+6] = output[i][j*6+6] + halfh*(force[1][j*3+2] + force[0][j*3+2]);
		}
		//force_i for next loop is force_i+1 for this loop. No need for double computation.
		for(j=0;j<bodies;j++){
			force[0][j*3] = force[1][j*3];
			force[0][j*3+1] = force[1][j*3+1];
			force[0][j*3+2] = force[1][j*3+2];
		}
	}
	//delete dross--keep output
	for(i=0;i<bodies;i++){
		delete[] relcoord[i];
	}
	delete[] relcoord;
	for(i=0;i<2;i++){
		delete[] force[i];
	}
	delete[] force;

}
void RK4(massivesystem& sol, double**& output, int steps, double tmax){
	int i, j, k, m, n, bodies;
	bodies = sol.massivebody_count;
	double h, halfh, halfhsq, rcubed, fourpisq;
	double* k1, *k2, *k3, *k4;
	h = (tmax-0.0)/((double) (steps));
	halfh = 0.5*h;
	halfhsq = halfh*h;	
	fourpisq = 4.0*M_PI*M_PI;
	//placeholders arrays
	k1 = new double[bodies*6];
	k2 = new double[bodies*6];
	k3 = new double[bodies*6];
	k4 = new double[bodies*6];
	for(i=0;i<bodies*6;i++){
		k1[i] = 0;
		k2[i] = 0;
		k3[i] = 0;
		k4[i] = 0;
	}			
	//output to provide t and x, y, z, vx, vy, vz for each body at each step i
	output = new double*[steps+1];	
	for(int i=0;i<steps+1;i++){
		output[i] = new double[bodies*6+1];
	}
	//relcoord to provide relative coordinates, x, y, z, r, of body i with respect to body j
	//technically double counting--consider revising
	double** relcoord = new double*[bodies];
	for(i=0;i<bodies;i++){
		relcoord[i] = new double[bodies*4];
	}
	//force provides fx,fy,fz for each body at step i and i+1
	double** force = new double*[2];
	for(i=0;i<2;i++){
		force[i] = new double[bodies*3];
	}
	//initialize all matrices to zero
	for(i=0;i<steps+1;i++){
		for(j=0;j<bodies*6+1;j++){
			output[i][j]=0.0;
		}
	}
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies*4;j++){
			relcoord[i][j]=0.0;
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<bodies*3;j++){
			force[i][j]=0.0;
		}
	}

	//reinitialize relcoord to initial conditions	
	for(i=0;i<bodies;i++){
		for(j=0;j<bodies;j++){
			if(i!=j){
				relcoord[i][j*4] = sol.system[i].position[0] - sol.system[j].position[0];
				relcoord[i][j*4+1] = sol.system[i].position[1] - sol.system[j].position[1];
				relcoord[i][j*4+2] = sol.system[i].position[2] - sol.system[j].position[2];
				relcoord[i][j*4+3] = sqrt(relcoord[i][j*4]*relcoord[i][j*4] + relcoord[i][j*4+1]*relcoord[i][j*4+1]+relcoord[i][j*4+2]*relcoord[i][j*4+2]);
			}
		}
	}
	//initialize forces
	for(j=0;j<bodies;j++){
		for(k=0;k<bodies;k++){
			if(j!=k){
				rcubed = pow(relcoord[j][k*4+3],3.0);
				force[0][j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
				force[0][j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
				force[0][j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
			}
		}
		force[0][j*3] *= fourpisq;
		force[0][j*3+1] *= fourpisq;
		force[0][j*3+2] *= fourpisq;
	}
	//initialize initial conditions for output
	for(j=0;j<bodies;j++){
		output[0][j*6+1] = sol.system[j].position[0];
		output[0][j*6+2] = sol.system[j].position[1];
		output[0][j*6+3] = sol.system[j].position[2];
		output[0][j*6+4] = sol.system[j].velocity[0];
		output[0][j*6+5] = sol.system[j].velocity[1];
		output[0][j*6+6] = sol.system[j].velocity[2];
	}

	//Begin RK4 loop
	for(i=0;i<steps;i++){
		//calculate first positions [1+1/2] (k1)
		for(j=0;j<bodies;j++){
			output[i+1][j*6+1] = output[i][j*6+1] + halfh*output[i][j*6+4];
			output[i+1][j*6+2] = output[i][j*6+2] + halfh*output[i][j*6+5];
			output[i+1][j*6+3] = output[i][j*6+3] + halfh*output[i][j*6+6];
			cout << setw(15) << output[i+1][j*6+1];
			cout << setw(15) << output[i+1][j*6+2];
			cout << setw(15) << output[i+1][j*6+3];
		}
cout << endl;
		//new relative positions at [i+1/2] (k2)
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*6+1] - output[i+1][n*6+1];
					relcoord[m][n*4+1] = output[i+1][m*6+2] - output[i+1][n*6+2];
					relcoord[m][n*4+2] = output[i+1][m*6+3] - output[i+1][n*6+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//new forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
			k2[j*6+3] = force[1][j*3];
			k2[j*6+4] = force[1][j*3+1];
			k2[j*6+5] = force[1][j*3+2];
		}
		//k2 using i+1/2 positions
		for(j=0;j<bodies;j++){
			k2[j*6] = output[i][j*6+4] + halfh*force[1][j*3];
			k2[j*6+1] = output[i][j*6+5] + halfh*force[1][j*3+1];
			k2[j*6+2] = output[i][j*6+6] + halfh*force[1][j*3+2];
		}
		//new temp positions at i+1/2 with k2
		for(j=0;j<bodies;j++){
			output[i+1][j*6+1] = output[i][j*6+1] + halfh*k2[j*6];
			output[i+1][j*6+2] = output[i][j*6+2] + halfh*k2[j*6+1];
			output[i+1][j*6+3] = output[i][j*6+3] + halfh*k2[j*6+2];
			cout << setw(15) << output[i+1][j*6+1];
			cout << setw(15) << output[i+1][j*6+2];
			cout << setw(15) << output[i+1][j*6+3];
		}
cout << endl;
		//new relative positions at [i+1/2] (k3)
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*6+1] - output[i+1][n*6+1];
					relcoord[m][n*4+1] = output[i+1][m*6+2] - output[i+1][n*6+2];
					relcoord[m][n*4+2] = output[i+1][m*6+3] - output[i+1][n*6+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//new forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
			k3[j*6+3] = force[1][j*3];
			k3[j*6+4] = force[1][j*3+1];
			k3[j*6+5] = force[1][j*3+2];
		}
		//k3 using i+1/2 positions
		for(j=0;j<bodies;j++){
			k3[j*6] = output[i][j*6+4] + halfh*force[1][j*3];
			k3[j*6+1] = output[i][j*6+5] + halfh*force[1][j*3+1];
			k3[j*6+2] = output[i][j*6+6] + halfh*force[1][j*3+2];
		}
		//new temp positions at i+1 k3
		for(j=0;j<bodies;j++){
			output[i+1][j*6+1] = output[i][j*6+1] + h*k3[j*6];
			output[i+1][j*6+2] = output[i][j*6+2] + h*k3[j*6+1];
			output[i+1][j*6+3] = output[i][j*6+3] + h*k3[j*6+2];
			cout << setw(15) << output[i+1][j*6+1];
			cout << setw(15) << output[i+1][j*6+2];
			cout << setw(15) << output[i+1][j*6+3];
		}
cout << endl;
		//new relative positions at [i+1] (k4)
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*6+1] - output[i+1][n*6+1];
					relcoord[m][n*4+1] = output[i+1][m*6+2] - output[i+1][n*6+2];
					relcoord[m][n*4+2] = output[i+1][m*6+3] - output[i+1][n*6+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//new forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
			k4[j*6+3] = force[1][j*3];
			k4[j*6+4] = force[1][j*3+1];
			k4[j*6+5] = force[1][j*3+2];
		}
		//k4 using i+1/2 positions
		for(j=0;j<bodies;j++){
			output[i+1][j*6+4] = output[i][j*6+4] + halfh*force[1][j*3];
			output[i+1][j*6+5] = output[i][j*6+5] + halfh*force[1][j*3+1];
			output[i+1][j*6+6] = output[i][j*6+6] + halfh*force[1][j*3+2];
		}
		//final positions at i+1
		for(j=0;j<bodies;j++){
			output[i+1][j*6+1] = output[i][j*6+1] + (h/6.0)*(output[i][j*6+4] + 2.0*k2[j*6] + 2.0*k3[j*6] + output[i+1][j*6+4]);
			output[i+1][j*6+2] = output[i][j*6+2] + (h/6.0)*(output[i][j*6+5] + 2.0*k2[j*6+1] + 2.0*k3[j*6+1] + output[i+1][j*6+5]);
			output[i+1][j*6+3] = output[i][j*6+3] + (h/6.0)*(output[i][j*6+6] + 2.0*k2[j*6+2] + 2.0*k3[j*6+2] + output[i+1][j*6+6]);
			cout << setw(15) << output[i+1][j*6+1];
			cout << setw(15) << output[i+1][j*6+2];
			cout << setw(15) << output[i+1][j*6+3];
		}
cout << endl;		
		//final relative positions at i+1
		for(m=0;m<bodies;m++){
			for(n=0;n<bodies;n++){
				if(m!=n){
					relcoord[m][n*4] = output[i+1][m*6+1] - output[i+1][n*6+1];
					relcoord[m][n*4+1] = output[i+1][m*6+2] - output[i+1][n*6+2];
					relcoord[m][n*4+2] = output[i+1][m*6+3] - output[i+1][n*6+3];
					relcoord[m][n*4+3] = sqrt(relcoord[m][n*4]*relcoord[m][n*4] + relcoord[m][n*4+1]*relcoord[m][n*4+1]+relcoord[m][n*4+2]*relcoord[m][n*4+2]);
				}
			}
		}
		//final forces with relative positions
		for(j=0;j<bodies*3;j++){
			force[1][j] = 0.0;
		}
		for(j=0;j<bodies;j++){
			for(k=0;k<bodies;k++){
				if(j!=k){
					rcubed = pow(relcoord[j][k*4+3],3.0);
					force[1][j*3] += -sol.system[k].mass*relcoord[j][k*4]/rcubed;
					force[1][j*3+1] += -sol.system[k].mass*relcoord[j][k*4+1]/rcubed;
					force[1][j*3+2] += -sol.system[k].mass*relcoord[j][k*4+2]/rcubed;
				}
			}
			force[1][j*3] *= fourpisq;
			force[1][j*3+1] *= fourpisq;
			force[1][j*3+2] *= fourpisq;
		}
		//final velocities at i+1
		for(j=0;j<bodies;j++){
			output[i+1][j*6+4] = output[i][j*6+4] + (h/6.0)*(force[0][j*3] + 2.0*k2[j*6+3] + 2.0*k3[j*6+3] + force[1][j*3]);
			output[i+1][j*6+5] = output[i][j*6+5] + (h/6.0)*(force[0][j*3+1] + 2.0*k2[j*6+4] + 2.0*k3[j*6+4] + force[1][j*3+1]);
			output[i+1][j*6+6] = output[i][j*6+6] + (h/6.0)*(force[0][j*3+2] + 2.0*k2[j*6+5] + 2.0*k3[j*6+5] + force[1][j*3+2]);
		}
		for(j=0;j<bodies;j++){
			force[0][j*3] = force[1][j*3];
			force[0][j*3+1] = force[1][j*3+1];
			force[0][j*3+2] = force[1][j*3+2];
		}

	}

	//delete dross--keep output
	for(i=0;i<bodies;i++){
		delete[] relcoord[i];
	}
	delete[] relcoord;
	for(i=0;i<2;i++){
		delete[] force[i];
	}
	delete[] force;
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
}

/*	double* x, *y, *r, *vx, *vy, *v;
	double x0, y0, rcubedinv, vx0, vy0, tmax, tmin;
	double h, halfh, hsq, halfhsq, pi, pisq, twopisqh, fourpisqh, twopisqhsq, vareu, varvv1, varvv2;
	double forcei, forcei1, halfx, halfy, k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y;
	char* outfilename;
	int n, i, j, k, method;
	pi = acos(-1);

	if(argc<5){
		cout << "Bad usage. Please enter also 'outfilename method step_count tmax." << endl;
		exit(1);
	}
	else{
		outfilename = argv[1];
		method = atoi(argv[2]);
		n = atoi(argv[3]);
		tmax = atof(argv[4]);
		tmin = 0.0;
		h = (tmax-tmin)/((double) (n));
		hsq = h*h;
		halfh=0.5*h;
	}

	twopisqh = 2.0*pi*pi*h;
	fourpisqh = 2.0*twopisqh;
	twopisqhsq = 2.0*pi*pi*h*h;
	array_alloc(x,n+1);
	array_alloc(y,n+1);
	array_alloc(r,n+1);
	array_alloc(vx,n+1);
	array_alloc(vy,n+1);
	array_alloc(v,n+1);

	x[0] = 1.0;
	y[0] = 0.0;
	r[0] = sqrt(x[0]*x[0] + y[0]*y[0]);
	vx[0] = 0.0;
	vy[0] = 2.0*pi;
	v[0] = sqrt(vx[0]*vx[0] + vy[0]*vy[0]);	

	if(method==0){
	//Euler
		for(i=0;i<n;i++){
			x[i+1] = x[i] + vx[i]*h;
			y[i+1] = y[i] + vy[i]*h;
			r[i+1] = sqrt(x[i+1]*x[i+1] + y[i+1]*y[i+1]);
			vareu = fourpisqh*pow(r[i+1], -3.0);
			vx[i+1] = vx[i] - vareu*x[i];
			vy[i+1] = vy[i] - vareu*y[i];
			v[i+1] = sqrt(vx[i+1]*vx[i+1] + vy[i+1]*vy[i+1]);
		}
	}
	else if(method==1){
	//Velocity Verlet
		rcubedinv = pow(r[0], -3.0);
		for(i=0;i<n;i++){
			forcei = twopisqh*rcubedinv;
			x[i+1] = x[i] + vx[i]*h - h*forcei*x[i];
			y[i+1] = y[i] + vy[i]*h - h*forcei*y[i];
			r[i+1] = sqrt(x[i+1]*x[i+1] + y[i+1]*y[i+1]);
			rcubedinv = pow(r[i+1],-3.0);
			forcei1 = twopisqh*rcubedinv;
			vx[i+1] = vx[i] - forcei1*x[i+1] - forcei*x[i];
			vy[i+1] = vy[i] - forcei1*y[i+1] - forcei*y[i];
			v[i+1] = sqrt(vx[i+1]*vx[i+1]+vy[i+1]*vy[i+1]);
		}
	}	
	else if(method==2){
		for(i=0;i<n;i++){
			x[i+1]=x[i]+halfh*vx[i];
			y[i+1]=y[i]+halfh*vy[i];
			rcubedinv = pow(x[i+1]*x[i+1] + y[i+1]*y[i+1], -1.5);
			forcei = twopisqh*rcubedinv;
			k2x = vx[i] - forcei*x[i+1];
			k2y = vy[i] - forcei*y[i+1];
			x[i+1] = x[i] + halfh*k2x;
			y[i+1] = y[i] + halfh*k2y;
			rcubedinv = pow(x[i+1]*x[i+1] + y[i+1]*y[i+1], -1.5);
			forcei = twopisqh*rcubedinv;
			k3x = vx[i] - forcei*x[i+1];
			k3y = vy[i] - forcei*y[i+1];
			x[i+1] = x[i] + h*k3x;
			y[i+1] = y[i] + h*k3y;
			rcubedinv = pow(x[i+1]*x[i+1] + y[i+1]*y[i+1], -1.5);
			forcei = fourpisqh*rcubedinv;
			vx[i+1] = vx[i] - forcei*x[i+1];
			vy[i+1] = vy[i] - forcei*y[i+1];
			x[i+1] = x[i] + (h/6.0)*(vx[i] + 2.0*k2x + 2.0*k3x + vx[i+1]);
			y[i+1] = y[i] + (h/6.0)*(vy[i] + 2.0*k2y + 2.0*k3y + vy[i+1]);
			r[i+1] = sqrt(x[i+1]*x[i+1] + y[i+1]*y[i+1]);		
			v[i+1] = sqrt(vx[i+1]*vx[i+1] + vy[i+1]*vy[i+1]);
		}
	}
			

	ofile.open(outfilename);
	ofile << "#t_i" << setw(15) << "x" << setw(15) << "y" << setw(15) << "vx" << setw(15) << "vy" << setw(15) << "r" << setw(15) << "v" << endl;
	for(i=0;i<n+1;i++){
		ofile << i*h << setw(15) << x[i] << setw(15) << y[i] << setw(15) << vx[i] << setw(15) << vy[i] << setw(15) << r[i] << setw(15) << v[i] << endl;
	}
	ofile.close();

	array_delete(x);
	array_delete(y);
	array_delete(r);
	array_delete(vx);
	array_delete(vy);
	array_delete(v);
*/
