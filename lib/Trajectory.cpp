#include <fstream>
#include <vector>
#include "Eigen/Core"
#include "Crystal.h"
#include "Trajectory.h"
#include "Constants.h"
#include "AuxFns.h"
#include "WaveFn.h"

using namespace std;

Trajectory::Trajectory() {}

Trajectory::~Trajectory(){}

void Trajectory::init_position(const double xg, const double yg){
	xcoord = xg;
	ycoord = yg;
	zcoord = 0.0;
	inside = true;
}

void Trajectory::init_velocity(const double vx_0, const double vy_0, const double vz_0){
	vx = vx_0;
	vy = vy_0;
	vz = vz_0;
}

/*Determine the nearest grid point to a point*/
int Trajectory::near_point(char coord,vector<double> xg){
	double ds;
	int n=xg.size();
	ds = abs(xg[0]-xg[1]);

	/*Determine box on grid where trajectory is located*/
	int n_min;
	if(coord == 'x')
		n_min = floor(xcoord/ds);
	else if(coord == 'y')
		n_min = floor(ycoord/ds);
	else{
		cout<<"Invalid coordinate\n";
		exit(0);
	}

	return n_min;
}

/* Second order Runge Kutta */
void Trajectory::RK2(const Eigen::VectorXcd& psig, const Beam& beam, const Eigen::MatrixXd& g, double dt){
	int g_size = g.rows();

	double xnew, ynew, znew;
	dcomp psi = 0, dpsix = 0, dpsiy = 0, dpsiz = 0;
	double kx1, ky1, kz1;
	double kx2, ky2, kz2;

	/*First coefficient k1 = dt*f*/
	/*Wave function*/
	for(int i=0;i<g_size;i++)
		psi = psi + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xcoord + (beam.k_in(1)+g(i,1))*ycoord +(beam.k_in(2)+g(i,2))*zcoord));

	/*Derivative of psi*/
	for(int i=0;i<g_size;i++){	
		dpsix = dpsix + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xcoord + (beam.k_in(1)+g(i,1))*ycoord + (beam.k_in(2)+g(i,2))*zcoord))*2.0*pi*ri*(beam.k_in(0)+g(i,0));
		dpsiy = dpsiy + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xcoord + (beam.k_in(1)+g(i,1))*ycoord + (beam.k_in(2)+g(i,2))*zcoord))*2.0*pi*ri*(beam.k_in(1)+g(i,1));
		dpsiz = dpsiz + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xcoord + (beam.k_in(1)+g(i,1))*ycoord + (beam.k_in(2)+g(i,2))*zcoord))*2.0*pi*ri*(beam.k_in(2)+g(i,2));
	}

	/*Velocity field*/
	kx1 = hbar/beam.get_mass()*imag(dpsix/psi)*1e10;
	ky1 = hbar/beam.get_mass()*imag(dpsiy/psi)*1e10;
	kz1 = hbar/beam.get_mass()*imag(dpsiz/psi)*1e10;

	/*Second coefficient k2 = dt*f(x+0.5*k1*dt)*/
	xnew = xcoord + 0.5*kx1*dt;
	ynew = ycoord + 0.5*ky1*dt;
	znew = zcoord + 0.5*kz1*dt;

	/*Wave function*/
	psi = 0; dpsix = 0; dpsiy = 0; dpsiz = 0;
	for(int i=0;i<g_size;i++)
		psi = psi + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xnew + (beam.k_in(1)+g(i,1))*ynew +(beam.k_in(2)+g(i,2))*znew));

	/*Derivative of psi*/
	for(int i=0;i<g_size;i++){	
		dpsix = dpsix + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xnew + (beam.k_in(1)+g(i,1))*ynew + (beam.k_in(2)+g(i,2))*znew))*2.0*pi*ri*(beam.k_in(0)+g(i,0));
		dpsiy = dpsiy + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xnew + (beam.k_in(1)+g(i,1))*ynew + (beam.k_in(2)+g(i,2))*znew))*2.0*pi*ri*(beam.k_in(1)+g(i,1));
		dpsiz = dpsiz + psig(i)*exp(2*pi*ri*((beam.k_in(0)+g(i,0))*xnew + (beam.k_in(1)+g(i,1))*ynew + (beam.k_in(2)+g(i,2))*znew))*2.0*pi*ri*(beam.k_in(2)+g(i,2));
	}

	kx2 = hbar/beam.get_mass()*imag(dpsix/psi)*1e10;
	ky2 = hbar/beam.get_mass()*imag(dpsiy/psi)*1e10;
	kz2 = hbar/beam.get_mass()*imag(dpsiz/psi)*1e10;

	xcoord += dt*kx2;
	ycoord += dt*ky2;
	zcoord += dt*kz2;
}

/*Overload for psig as an array*/
void Trajectory::RK2(const Eigen::ArrayXXcd& psig, const Beam& beam, const Eigen::MatrixXd& g, const Eigen::MatrixXd& q, double dt){
	int g_size = g.rows();

	double xnew, ynew, znew;
	dcomp psi = 0, dpsix = 0, dpsiy = 0, dpsiz = 0;
	double kx1, ky1, kz1;
	double kx2, ky2, kz2;

	/*First coefficient k1 = dt*f*/
	/*Wave function*/
	for(int j=0;j<q.rows();j++)
		for(int i=0;i<g_size;i++)
			psi = psi + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xcoord + (beam.k_in(1)+g(i,1)+q(j,1))*ycoord +(beam.k_in(2)+g(i,2)+q(j,2))*zcoord));

	/*Derivative of psi*/
	for(int j=0;j<q.rows();j++){
		for(int i=0;i<g_size;i++){	
			dpsix = dpsix + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xcoord + (beam.k_in(1)+g(i,1)+q(j,1))*ycoord + (beam.k_in(2)+g(i,2)+q(j,2))*zcoord))*2.0*pi*ri*(beam.k_in(0)+g(i,0)+q(j,0));
			dpsiy = dpsiy + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xcoord + (beam.k_in(1)+g(i,1)+q(j,1))*ycoord + (beam.k_in(2)+g(i,2)+q(j,2))*zcoord))*2.0*pi*ri*(beam.k_in(1)+g(i,1)+q(j,1));
			dpsiz = dpsiz + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xcoord + (beam.k_in(1)+g(i,1)+q(j,1))*ycoord + (beam.k_in(2)+g(i,2)+q(j,2))*zcoord))*2.0*pi*ri*(beam.k_in(2)+g(i,2)+q(j,2));
		}
	}

	/*Velocity field*/
	kx1 = hbar/beam.get_mass()*imag(dpsix/psi)*1e10;
	ky1 = hbar/beam.get_mass()*imag(dpsiy/psi)*1e10;
	kz1 = hbar/beam.get_mass()*imag(dpsiz/psi)*1e10;

	/*Second coefficient k2 = dt*f(x+0.5*k1*dt)*/
	xnew = xcoord + 0.5*kx1*dt;
	ynew = ycoord + 0.5*ky1*dt;
	znew = zcoord + 0.5*kz1*dt;

	/*Wave function*/
	psi = 0; dpsix = 0; dpsiy = 0; dpsiz = 0;
	for(int j=0;j<q.rows();j++)
		for(int i=0;i<g_size;i++)
			psi = psi + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xnew + (beam.k_in(1)+g(i,1)+q(j,1))*ynew +(beam.k_in(2)+g(i,2)+q(j,2))*znew));

	/*Derivative of psi*/
	for(int j=0;j<q.rows();j++){
		for(int i=0;i<g_size;i++){	
			dpsix = dpsix + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xnew + (beam.k_in(1)+g(i,1)+q(j,1))*ynew + (beam.k_in(2)+g(i,2)+q(j,2))*znew))*2.0*pi*ri*(beam.k_in(0)+g(i,0)+q(j,0));
			dpsiy = dpsiy + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xnew + (beam.k_in(1)+g(i,1)+q(j,1))*ynew + (beam.k_in(2)+g(i,2)+q(j,2))*znew))*2.0*pi*ri*(beam.k_in(1)+g(i,1)+q(j,1));
			dpsiz = dpsiz + psig(i,j)*exp(2*pi*ri*((beam.k_in(0)+g(i,0)+q(j,0))*xnew + (beam.k_in(1)+g(i,1)+q(j,1))*ynew + (beam.k_in(2)+g(i,2)+q(j,2))*znew))*2.0*pi*ri*(beam.k_in(2)+g(i,2)+q(j,2));
		}
	}

	kx2 = hbar/beam.get_mass()*imag(dpsix/psi)*1e10;
	ky2 = hbar/beam.get_mass()*imag(dpsiy/psi)*1e10;
	kz2 = hbar/beam.get_mass()*imag(dpsiz/psi)*1e10;

	xcoord += dt*kx2;
	ycoord += dt*ky2;
	zcoord += dt*kz2;
}

/*Runge Kutta for f=ma (second derivative ODE)*/
void Trajectory::RK2_fe(const Crystal& crystal, const Eigen::MatrixXd& g, bool absorb, double dt, const Beam& beam){

	double k1v[3], k1x[3], k2v[3], k2x[3];
	Eigen::Vector3d r,dr,force, v;
	
	r << xcoord, ycoord, zcoord;
	v << vx, vy, vz;
	dr << 0.001, 0.001, 0.001;

	crystal.fe_CD(r,dr,force,g,absorb, beam);

	for(int i=0;i<3;i++){
		k1v[i] = force(i)/beam.get_mass()*dt;
		k1x[i] = v(i)*dt;
	}

	r << xcoord + 0.5*k1x[0]*dt, ycoord + 0.5*k1x[1]*dt, zcoord + 0.5*k1x[2]*dt;
	crystal.fe_CD(r,dr,force,g,absorb, beam);

	for(int i=0;i<3;i++){
		k2v[i] = force(i)/beam.get_mass();
		k2x[i] = v(i)+0.5*k1v[i];
	}
	

	xcoord += k2x[0]*dt;
	ycoord += k2x[1]*dt;
	zcoord += k2x[2]*dt;

	v(0) += k2v[0]*dt;
	v(1) += k2v[1]*dt;
	v(2) += k2v[2]*dt;
}

void Trajectory::RK2_fq(const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Beam& beam, const Eigen::MatrixXd& g, double dt, const WaveFn& wf){
	double k1v[3], k1x[3], k2v[3], k2x[3];
	Eigen::Vector3d r,force, v;
	
	r << xcoord, ycoord, zcoord;
	v << vx, vy, vz;
	double dr = 0.01;

	wf.fq_CD(r,dr,force, psi0, gamma, Coeff, beam, g);

	for(int i=0;i<3;i++){
		k1v[i] = force(i)/beam.get_mass()*dt;
		k1x[i] = v(i)*dt;
	}

	r << xcoord + 0.5*k1x[0]*dt, ycoord + 0.5*k1x[1]*dt, zcoord + 0.5*k1x[2]*dt;
	wf.fq_CD(r,dr,force, psi0, gamma, Coeff, beam, g);

	for(int i=0;i<3;i++){
		k2v[i] = force(i)/beam.get_mass();
		k2x[i] = v(i)+0.5*k1v[i];
	}
		

	xcoord += k2x[0]*dt;
	ycoord += k2x[1]*dt;
	zcoord += k2x[2]*dt;

	v(0) += k2v[0]*dt;
	v(1) += k2v[1]*dt;
	v(2) += k2v[2]*dt;
}
