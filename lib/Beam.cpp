#include <complex>
#include <fstream>
#include <iostream>
#include "Eigen/Core"
#include "Beam.h"
#include "Crystal.h"
#include "Constants.h"
#include "AuxFns.h"

Beam::Beam(){}

Beam::~Beam(){}

void Beam::set_E0(double EeV){
	E0 = EeV/toJ; //Convert eV to joules
	double lorentz;
	lorentz = 1. + E0/(rest_mass*c*c);
	mass = rest_mass*lorentz;

	lambda = h*c/sqrt(E0*(2*rest_mass*c*c + E0));
	lambda = lambda*1e10;
	std::cout<<"Electron wavelength in Angstrom: "<<lambda<<"\n";

}

void Beam::set_microscope_settings(double ac, double df, double cs){
	aconv = ac*0.001; //convergence angle in rad
	dF = df;
	Cs = cs;
}

/*************************************************
 * Zone axis is in real space, kt is in reciprocal space
 * So zone axis gets converted to reciprocal space before calculation of kin
 * *************************************************/
void Beam::calc_kin(const Eigen::Vector3d& kt){
	Eigen::Vector3d kza;
	kza << 0,0,1;

	k_in = kza + kt;
	k_in.normalize();
	k_in /= lambda;
	
	v_in = kza + kt;
	v_in.normalize();
	v_in *= speed;
}

double Beam::get_lambda() const{
	return lambda;
}

double Beam::get_aconv() const{
	return aconv;
}

double Beam::get_dF() const{
	return dF;
}

double Beam::get_E0() const{
	return E0;
}

double Beam::get_Cs() const{
	return Cs;
}

double Beam::get_v() const{
	return speed;
}

double Beam::get_mass() const{
	return mass;
}

void Beam::build_planewave(Eigen::VectorXcd& psi0){
	psi0(0) = 1;
}

void Beam::build_probe(const Eigen::Vector3d& R0, const Eigen::MatrixXd& g,Eigen::VectorXcd& psi0){
	int psize = psi0.rows();
	double tophat,chi;
	double Kmax = lambda/aconv;
	cout<<"Kmax:"<<Kmax<<"\n";
	for(int i=0;i<psize;i++){
		if(g.row(i).norm()<=Kmax){
			tophat = 1;
		}
		else
			tophat = 0;
		chi = pi*dF*lambda*g.row(i).norm()*g.row(i).norm()+0.5*pi*Cs*pow(lambda,3)*pow(g.row(i).norm(),4);
		chi = 0;
		psi0(i) = tophat*exp(-ri*chi)*exp(-2*pi*ri*(g.row(i).dot(R0)));
	}

}
