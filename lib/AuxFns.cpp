#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include <ctime>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include <fftw3.h>
#include "Constants.h"
#include "AuxFns.h"
#include "Crystal.h"
#include "Beam.h"
#include "Trajectory.h"

//Insertion sort
void AuxFns::sort(Eigen::Vector3d& a){
	signed int temp;
	for(int j=1;j<a.size();j++){ 
		temp=a(j);
		for(int i=j;i>0;i--){
		        if(a(i-1)<temp){
				a(i)=temp;
				break;
			}
			a(i)=a(i-1);
			a(i-1) = temp;
		}
	}

}

/*Convert between from direct space to reciprocal space
 * *******************************************/
void AuxFns::convertKSpace(const Eigen::Vector3d& x, Eigen::Vector3d& x_new, const Crystal& crystal){

	double ax = crystal.get_aax();
	double ay = crystal.get_aay();
	double az = crystal.get_aaz();
	Eigen::Matrix3d mt;

	mt << ax*ax, 0, 0,
	    	0, ay*ay, 0,
		0, 0, az*az;

	x_new = x.transpose()*mt;
}
/************************************************************************
 * Build realspace grid
 * aa = {aax,aay} and npoints = {npointsx,npointsy}
 * **********************************************************************/
double AuxFns::RSgrid(double aa[2], int npoints[2], vector <double>& xg, vector<double>& yg){
	double dx, dy;
	dx = aa[0]/npoints[0];
	dy = aa[1]/npoints[1];
	
	for(int i=0;i<npoints[0];i++)
		xg.push_back(i*dx);
	for(int j=0;j<npoints[1];j++)
		yg.push_back(j*dy);
	cout<<"x axis from "<<xg[0]<<" to "<<xg[npoints[0]-1]<<" Angstrom\n";
	cout<<"y axis from "<<yg[0]<<" to "<<yg[npoints[1]-1]<<" Angstrom\n";

	return npoints[0]/(2*aa[0]);
}

/*********************************************************************
 * Build kspace grid
 * ******************************************************************/
void AuxFns::KSgrid(double aa[2], int npoints[2], vector<double>& kxg, vector<double>& kyg,vector<double>& ckx, vector<double>& cky){
	double dx, dy;
	ckx.resize(npoints[0]); cky.resize(npoints[1]);
	kxg.resize(npoints[0]); kyg.resize(npoints[1]);
	dx = aa[0]/npoints[0];
	dy = aa[1]/npoints[1];

	double dkx = 1.0/dx;
	double dky = 1.0/dy;

	ckx[0] = 0; cky[0] = 0;

	/*First derivative*/
	for(int i=1;i<npoints[0]/2;i++){
		ckx[i] = i*dkx;
		ckx[npoints[0]-i] = -ckx[i];
	}
	for(int i=1;i<npoints[1]/2;i++){
		cky[i] = i*dky;
		cky[npoints[1]-i] = -cky[i];
	}

	/*Grid k-values*/
	for(int i=0;i<npoints[0]/2;i++){
		kxg[i] = ckx[npoints[0]/2+i];
		kxg[i+npoints[0]/2] = ckx[i];
	}
	for(int i=0;i<npoints[1]/2;i++){
		kyg[i] = cky[npoints[1]/2+i];
		kyg[i+npoints[1]/2] = cky[i];
	}

}

/******************************************************************
 * Compute list of reflexions for Bloch wave computation in zone axis orientation
 * Taken from De Graef
 * **************************************************************/
void AuxFns::gZA(const Beam& beam, const Crystal& crystal, double gmax, double Climit[2], bool absorb, Eigen::MatrixXd& g, Eigen::MatrixXd& gpert, vector <double>& sg){
	g.resize(1,3);
	gpert.resize(1,3);
	g << 0,0,0;
	sg.push_back(0);

	Eigen::Vector3d gtemp;
	double sgtemp;
	complex<double> Vg,Ug;
	bool first_time = true;
	ofstream data_o("Input/hkldata.txt");
	data_o<<"h\tk\tl\t|g|\tsg\treal(Ug)\timag(Ug)\t|Ug|\n";
	cout<<crystal.ctype<<"\n";

	int max_hkl = 30;
	for(int i=-max_hkl;i<=max_hkl;i++){
		for(int j=-max_hkl;j<=max_hkl;j++){
			for(int k=-max_hkl;k<=max_hkl;k++){
				if(i == 0 && j == 0 && k == 0)
					continue;
				gtemp = i*crystal.b.col(0)+j*crystal.b.col(1)+k*crystal.b.col(2);
				/*Excitation error*/
				sgtemp = sgcalc(gtemp,beam.k_in);
				/*Fourier coefficient*/
				Vg = crystal.calcVg(gtemp, beam);
				Ug = 2*beam.get_mass()*charge/(h*h*1e20)*Vg;
					
				/*************************/
				if(crystal.allowReflection(i,j,k) && gtemp.norm()<gmax){
					if(abs(sgtemp)<=beam.get_lambda()*abs(Ug)*Climit[0]){
						g.conservativeResize(g.rows()+1,g.cols());
						g.row(g.rows()-1) = gtemp;
						sg.push_back(sgtemp);
						if(absorb == true){
							dcomp V_imag = crystal.calcVgAbsorb(gtemp,beam);
							Vg = Vg + ri*V_imag;
							Ug = 2*beam.get_mass()*charge/(h*h*1e20)*Vg;
						}
						data_o<<i<<"\t"<<j<<"\t"<<k<<"\t"<<gtemp.norm()<<"\t"<<sgtemp<<"\t"<<real(Ug)<<"\t"<<imag(Ug)<<"\t"<<abs(Ug)<<"\n";
					}
					else if(abs(sgtemp)<=beam.get_lambda()*abs(Ug)*Climit[1]){
						if(first_time == true){
							gpert.row(gpert.rows()-1) = gtemp;
							first_time = false;
						}
						else{
							gpert.conservativeResize(gpert.rows()+1,gpert.cols());
							gpert.row(gpert.rows()-1) = gtemp;
						}
					}
				}
			}
		}
	}
	data_o.close();
}

/************************************************************
 * Calculate extinction distance of point g
 * **********************************************************/
double AuxFns::calcExtinctionDistance(const Eigen::Vector3d& g, const Crystal& crystal, const Beam& beam){

	double xi_g;
	dcomp Vg, Ug;
	Vg = crystal.calcVg(g, beam);
	Ug = 2*beam.get_mass()*abs(charge)*Vg/(h*h)*1e-20;
	cout<<"Ug = "<<abs(Vg)<<"\n";

	Eigen::Vector3d un;
	un << 0,0,1;
	Eigen::Vector3d kg;
	kg = beam.k_in+g;

	double calpha = kg.dot(un)/(kg.norm()*un.norm());
	xi_g = abs(Ug)/(kg.norm()*calpha);

	xi_g = 1.0/xi_g;

	return xi_g;
}

/******************************************************************
 * Compute list of reflexions for Bloch wave computation in 
 * systematic row orientation (not sure this works...)
 * **************************************************************/
void AuxFns::gSR(const Beam& beam, const Crystal& crystal, int maxg, bool absorb, Eigen::MatrixXd& g, vector <double>& sg, const Eigen::Vector3d& gSR){
	g.resize(1,3);
	g << 0,0,0;
	sg.push_back(0);

	Eigen::Vector3d gtemp;
	double sgtemp;
	complex<double> Vg,Ug;
	bool first_time = true;
	ofstream data_o("Input/hkldata.txt");
	data_o<<"h\tk\tl\t|g|\tsg\treal(Ug)\timag(Ug)\t|Ug|\n";
	cout<<crystal.ctype<<"\n";

	for(int i=-maxg;i<=maxg;i++){
		if(i == 0)
			continue;
		gtemp = i*(gSR(0)*crystal.b.col(0)+gSR(1)*crystal.b.col(1)+gSR(2)*crystal.b.col(2));
		/*Excitation error*/
		sgtemp = sgcalc(gtemp,beam.k_in);
		/*Fourier coefficient*/
		Vg = crystal.calcVg(gtemp, beam);
		Ug = 2*beam.get_mass()*abs(charge)/(h*h*1e20)*Vg;
			
		/*************************/
		g.conservativeResize(g.rows()+1,g.cols());
		g.row(g.rows()-1) = gtemp;
		sg.push_back(sgtemp);
		if(absorb == true){
			dcomp V_imag = crystal.calcVgAbsorb(gtemp,beam);
			Vg = Vg + ri*V_imag;
			Ug = 2*beam.get_mass()*abs(charge)/(h*h*1e20)*Vg;
		}
		data_o<<i*gSR(0)<<"\t"<<i*gSR(1)<<"\t"<<i*gSR(2)<<"\t"<<gtemp.norm()<<"\t"<<sgtemp<<"\t"<<real(Ug)<<"\t"<<imag(Ug)<<"\t"<<abs(Ug)<<"\n";
	}
	data_o.close();
}

/*****************************************************************
 * Compute excitation error for reflexion g and incident wave vector k
 * ***************************************************************/
double AuxFns::sgcalc(const Eigen::Vector3d& g, const Eigen::Vector3d& k){
	double sg, calpha;
	Eigen::Vector3d un;
	un << 0,0,1;
	Eigen::Vector3d kg;
	kg = k+g;

	calpha = kg.dot(un)/(kg.norm()*un.norm());

	sg = -g.dot(2*k+g)/(2*(k+g).norm()*calpha);

	return sg;
}

/*Functions for writing to file*/
void AuxFns::write_psif(string iterator, const Eigen::ArrayXXcd& psif){
	string filename = "Results/I0_z"+iterator+".txt";
	ofstream output(filename.c_str());
	output<<(psif*psif.conjugate()).real();
	output.close();

	filename = "Results/Amplitude_z"+iterator+".txt";
	ofstream realpsi_o(filename.c_str());
	realpsi_o<<psif.real();
	realpsi_o.close();

	filename = "Results/Phase_z"+iterator+".txt";
	ofstream imagpsi_o(filename.c_str());
	imagpsi_o<<psif.imag();
	imagpsi_o.close();
}

void AuxFns::write_psig(string iterator, const Eigen::MatrixXd& g, const Eigen::VectorXcd& psig){ 
	int g_size = g.rows();
	string filename = "Results/Psig_z"+iterator+".txt";
	ofstream psig_o(filename.c_str());
	for(int i=0;i<g_size;i++)
		psig_o<<g(i,0)<<"\t"<<g(i,1)<<"\t"<<g(i,2)<<"\t"<<real(psig(i))<<"\t"<<imag(psig(i))<<"\n";
	psig_o.close();
}

void AuxFns::write_fpsi(string iterator, const Eigen::ArrayXXcd& fpsi){
	ofstream fpsi_o("Results/Diff"+iterator+".txt");
	fpsi_o<<real(fpsi*fpsi.conjugate());
	fpsi_o.close();
}

void AuxFns::write_vfield(string iterator, const Eigen::ArrayXXd& vfieldx,const Eigen::ArrayXXd& vfieldy){
	ofstream vfield_ox("Results/Vfield"+iterator+"_X.txt");
	ofstream vfield_oy("Results/Vfield"+iterator+"_Y.txt");
	vfield_ox<<vfieldx;
	vfield_oy<<vfieldy;
	vfield_ox.close();
	vfield_oy.close();
}

void AuxFns::write_qp(string iterator, const Eigen::ArrayXXd& qpot){
	ofstream qpot_o("Results/Qpotential"+iterator+".txt");
	qpot_o<<qpot;
	qpot_o.close();
}

void AuxFns::write_traj(ofstream &out, Trajectory traj){
	out<<traj.xcoord<<"\t"<<traj.ycoord<<"\t"<<traj.zcoord<<"\n";
}
void AuxFns::write_potential(string iterator, const Eigen::ArrayXXcd& pot){
	string filename = "Results/Potential_z"+iterator+".txt";
	ofstream output(filename.c_str());
	output<<pot.real();
	output.close();
}	

void AuxFns::write_fe(string iterator, const Eigen::ArrayXXd& fex, const Eigen::ArrayXXd& fey){
	string filenamex = "Results/fex_z"+iterator+".txt";
	string filenamey = "Results/fey_z"+iterator+".txt";
	ofstream outputx(filenamex.c_str());
	outputx<<fex;
	outputx.close();
	ofstream outputy(filenamey.c_str());
	outputy<<fey.real();
	outputy.close();
}

void AuxFns::write_fq(string iterator, const Eigen::ArrayXXd& fqx, const Eigen::ArrayXXd& fqy){
	string filenamex = "Results/fqx_z"+iterator+".txt";
	string filenamey = "Results/fqy_z"+iterator+".txt";
	ofstream outputx(filenamex.c_str());
	outputx<<fqx;
	outputx.close();
	ofstream outputy(filenamey.c_str());
	outputy<<fqy;
	outputy.close();
}

/******************************************************
 * Random number generator for Gaussian distribution
 * ****************************************************/
double AuxFns::gasdev(){
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	time_t timebuffer;
	const double seconds = time(&timebuffer);
	int idum = -(seconds);
			
        if (idum < 0) iset=0;
	if (iset == 0) {
		do {
			v1=2.0*ran1()-1;
			v2=2.0*ran1()-1;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

/*****************************************************************
 * Random number generator for uniform deviates (Press et al.)   *
 * ***************************************************************/
double AuxFns::ran1(){
	const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
	const int NDIV=(1+(IM-1)/NTAB);
	const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
	static int iy=0;
	static vector<int> iv(NTAB);
	int j,k;
	double temp;
	time_t timebuffer;
	const double seconds = time(&timebuffer);
	int idum = -(seconds);

        if (idum <= 0 || !iy) {
	        if (-idum < 1) idum=1;
	        else idum = -idum;
	        for (j=NTAB+7;j>=0;j--) {
		        k=idum/IQ;
		        idum=IA*(idum-k*IQ)-IR*k;
		        if (idum < 0) idum += IM;
		        if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
