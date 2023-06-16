#include <complex>
#include <vector>
#include "Eigen/Core"
#include "Constants.h"
#include "Beam.h"
#include "WaveFn.h"

WaveFn::WaveFn(vector<double> yg, vector<double> xg){

	xlength = xg.size(); ylength = yg.size();	
	psif = Eigen::ArrayXXcd::Zero(ylength,xlength);
	dpsix = Eigen::ArrayXXcd::Zero(ylength,xlength);
	dpsiy = Eigen::ArrayXXcd::Zero(ylength,xlength);
}

WaveFn::~WaveFn(){}

void WaveFn::set_thickness(double z){
	zf = z;
}

double WaveFn::get_zf(){
	return zf;
}

void WaveFn::build_psig(const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Eigen::VectorXcd& psi0){

	int g_size = gamma.rows();
	psig.resize(g_size);
	Eigen::VectorXcd alpha = Coeff.adjoint()*psi0;
	Eigen::VectorXcd gamma_new(g_size);

	for(int k=0;k<g_size;k++)
		gamma_new(k) = exp(pi*ri*gamma(k)*zf); 

	Eigen::MatrixXcd phase = gamma_new.asDiagonal();
	psig = Coeff*phase*alpha;

}

void WaveFn::build_TwoBeam_psig(const double xig, const double sg){
	
	psig.resize(2);
	double w = sg*xig;

	psig(0) = cos(pi*sqrt(1+w*w)*zf/xig)-ri*w/sqrt(1+w*w)*sin(pi*sqrt(1+w*w)*zf/xig);
	psig(1) = ri/sqrt(1+w*w)*sin(pi*sqrt(1+w*w)*zf/xig);

}

void WaveFn::build_EW(const Beam& beam, const Eigen::MatrixXd& g, const vector<double>& xg, const vector<double>& yg){

	int g_size = g.rows();

	/*Wave function at plane z*/
	for(int i=0;i<ylength;i++)
		for(int j=0;j<xlength;j++)
		 	for(int k=0;k<g_size;k++)
				psif(i,j) = psif(i,j) + psig(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*xg[j]+(beam.k_in(1)+g(k,1))*yg[i]+(beam.k_in(2)+g(k,2))*zf));

}

void WaveFn::build_BW(uint64_t bw, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Eigen::VectorXcd& psi0){ 

	int g_size = gamma.rows();
	psig.resize(g_size);
	Eigen::VectorXcd alpha = Coeff.adjoint()*psi0;
	Eigen::VectorXcd gamma_new(g_size);

	for(int k=0;k<g_size;k++)
		gamma_new(k) = exp(pi*ri*gamma(k)*zf); 

	Eigen::MatrixXcd phase = gamma_new.asDiagonal();
	for(int i=0;i<g_size;i++)
		psig(i) = Coeff(i,bw)*phase(bw,bw)*alpha(bw);
	
}

void WaveFn::build_vfield(const vector<double>& xg, const vector<double>& yg, const Eigen::MatrixXd& g, const Beam& beam) {

	int g_size = g.rows();

	for(int i=0;i<ylength;i++){
		for(int j=0;j<xlength;j++){
			for(int k=0;k<g_size;k++){
				dpsix(i,j) = dpsix(i,j) + psig(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*xg[j]+(beam.k_in(1)+g(k,1))*yg[i]+(beam.k_in(2)+g(k,2))*zf))*2.0*pi*ri*(beam.k_in(0)+g(k,0));
				dpsiy(i,j) = dpsiy(i,j) + psig(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*xg[j]+(beam.k_in(1)+g(k,1))*yg[i]+(beam.k_in(2)+g(k,2))*zf))*2.0*pi*ri*(beam.k_in(1)+g(k,1));
			}
		}
	}
	vfieldx = (hbar)/beam.get_mass()*imag(dpsix/psif)*(1e10);
	vfieldy = (hbar)/beam.get_mass()*imag(dpsiy/psif)*(1e10);
}

void WaveFn::build_qpot(double dx, double dy, const Beam& beam){
	Eigen::ArrayXXd R(ylength,xlength);
	Eigen::ArrayXXd grad2Rx(ylength,xlength);
	Eigen::ArrayXXd grad2Ry(ylength,xlength);
	Eigen::ArrayXXd grad2R(ylength,xlength);

	R = real(psif*psif.conjugate()).pow(0.5);

	/*grad^2R inside grid points (centered difference)*/
	for(int i=1;i<ylength-1;i++){
		for(int j=1;j<xlength-1;j++){
			grad2Rx(i,j) = (R(i,j-1) - 2*R(i,j) + R(i,j+1))/(dx*dx);
			grad2Ry(i,j) =(R(i-1,j) - 2*R(i,j) + R(i+1,j))/(dy*dy);
		}
	}

	/*grad^2R at boundaries*/
	for(int i=1;i<ylength-1;i++){
		grad2Rx(i,0) = (R(i,xlength-1) - 2*R(i,0) + R(i,1))/(dx*dx);
		grad2Rx(i,xlength-1) = (R(i,xlength-2) - 2*R(i,xlength-1) + R(i,0))/(dx*dx);
		grad2Ry(i,0) = (R(i-1,0) - 2*R(i,0) + R(i+1,0))/(dy*dy);
		grad2Ry(i,xlength-1) = (R(i-1,xlength-1) - 2*R(i,xlength-1) + R(i+1,xlength-1))/(dy*dy);
	}
	for(int j=1;j<xlength-1;j++){
		grad2Rx(0,j) = (R(0,j-1) - 2*R(0,j) + R(0,j+1))/(dx*dx);
		grad2Rx(ylength-1,j) = (R(ylength-1,j-1) - 2*R(ylength-1,j) + R(ylength-1,j+1))/(dx*dx);
		grad2Ry(0,j) = (R(ylength-1,j) - 2*R(0,j) + R(1,j))/(dy*dy);
		grad2Ry(ylength-1,j) = (R(ylength-2,j) - 2*R(ylength-1,j) + R(0,j))/(dy*dy);
	}

	/*Last is at the corners*/
	grad2Rx(0,0) = (R(0,xlength-1) - 2*R(0,0) + R(0,1))/(dx*dx);
	grad2Ry(0,0) = (R(ylength-1,0) - 2*R(0,0) + R(1,0))/(dy*dy);

	grad2Rx(0,xlength-1) = (R(0,xlength-2) - 2*R(0,xlength-1) + R(0,0))/(dx*dx);
	grad2Ry(0,xlength-1) = (R(ylength-1,xlength-1) - 2*R(0,xlength-1) + R(1,xlength-1))/(dy*dy);

	grad2Rx(ylength-1,0) = (R(ylength-1,xlength-1) - 2*R(ylength-1,0) + R(ylength-1,1))/(dx*dx);
	grad2Ry(ylength-1,0) = (R(ylength-2,0) - 2*R(ylength-1,0) + R(0,0))/(dy*dy);

	grad2Rx(ylength-1,xlength-1) = (R(ylength-1,xlength-2) - 2*R(ylength-1,xlength-1) + R(ylength-1,0))/(dx*dx);
	grad2Ry(ylength-1,xlength-1) = (R(ylength-2,xlength-1) - 2*R(ylength-1,xlength-1) + R(0,xlength-1))/(dy*dy);

	grad2R = grad2Rx + grad2Ry;
	Qpot = -(hbar*hbar)/(2*beam.get_mass())*R.inverse()*grad2R*1e20*toJ; //QP in eV

	/*Quantum force*/
	fqx.resize(ylength,xlength);
	fqy.resize(ylength,xlength);

	for(int i=0;i<ylength-1;i++){
		for(int j=0;j<xlength-1;j++){
			fqx(i,j) = -(Qpot(i,j+1) - Qpot(i,j))/dx;
			fqy(i,j) = -(Qpot(i+1,j) - Qpot(i,j))/dy;
		}
	}

	for(int i=0;i<ylength-1;i++){
		fqx(i,xlength-1) = -(Qpot(i,0) - Qpot(i,xlength-1))/dx;
		fqy(i,xlength-1) = -(Qpot(i+1,xlength-1) - Qpot(i,xlength-1))/dy;
	}

	for(int j=0;j<xlength-1;j++){
		fqx(ylength-1,j) = -(Qpot(ylength-1,j+1) - Qpot(ylength-1,j))/dx;
		fqy(ylength-1,j) = -(Qpot(0,j) - Qpot(ylength-1,j))/dy;
	}
	fqx(ylength-1,xlength-1) = -(Qpot(ylength-1,xlength-1) - Qpot(ylength-1,0))/dx;
	fqy(ylength-1,xlength-1) = -(Qpot(0,xlength-1) - Qpot(ylength-1,xlength-1))/dy;

}

/*Compute Q at a point, centered difference for second order derivative (R)*/
double WaveFn::Qcalc_pt(const Eigen::Vector3d& r, const double dr, const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Beam& beam, const Eigen::MatrixXd& g) const{
	int g_size = g.rows();
	Eigen::VectorXcd alpha = Coeff.adjoint()*psi0;
	Eigen::VectorXcd gamma_new(g_size);
	dcomp psi = 0, psi_m = 0, psi_p = 0;
	double R, Rm, Rp, grad2x, grad2y, grad2z;

	for(int i=0;i<g_size;i++)
		gamma_new(i) = exp(pi*ri*gamma(i)*r(2)); 

	Eigen::MatrixXcd phase = gamma_new.asDiagonal();
	Eigen::VectorXcd pg = Coeff*phase*alpha;

	/*gradx*/
	for(int k=0;k<g_size;k++){
		psi = psi + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*r(0) + (beam.k_in(1)+g(k,1))*r(1) + (beam.k_in(2)+g(k,2))*r(2)));	
		psi_m = psi_m + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*(r(0)-dr) + (beam.k_in(1)+g(k,1))*r(1) + (beam.k_in(2)+g(k,2))*r(2)));	
		psi_p = psi_p + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*(r(0)+dr) + (beam.k_in(1)+g(k,1))*r(1) + (beam.k_in(2)+g(k,2))*r(2)));	
	}

	R = sqrt(real(psi*conj(psi)));
	Rm = sqrt(real(psi_m*conj(psi_m)));
	Rp = sqrt(real(psi_p*conj(psi_p)));

	grad2x = (Rm - 2*R + Rp)/(dr*dr);

	/*grady*/
	for(int k=0;k<g_size;k++){
		psi = psi + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*r(0) + (beam.k_in(1)+g(k,1))*r(1) + (beam.k_in(2)+g(k,2))*r(2)));	
		psi_m = psi_m + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*r(0) + (beam.k_in(1)+g(k,1))*(r(1)-dr) + (beam.k_in(2)+g(k,2))*r(2)));	
		psi_p = psi_p + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*r(0) + (beam.k_in(1)+g(k,1))*(r(1)+dr) + (beam.k_in(2)+g(k,2))*r(2)));	
	}

	R = sqrt(real(psi*conj(psi)));
	Rm = sqrt(real(psi_m*conj(psi_m)));
	Rp = sqrt(real(psi_p*conj(psi_p)));

	grad2y = (Rm - 2*R + Rp)/(dr*dr);

	/*gradz*/
	for(int k=0;k<g_size;k++){
		psi = psi + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*r(0) + (beam.k_in(1)+g(k,1))*r(1) + (beam.k_in(2)+g(k,2))*r(2)));	
		psi_m = psi_m + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*r(0) + (beam.k_in(1)+g(k,1))*r(1) + (beam.k_in(2)+g(k,2))*(r(2)-dr)));	
		psi_p = psi_p + pg(k)*exp(2*pi*ri*((beam.k_in(0)+g(k,0))*r(0) + (beam.k_in(1)+g(k,1))*r(1) + (beam.k_in(2)+g(k,2))*(r(2)+dr)));	
	}

	R = sqrt(real(psi*conj(psi)));
	Rm = sqrt(real(psi_m*conj(psi_m)));
	Rp = sqrt(real(psi_p*conj(psi_p)));

	grad2z = (Rm - 2*R + Rp)/(dr*dr);

	double Q = -(hbar*hbar)/(2*beam.get_mass())*pow(R,-1)*(grad2x+grad2y+grad2z)*1e20*toJ;
	
	return Q;
}
/*******************************************************
 * Centered difference to calculate electrostatic force
 * at single point
 * *****************************************************/
void WaveFn::fq_CD(const Eigen::Vector3d& r, const double dr, Eigen::Vector3d& f, const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Beam& beam, const Eigen::MatrixXd& g) const{
	
	double Qpx, Qpy, Qpz, Qmx, Qmy, Qmz;
	Eigen::Vector3d trans;

	trans << 0.5, 0, 0;
	Qmx = Qcalc_pt(r - dr*trans, dr, psi0, gamma, Coeff, beam, g);

	trans << 0, 0.5, 0;
	Qmy = Qcalc_pt(r - dr*trans, dr, psi0, gamma, Coeff, beam, g);

	trans << 0, 0, 0.5;
	Qmz = Qcalc_pt(r - dr*trans, dr, psi0, gamma, Coeff, beam, g);

	trans << 0.5, 0, 0;
	Qpx = Qcalc_pt(r + dr*trans, dr, psi0, gamma, Coeff, beam, g);

	trans << 0, 0.5, 0;
	Qpy = Qcalc_pt(r + dr*trans, dr, psi0, gamma, Coeff, beam, g);

	trans << 0, 0, 0.5;
	Qpz = Qcalc_pt(r + dr*trans, dr, psi0, gamma, Coeff, beam, g);

	f << -(Qpx-Qmx)/dr, -(Qpy-Qmy)/dr, -(Qpz-Qmz)/dr; 
}
