#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include <ctime>
#include <array>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Core"
#include <map>
#include "Crystal.h"
#include "AuxFns.h"
#include "Constants.h"
#include <boost/math/special_functions/expint.hpp>

Crystal::Crystal(){}
Crystal::~Crystal(){}

void Crystal::read_cfg(std::ifstream& datafile){

    //Read data from file in to strings, map used for lines with '='
    std::map<std::string, std::string> datamap;
    std::vector<std::string> elementData;
    std::string line;
    if(datafile.is_open()){
        while(getline(datafile, line)){
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), 
			    line.end());
            std::istringstream is_line(line);
            std::string key;
            if(line[0] == '#' || line.empty())
                continue;
            if(line.find("=") != std::string::npos){
                if(std::getline(is_line, key, '=')){
                    std::string value;
                    if(std::getline(is_line, value))
                        datamap.insert(std::pair<std::string,std::string>(key,value));
                }
            }
            else{
                std::string value;
                std::getline(is_line, value);
                elementData.push_back(value);
            }
        }
    }
    else{
        std::cerr << "Couldn't open material config file\n";
	exit(1);
    }

    //Assign datamap to values
    ctype = datamap["CrystalType"];
    unitCellvx << stod(datamap["H0(1,1)"]), stod(datamap["H0(1,2)"]), stod(datamap["H0(1,3)"]);
    unitCellvy << stod(datamap["H0(2,1)"]), stod(datamap["H0(2,2)"]), stod(datamap["H0(2,3)"]);
    unitCellvz << stod(datamap["H0(3,1)"]), stod(datamap["H0(3,2)"]), stod(datamap["H0(3,3)"]);
    nElements = stoi(datamap["NumberOfElements"]);
    uint64_t nAtoms = stoi(datamap["NumberOfAtoms"]);

    //Assign element data
    for(uint64_t i=0; i<nElements;i++){
        for(uint64_t st=0;st<41;st++)
            if(!elementData[0+2*i].compare(elts[st]))
                Z.push_back(st);

        DWfactors.push_back(stod(elementData[1 + 2*i]));
    }

    //Build fractional crystal and assign types
    std::string delimiter = ",";
    std::string strOrig, temp;
    for(uint64_t i=elementData.size() - nAtoms;i<elementData.size();i++){
        s.conservativeResize(s.rows()+1, 4);
        strOrig = elementData[i];

        for(uint64_t j=0;j<4;j++){
            temp = strOrig.substr(0,strOrig.find(delimiter));
            s(s.rows()-1,j) = stod(temp);
            strOrig.erase(0, strOrig.find(delimiter) + delimiter.length());
        }
    }
    for(uint64_t i=0; i<s.rows();i++)
        atomTypes.push_back(s(i,3)-1);

    s.conservativeResize(s.rows(),3);
}

/*Lattice parameters of super cell*/
double Crystal::set_aax(){
    return unitCellvx.norm();
}
double Crystal::set_aay(){
    return unitCellvy.norm();
}
double Crystal::set_aaz(){
    return unitCellvz.norm();
}
double Crystal::get_aax() const{
	return aax;
}
double Crystal::get_aay() const{
	return aay;
}
double Crystal::get_aaz() const{
	return aaz;
}

/*****************************************************
*Build complete crystal and reciprocal lattice vectors
******************************************************/
void Crystal::buildCrystal(){
	Eigen::Matrix3d G;
	G <<unitCellvx, unitCellvy, unitCellvz;
	s = s*G;
	omega_inv = 1./(unitCellvx.dot(unitCellvy.cross(unitCellvz)));
	aax = set_aax();
	aay = set_aay();
	aaz = set_aaz();
	
	/*Reciprocal vectors*/
	Eigen::Vector3d b1, b2, b3;
	b1 = unitCellvy.cross(unitCellvz)/(unitCellvx.dot(unitCellvy.cross(unitCellvz)));
	b2 = unitCellvz.cross(unitCellvx)/(unitCellvx.dot(unitCellvy.cross(unitCellvz)));
	b3 = unitCellvx.cross(unitCellvy)/(unitCellvx.dot(unitCellvy.cross(unitCellvz)));
	b.resize(3,3);
	b << b1,b2,b3;
}

/***************************************************
 * Compute real part of fourier coefficient of potential reflexion g
 * *************************************************/
dcomp Crystal::calcVg(const Eigen::Vector3d& gvector, const Beam& beam) const{

	dcomp Vg = 0;
	for(int atom=0;atom<s.rows();atom++){
		double fg = 0;
		for(int l=0;l<3;l++)
		 	fg = fg + scatF_A[Z[atomTypes[atom]]][l]/(gvector.norm()*gvector.norm() + scatF_B[Z[atomTypes[atom]]][l]) + scatF_C[Z[atomTypes[atom]]][l] * exp(-scatF_D[Z[atomTypes[atom]]][l]*gvector.norm()*gvector.norm());
		Vg = Vg + fg*exp(-2*pi*ri*gvector.dot(s.row(atom)));
	}
	Vg = Vg*h*h*omega_inv/(2*pi*beam.get_mass()*abs(charge))*1e20;

	return Vg;
}
/**************************************************
 * Compute absorption part of fourier coefficient of potential
 * ************************************************/
dcomp Crystal::calcVgAbsorb(const Eigen::Vector3d& gvector, const Beam& beam) const{
	double I1ij, I2ij;
	double A[6], B[6];
	double fg = 0;
	dcomp VgAbsorb = 0;
	double gn = gvector.norm();
	double DW;
	
	/*For each atom*/
	for(int atom=0;atom<s.rows();atom++){
		
		fg = 0;
		DW = DWfactors[atomTypes[atom]];
		/*Parameters first*/
		A[0] = 0.02395*Z[atomTypes[atom]]/(3*(1 + WK_V[Z[atomTypes[atom]]]));
		A[1] = A[0]; A[2] = A[0];
		A[3] = WK_V[Z[atomTypes[atom]]]*A[0];
		A[4] = A[3]; A[5] = A[3];
		for(uint64_t i=0;i<6;i++)
			B[i] = WK_B[Z[atomTypes[atom]]][i];
		
		/*Now cross over all fitting parameters*/
		for(uint64_t i=0;i<6;i++){
			for(uint64_t j=0;j<6;j++){

				if(gn == 0){
					I1ij = pi*(B[i]*log((B[i] + B[j])/B[i]) + B[j]*log((B[i] + B[j])/B[j]));
					I2ij = pi*((B[i] + 2*DW)*log((B[i] + B[j] + 2*DW)/(B[i] + 2*DW)) + B[j]*log((B[i] + B[j] + 2*DW)/(B[i] + 2*DW)) + 2*DW*log(2*DW/(B[j] + 2*DW)));
				}
				else{
					I1ij = (pi/(gn*gn)) * (2*exp(1) + log(B[i]*gn*gn) + log(B[j]*gn*gn) - 2*boost::math::expint(-(B[i]*B[j])/(B[i]+B[j])*gn*gn) + exp(-B[i]*gn*gn) * (boost::math::expint(B[i]*B[i]/(B[i] + B[j])*gn*gn)) - boost::math::expint(B[i]*gn*gn) + exp(-B[j]*gn*gn) * (boost::math::expint((B[j]*B[j])/(B[i]+B[j])*gn*gn) - boost::math::expint(B[j]*gn*gn)));
					I2ij = (pi/(gn*gn)) * (2*boost::math::expint(-DW*(B[i]+DW)/(B[i]+2*DW)*gn*gn) + 2*boost::math::expint(-DW*(B[j]+DW)/(B[j]+2*DW)*gn*gn) - 2*boost::math::expint(-(B[i]+DW)*(B[j]+DW)/(B[i]+B[j]+2*DW)*gn*gn) - 2*boost::math::expint(-0.5*DW*gn*gn) + exp(-DW*gn*gn) * (2*boost::math::expint(0.5*DW*gn*gn) - boost::math::expint(DW*DW/(B[i]+2*DW)*gn*gn) - boost::math::expint(DW*DW/(B[j]+2*DW)*gn*gn)) + exp(-(B[i]+DW)*gn*gn) * (boost::math::expint(pow(B[i]+DW,2)/(B[i]+B[j]+2*DW)*gn*gn) - boost::math::expint(pow(B[i]+DW,2)/(B[i]+2*DW)*gn*gn)) + exp(-(B[j]+DW)*gn*gn) * (boost::math::expint(pow(B[j]+DW,2)/(B[i]+B[j]+2*DW)*gn*gn) - boost::math::expint(pow(B[j]+DW,2)/(B[j]+2*DW)*gn*gn)));
				}
				fg = fg + A[i]*A[j]*(exp(-DW*gn*gn)*I1ij + I2ij);
			}
		}
		fg = fg*(1./beam.k_in.norm());
		VgAbsorb = VgAbsorb + exp(-ri*gvector.dot(s.row(atom)))*exp(-DW*gn*gn)*fg;
	}
	VgAbsorb = VgAbsorb*4.*pi*hbar*hbar*omega_inv/(2*beam.get_mass())*1e20;
	
	return VgAbsorb;
	
}
/****************************************
 * Determine whether structure factor of given reflection is non-zero
 * ****************************************/
bool Crystal::allowReflection(int h, int k, int l) const{
	if(!ctype.compare("FCC")){
		if((h%2 == 0 && k%2 == 0 && l%2 == 0) || (abs(h)%2 == 1 && abs(k)%2 == 1 && abs(l)%2 == 1))
			return true;
		return false;

	}
	if(!ctype.compare("BCC")){
		if(abs((h+k+l)%2) == 0)
			return true;
		return false;
	}
	if(!ctype.compare("DC")){
		if((h%2 == 0 && k%2 == 0 && l%2 == 0) || (abs(h)%2 == 1 && abs(k)%2 == 1 && abs(l)%2 == 1)){
			if(abs((h+k+l)%4) == 0 || abs((h+k+l)%4) == 1 || abs((h+k+l)%4) == 3)
				return true;
			return false;
		}
	}
	return false;
}
/*****************************************
 * Build off diagonal elements of eigenvalue matrix
 * *****************************************/
void Crystal::buildXmatrix(const Eigen::MatrixXd& g, bool absorb, Eigen::MatrixXcd& X, const Beam& beam){
	int gsize = g.rows();
	int s_size = s.rows();
	X.resize(gsize,gsize);
	dcomp Vg, Ug;
	Eigen::Vector3d gdiff;

	/*Fill the off diagonal elements*/
	
	clock_t start_time = clock();

	for(int i=0;i<gsize;i++){
		for(int j=0;j<i;j++){
			gdiff = g.row(i)-g.row(j);
			Vg = calcVg(gdiff, beam);
			if(absorb == true){
				dcomp V_imag = calcVgAbsorb(gdiff,beam);
				Vg = Vg + ri*V_imag;
			}
			Ug = 2*beam.get_mass()*abs(charge)/(h*h*1e20)*Vg;
			X(i,j) = Ug;

			/*For the upper triangle*/
			if(absorb == false)
				X(j,i) = conj(X(i,j));
			else if(absorb == true){
				gdiff = g.row(j)-g.row(i);
				Vg = calcVg(gdiff, beam);
				dcomp V_imag = calcVgAbsorb(gdiff, beam);
				Vg = Vg + ri*V_imag;
				Ug = 2*beam.get_mass()*abs(charge)/(h*h*1e20)*Vg;
				X(j,i) = Ug;

			}
		}
	}
	clock_t end_time = clock();
	clock_t nbrticks = end_time-start_time;
	double nbrseconds = nbrticks/(double) CLOCKS_PER_SEC;
//	cout<<"Time to build matrix (in s): "<<nbrseconds<<"\n";

}

/****************************************
 * Build electrostatic potential
 * **************************************/
void Crystal::build_pot(const Eigen::MatrixXd& g, bool absorb, vector<double> xg, vector<double> yg, double zf, const Beam& beam){
	int xlength = xg.size(); int ylength = yg.size();
	potential = Eigen::ArrayXXcd::Zero(ylength,xlength);
	int g_size = g.rows();
	int s_size = s.rows();
	dcomp Vg;

	for(int i=0;i<ylength;i++){
		for(int j=0;j<xlength;j++){
			for(int k=0;k<g_size;k++){
				Vg = calcVg(g.row(k), beam);
				if(absorb == true){
					dcomp V_imag = calcVgAbsorb(g.row(k), beam);
					Vg = Vg + ri*V_imag;
				}
				potential(i,j) = potential(i,j) + Vg*exp(2*pi*ri*(g(k,0)*xg[j] + g(k,1)*yg[i] + g(k,2)*zf));
			}
		}
	}

	potential = charge*potential*toJ;

	/*Electrostatic force*/
	double dx = abs(xg[0]-xg[1]);
	double dy = abs(yg[0]-yg[1]);
	dcomp partialx, partialy;
	fex.resize(ylength,xlength);
	fey.resize(ylength,xlength);

	for(int i=0;i<ylength-1;i++){
		for(int j=0;j<xlength-1;j++){
			fex(i,j) = -real(potential(i,j+1) - potential(i,j))/dx;
			fey(i,j) = -real(potential(i+1,j) - potential(i,j))/dy;
		}
	}

	for(int i=0;i<ylength-1;i++){
		fex(i,xlength-1) = -real(potential(i,0) - potential(i,xlength-1))/dx;
		fey(i,xlength-1) = -real(potential(i+1,xlength-1) - potential(i,xlength-1))/dy;
	}

	for(int j=0;j<xlength-1;j++){
		fex(ylength-1,j) = -real(potential(ylength-1,j+1) - potential(ylength-1,j))/dx;
		fey(ylength-1,j) = -real(potential(0,j) - potential(ylength-1,j))/dy;
	}
	fex(ylength-1,xlength-1) = -real(potential(ylength-1,xlength-1) - potential(ylength-1,0))/dx;
	fey(ylength-1,xlength-1) = -real(potential(0,xlength-1) - potential(ylength-1,xlength-1))/dy;
}

/*******************************************************
 * Centered difference to calculate electrostatic force
 * at single point
 * *****************************************************/
void Crystal::fe_CD(const Eigen::Vector3d& r, const Eigen::Vector3d& dr, Eigen::Vector3d& f, const Eigen::MatrixXd& g, bool absorb, const Beam& beam) const{
	
	int g_size = g.rows();
	int s_size = s.rows();
	dcomp Vg;
	Eigen::Vector3cd potnm1, potnp1;
	potnm1 << 0,0,0;
	potnp1 << 0,0,0;
	Eigen::Matrix3d dr_m;
	dr_m << 0.5, 0, 0,
		0, 0.5, 0,
		0, 0, 0.5;

	/*f_(n-1)*/
	for(int i=0;i<3;i++){
		for(int k=0;k<g_size;k++){
			Vg = calcVg(g.row(k), beam);
			if(absorb == true){
				dcomp V_imag = calcVgAbsorb(g.row(k),beam);
				Vg = Vg + ri*V_imag;
			}
			potnm1(i) = potnm1(i) + Vg*exp(2*pi*ri*(g(k,0)*(r(0)-dr_m(i,0)*dr(0)) + g(k,1)*(r(1)-dr_m(i,1)*dr(1)) + g(k,2)*(r(2)-dr_m(i,2)*dr(2))));
		}
	}

	/*f_(n+1)*/
	for(int i=0;i<3;i++){
		for(int k=0;k<g_size;k++){
			Vg = calcVg(g.row(k), beam);
			if(absorb == true){
				dcomp V_imag = calcVgAbsorb(g.row(k), beam);
				Vg = Vg + ri*V_imag;
			}
			potnp1(i) = potnp1(i) + Vg*exp(2*pi*ri*(g(k,0)*(r(0)+dr_m(i,0)*dr(0)) + g(k,1)*(r(1)+dr_m(i,1)*dr(1)) + g(k,2)*(r(2)+dr_m(i,2)*dr(2))));
		}
	}
	potnp1 *= charge*toJ;
	potnm1 *= charge*toJ;
	Eigen::Vector3d dr_inv;
	dr_inv = dr.cwiseInverse();
	f = -(potnp1.real() - potnm1.real()).cwiseProduct(dr_inv);
}
/****************************************
 * Displace atoms for TDS calculation
 * **************************************/
void Crystal::TDSdisplacement(double sigmaTDS){

		double gauss1,gauss2,gauss3;
		int s_size = s.rows();
		for(int i=0;i<s_size;i++){
			gauss1=AuxFns::gasdev()*sigmaTDS; 
			gauss2=AuxFns::gasdev()*sigmaTDS; 
			gauss3=AuxFns::gasdev()*sigmaTDS;
			s(i,0) = gauss1 + s(i,0); s(i,1) = gauss2 + s(i,1); s(i,2) = gauss3 + s(i,2);
		}
}
void Crystal::addPerturbation(const Eigen::MatrixXd& g, const Eigen::MatrixXd& gpert, const Beam& beam, bool absorb, Eigen::MatrixXcd& X){
	int gsize = g.rows();
	int s_size = s.rows();
	int gpert_size = gpert.rows();
	Eigen::Vector3d gdiff1, gdiff2, gdiff;
	dcomp Vg, Vg1, Vg2, Ug, Ug1, Ug2, Ugsum; 	
	double fg, fg1, fg2, sh;

	/*Compute perturbation contributions*/
	for(int i=0;i<gsize;i++){
		for(int j=0;j<gsize;j++){
			if(i==j){
				Ugsum = 0;
				for(int m=0;m<gpert_size;m++){
					gdiff = g.row(i)-gpert.row(m);
					sh = AuxFns::sgcalc(gpert.row(m),beam.k_in);
					Vg = calcVg(gdiff, beam);
					if(absorb == true){
						dcomp V_imag = calcVgAbsorb(gdiff,beam);
						Vg = Vg + ri*V_imag;
					}
					Ug = 2*beam.get_mass()*abs(charge)/(h*h*1e20)*Vg;	
					Ugsum = Ugsum - Ug*Ug/(2*beam.k_in.norm()*sh);
				}


			}
			else{
				Ugsum = 0;
				for(int m=0;m<gpert_size;m++){
					gdiff1 = g.row(i)-gpert.row(m);
					gdiff2 = gpert.row(m)-g.row(j);
					sh = AuxFns::sgcalc(gpert.row(m),beam.k_in);
	
					Vg1 = calcVg(gdiff1, beam);
					Vg2 = calcVg(gdiff2, beam);
					if(absorb == true){
						dcomp V_imag1 = calcVgAbsorb(gdiff1,beam);
						dcomp V_imag2 = calcVgAbsorb(gdiff2,beam);
						Vg1 = Vg1 + ri*V_imag1;
						Vg2 = Vg2 + ri*V_imag2;
					}
					
					Ug1 = 2*beam.get_mass()*abs(charge)/(h*h*1e20) * Vg1; Ug2 = 2*beam.get_mass()*abs(charge)/(h*h*1e20) * Vg2;	
					Ugsum = Ugsum - Ug1*Ug2/(2*beam.k_in.norm()*sh);

				}
			}
			X(i,j) = X(i,j) - Ugsum;
		}
	}

}

void Crystal::atomwrite(){
	/*Print atom coordinates*/
	ofstream s_out("Input/atoms.txt");
	s_out << s;
	s_out.close();
}

/*****************************************************
 * Build rotation matrix for given direction zone_axis
 * ***************************************************/
void Crystal::buildRot(Eigen::Matrix3d& Rot, const Eigen::Vector3d& zone_axis){
	Eigen::Vector3d u,nu;
	u << 0,0,1;
	double alpha;
	alpha = acos(u.dot(zone_axis)/(u.norm()*zone_axis.norm()));
	nu = u.cross(zone_axis);
	nu.normalize();

	Rot(0,0) = cos(alpha)+nu(0)*nu(0)*(1-cos(alpha));
	Rot(0,1) = nu(0)*nu(1)*(1-cos(alpha))-nu(2)*sin(alpha);
	Rot(0,2) = nu(0)*nu(2)*(1-cos(alpha))+nu(1)*sin(alpha);
	Rot(1,0) = nu(1)*nu(0)*(1-cos(alpha))+nu(2)*sin(alpha);
	Rot(1,1) = cos(alpha)+nu(1)*nu(1)*(1-cos(alpha));
	Rot(1,2) = nu(1)*nu(2)*(1-cos(alpha))-nu(0)*sin(alpha);
	Rot(2,0) = nu(2)*nu(0)*(1-cos(alpha))-nu(1)*sin(alpha);
	Rot(2,1) = nu(2)*nu(1)*(1-cos(alpha))+nu(0)*sin(alpha);
	Rot(2,2) = cos(alpha)+nu(2)*nu(2)*(1-cos(alpha));
	
}
