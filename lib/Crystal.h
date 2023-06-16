#ifndef CRYSTAL_H
#define CRYSTAL_H
#include "Eigen/Core"
#include <vector>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include "Beam.h"

/*****************************************************************************
 * Class definition for a crystal 
 *****************************************************************************/ 	   

class Crystal {
	typedef std::complex<double> dcomp;
	private:
	        Eigen::Vector3d unitCellvx, unitCellvy, unitCellvz;
 	       	double aa,aax,aay,aaz;
		/*Todo: use Weickenmeier and Kohl for absorption*/
		void buildRot(Eigen::Matrix3d& Rot,const Eigen::Vector3d& zone_axis);
		double set_aax();
		double set_aay();
		double set_aaz();
	public:	
		Crystal();
		~Crystal();
		double omega_inv;
		std::vector<uint64_t> Z, atomTypes;
		std::vector<double> DWfactors;
		Eigen::ArrayXXcd potential;
		std::string ctype;
		Eigen::MatrixXd s,b;
        	uint64_t nElements;
	        void read_cfg(std::ifstream& datafile);
		double get_aax() const;
		double get_aay() const;
		double get_aaz() const;

		/*Crystal and g-vectors*/
		void buildCrystal();
		/*TODO: make general allow reflections*/
		dcomp calcVg(const Eigen::Vector3d& gvector, const Beam& beam) const;
		dcomp calcVgAbsorb(const Eigen::Vector3d& gvector, const Beam& beam) const;
		bool allowReflection(int h, int k, int l) const;
		void addPerturbation(const Eigen::MatrixXd& g, const Eigen::MatrixXd& gpert, const Beam& beam, bool absorb, Eigen::MatrixXcd& X);
		void buildXmatrix(const Eigen::MatrixXd& g, bool absorb, Eigen::MatrixXcd& X, const Beam& beam);
		void atomwrite();

		/*Electrostatic potential and force*/
		Eigen::ArrayXXd fex,fey;
		void fe_CD(const Eigen::Vector3d& r, const Eigen::Vector3d& dr, Eigen::Vector3d& f, const Eigen::MatrixXd& g, bool absorb, const Beam& beam) const;
		void build_pot(const Eigen::MatrixXd& g, bool absorb, std::vector<double> xg, std::vector<double> yg, double zf, const Beam& beam);
		void TDSdisplacement(double sigmaTDS);
};
#endif
