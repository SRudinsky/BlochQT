#ifndef WAVEFN_H
#define WAVEFN_H
#include <complex>
#include <vector>
#include "Eigen/Core"
#include "Beam.h"

using namespace std;

/*****************************************************
 * Class definition for the exit wave function and all 
 * properties calculated from it
 * ***************************************************/

class WaveFn {
	private:
		double zf;
		int xlength, ylength;
		Eigen::ArrayXXcd dpsix;
		Eigen::ArrayXXcd dpsiy;

	public:
		WaveFn(vector<double> yg,vector<double> xg);
		~WaveFn();
		Eigen::VectorXcd psig;
		Eigen::ArrayXXcd psif;
		Eigen::ArrayXXd vfieldx, vfieldy, Qpot, fqx, fqy;

		void set_thickness(double z);
		double get_zf();

		void build_psig(const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Eigen::VectorXcd& psi0);
		void build_TwoBeam_psig(const double xig, const double sg);
		void build_EW(const Beam& beam, const Eigen::MatrixXd& g, const vector<double>& xg, const vector<double>& yg);
		void build_BW(uint64_t bw, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Eigen::VectorXcd& psi0); 
		void build_vfield(const vector<double>& xg, const vector<double>& yg, const Eigen::MatrixXd& g, const Beam& beam);
		void build_qpot(double dx, double dy, const Beam& beam);
		void build_qforce();
		double Qcalc_pt(const Eigen::Vector3d& r, const double dr, const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Beam& beam, const Eigen::MatrixXd& g) const;
		void fq_CD(const Eigen::Vector3d& r, const double dr, Eigen::Vector3d& f, const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Beam& beam, const Eigen::MatrixXd& g) const;
};
#endif
