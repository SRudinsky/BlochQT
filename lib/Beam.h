#ifndef BEAM_H
#define BEAM_H
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include "Eigen/Core"

using namespace std;

/**********************************************
 * Class definition for the indicent beam and microscope
 * settings
 * *******************************************/

class Beam {
	private:
		/*Lambda in Angstrom, E0 in joules*/
		double E0,lambda,aconv,dF,Cs,speed,mass;

	public:
		Beam();
		~Beam();
		Eigen::Vector3d k_in, v_in;
		void calc_kin(const Eigen::Vector3d& kt);
		void set_E0(double EeV);
		void set_microscope_settings(double ac,double df,double cs);
		double get_lambda() const;
		double get_aconv() const;
		double get_dF() const;
		double get_E0() const;
		double get_Cs() const;
		double get_v() const;
		double get_mass() const;
		void build_planewave(Eigen::VectorXcd& psi0);
		void build_probe(const Eigen::Vector3d& R0, const Eigen::MatrixXd& g,Eigen::VectorXcd& psi0);

};
#endif
