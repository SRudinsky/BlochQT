#ifndef TRAJECTORY_H
#define TRAJECTORY_H
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include "Eigen/Core"
#include <fstream>
#include "Crystal.h"
#include "Beam.h"
#include "WaveFn.h"

using namespace std;

/****************************
 * Class definition for electron trajectories
 * **************************/

class Trajectory {
	private:
		bool inside;

	public:
		Trajectory();
		~Trajectory();
		double xcoord, ycoord, zcoord;
		double vx, vy, vz;
		void init_position(const double xg, const double yg);
		void init_velocity(const double vx_0, const double vy_0, const double vz_0);
		int near_point(char coord,vector<double> xg);
		void RK2(const Eigen::ArrayXXcd& psig, const Beam& beam, const Eigen::MatrixXd& g, const Eigen::MatrixXd& q, double dt);
		void RK2(const Eigen::VectorXcd& psig, const Beam& beam, const Eigen::MatrixXd& g, double dt);

		void RK2_fe(const Crystal& crystal, const Eigen::MatrixXd& g, bool absorb, double dt, const Beam& beam);
		void RK2_fq(const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& gamma, const Eigen::MatrixXcd& Coeff, const Beam& beam, const Eigen::MatrixXd& g, double dt, const WaveFn& wf);
};
#endif
