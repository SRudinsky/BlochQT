#include "Eigen/Core"
#include <vector>
#include <complex>
#include "Crystal.h"
#include "Beam.h"
#include "Trajectory.h"

using namespace std;

namespace AuxFns {

	/*Functions pertaining to the g-vector*/
	double sgcalc(const Eigen::Vector3d& g, const Eigen::Vector3d& k);
	void gZA(const Beam& beam, const Crystal& crystal, double gmax, double Climit[2], bool absorb, Eigen::MatrixXd& g, Eigen::MatrixXd& gpert, vector<double>& sg);
	void gSR(const Beam& beam, const Crystal& crystal, int maxg, bool absorb, Eigen::MatrixXd& g, vector <double>& sg, const Eigen::Vector3d& gSR);
	double calcExtinctionDistance(const Eigen::Vector3d& g, const Crystal& crystal, const Beam& beam);

	/*Grid functions*/
	void convertKSpace(const Eigen::Vector3d& x, Eigen::Vector3d& x_new, const Crystal& crystal);
	double RSgrid(double aa[2], int npoints[2], vector<double>& xg, vector<double>& yg);
	void KSgrid(double aa[2], int npoints[2], vector<double>& kxg, vector<double>& kyg, vector<double>& ckx, vector<double>& cky);

	/*Write functions*/
	void write_psif(string iterator, const Eigen::ArrayXXcd& psif);
	void write_psig(string iterator, const Eigen::MatrixXd& g, const Eigen::VectorXcd& psig);
	void write_fpsi(string iterator, const Eigen::ArrayXXcd& fpsi);
	void write_vfield(string iterator, const Eigen::ArrayXXd& vfieldx, const Eigen::ArrayXXd& vfieldy);
	void write_traj(ofstream &out, Trajectory traj);
	void write_qp(string iterator, const Eigen::ArrayXXd& qpot);
	void write_potential(string iterator, const Eigen::ArrayXXcd& pot);
	void write_fe(string iterator, const Eigen::ArrayXXd& fex, const Eigen::ArrayXXd& fey);
	void write_fq(string iterator, const Eigen::ArrayXXd& fqx, const Eigen::ArrayXXd& fqy);

	/*Random functions*/
	double ran1();
	double gasdev();
	void sort(Eigen::Vector3d& a);

}
