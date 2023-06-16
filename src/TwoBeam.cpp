#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <complex>
#include <vector>
#include <ctime>
#include <cstdio>
#include <array>	
#include <algorithm>
#include "../include/Eigen/Dense"
#include "../include/Eigen/Core"
#include "../include/Eigen/Eigenvalues"
#ifdef __linux__
#include <bits/stdc++.h>
#elif _WIN32
#include <Windows.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include "../lib/Crystal.h"
#include "../lib/Beam.h"
#include "../lib/WaveFn.h"
#include "../lib/Trajectory.h"
#include "../lib/AuxFns.h"
#include "../lib/Constants.h"

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;

/* All calculations done in Angstrom! */

int main(){
 	system("exec rm Results/*.txt");

/*Create output directories*/
#ifdef __linux__
	if(mkdir("Results", 0777) == -1)
		cerr << "Error : "<< strerror(errno) << "\n";
	else
		cout << "Results directory created\n";

	if(mkdir("Input",0777) == -1)
		cerr << "Error : "<<strerror(errno)<<"\n";
	else
		cout<<"Input directory created\n";
#elif _WIN32
	if (!(CreateDirectoryW(L"Results", NULL)))
	{
		std::cerr << "Results folder exists.\n";
	}
	if (!(CreateDirectoryW(L"Input", NULL)))
	{
		std::cerr << "Input folder exists.\n";
	}
#endif

	
/*Input parameters for beam*/
	Beam beam;
	double E0_in;
	cout<<"Energy (eV): \n";
	cin>>E0_in;
	beam.set_E0(E0_in); //in eV
	beam.set_microscope_settings(10,0,1.1e8); //convergence angle in mrad and dF and Cs in Angstrom

/*The Crystal!*/
	Crystal cryst;
	
/*Input parameters for material*/
	string mat_in;
	cout << "File name including location w.r.t. executable\n";
	cin >> mat_in;
	ifstream cinput(mat_in);
	cryst.read_cfg(cinput);

/*Build crystal*/	
	cryst.buildCrystal();
	cryst.atomwrite();

	double aax = cryst.get_aax();
	double aay = cryst.get_aay();
	double aaz = cryst.get_aaz();

/*Compute incident wave vector*/
	Vector3d kt;
	char tilt;
	bool SysR = false;
	Vector3d gSR;

	cout<<"Tilt (y/n)?\n";
	cin>>tilt;

	if(tilt == 'n')
		kt << 0,0,0;
	else if(tilt == 'y'){
		char BorA;
		cout<<"Bragg (B) or choose angle (A) or systematic row (S)?\n";
		cin>>BorA;
		if(BorA == 'B'){
			Vector3d bSpot;
			cout<<"Enter Bragg reflection:\n";
			cin>>bSpot(0)>>bSpot(1)>>bSpot(2);
			kt = -(0.5)*(bSpot(0)*cryst.b.col(0) + bSpot(1)*cryst.b.col(1) + bSpot(2)*cryst.b.col(2));
			cout<<"Bragg angle: "<<asin(beam.get_lambda()*(bSpot(0)*cryst.b.col(0) + bSpot(1)*cryst.b.col(1) + bSpot(2)*cryst.b.col(2)).norm())*0.5<<"\n";
		}
		else if(BorA == 'A'){
			double tiltangle;
			cout<<"Enter tilt angle (single direction in x):\n";
			cin>>tiltangle;
			double ktx_B = (1/beam.get_lambda())*sin(tiltangle);
			kt << ktx_B,0,0;
		}
		else if(BorA == 'S'){
			SysR = true;
			double tiltSR;
			cout<<"Enter beam for systematic row: \n";
			cin>>gSR(0)>>gSR(1)>>gSR(2);
			cout<<"Enter integer factor for g:\n";
			cin>>tiltSR;

			kt = tiltSR*(gSR(0)*cryst.b.col(0) + gSR(1)*cryst.b.col(1) + gSR(2)*cryst.b.col(2));
		}
		else{
			cout<<"Invalid input (B or A)\n";
			exit(0);
		}
	}
	else{
		cout<<"Invalid input (y or n)\n";
		exit(0);
	}
	cout<<"kt = "<<kt.norm()<<"\n";
	beam.calc_kin(kt);

/*Compute Two beam */
	int hh,kk,ll;
	cout<<"Enter choice of beam (h k l)\n";
	cin>>hh>>kk>>ll;
	MatrixXd g;
	g.resize(2,3);
	g.row(0) << 0,0,0;
	g.row(1) = hh*cryst.b.col(0)+kk*cryst.b.col(1)+ll*cryst.b.col(2);
	double sg = AuxFns::sgcalc(g.row(1),beam.k_in);
	double xi_g = AuxFns::calcExtinctionDistance(g.row(1),cryst,beam);
	cout<<"Extinction distance = "<<xi_g<<"\n";

/*Determine thicknesses for calculation*/
	double zf;
	cout<<"Enter max thickness (in Angstrom): \n";
	cin>>zf;
	int nslices;
	cout<<"Enter number of slices\n";
	cin >> nslices;
	unique_ptr<double[]> z(new double[nslices]);
	double dz = zf/(double)nslices;
	for(int i=1;i<=nslices;i++)
		z[i-1] = i*dz;

/*The grid and real space representation*/
	vector<double> xg,yg,kxg,kyg,ckx,cky;
	int npointsx,npointsy;
	cout<<"Number of real space grid points in x and y?\n";
	cin>> npointsx >> npointsy;
	int npoints[2] = {npointsx,npointsy};
	double grid_limit[2] = {aax,aay};
	double gmax = AuxFns::RSgrid(grid_limit,npoints,xg,yg);
	double dx = grid_limit[0]/npointsx;
	double dy = grid_limit[1]/npointsy;
	AuxFns::KSgrid(grid_limit,npoints,kxg,kyg,ckx,cky);

/*Compute zeries of wave functions for each z*/
	WaveFn wf(yg,xg);
	bool absorb = false;
	ofstream z_out("Results/Depth.txt");
	z_out << "Depth in angstrom\n";
	for(int slice=0;slice<nslices;slice++){
		cout<<"Slice: "<<z[slice]<<" Angstrom\n";
		z_out << z[slice]<<"\n";
		ostringstream ss;
		ss<<slice;
		string zstr = ss.str();

		wf.set_thickness(z[slice]);
		wf.build_TwoBeam_psig(xi_g,sg);
		wf.build_EW(beam, g, xg, yg);
		wf.build_vfield(xg, yg, g, beam);
		wf.build_qpot(dx, dy, beam);
		cryst.build_pot(g, absorb, xg, yg, z[slice], beam);

		AuxFns::write_psif(zstr, wf.psif);
		AuxFns::write_psig(zstr, g, wf.psig);
		AuxFns::write_vfield(zstr, wf.vfieldx, wf.vfieldy);
		AuxFns::write_qp(zstr, wf.Qpot);
		AuxFns::write_potential(zstr, cryst.potential);
		AuxFns::write_fe(zstr, cryst.fex, cryst.fey);
		AuxFns::write_fq(zstr, wf.fqx, wf.fqy);
	}

/*Compute trajectories in time*/
	char traj_choice;
	cout<<"Compute trajectories?\n";
	cin>>traj_choice;
	
	if(traj_choice == 'y'){
/*Initialize trajectories*/	
		int ntrajx,ntrajy;
		cout<<"Number of trajectories in x and y?\n";
		cin>>ntrajx>>ntrajy;
		int ntraj_total = ntrajx*ntrajy;
		vector<Trajectory> traj(ntraj_total);
		vector<ofstream> traj_out;
		cout<<traj.size()<<"\n";
		int cnt = 0;
		for(int i=0;i<ntrajx;i++){
			for(int j=0;j<ntrajy;j++){
				traj[cnt].init_position(aax*i/(double)ntrajx,j*aay/(double)ntrajy);
				ostringstream ss;
				ss<<cnt;
				string zstr = ss.str();
				ofstream out;
				out.open("Results/Traj"+zstr+".txt");
				traj_out.push_back(move(out));
				AuxFns::write_traj(traj_out[cnt],traj[cnt]);
				cnt++;
			}
		}
		double dt = 1e-10; //Time invariant system, dt is just an increment in space, unless using force equation
	
		double avg_z = 0;
		while(avg_z<zf){
			if(abs(avg_z/zf-0.25) < 0.00001)
				cout<<"25% done\n";
			if(abs(avg_z/zf-0.5) < 0.00001)
				cout<<"50% done\n";
			if(abs(avg_z/zf-0.75) < 0.00001)
				cout<<"75% done\n";
	
			avg_z = 0;
			for(int j=0;j<ntraj_total;j++){
				wf.set_thickness(traj[j].zcoord);
				wf.build_TwoBeam_psig(xi_g, sg);
				traj[j].RK2(wf.psig, beam, g, dt);
				AuxFns::write_traj(traj_out[j],traj[j]);
				avg_z += traj[j].zcoord;
			}
			avg_z /= ntraj_total;
		}
	
		for(int j=0;j<ntraj_total;j++)
			traj_out[j].close();
	}

	return 0;

}
