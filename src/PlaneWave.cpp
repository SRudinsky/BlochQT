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
#include <sys/stat.h>
#include <sys/types.h>
#ifdef __linux__
#include <bits/stdc++.h>
#elif _WIN32
#include <Windows.h>
#endif
#include "../lib/WaveFn.h"
#include "../lib/Crystal.h"
#include "../lib/Beam.h"
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
	
	/*****************************/
	/*The Crystal!*/
	Crystal cryst;
	
/*Input parameters for material*/
	string mat_in;
	int isCu;
	cout << "Select material: (0) Al, (1) Cu\n";
	cin >> isCu;

#ifdef __linux__
	if(isCu)
		mat_in = "../Materials/Cu.cfg";
	else
		mat_in = "../Materials/Al.cfg";

#elif _WIN32
	if(isCu)
		mat_in = "..\\..\\Materials\\Cu.cfg";
	else
		mat_in = "..\\..\\Materials\\Al.cfg";
#endif
		
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
			cout<<"Bragg angle: "<<asin(beam.get_lambda()*(bSpot(0)*cryst.b.col(0) + bSpot(1)*cryst.b.col(1) + bSpot(2)*cryst.b.col(2)).norm()*0.5)<<"\n";
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
	cout<<"k0 = "<<beam.k_in.norm()<<"\n";

/*The grid and real space representation*/
	vector<double> xg,yg;
	int npointsx,npointsy;
	cout<<"Number of real space grid points in x and y?\n";
	cin>> npointsx >> npointsy;
	int npoints[2] = {npointsx,npointsy};
	double grid_limit[2] = {aax,aay};
	double dx = grid_limit[0]/npointsx;
	double dy = grid_limit[1]/npointsy;
	double gmax = AuxFns::RSgrid(grid_limit,npoints,xg,yg);

/*Absorption?*/
	char ans;
	bool absorb = false;
	cout<<"Include absorption (y/n)?\n";
	cin>>ans;

	if(ans=='n')
		absorb = false;
	else if(ans=='y')
		absorb = true;
	else{
		cout<<"Not an option!\n";
		exit(0);
	}
	/*Haven't implemented absorption properly yet*/
	cout << "Absorption not implemented yet\n";
	absorb = false;

/*Compute g vectors to be used*/
	gmax = 5;
	MatrixXd g,gpert;
	Vector3d g1,g2;
	vector<double> sg;
	if(SysR){
		int max_int;
		cout<<"Enter max reflexion integer multiple:\n";
		cin >> max_int;
		AuxFns::gSR(beam, cryst, max_int, absorb, g, sg, gSR);
	}
	else{
		char TB;
		cout<<"Two beam (y/n)?\n";
		cin >> TB;
		if(TB == 'y'){
			int TB_i, TB_j, TB_k;
			cout<<"Enter beam:\n";
			cin >> TB_i >> TB_j >> TB_k;
			g.resize(2,3);
			g.row(0) << 0,0,0;
			g.row(1) = TB_i*cryst.b.col(0) + TB_j*cryst.b.col(1) + TB_k*cryst.b.col(2);
		}
		else{
			double c1_in, c2_in;
			cout<<"Input limit c1 (strong beam) and c2 (weak beam) in 1/V/Angstrom^2 \n";
			cin>>c1_in>>c2_in;
			double Climit[2] = {c1_in,c2_in};
			cout<<"gmax: "<<gmax<<" Climit: ["<<Climit[0]<<","<<Climit[1]<<"]"<<"\n";
			AuxFns::gZA(beam,cryst,gmax,Climit,absorb,g,gpert,sg);
		}
	}
	const int g_size = g.rows();
	cout<<"Number of g-vectors in diagonalization: "<<g_size<<"\n";
	ofstream g_out("Input/gvectors.txt");
	g_out<<g;
	g_out.close();

	cout<<"Number of g-vectors in perturbation: "<<gpert.rows()<<"\n";
	ofstream gpert_out("Input/gpertvectors.txt");
	gpert_out<<gpert;
	gpert_out.close();

/*Build bloch wave matrix*/
	MatrixXcd X;

	cryst.buildXmatrix(g,absorb,X, beam);

	/*Compute diagonal of matrix*/
	Vector3d gX;
	for(int i=0;i<g_size;i++){
		gX = g.row(i);
		if(absorb == true){
			dcomp V_imag = cryst.calcVgAbsorb(gX,beam);
			X(i,i) = -(kt+gX).dot(kt+gX) + ri*(2*beam.get_mass()*abs(charge)*V_imag/(h*h)*1e-20); 
		}
		else
			X(i,i) = -(kt+gX).dot(kt+gX);
	}

	/*Add Bethe potentials*/
	cryst.addPerturbation(g,gpert,beam,absorb,X);

	ofstream x_o("Input/Xmatrix.txt");
	x_o << X;
	x_o.close();
	
/*Compute eigenvalues and eigenvectors of Bloch matrix*/
	VectorXcd gamma;
	MatrixXcd Coeff;
	Vector3d un;
	un << 0,0,1;

	clock_t start_time = clock();
	if(absorb == true){
		ComplexEigenSolver<MatrixXcd> es(X);
		gamma = es.eigenvalues()/beam.k_in.dot(un);
		Coeff = es.eigenvectors();
		ofstream g_o("Input/Gamma.txt");
		g_o<<gamma<<"\n";
		g_o.close();
	}
	else{ 
		SelfAdjointEigenSolver<MatrixXcd> es(X);
		gamma = es.eigenvalues()/beam.k_in.dot(un);
		Coeff = es.eigenvectors();
		ofstream g_o("Input/Gamma.txt");
		g_o<<gamma<<"\n";
		g_o.close();
	}
	clock_t end_time = clock();
	clock_t nbrticks = end_time-start_time;
	double nbrseconds = nbrticks/(double) CLOCKS_PER_SEC;
	cout<<"Time to compute eigenvalues (in s): "<<nbrseconds<<"\n";

/*Initial boundary conditions of wave function*/
	VectorXcd psi0 = VectorXcd::Zero(g_size);
	beam.build_planewave(psi0);

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

/*Compute zeries of wave functions for each z*/
	WaveFn wf(yg,xg);
	ofstream z_out("Results/Depth.txt");
	z_out << "Depth in angstrom\n";
	for(int slice=0;slice<nslices;slice++){
		cout<<"Slice: "<<z[slice]<<" Angstrom\n";
		z_out << z[slice]<<"\n";
		ostringstream ss;
		ss<<slice;
		string zstr = ss.str();

		wf.set_thickness(z[slice]);
		wf.build_psig(gamma, Coeff, psi0);
	//	wf.build_BW(0, gamma, Coeff, psi0);
		wf.build_EW(beam, g, xg, yg);
//		wf.build_vfield(xg, yg, g, beam);
//		wf.build_qpot(dx, dy, beam);
		//cryst.build_pot(g, absorb, xg, yg, z[slice], beam);
		

		AuxFns::write_psif(zstr, wf.psif);
		AuxFns::write_psig(zstr, g, wf.psig);
		AuxFns::write_vfield(zstr, wf.vfieldx, wf.vfieldy);
		//AuxFns::write_qp(zstr, wf.Qpot);
		//AuxFns::write_potential(zstr, cryst.potential);
		//AuxFns::write_fe(zstr, cryst.fex, cryst.fey);
		//AuxFns::write_fq(zstr, wf.fqx, wf.fqy);
	}

/*Compute trajectories in time*/
	char traj_choice;
	cout<<"Compute trajectories?\n";
	cin>>traj_choice;
	
	if(traj_choice == 'y'){
		int ntrajx,ntrajy;
		cout<<"Number of trajectories in x and y?\n";
		cin>>ntrajx>>ntrajy;
		vector<Trajectory> traj(ntrajx*ntrajy);
		vector<ofstream> traj_out;
		cout<<traj.size()<<"\n";
		int cnt = 0;	
		for(int i=0;i<ntrajx;i++){
			for(int j=0;j<ntrajy;j++){
				traj[cnt].init_position(aax*i/(double)ntrajx,aay*j/(double)ntrajy);
				traj[cnt].init_velocity(beam.v_in(0),beam.v_in(1),beam.v_in(2));
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
		int ntraj_total = ntrajx*ntrajy;
		double dt = 5e-10; //Time invariant system, dt is just an increment in space, unless using force equation
	
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
	//			traj[j].RK2_fq(psi0, gamma, Coeff, beam.k_in, g, dt, wf);
	//			traj[j].RK2_fe(cryst, g, absorb, dt);
				wf.set_thickness(traj[j].zcoord);
				wf.build_psig(gamma, Coeff, psi0);
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

