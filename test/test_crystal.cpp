#define BOOST_TEST_MODULE TEST_Crystal
#include <fstream>
#include <complex>
#include <boost/test/included/unit_test.hpp>
#include "../include/Eigen/Core"
#include "../include/Eigen/Dense"
#include "../lib/Crystal.h"
#include "../lib/Beam.h"

BOOST_AUTO_TEST_CASE(test_is_working){
    //BOOST_FAIL("Nothing to test")
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(readcfg){
    Crystal crystal;
    std::ifstream file("Cu.cfg");

    crystal.read_cfg(file);
    Eigen::MatrixXd stest(4,3);
    stest << 0, 0, 0,
            0.5, 0.5, 0,
            0, 0.5, 0.5,
            0.5, 0, 0.5;

    BOOST_CHECK_EQUAL(crystal.s,stest);

    crystal.buildCrystal();
    BOOST_CHECK_EQUAL(crystal.atomTypes.size(),4);
    BOOST_CHECK_EQUAL(crystal.s.rows(),4);
}

BOOST_AUTO_TEST_CASE(calc_absorb){
	Crystal crystal;
	std::ifstream file("Cu.cfg");
	crystal.read_cfg(file);
	crystal.buildCrystal();
	Beam beam;
	beam.set_E0(200000);
	Eigen::Vector3d kt,g1,g2;
	kt << 0,0,0;
	beam.calc_kin(kt);
	
	g1 << 0,0,0;
	g2 << -8.2989847575313291, -8.2989847575313291, -8.2989847575313291;
//	g2 = 2*crystal.b.col(0) + 4*crystal.b.col(1) + 2*crystal.b.col(2);

	complex<double> Vg1 = crystal.calcVgAbsorb(g1,beam);
	complex<double> Vg2 = crystal.calcVgAbsorb(g2, beam);
	std::cout<<Vg1<<"\t"<<Vg2<<"\n";
}	
