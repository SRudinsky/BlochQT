#define BOOST_TEST_MODULE TEST_AuxFns
#include <vector>
#include "../include/Eigen/Core"
#include <boost/test/included/unit_test.hpp>
#include "../lib/AuxFns.h"

using namespace std;

//BOOST_AUTO_TEST_SUITE(test_traj)

BOOST_AUTO_TEST_CASE(test_is_working){
    //BOOST_FAIL("Nothing to test");
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(extinctionDistance){
	Crystal crystal;
	ifstream cin("../bin/Materials/Cu.cfg");
	crystal.read_cfg(cin);
	crystal.buildCrystal();
	
	Beam beam;
	beam.set_E0(100000);

	Eigen::Vector3d g1; g1 << 1,1,1;
	double xig = AuxFns::calcExtinctionDistance(g1,crystal,beam);

	BOOST_CHECK_CLOSE(xig,286,10.0);
}
