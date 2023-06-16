#define BOOST_TEST_MODULE TEST_Traj
#include <vector>
#include "../include/Eigen/Core"
#include <boost/test/included/unit_test.hpp>
#include "../lib/Trajectory.h"

using namespace std;

//BOOST_AUTO_TEST_SUITE(test_traj)

BOOST_AUTO_TEST_CASE(test_is_working){
    //BOOST_FAIL("Nothing to test");
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(near_point){
	/*Grid*/
	vector<double> xg,yg;
	double xf = 5.0, yf = 5.0;
	int ngrid = 10;
	double dx = xf/(double)ngrid, dy = yf/(double)ngrid;
	
	for(int i=0;i<ngrid;i++){
		xg.push_back(i*dx);
		yg.push_back(i*dy);
	}

	Trajectory t1;
	t1.init_position(2.75,3.25);

	int xp=t1.near_point('x',xg);
	int yp=t1.near_point('y',yg);

	BOOST_CHECK_EQUAL(xp,5);
	BOOST_CHECK_EQUAL(yp,6);

	Trajectory t2;
	t2.init_position(0.25,0.25);

	xp=t2.near_point('x',xg);
	yp=t2.near_point('y',yg);

	BOOST_CHECK_EQUAL(xp,0);
	BOOST_CHECK_EQUAL(yp,0);
}

