#include "Planets.h"
#include "Physics.h"
#include <iostream>
#include <cmath>
#include <cassert>

//Test case 1 Stationary bodies
void testStationary() {

	Planets sun(1.989e30, 0, 0, 0, 0, 0, 0);
	Planets rock(1e10, 1e7, 0, 0, 0, 0, 0);

	// Time skip (not fully accurate due to non constant acceleration)
	// (However 1s constant acceleration is a decent approximation for astronomical time periods)
	double dt = 1;

	applyGravity(rock, sun);
	rock.update(dt);

	std::cout << "[Stationary test]: Rock Postion:"
		      << rock.getX() << ", " << rock.getY() << ", " << rock.getZ() << "\n";

	// Manually calculated test case values
	double expectedX = 1e7 - 1.32751827e6 * dt;
	double expectedY = 0;

	//Tolerance for floting point error
	double tol = 1e-3;

	//Test Case verification, temporary will recode
	assert(std::abs(rock.getX() - expectedX) < tol);
	assert(std::abs(rock.getY() - expectedY) < tol);

	std::cout << "[Verification Passed Test Case 1]\n";
};
//Test case 2 Circular orbit (Work in progress)
void testCircular() {

	// Determine planet objects with values

	// Manually calculated test case values

	// Simulate for one full expected orbit

	// Calculate excentricity of simulated orbit

	// Compare to expected excentricity

	// Pass/Fail test case

}
//Test case 3 Post-Newtonian Precession
void testPrecession() {

	//test case for Post-Newtonian Precession of the orbits

};
//Test case 4 3D Orbit
void test3DimOrbit() {

	//Test Case with a plane not perpendicular to the z-x and z-y planes.

};

//Test
int main() {

	std::cout << "This is the Test cases...\n\n";

	testStationary();
	testCircular();

	std::cout << "All test run.\n";

	return 0;

}