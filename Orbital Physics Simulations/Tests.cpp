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
	Planets sun(1.989e30, 0, 0, 0, 0, 0, 0);
	Planets planet(1e10, 1e7, 0, 0, 0, 0, 0);
	double G = 6.67430e-11;

	// Manually calculated test case values
	double expected_rmin = 1e7;
	double expected_rmax = 1e7;

	double tol = 1e-3;

	// Simulate for one full expected orbit

	double rMin = find_rMin(planet, sun);
	double rMax = find_rMax(planet, sun);

	// Compare to expected excentricity
	assert(std::abs(rMin - expected_rmin) < tol);
	assert(std::abs(rMax - expected_rmax) < tol);

	// Pass/Fail test case
	std::cout << "[Verification Passed Test Case 2]\n";
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