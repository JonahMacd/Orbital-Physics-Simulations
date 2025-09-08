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
//Test case 2 Circular orbit
void testCircular() {

	//test case for circular orbits
	Planets sun(1.989e30, 0, 0, 0, 0, 0, 0);
	Planets planet(1e10, 1e7, 0, 0, 0, 0, 0);

	double G = 6.67430e-11;
	double c = 3e8;
	double dt = 0.001;
	int steps = 1000;

	double r0 = 1e7;
	double v_circ = std::sqrt(G * sun.getMass() / r0);
	planet.setVelocity(0, v_circ, 0);

	double tol = 1e-2; // meters
	bool pass = true;

	for (int i = 0; i < steps; ++i) {
		applyGravity(planet, sun, G, c);
		planet.update(dt);

		double dx = planet.getX() - sun.getX();
		double dy = planet.getY() - sun.getY();
		double dz = planet.getZ() - sun.getZ();
		double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

		// For early timesteps, distance should stay roughly near r0
		if (std::abs(distance - r0) > tol) {
			pass = false;
			break;
		}
	}

	if (pass)
		std::cout << "[PASS] Elliptical orbit stays within expected bounds\n";
	else
		std::cout << "[FAIL] Orbit deviates beyond expected tolerance\n";
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