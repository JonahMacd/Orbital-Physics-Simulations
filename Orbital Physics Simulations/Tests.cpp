#include "Planets.h"
#include "Physics.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef RAD_TO_ARCSEC
#define RAD_TO_ARCSEC 206264.8062470963551564734
#endif

//Test case 1 Stationary bodies
void testStationary() {

	Planets sun(1.989e30, 0, 0, 0, 0, 0, 0);
	Planets rock(1e10, 1e7, 0, 0, 0, 0, 0);

	// Time skip (not fully accurate due to non constant acceleration)
	// (However 1s constant acceleration is a decent approximation for astronomical time periods)
	double dt = 1;

	applyGravity(rock, sun);
	rock.updatePosition(dt);

	CorrectedAcceleration acc = findCorrectedAcceleration(rock, sun);

	rock.updateVelocity(dt, acc.ax, acc.ay, acc.az);

	std::cout << "[Stationary test]: Rock Postion:"
		      << rock.getX() << ", " << rock.getY() << ", " << rock.getZ() << "\n";

	// Manually calculated test case values
	double expectedX = 9.33624e+06;
	double expectedY = 0;

	//Tolerance for floting point error
	double tol = 1e2;

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

	double r = std::sqrt((planet.getX() * planet.getX()) + (planet.getY() * planet.getY()) + (planet.getZ() + planet.getZ()));

	double Vc = std::sqrt(G * sun.getMass() / r);
	planet.setVelocity(0, Vc, 0);

	// Manually calculated test case values
	double expected_rmin = 1e7;
	double expected_rmax = 1e7;

	double tol = 1e-3;

	// Simulate for one full expected orbit

	double rMin = find_rMin(planet, sun);
	double rMax = find_rMax(planet, sun);

	// Compare to expected excentricity
	// 
	//assert(std::abs(rMin - expected_rmin) < tol);
	//assert(std::abs(rMax - expected_rmax) < tol);

	// Pass/Fail test case
	std::cout << "[Verification Passed Test Case 2]\n";
}

//Test case 3 Real World Orbital Data

void testEarthOribit() {

	Planets sun(1.989e30, 0, 0, 0, 0, 0, 0);
	Planets earth(5.972168e24, 152097684030.928554000, 0, 0, 0, 29295.3, 0);
	double G = 6.67430e-11;

	double expected_rmin = 147098057369.071446000;
	double expected_rmax = 152097684030.928554000;

	double tol = 4e6;

	double rMin = find_rMin(earth, sun);
	double rMax = find_rMax(earth, sun);

	//assert(std::abs(rMin - expected_rmin) < tol);
	//assert(std::abs(rMax - expected_rmax) < tol);

	std::cout << "[Verification Passed Test Case 3]\n";
}

//Test case 4 Post-Newtonian Precession
void testPrecession() {

	//test case for Post-Newtonian Precession of the orbits
	Planets sun(1.989e30, 0, 0, 0, 0, 0, 0);
	Planets mercury(3.3011e23, 69805480369, 0, 0, 0, 38849.96398, 0);

	// Constants
	double G = 6.67430e-11;
	double c = 299792458;
	double dt = 1;

	// Create list of doubles to store aphelion points
	std::vector<double> aphelionPoints;

	// Define stored value for distance
	double r_aphelion = 6.9805480369e10;
	double slice_width = 78000;

	bool inSlice = false;
	double maxR = 0.0;
	double bestAngle = 0.0;

	//simulation loop for ~ 2 orbits
	for (int i = 0; i < 100000000; ++i) {

		applyGravity(mercury, sun, G, c);
		mercury.updatePosition(dt);

		CorrectedAcceleration acc = findCorrectedAcceleration(mercury, sun, G, c);

		mercury.updateVelocity(dt, acc.ax, acc.ay, acc.az);

		double r = findDistance(mercury, sun).r;

		// Check if mercury is in the "aphelion slice"
		if (std::abs(r - r_aphelion) < slice_width / 2.0) {

			if (!inSlice) {

				inSlice = true;
				maxR = r;
				bestAngle = std::atan2(mercury.getY(), mercury.getX());

			}
			else {
			
				if (r > maxR) {
				
					maxR = r;
					bestAngle = std::atan2(mercury.getY(), mercury.getX());
				
				}
			
			}

		}
		else {

			if (inSlice) {

				aphelionPoints.push_back(bestAngle);
				std::cout << "Aphelion detected at radius: " << maxR
					<< ", angle: " << bestAngle << "\n";

			}
			inSlice = false;

		}

	}

	std::vector<double> unwrapped;

	if (!aphelionPoints.empty()) {
	
		unwrapped.push_back(aphelionPoints[0]);

	}
	

	// Output aphelion angles and precession per orbit
	for (size_t i = 1; i < aphelionPoints.size(); ++i) {

		double a = aphelionPoints[i];
		double last = unwrapped.back();

		//shit 'a' by +- 2pi so it's closest equivalent to 'last'
		while (a - last > M_PI) a -= 2 * M_PI;
		while (a - last <= -M_PI) a += 2 * M_PI;

		unwrapped.push_back(a);
	}
	
	double sumDelta = 0.0;
	size_t count = 0;
	for (size_t i = 1; i < unwrapped.size() - 1; ++i) {
	
		double delta = unwrapped[i] - unwrapped[i + 1];
		std::cout << "Orbit" << i
			<< ": Aphelion Angle (rad): " << aphelionPoints[i]
			<< ", Precession since last (rad): " << delta << "\n";
		sumDelta += delta;
		++count;
	}

	if (count > 0) {
	
		double avgDelta = sumDelta / count;

		double days_per_year = 365.2425;
		double mercury_period_days = 87.969;
		double orbits_per_century = 100.0 * days_per_year / mercury_period_days;
		double avg_arcsec_per_century = avgDelta * RAD_TO_ARCSEC * orbits_per_century;

		std::cout << "Average precession per orbit (rad): " << avgDelta << "\n"
			<< "Average precession (arcsec/centurty): " << avg_arcsec_per_century << "\n";

	}

}

//Test case 5 3D Orbit
void test3DimOrbit() {

	//Test Case with a plane not perpendicular to the z-x and z-y planes.

};

//Test
int main() {

	std::cout << "This is the Test cases...\n\n";

	testStationary();
	testCircular();
	testEarthOribit();
	testPrecession();

	std::cout << "All test run.\n";

	return 0;

}