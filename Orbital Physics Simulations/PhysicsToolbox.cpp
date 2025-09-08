#include "Physics.h"
#include "Planets.h"
#include <cmath>

Planets dummy(1, 0, 0, 0, 0, 0, 0);

double findDistance(Planets& body1, Planets& body2) {

	//Compute radius vector from body 1 to body 2
	double dx = body2.getX() - body1.getX();
	double dy = body2.getY() - body1.getY();
	double dz = body2.getZ() - body1.getZ();
	double r = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));

	// return distance values
	return r, dx, dy, dz;
};

double findAngularMomentum(Planets& body) {

	// Initilize (oooh fancy word) variables (trying to be object oriented :/ )
	double r, dx, dy, dz = findDistance(body, dummy);
	double v, vx, vy, vz = findVelocity(body);

	// calculate angular momentum paramaters
	double Lx = dy * vz - dz * vy;
	double Ly = dz * vx - dx * vz;
	double Lz = dx * vy - dy * vx;

	double L = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	// return angular momentum values
	return L, Lx, Ly, Lz;
};

double findCorrectedAcceleration(Planets& body1, Planets& body2, double G, double c) {

	// Initilize variables again
	double r, dx, dy, dz = findDistance(body1, body2);
	double L, Lx, Ly, Lz = findAngularMomentum(body1);

	// General relativistic correctoin factor
	double correctionFactor = 1 + (3 * L * L) / ((r * r) * (c * c));

	// calculate acceleration values
	double ax = G * body2.getMass() * dx / (r * r * r) * correctionFactor;
	double ay = G * body2.getMass() * dy / (r * r * r) * correctionFactor;
	double az = G * body2.getMass() * dz / (r * r * r) * correctionFactor;

	// return acceleration values
	return ax, ay, az;
};

double findVelocity(Planets& body1) {

	// same as all the other fucntions but is velocity instead

	double vx = body1.getVx();
	double vy = body1.getVy();
	double vz = body1.getVz();

	double v = std::sqrt((vx * vx) + (vy * vy) + (vz * vz));

	return v, vx, vy, vz;
};

double findSpecificAngularMomentum(Planets& body1, Planets& body2) {

	//specific angluar momentum (I wont lie, I dont know how this is different from normal angular momentum but okay)
	double r, dx, dy, dz = findDistance(body1, body2);
	double v, vx, vy, vz = findVelocity(body1);

	double hx = dy * vz - dz * vy;
	double hy = dz * vx - dx * vz;
	double hz = dx * vy - dy * vx;

	double h = std::sqrt((hx * hx) + (hy * hy) + (hz * hz));

	return h, hx, hy, hz;
};

double findMu(Planets& body1, Planets& body2, double G) {

	// this is another weird thing of orbital physics
	double m1 = body1.getMass();
	double m2 = body2.getMass();

	double mu = G * (m1 + m2);

	return mu;
};

double findSpecificEnergy(Planets& body1, Planets& body2) {

	// specific energy also confuses me but this is a function to find it, I can calculate it but can't explain it
	double r, dx, dy, dz = findDistance(body1, body2);
	double v, vx, vy, vz = findVelocity(body1);
	double mu = findMu(body1, body2);

	double epsilon = ((v * v) / 2) - (mu / r);

	return epsilon;
};

double findEccentricity(Planets& body1, Planets& body2, double G, double c) {

	// calculates how non-circle your circles are
	double h, hx, hy, hz = findSpecificAngularMomentum(body1, body2);
	double mu = findMu(body1, body2, G);
	double epsilon = findSpecificEnergy(body1, body2);

	double eccentricity = std::sqrt(1 + (2 * epsilon * h * h) / (mu * mu));

	return eccentricity;

};

double findSemiMajorAxis(Planets& body1, Planets& body2) {

	//semi majour axis or your "a" value in your orbit
	double epsilon = findSpecificEnergy(body1, body2);
	double mu = findMu(body1, body2);

	double a = -mu / (2 * epsilon);

	return a;

};

double find_rMax(Planets& body1, Planets& body2) {
	
	//maximum radius of ellipse
	double a = findSemiMajorAxis(body1, body2);
	double e = findEccentricity(body1, body2);

	double rMax = a * (1 + e);

	return rMax;
};

double find_rMin(Planets& body1, Planets& body2) {

	// minimum radius of ellipse
	double a = findSemiMajorAxis(body1, body2);
	double e = findEccentricity(body1, body2);

	double rMin = a * (1 - e);

	return rMin;
};

void applyGravity(Planets& body1, Planets& body2, double G, double c) {

	//Compute radius vectors

		//compute radius vector from body1 to body2
	double r, dx, dy, dz = findDistance(body1, body2);

	//Compute angular momentum per unit mass

	double L, Lx, Ly, Lz = findAngularMomentum(body1);

	//Compute Newtonian and Post-Newtonian factor

	double correctionFactor = 1 + (3 * L * L) / ((r * r) * (c * c));

	//Compute accelerations

	double ax, ay, az = findCorrectedAcceleration(body1, body2, G, c);

	//Apply acceleration as a force

	body1.applyForce(ax * body1.getMass(), ay * body1.getMass(), az * body1.getMass());
};