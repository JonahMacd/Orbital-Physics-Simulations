#pragma once
#include "Planets.h"

// Define Structures to store data, makes it easier to return multiple values (like classes but less complex)
struct Distance {
	double r;
	double dx;
	double dy;
	double dz;
};

Distance findDistance(Planets& body1, Planets& body2);

struct AngularMomentum {
	double L;
	double Lx;
	double Ly;
	double Lz;
};

AngularMomentum findAngularMomentum(Planets& body, Planets& body2);

struct CorrectedAcceleration {
	double ax;
	double ay;
	double az;
};

CorrectedAcceleration findCorrectedAcceleration(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 299792458);

struct Velocity {
	double v;
	double vx;
	double vy;
	double vz;
};

Velocity findVelocity(Planets& body1, Planets& body2);

struct SpecificAngularMomentum {
	double h;
	double hx;
	double hy;
	double hz;
};

SpecificAngularMomentum findSpecificAngularMomentum(Planets& body1, Planets& body2);

double findMu(Planets& body1, Planets& body2, double G = 6.67430e-11);

double findSpecificEnergy(Planets& body1, Planets& body2, double G = 6.67430e-11);

void applyGravity(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 299792458);

double findEccentricity(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 299792458);

double findSemiMajorAxis(Planets& body1, Planets& body2, double G = 6.67430e-11);

double find_rMax(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 299792458);

double find_rMin(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 299792458);