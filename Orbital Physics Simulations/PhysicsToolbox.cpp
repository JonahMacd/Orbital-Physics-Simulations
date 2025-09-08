#include "Physics.h"
#include <cmath>

void applyGravity(Planets& body1, Planets& body2, double G, double c) {

	//Compute radius vectors

		//compute radius vector for origin to body1
		double dx1 = body1.getX() - 0.0;
		double dy1 = body1.getY() - 0.0;
		double dz1 = body1.getZ() - 0.0;

		//compute radius vector for origin to body2
		double dx2 = body2.getX() - 0.0;
		double dy2 = body2.getY() - 0.0;
		double dz2 = body2.getZ() - 0.0;

		//compute radius vector from body1 to body2
		double dx = body2.getX() - body1.getX();
		double dy = body2.getY() - body1.getY();
		double dz = body2.getZ() - body1.getZ();
		double r = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));

	//Compute angular momentum per unit mass
		double vx1 = body1.getVx();
		double vy1 = body1.getVy();
		double vz1 = body1.getVz();

		double Lx = dy * vz1 - dz * vy1;
		double Ly = dz * vx1 - dx * vz1;
		double Lz = dx * vy1 - dy * vx1;

		double L = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	//Compute Newtonian and Post-Newtonian factor

		double correctionFactor = 1 + (3 * L * L) / ((r * r) * (c * c));

	//Compute accelerations

		double ax = G * body2.getMass() * dx / (r * r * r) * correctionFactor;
		double ay = G * body2.getMass() * dy / (r * r * r) * correctionFactor;
		double az = G * body2.getMass() * dz / (r * r * r) * correctionFactor;

	//Apply acceleration as a force

		body1.applyForce(ax * body1.getMass(), ay * body1.getMass(), az * body1.getMass());
}