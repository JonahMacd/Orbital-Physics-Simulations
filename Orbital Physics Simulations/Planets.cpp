#include "Planets.h"

Planets::Planets(double m, double xPOS, double yPOS, double zPOS, double Vx, double Vy, double Vz) 
	: mass(m), xPos(xPOS), yPos(yPOS), zPos(zPOS), vx(Vx), vy(Vy), vz(Vz), ax(0), ay(0), az(0) {}

void Planets::applyForce(double fx, double fy, double fz) {

	//Determine acceleration from force

	ax += fx / mass;
	ay += fy / mass;
	az += fz / mass;

}

void Planets::updatePosition(double dt) {

	//Redefine all time dependent variables

	//Define position based on time passed dt (Remake later to a for loop for simulation animation for graphics)
	xPos += vx * dt + (0.5 * ax * dt * dt);
	yPos += vy * dt + (0.5 * ay * dt * dt);
	zPos += vz * dt + (0.5 * az * dt * dt);

}

void Planets::updateVelocity(double dt, double ax_new, double ay_new, double az_new) {
	//Update velocity based on time passed dt (Remake later to a for loop for simulation animation for graphics)
	vx += 0.5 * (ax + ax_new) * dt;
	vy += 0.5 * (ay + ay_new) * dt;
	vz += 0.5 * (az + az_new) * dt;

	ax = ax_new;
	ay = ay_new;
	az = az_new;
}

double Planets::getX() const { return xPos; }
double Planets::getY() const { return yPos; }
double Planets::getZ() const { return zPos; }
double Planets::getVx() const { return vx; }
double Planets::getVy() const { return vy; }
double Planets::getVz() const { return vz; }
double Planets::getAx() const { return ax; }
double Planets::getAy() const { return ay; }
double Planets::getAz() const { return az; }
double Planets::getMass() const { return mass; }