#include "Planets.h"

Planets::Planets(double m, double xPOS, double yPOS, double zPOS, double Vx, double Vy, double Vz) 
	: mass(m), xPos(xPOS), yPos(yPOS), zPos(zPOS), vx(Vx), vy(Vy), vz(Vz), ax(0), ay(0), az(0) {}

void Planets::applyForce(double fx, double fy, double fz) {

	//Determine acceleration from force

	ax += fx / mass;
	ay += fy / mass;
	ax += fz / mass;

}

void Planets::update(double dt) {

	//Redefine all time dependent variables
	//Define velocity based on time passed dt (Remake later to a for loop for simulation animation for graphics)
	vx += ax * dt;
	vy += ay * dt;
	vz += az * dt;

	//Define position based on time passed dt (Remake later to a for loop for simulation animation for graphics)
	xPos += vx * dt;
	yPos += vy * dt;
	zPos += vz * dt;

	//Reset acceleration for next force application as it might change depending on multiple planets
	ax = ay = az = 0;

}

double Planets::getX() const { return xPos; }
double Planets::getY() const { return yPos; }
double Planets::getZ() const { return zPos; }
double Planets::getVx() const { return vx; }
double Planets::getVy() const { return vy; }
double Planets::getVz() const { return vz; }
double Planets::getMass() const { return mass; }