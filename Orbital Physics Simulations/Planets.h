#pragma once

class Planets {

public:
	Planets(double Mass, double xPOS, double yPOS, double zPOS, double Vx, double Vy, double Vz);
	
	// Force and orbital values will be coded here, later on in a cpp file.    

	void setVelocity(double vx_new, double vy_new, double vz_new) {

		vx = vx_new;
		vy = vy_new;
		vz = vz_new;

	};
	void setPostition(double x_new, double y_new, double z_new) {
	
		xPos = x_new;
		yPos = y_new;
		zPos = z_new;

	};

	void updatePosition(double dt);
	void updateVelocity(double dt, double ax_new, double ay_new, double az_new);
	void applyForce(double fx, double fy, double fz);

	double getX() const;
	double getY() const;
	double getZ() const;
	double getVx() const;
	double getVy() const;
	double getVz() const;
	double getAx() const;
	double getAy() const;
	double getAz() const;
	double getMass() const;

private:
	double mass;
	double xPos, yPos, zPos;
	double vx, vy, vz;
	double ax, ay, az;

};
