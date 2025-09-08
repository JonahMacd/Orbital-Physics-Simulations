#pragma once
#include "Planets.h"

double findDistance(Planets& body1, Planets& body2);

double findAngularMomentum(Planets& body);

double findCorrectedAcceleration(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 3e8);

double findVelocity(Planets& body1);

double findSpecificAngularMomentum(Planets& body1, Planets& body2);

double findMu(Planets& body1, Planets& body2, double G = 6.67430e-11);

double findSpecificEnergy(Planets& body1, Planets& body2);

void applyGravity(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 3e8);

double findEccentricity(Planets& body1, Planets& body2, double G = 6.67430e-11, double c = 3e8);

double findSemiMajorAxis(Planets& body1, Planets& body2);

double find_rMax(Planets& body1, Planets& body2);

double find_rMin(Planets& body1, Planets& body2);