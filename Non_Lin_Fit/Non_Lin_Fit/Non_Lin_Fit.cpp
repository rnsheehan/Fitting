// Add an implementation to the DLL
// MathLibrary.cpp : Defines the exported functions for the DLL.
#include "pch.h" // use stdafx.h in Visual Studio 2017 and earlier
#include <utility>
#include <limits.h>
#include "Non_Lin_Fit.h"

// Function Definitions for the Functions declared in "mathLibrary.h"
double DSQR(double a)
{
	// Efficient squaring operator
	// Write injuries in dust, benefits in marble

	double darg;
	return ((darg = (a)) == static_cast<double>(0) ? static_cast<double>(0) : darg * darg);
}

void Lorentzian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Lorentzian function to be fitted
	// a stores Lorentzian parameters a = { A, x_{centre}, G/2}
	// a[0] = A, a[1] = x_{centre}, a[2] = G / 2
	// G/2 is the Lorentzian half-width at half-maximum (i.e. linewidth)
	// Lorentzian value is given by *y
	// dyda is array that stores value of derivative of Lorentzian function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 21 - 10 - 2021

	double t1 = x - a[1]; // ( x - x_{centre} )
	double t2 = (a[0] * a[2]); // A (G/2)
	double denom = (DSQR(t1) + DSQR(a[2])); // ( x - x_{centre} )^{2} + (G/2)^{2}
	*y = t2 / denom; // A (G/2) / [ ( x - x_{centre} )^{2} + (G/2)^{2} ]
	dyda[0] = (*y) / a[0]; // \partial L / \partial A
	dyda[1] = (2.0 * t1 * DSQR(*y)) / t1; // \partial L / \partial x_{centre}
	dyda[2] = (*y) * ((1.0 / a[2]) - ((2.0 * (*y)) / a[0])); // \partial L / \partial G
}

void Gaussian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Gaussian function to be fitted
	// a stores Gaussian parameters a = { amplitude, mean, standard deviation} = {A, b, c}
	// a[0] = A, a[1] = b, a[2] = c
	// c is related the Gaussian half-width at half-maximum (i.e. linewidth) HWHM = sqrt{ 2 log(2) } c ~ 1.17741 c
	// Gaussian value is given by *y
	// dyda is array that stores value of derivative of Lorentzian function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 27 - 10 - 2021
	
	double t1 = x - a[1]; // ( x - b )
	double t1sqr = DSQR(t1); // ( x - b )^{2}
	double csqr = DSQR(a[2]); // c^{2}
	double arg = (-1.0 * t1sqr) / (2.0 * csqr); // -( x - b )^{2} / 2 c^{2}
	*y = a[0] * exp(arg); // A exp( -( x - b )^{2} / 2 c^{2} )
	dyda[0] = (*y) / a[0]; // \partial G / \partial A
	dyda[1] = (t1 / csqr) * (*y); // \partial G / \partial b
	dyda[2] = (t1sqr / (a[2] * csqr)) * (*y); // \partial G / \partial c	
}