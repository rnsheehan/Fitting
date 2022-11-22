#ifndef PROBABILITY_FUNCTIONS_H
#define PROBABILITY_FUNCTIONS_H

// Declaration of a namespace that contains functions used to compute probabilities in Statistics
// R. Sheehan 5 - 9 - 2017

namespace probability{

	// Gamma Function
	double gammln(double x); // Computes the value of ln[gamma(xx)] for xx>0

	// Incomplete Gamma Function
	double gammp(double a, double x); // Computes the incomplete gamma function P(a,x) from the functions gser and gcf
	double gammq(double a, double x); // Computes the incomplete gamma function Q(a,x)=1-P(a,x) from the functions gser and gcf

	void gser(double *gamser,double a,double x,double *gln); // Computes the incomplete gamma function P(a,x), calculated by its series representation
	void gcf(double *gammcf,double a,double x,double *gln); // Computes the incomplete gamma function Q(a,x), calculated by its continued fraction representation
	
	// Error Function
	double erff(double x); // erf(x)
	double erffc(double x); // 1 - erf(x)
	double erfcc(double x); // concise routine for 1-erf(x)

	// Factorial Function
	double factorial(int n);

	// Incomplete Beta Function
	double betacf(double a, double b, double x);
	double betai(double a, double b, double x); 

	// Kolmorgorov Smirnov Probability Function
	double probks(double alam); 
}

#endif