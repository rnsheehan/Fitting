#ifndef ATTACH_H
#include "Attach.h"
#endif

void Voigt(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na); 

int main(int argc, char* argv[])
{
	// run program from command line
	// assuming that inputs are going to be Voigt fval fh fc fs fg
	// program should return the value of the Voigt function
	try {
		if (argc > 1) {
			// List off the input parameters
			// Program needs 2 or more parameters to run, remember that the name of the program is also considered a parameter
			// argv[0] = program name

			std::cout << "Name of the program is " << argv[0] << ".exe\n";
			std::cout << argc - 1 << " parameters were input into the program\n";
			for (int count = 1; count < argc; count++) {
				std::cout << "argv[" << count << "] = " << argv[count] << "\n";
			}
			std::cout << "\n";

			int npars = 4; 
			double frq_val = atof( argv[1] ); 
			double Voigt_val = 0.0; 
			std::vector<double> Voigt_pars(npars, 0.0); // vector to hold the Voigt function parameter values
			std::vector<double> Voigt_ders(npars, 0.0); // vector to hold the computed function derivatives

			for (int i = 0; i < npars; i++) { 
				Voigt_pars[i] = atof( argv[i + 2] );  
			}

			Voigt(frq_val, Voigt_pars, &Voigt_val, Voigt_ders, npars);

			return Voigt_val; 
		}
		else {
			throw std::invalid_argument("Not enough arguments were input\n");
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
		return 0;
	}
}

void Voigt(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Voigt function to be fitted
	// a stores Voigt parameters a = { h, x_{centre}, g, sigma}
	// a[0] = h, a[1] = x_{centre}, a[2] = g, a[3] = sigma
	// h is an amplitude fitting factor
	// x_centre is the centre of the Voigt frunction
	// g is the half-width at half-maximum of the Lorentzian portion of Voigt
	// sigma is the std. dev. of the Gaussian portion of Voigt HWHM_{Gauss} = sqrt( 2 log(2) ) c
	// It should be possible to express HWHM_{Voigt} in terms of g and sigma
	// Voigt value is given by *y
	// dyda is array that stores value of derivative of Voigt function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 30 - 11 - 2021

	try {
		std::complex<double> z = ((x - a[1] + eye * a[2]) / a[3]); // z = (x - x_{0} + i g) / sigma
		std::complex<double> W = Faddeeva::w(z); // w(z) = exp(-z^{2}) Erfc(-i z)
		std::complex<double> dW = ((2.0 * eye) / SQRT_PI) - (2.0 * z * W); // dw / dz by definition
		double h_sig = a[0] / a[3]; // h / sigma

		*y = real(a[0] * W); // V = re( h w(z) )
		dyda[0] = real(W); // \partial V / \partial h
		dyda[1] = real(-1.0 * h_sig * dW); // \partial V / \partial x_{0} = - \partial V / \partial x
		dyda[2] = real(eye * h_sig * dW); // \partial V / \partial g
		dyda[3] = real(-1.0 * h_sig * z * dW); // \partial V / \partial sigma
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}