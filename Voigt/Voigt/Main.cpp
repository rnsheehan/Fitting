#ifndef ATTACH_H
#include "Attach.h"
#endif

// Command line approach for computing the Voigt spectrum and the Voigt HWHM
// There are two options 
// 1. Enter command Voigt fh fc fs fg flow fhigh Nsteps filename to compute a VOIGT spectrum based on Voigt parameters fh fc fs fg 
// in the range [flow, fhigh] with results stored in the file filename
// 2. Enter command Voigt fh fc fs fg filename to determine the 3dB HWHM and the 20dB down HW of the VOigt spectrum based on parameters fh fc fs fg
// Computed results are output to a file
// R. Sheehan 21 - 11 - 2022

// forward declarations
void Voigt(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na); 
void Voigt_HWHM(double xlow, double xhigh, std::vector<double>& a, int& na, double* HWHM, bool loud = false);
void Compute_HWHM(std::string& res_file, std::vector<double>& spctr_pars); 
void Compute_Spectrum(std::string& res_file, interval& spctr_rng, std::vector<double>& spctr_pars, void (*f)(double, std::vector<double>&, double*, std::vector<double>&, int&));
void parse_inputs(int argc, char* argv[], int &npars, std::vector<double> &spctr_pars, int &calc_type, interval &spctr_rng, std::string& res_file);

// main
int main(int argc, char* argv[])
{
	// run program from command line
	// there will be two options
	// 1. Enter command Voigt fh fc fs fg flow fhigh Nsteps filename to compute a VOIGT spectrum based on Voigt parameters fh fc fs fg 
	// in the range [flow, fhigh] with results stored in the file filename
	// 2. Enter command Voigt fh fc fs fg filename to determine the 3dB HWHM and the 20dB down HW of the VOigt spectrum based on parameters fh fc fs fg
	try {
		if (argc > 1) {
			int npars = 4; 
			int CALC = 0; 
			std::vector<double> Voigt_pars(npars, 0.0); // vector to hold the Voigt function parameter values
			interval plt_rng; 
			std::string filename = empty_str; 

			parse_inputs(argc, argv, npars, Voigt_pars, CALC, plt_rng, filename);

			if (CALC == HWHM_CALC) {
				std::cout << "Computing HWHM\n";
				Compute_HWHM(filename, Voigt_pars); 
			}
			else if (CALC == SPCTR_CALC) {
				std::cout << "Computing Spectrum\n";
				Compute_Spectrum(filename, plt_rng, Voigt_pars, Voigt);
			}
			else {
				std::cout << "Doing Nothing\n"; 
			}
		}
		else {
			std::string reason = "Error: Voigt-main\n";
			reason += "Not enough arguments were input\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}

	return 0; 
}

// function definitions
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

void parse_inputs(int argc, char* argv[], int& npars, std::vector<double>& spctr_pars, int &calc_type, interval& spctr_rng, std::string &res_file)
{
	// Assess and convert the inputs into a useful form
	// Decide if inputs are correct and which calculation to perform
	// R. Sheehan 22 - 11 - 2022

	// 1. Enter command Voigt fh fc fs fg flow fhigh Nsteps filename (9 inputs) to compute a VOIGT spectrum
	// 2. Enter command Voigt fh fc fs fg filename (6 inputs) to determine the 3dB HWHM and the 20dB down HW

	try {
		if (argc >= 6) {
			// List off the input parameters
			// Program needs 6 or more parameters to run, remember that the name of the program is also considered a parameter
			// argv[0] = program name
			std::cout << "Name of the program is " << argv[0] << ".exe\n";
			std::cout << argc - 1 << " parameters were input into the program\n";
			for (int count = 1; count < argc; count++) {
				std::cout << "argv[" << count << "] = " << argv[count] << "\n";
			}
			std::cout << "\n";

			// extract the spectrum parameters
			for (int i = 0; i < npars; i++) {
				spctr_pars[i] = atof( argv[i + 1] ); 
			}

			// Decide on the calculation type
			calc_type = ( argc == 9 ) ? SPCTR_CALC : HWHM_CALC;

			// assign the filename
			res_file = ( argc == 9 ) ? argv[8] : argv[5];

			// Set the plot range for the spectral plot
			if (calc_type == SPCTR_CALC) {
				spctr_rng.set_xl_xu( atoi(argv[7]), atof(argv[5]), atof(argv[6]) );
				std::cout << "Interval Properties: [" << spctr_rng.get_x_lower() << " , " << spctr_rng.get_x_upper() << "], dx = " << spctr_rng.get_dx() << "\n"; 
				std::cout << "Filename: " << res_file << "\n"; 
			}
			else {
				std::cout << "Computing HWHM\n"; 
				std::cout << "Filename: " << res_file << "\n";
			}			 
		}
		else {
			std::string reason = "Error: void assess_inputs(int argc, char* argv[], int& npars, std::vector<double>& spctr_pars, double& delta_f, int &calc_type)\n"; 
			reason += "Not enough arguments were input\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	 
}

void Compute_Spectrum(std::string& res_file, interval& spctr_rng, std::vector<double>& spctr_pars, void (*f)(double, std::vector<double>&, double*, std::vector<double>&, int&))
{
	// Compute the spectrum for the function defined in f
	// res_file is the name of the file that will contain the computed data
	// spctr_rng is the interval over which the spectrum will be computed
	// *f is a pointer to the function whose spectrum to be computed
	// R. Sheehan 22 - 11 - 2022

	try {
		bool c1 = res_file != empty_str ? true : false; 
		bool c2 = spctr_rng.has_bounds(); 
		bool c3 = spctr_pars[0] > 0.0 ? true : false;
		bool c4 = useful_funcs::valid_filename_length(res_file); 
		bool c10 = c1 && c2 && c3 && c4; 

		if (c10) {
			std::ofstream write(res_file, std::ios_base::out, std::ios_base::trunc);
			if (write.is_open()) {
				int npars = static_cast<int>(spctr_pars.size());
				std::vector<double> dyda(npars, 0.0);
				double fval = spctr_rng.get_x_lower();
				double delta_f = spctr_rng.get_dx(); 
				double spctr_val = 0.0; 
				for (int i = 0; i < spctr_rng.get_Nsteps(); i++) {
					Voigt(fval, spctr_pars, &spctr_val, dyda, npars); 
					write << std::setprecision(10) << fval << " , " << spctr_val << " , " << -dyda[1] << "\n"; 
					fval += delta_f; 
				}
				write.close(); 
			}
			else {
				std::string reason = "Error: void Compute_Spectrum(std::string& res_file, interval& spctr_rng, std::vector<double>& spctr_pars, void (*f)(double, std::vector<double>&, double*, std::vector<double>&, int&))\n";
				reason += "Cannot open file: " + res_file + "\n"; 
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason = "Error: void Compute_Spectrum(std::string& res_file, interval& spctr_rng, std::vector<double>& spctr_pars, void (*f)(double, std::vector<double>&, double*, std::vector<double>&, int&))\n"; 
			reason += "Incorrect input arguments\n"; 
			if (!c1) reason += "res_file is not defined\n"; 
			if (!c2) reason += "sweep interval is not defined\n"; 
			if (!c3) reason += "spectrum parameters is not defined\n"; 
			if (!c4) reason += "filename length is too long\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Voigt_HWHM(double xlow, double xhigh, std::vector<double>& a, int& na, double* HWHM, bool loud)
{
	// determine the HWHM of a Voigt function with the parameter set defined in a
	// use bisection method algorithm to look for the HWHM on the interval [xlow, xhigh]
	// a stores Voigt parameters a = { h, x_{centre}, g, sigma}
	// a[0] = h, a[1] = x_{centre}, a[2] = g, a[3] = sigma
	// h is an amplitude fitting factor
	// x_centre is the centre of the Voigt frunction
	// g is the half-width at half-maximum of the Lorentzian portion of Voigt
	// sigma is the std. dev. of the Gaussian portion of Voigt HWHM_{Gauss} = sqrt( 2 log(2) ) c
	// parameters in a will have been determined by Levenberg-Marquardt non-linear fit procedure
	// result will be stored in HWHM
	// R. Sheehan 1 - 12 - 2021

	try {

		bool c1 = xhigh > xlow ? true : false;
		bool c2 = xlow > 0 ? true : false;
		bool c3 = a[3] > 0 ? true : false;
		bool c10 = c1 && c2;

		if (c10) {
			int count = 0, MAXIT = 50; // max. no iterations of bisection method
			bool converged = false;
			double TOL = 1.0e-3; // desired accuracy of root computation
			double root = xlow, left = xlow, right = xhigh, dx = 0.0, fl = 0.0, fr = 0.0;
			double lor_gau = a[2] / a[3]; // g_{lor} / sigma_{gau}
			double Vmax_half = 0.5 * (a[0] * exp(template_funcs::DSQR(lor_gau)) * probability::erffc(lor_gau)); // half the peak value of the Voigt function

			if (loud)std::cout << "Peak value: " << 2.0 * Vmax_half << "\n";

			std::vector<double> dyda(na, 0.0);

			// subtract the peak value since you're looking for point where Voigt = 0.5*Vmax

			Voigt(left, a, &fl, dyda, na); fl -= Vmax_half; // compute the value of the function at the interval endpoints
			Voigt(right, a, &fr, dyda, na); fr -= Vmax_half; // compute the value of the function at the interval endpoints 
			fl = template_funcs::Signum(fl); fr = template_funcs::Signum(fr);

			// test the interval to ensure it contains a root ivttest == -1 => interval has root
			if ((fl * fr) < 0.0) {

				// the interval contains a root, search can proceed
				if (loud)std::cout << "Initial approximation to the root is " << root << ", fl = " << fl << " , fr = " << fr << "\n";

				// Compute MAXIT iterations of BisectRoot, stop if the desired tolerance is reached
				count = 0;
				while (count < MAXIT) {

					// Compute the amount by which the root position must be updated
					dx = 0.5 * (right - left);

					// Update the position of the root
					root = left + dx;

					if (loud)std::cout << "Iteration: " << count << ", root = " << root << ", fl = " << fl << " , fr = " << fr << "\n";

					// Test for convergence
					if (fabs(dx) < TOL) {
						if (loud)std::cout << "Bisection has converged to a root within tolerance " << TOL << " after " << count << " iterations\n";
						converged = true;
						break;
					}
					else {
						// Update the endpoints of the interval containing the root
						Voigt(root, a, &fr, dyda, na); fr -= Vmax_half; fr = template_funcs::Signum(fr);

						if ((fl * fr) > 0.0) {
							left = root;
							fl = fr;
						}
						else {
							right = root;
						}
					}

					count++;
				}

				if (converged) {
					if (loud)std::cout << "Root is located at " << root << "\n";
					*HWHM = fabs(root - a[1]);
				}
				else {
					if (loud)std::cout << "Bisection has not converged to a root within tolerance " << TOL << " after " << count << " iterations\n";
					*HWHM = fabs(root - a[1]);
				}
			}
			else {
				// the interval contains no root, search cannot proceed
				std::string reason;
				reason = "Error: Voigt_FWHM()\n";
				reason += "Search interval improperly defined\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: Voigt_FWHM()\n";
			if (!c1 || !c2) reason += "Search interval improperly defined\n";
			if (!c3) reason += "Voigt parameters improperly defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void Compute_HWHM(std::string& res_file, std::vector<double>& spctr_pars)
{
	// Compute the HWHM of the spectrum
	// R. Sheehan 22 - 11 - 2022

	try {
		bool c1 = res_file != empty_str ? true : false;
		bool c2 = spctr_pars[0] > 0.0 ? true : false;
		bool c3 = useful_funcs::valid_filename_length(res_file);
		bool c4 = spctr_pars[1] > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4;

		if (c10) {
			std::ofstream write(res_file, std::ios_base::out, std::ios_base::trunc);
			if (write.is_open()) {
				int npars = static_cast<int>(spctr_pars.size());
				double HWHM = 0.0;
				double f_span = 20; 
				Voigt_HWHM(spctr_pars[1], spctr_pars[1] + f_span, spctr_pars, npars, &HWHM);
				write << "Spectrum Parameters: "; 
				for (int i = 0; i < npars; i++) {
					write << std::setprecision(5) << spctr_pars[i] << " , ";
				}
				write << std::setprecision(5) << "HWHM: " << HWHM << "\n"; 				
				write.close(); 
			}
			else {
				std::string reason = "Error: void void Compute_HWHM(std::string& res_file, std::vector<double>& spctr_pars)\n";
				reason += "Cannot open file: " + res_file + "\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason = "Error: void Compute_HWHM(std::string& res_file, std::vector<double>& spctr_pars)\n";
			reason += "Incorrect input arguments\n";
			if (!c1) reason += "res_file is not defined\n";
			if (!c2 || !c3) reason += "spectrum parameters is not defined\n";
			if (!c4) reason += "filename length is too long\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}
