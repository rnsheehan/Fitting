#ifndef ATTACH_H
#include "Attach.h"
#endif

// Command line approach for computing the Lorentz spectrum and the Lorentz HWHM
// There are two options 
// 1. Enter command Lorentz fh fc fs fg flow fhigh Nsteps filename to compute a Lorentz spectrum based on Voigt parameters fh fc fs fg 
// in the range [flow, fhigh] with results stored in the file filename
// 2. Enter command Lorentz fh fc fs fg filename to determine the 3dB HWHM and the 20dB down HW of the Lorentz spectrum based on parameters fh fc fs fg
// Computed results are output to a file
// R. Sheehan 22 - 11 - 2022

// forward declarations
void Lorentzian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na);
void Compute_HWHM(std::string& res_file, std::vector<double>& spctr_pars); 
void Compute_Spectrum(std::string& res_file, interval& spctr_rng, std::vector<double>& spctr_pars, void (*f)(double, std::vector<double>&, double*, std::vector<double>&, int&));
void parse_inputs(int argc, char* argv[], int &npars, std::vector<double> &spctr_pars, int &calc_type, interval &spctr_rng, std::string& res_file);

// main
int main(int argc, char* argv[])
{
	// run program from command line
	// there will be two options
	// 1. Enter command Lorentz fh fc fs flow fhigh Nsteps filename to compute a VOIGT spectrum based on Voigt parameters fh fc fs
	// in the range [flow, fhigh] with results stored in the file filename
	// 2. Enter command Lorentz fh fc fs filename to determine the 3dB HWHM and the 20dB down HW of the VOigt spectrum based on parameters fh fc fs
	try {
		if (argc > 1) {
			int npars = 3; 
			int CALC = 0; 
			std::vector<double> Lorentz_pars(npars, 0.0); // vector to hold the Voigt function parameter values
			interval plt_rng; 
			std::string filename = empty_str; 

			parse_inputs(argc, argv, npars, Lorentz_pars, CALC, plt_rng, filename);

			if (CALC == HWHM_CALC) {
				std::cout << "Computing HWHM\n";
				Compute_HWHM(filename, Lorentz_pars);
			}
			else if (CALC == SPCTR_CALC) {
				std::cout << "Computing Spectrum\n";
				Compute_Spectrum(filename, plt_rng, Lorentz_pars, Lorentzian);
			}
			else {
				std::cout << "Doing Nothing\n"; 
			}
		}
		else {
			std::string reason = "Error: Lorentz-main\n";
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
void Lorentzian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Lorentzian function to be fitted
	// a stores Lorentzian parameters a = { A, x_{centre}, G/2}
	// a[0] = A, a[1] = x_{centre}, a[2] = G / 2
	// G/2 is the Lorentzian half-width at half-maximum (i.e. linewidth)
	// Lorentzian value is given by *y
	// This Lorentzian is not normalised, to normalise multiply *y by 1/pi
	// Normalisation not required for my purposes
	// dyda is array that stores value of derivative of Lorentzian function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 21 - 10 - 2021

	try {
		double t1 = x - a[1]; // ( x - x_{centre} )
		double t2 = (a[0] * a[2]); // A (G/2)
		double denom = (template_funcs::DSQR(t1) + template_funcs::DSQR(a[2])); // ( x - x_{centre} )^{2} + (G/2)^{2}
		*y = t2 / denom; // A (G/2) / [ ( x - x_{centre} )^{2} + (G/2)^{2} ]
		dyda[0] = (*y) / a[0]; // \partial L / \partial A
		dyda[1] = (2.0 * t1 * template_funcs::DSQR(*y)) / t1; // \partial L / \partial x_{centre}
		dyda[2] = (*y) * ((1.0 / a[2]) - ((2.0 * (*y)) / a[0])); // \partial L / \partial G
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

	// 1. Enter command Lorentz fh fc fs flow fhigh Nsteps filename (8 inputs) to compute a Lorentz spectrum
	// 2. Enter command Lorentz fh fc fs filename (5 inputs) to determine the 3dB HWHM and the 20dB down HW

	try {
		if (argc >= 5) {
			// List off the input parameters
			// Program needs 5 or more parameters to run, remember that the name of the program is also considered a parameter
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
			calc_type = ( argc == 8 ) ? SPCTR_CALC : HWHM_CALC;

			// assign the filename
			res_file = ( argc == 8 ) ? argv[7] : argv[4];

			// Set the plot range for the spectral plot
			if (calc_type == SPCTR_CALC) {
				spctr_rng.set_xl_xu( atoi(argv[6]), atof(argv[4]), atof(argv[5]) );
				std::cout << "Interval Properties: [" << spctr_rng.get_x_lower() << " , " << spctr_rng.get_x_upper() << "], dx = " << spctr_rng.get_dx() << "\n"; 
				std::cout << "Filename: " << res_file << "\n"; 
			}
			else {
				std::cout << "Computing HWHM\n"; 
				std::cout << "Filename: " << res_file << "\n";
			}			 
		}
		else {
			std::string reason = "Error: void parse_inputs(int argc, char* argv[], int& npars, std::vector<double>& spctr_pars, double& delta_f, int &calc_type)\n"; 
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
					Lorentzian(fval, spctr_pars, &spctr_val, dyda, npars);
					//write << std::setprecision(10) << fval << " , " << spctr_val << " , " << -dyda[1] << "\n"; 
					write << std::setprecision(10) << fval << " , " << spctr_val << "\n"; 
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
				double HWHM = spctr_pars[2];
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
