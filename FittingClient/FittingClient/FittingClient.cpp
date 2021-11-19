// FittingClient.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Non_Lin_Fit.h"

void read_into_matrix(std::string& filename, std::vector<std::vector<double>>& data, int& n_rows, int& n_cols, bool loud); 
void write_into_file(std::string& filename, std::vector<double>& data, bool loud = false); 
std::vector< std::vector< double > > array_2D(int& nrows, int& ncols); // create a 2D array of given size
void Lorentzian_data_fit_test();

int main()
{
	Lorentzian_data_fit_test(); 
}

void read_into_matrix(std::string& filename, std::vector<std::vector<double>>& data, int& n_rows, int& n_cols, bool loud)
{
	// read an array of data from a file
	// store the data in a matrix of size n_rows * n_cols
	// R. Sheehan 18 - 12 - 2018

	try {
		std::ifstream the_file;
		the_file.open(filename, std::ios_base::in);

		if (the_file.is_open()) {

			if (loud) std::cout << filename << " opened for reading\n";

			if (!data.empty()) data.clear(); // make sure that data is empty when values are being read in

			bool has_header = false;
			std::string line, item;

			char endline = '\n';
			char tab_token = '\t';
			char comma_token = ',';
			char sep_token;

			// Initialise the values to zero
			n_rows = 0; n_cols = 0;

			// Determine which token separates data in the file
			the_file.seekg(1, std::ios::beg); // move to the start of the file
			std::getline(the_file, line, endline);
			if (line.find(tab_token) != std::string::npos) sep_token = tab_token;
			if (line.find(comma_token) != std::string::npos) sep_token = comma_token;

			if (loud) std::cout << filename << " uses " << sep_token << " as a token\n";

			// Determine if first line is a file header
			if (isalpha(line[2])) has_header = true;

			// Count the number of rows and columns
			// This only seems to work when the data are separated by ',' also works for '\t' and ' '
			// http://www.cplusplus.com/reference/string/string/getline/
			// getline (istream& is, string& str, char delim)
			// Extracts characters from is and stores them into str until the delimitation character delim is found
			the_file.seekg(has_header ? 1 : 0, std::ios::beg); // move to the start of the file
			while (std::getline(the_file, line, endline)) {
				n_rows++;
				std::istringstream linestream(line);
				if (n_rows == 1) {
					while (std::getline(linestream, item, sep_token)) {
						n_cols++;
					}
				}
			}

			if (loud) std::cout << filename << " contains " << n_rows << " rows and " << n_cols << " columns\n";

			if (n_rows > 1 && n_cols > 0) {
				// Allocate the memory required to store the data
				data.resize(n_rows);
				for (size_t k = 0; k < data.size(); k++) {
					data[k].resize(n_cols, 0.0);
				}

				the_file.clear(); // empty a buffer, needed to ensure data can be read from the file
				the_file.seekg(has_header ? 1 : 0, std::ios::beg); // move to the start of the file

				int i, j;

				i = 0;
				while (std::getline(the_file, line, endline)) {
					std::istringstream linestream(line);
					j = 0;
					while (std::getline(linestream, item, sep_token)) {
						data[i][j] = atof(item.c_str());
						j++;
					}
					i++;
				}

				the_file.close();
			}
			else {
				std::string reason;
				reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
				reason = filename + " contains no data\n";
				//reason += "n_rows: " + template_funcs::toString(n_rows) + ", n_cols: " + template_funcs::toString(n_cols) + "\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
			reason += "Cannot open: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}

void write_into_file(std::string& filename, std::vector<double>& data, bool loud)
{
	// write a single column of data to a new file

	try {

		if (!data.empty()) {

			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {

				for (size_t k = 0; k < data.size(); k++) {
					write << std::setprecision(10) << data[k] << "\n";
				}

				write.close();
			}
			else {
				std::string reason;
				reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
				reason += "Could not open file: " + filename + "\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
			reason += "Filename: " + filename + " is not valid\n";
			throw std::runtime_error(reason);
		}
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}

std::vector< std::vector< double > > array_2D(int& nrows, int& ncols)
{
	// create an array of zeroes of size nrows*ncols
	// R. Sheehan 22 - 6 - 2018

	try {
		if (nrows > 0 && ncols > 0) {
			std::vector< std::vector< double > > name;
			name.resize(nrows);
			for (int i = 0; i < nrows; i++) name[i].resize(ncols, 0.0);
			return name;
		}
		else {
			std::vector< std::vector< double > > temp;
			return temp;
			std::string reason = "Error: std::vector< std::vector< double > > lin_alg::array_2D\n";
			//if (nrows <= 1) reason += "nrows = " + template_funcs::toString(nrows) + " too small\n";
			//if (ncols <= 1) reason += "ncols = " + template_funcs::toString(ncols) + " too small\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what(); 
	}
}

void Lorentzian_data_fit_test()
{
	// Apply LM method to measured linewidth spectrum data
	// R. Sheehan 26 - 10 - 2021

	// Read in the measured spectral data
	//std::string filename = "Sample_LLM.csv"; 
	//std::string filename = "Lorentz_iodeal.csv"; // this is the same data set as Sample_LLM.csv
	std::string filename = "Smpl_LLM_8.txt";

	int npts, n_rows, npars = 3, n_cols, indx_max = 0;
	long idum = (-1011);
	double spread = 0.1, spctr_max = -500.0, f_max = 0, f_start, f_end, scale_fac;

	std::vector<std::vector<double>> the_data;

	read_into_matrix(filename, the_data, n_rows, n_cols, true);

	std::vector<double> xdata;
	std::vector<double> ydata;

	scale_fac = 1.0e+6;  f_start = 70.0; f_end = 90.0;
	for (int i = 0; i < n_rows; i++) {
		if (the_data[i][0] > f_start && the_data[i][0] < f_end) {
			xdata.push_back(the_data[i][0]);
			ydata.push_back(the_data[i][1]); 
		}
	}

	npts = static_cast<int>(xdata.size());

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	a_guess[0] = 10;
	a_guess[1] = 80;
	a_guess[2] = 1; // initial guesses for the parameters

	// run the fitting algorithm

	// declarate the arrays
	int n_stats = 4; 
	double *freq_data = NULL;
	double * spctrm_data = NULL;
	double * fit_data = NULL;
	double a_pars[] = {a_guess[0], a_guess[1], a_guess[2]};
	double gof_stats[] = {0, 0, 0, 0};

	freq_data = new double[npts]; 
	spctrm_data = new double[npts];
	fit_data = new double[npts];

	for (int i = 0; i < npts; i++) {
		freq_data[i] = xdata[i];
		spctrm_data[i] = ydata[i];
		fit_data[i] = 0.0; 
	}

	Lorentz_Fit(npts, freq_data, spctrm_data, fit_data, npars, a_pars, n_stats, gof_stats);
	//Gauss_Fit(npts, freq_data, spctrm_data, fit_data, npars, a_pars, n_stats, gof_stats);

	std::cout << "\nFitted centre freq: " << a_pars[1] << " MHz\n";
	std::cout << "Computed peak val: " << a_pars[0] / a_pars[2] << " uW\n";
	std::cout << "Computed HWHM: " << a_pars[2] << " MHz\n";
	std::cout << "\nGoodness-of-fit statistics\n";
	std::cout << "The chi-sq value for the fit is " << gof_stats[0] << "\n";
	std::cout << "nu for the fit is " << gof_stats[0] / gof_stats[1] << "\n";
	std::cout << "chi-sq / nu = " << gof_stats[1] << "\n";
	std::cout << "goodness of fit = " << gof_stats[3] << "\n";
	std::cout << "coefficient of determination = " << gof_stats[2] << "\n\n";

	// see computed fit values
	filename = "Computed_Fit_Values.txt"; 

	std::vector<double> fit_vals; 
	for (int i = 0; i < npts; i++)fit_vals.push_back(fit_data[i]); 

	write_into_file(filename, fit_vals, true); 

	// deallocate memory

	the_data.clear(); xdata.clear(); ydata.clear(); 

	delete[] freq_data; delete[] spctrm_data; delete[] fit_data;
}