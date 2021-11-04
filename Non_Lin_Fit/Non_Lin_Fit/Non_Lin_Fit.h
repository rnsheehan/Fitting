// Non_Lin_Fit.h - Contains declarations of non-linear fitting functions
// This header file declares some functions perform various non-linear fits
#pragma once

#include <cstdlib>
#include <iostream> // cout, cin, cerr
#include <iomanip> // setw, setprecision, time
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

// preprocessor statements
// The new project template for a DLL project adds PROJECTNAME_EXPORTS to the defined preprocessor macros. 
// In this example, Visual Studio defines NON_LIN_FIT_EXPORTS when your NON_LIN_FIT DLL project is built.
// When the NON_LIN_FIT_EXPORTS macro is defined, the NON_LIN_FIT_API macro sets the __declspec(dllexport) modifier on the function declarations. 
// This modifier tells the compiler and linker to export a function or variable from the DLL for use by other applications.
// When NON_LIN_FIT_EXPORTS is undefined, for example, when the header file is included by a client application, NON_LIN_FIT_API applies the __declspec(dllimport) modifier to the declarations.
// This modifier optimizes the import of the function or variable in an application.

//#ifdef NON_LIN_FIT_EXPORTS
//#define NON_LIN_FIT_API __declspec(dllexport)
//#else
//#define NON_LIN_FIT_API __declspec(dllimport)
//#endif

//#define DllExport __declspec(dllexport)
//
//#define DllImport __declspec(dllimport)

// Function Declarations

// Useful Functions

//extern "C" NON_LIN_FIT_API 
double DSQR(double a);

double convert_dBm_to_mW(double dBm_val); 

double convert_mW_to_dBm(double mW_val); 

//extern "C" NON_LIN_FIT_API 
std::string toStringInt(const int& t);

//extern "C" NON_LIN_FIT_API 
std::string toStringDouble(const double& t, int places);

//extern "C" NON_LIN_FIT_API 
void SWAP(double& a, double& b);

// Vector utilities

//extern "C" NON_LIN_FIT_API 
bool array_2D_square(std::vector< std::vector< double > >& name); // test to see if an array is square

//extern "C" NON_LIN_FIT_API 
std::pair<int, int> array_2D_size(std::vector< std::vector< double > >& name); // return the size of an array as a pair

//extern "C" NON_LIN_FIT_API 
std::vector< std::vector< double > > array_2D(int& nrows, int& ncols); // create a 2D array of given size

// Linear Algebra
// Implementation of Gaussian elimination with partial pivoting 
//extern "C" NON_LIN_FIT_API 
void gaussj(std::vector< std::vector< double > >& a, int& n, std::vector< std::vector< double > >& b, int& m);

// Fitting Models

//extern "C" NON_LIN_FIT_API 
void Lorentzian(double x, std::vector<double>&a, double* y, std::vector<double>&dyda, int& na);

//extern "C" NON_LIN_FIT_API 
void Gaussian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na);

// Probability Functions

//extern "C" NON_LIN_FIT_API 
double gammln(double x); // Computes the value of ln[gamma(xx)] for xx>0

//extern "C" NON_LIN_FIT_API 
double gammq(double a, double x); // Computes the incomplete gamma function Q(a,x)=1-P(a,x) from the functions gser and gcf

//extern "C" NON_LIN_FIT_API 
void gser(double* gamser, double a, double x, double* gln); // Computes the incomplete gamma function P(a,x), calculated by its series representation

//extern "C" NON_LIN_FIT_API 
void gcf(double* gammcf, double a, double x, double* gln); // Computes the incomplete gamma function Q(a,x), calculated by its continued fraction representation

// Fitting Functions

// Function for expanding the storage of an array
//extern "C" NON_LIN_FIT_API
void covsrt(std::vector< std::vector< double > >& covar, int& ma, std::vector<int>& ia, int& mfit);

// Perform fit to non-linear function
//extern "C" NON_LIN_FIT_API 
void non_lin_fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& covar, std::vector<std::vector<double>>& alpha, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), int& itmax, double& toler, bool loud = false);

// Levenberg-Marquart Algorithm
//extern "C" NON_LIN_FIT_API 
void mrqmin(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& covar, std::vector<std::vector<double>>& alpha, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), double* alamda);

// Function for computing systems of equations used in LMA
//extern "C" NON_LIN_FIT_API 
void mrqcof(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& alpha, std::vector<double>& beta, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&));

// Function for evaluating goodness-of-fit statistics
//extern "C" NON_LIN_FIT_API 
void goodness_of_fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, int& ma, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), double* chisq, double* nu, double* rsqr, double* gof);

// Function for computing the residual of a fit
//extern "C" NON_LIN_FIT_API 
void residuals(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, int& ma, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), std::vector<std::vector<double>>& data);

// Actual functions that will be called by external user of DLL

// Fit a Lorentzian to a data set
extern "C" _declspec(dllexport) void Lorentz_Fit(int n_data, double freq_data[], double spctrm_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[]);

extern "C" _declspec(dllexport) void Testing();