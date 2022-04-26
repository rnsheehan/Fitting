// Non_Lin_Fit.h - Contains declarations of non-linear fitting functions
// This header file declares some functions perform various non-linear fits
#pragma once

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

double Signum(double a); 

double convert_dBm_to_mW(double dBm_val); 

double convert_mW_to_dBm(double mW_val); 

std::string toStringInt(const int& t);

std::string toStringDouble(const double& t, int places);

void SWAP(double& a, double& b);

// Vector utilities
bool array_2D_square(std::vector< std::vector< double > >& name); // test to see if an array is square

std::pair<int, int> array_2D_size(std::vector< std::vector< double > >& name); // return the size of an array as a pair

std::vector< std::vector< double > > array_2D(int& nrows, int& ncols); // create a 2D array of given size

// Linear Algebra
// Implementation of Gaussian elimination with partial pivoting 
void gaussj(std::vector< std::vector< double > >& a, int& n, std::vector< std::vector< double > >& b, int& m);

// Fitting Models
void Lorentzian(double x, std::vector<double>&a, double* y, std::vector<double>&dyda, int& na);

void Gaussian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na);

void Voigt(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na); 

void diode_voltage(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na);

void resistive_diode_voltage(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na); 

void Voigt_HWHM(double xlow, double xhigh, std::vector<double>& a, int& na, double* HWHM, bool loud = false); 

void Ring_Down(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na);

// Probability Functions
double gammln(double x); // Computes the value of ln[gamma(xx)] for xx>0

double gammp(double a, double x); // Computes the incomplete gamma function P(a,x) from the functions gser and gcf

double gammq(double a, double x); // Computes the incomplete gamma function Q(a,x)=1-P(a,x) from the functions gser and gcf

void gser(double* gamser, double a, double x, double* gln); // Computes the incomplete gamma function P(a,x), calculated by its series representation

void gcf(double* gammcf, double a, double x, double* gln); // Computes the incomplete gamma function Q(a,x), calculated by its continued fraction representation

double erffc(double x); // 1 - erf(x)

// Fitting Functions

// Function for expanding the storage of an array
void covsrt(std::vector< std::vector< double > >& covar, int& ma, std::vector<int>& ia, int& mfit);

// Perform fit to non-linear function
void non_lin_fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& covar, std::vector<std::vector<double>>& alpha, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), int& itmax, double& toler, bool loud = false);

// Levenberg-Marquart Algorithm
void mrqmin(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& covar, std::vector<std::vector<double>>& alpha, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), double* alamda);

// Function for computing systems of equations used in LMA
void mrqcof(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& alpha, std::vector<double>& beta, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&));

// Function for evaluating goodness-of-fit statistics
void goodness_of_fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, int& ma, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), double* chisq, double* nu, double* rsqr, double* gof, bool loud = false);

// Function for computing the residual of a fit
void residuals(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, int& ma, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), std::vector<std::vector<double>>& data);

// Actual functions that will be called by external user of DLL

// Fit a Lorentzian to a data set
extern "C" _declspec(dllexport) void Lorentz_Fit(int n_data, double freq_data[], double spctrm_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[]);

// Fit a Gaussian to a data set
extern "C" _declspec(dllexport) void Gauss_Fit(int n_data, double freq_data[], double spctrm_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[]);

// Fit a Voigt profile to a data set
extern "C" _declspec(dllexport) void Voigt_Fit(int n_data, double freq_data[], double spctrm_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[], double *HWHM);

extern "C" _declspec(dllexport) void Diode_Fit(int n_data, double current_data[], double voltage_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[]);

extern "C" _declspec(dllexport) void Resistive_Diode_Fit(int n_data, double current_data[], double voltage_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[]);

extern "C" _declspec(dllexport) void Ring_Down_Fit(int n_data, double time_data[], double voltage_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[]);

extern "C" _declspec(dllexport) void Testing();