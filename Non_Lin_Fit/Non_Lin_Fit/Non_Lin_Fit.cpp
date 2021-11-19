// Add an implementation to the DLL
// Non_Lin_Fit.cpp : Defines the exported functions for the DLL.
#include "pch.h" // use stdafx.h in Visual Studio 2017 and earlier


// Function Definitions for the Functions declared in "Non_Lin_Fit.h"

// Useful Functions

double DSQR(double a)
{
	// Efficient squaring operator
	// Write injuries in dust, benefits in marble

	double darg;
	return ((darg = (a)) == static_cast<double>(0) ? static_cast<double>(0) : darg * darg);
}

double convert_dBm_to_mW(double dBm_val)
{
	// convert a dBm power reading to a mW power reading
	return pow( 10.0, (dBm_val / 10.0) ); 
}

double convert_mW_to_dBm(double mW_val)
{
	// convert a mW power reading to dBm power reading
	if (mW_val > 1.0e-9) {
		return 10.0 * log10(mW_val); 
	}
	else {
		return -90.0; 
	}
}   

std::string toStringInt(const int& t)
{
	// This is just too convenient not to use
	// Is there a version that can include something similar to %0.5d ? 
	// There has to be, look into the setw method for strings
	// Requires the string-stream (sstream) standard library
	// R. Sheehan 16 - 5 - 2011

	std::ostringstream oss; // create a stream
	oss << t;				// insert value to stream
	return oss.str();		// return as a string
}

std::string toStringDouble(const double& t, int places)
{
	// toString function that allows for the
	// number of decimal places to be specified 
	// far too convenient
	// R. Sheehan 17 - 5 - 2011

	std::ostringstream oss; // create a stream

	oss << std::fixed << std::setprecision(places) << t; // insert value to stream

	return oss.str(); // return as a string
}

void SWAP(double& a, double& b)
{
	// Updated version of SWAP Macro

	double itemp = (a); (a) = (b); (b) = itemp;
}

// Vector utilities

bool array_2D_square(std::vector< std::vector< double > >& name)
{
	// test the 2D array to see if it is square
	// R. Sheehan 22 - 6 - 2018

	std::pair<int, int> sizes = array_2D_size(name);

	return sizes.first == sizes.second ? true : false;
}

std::pair<int, int> array_2D_size(std::vector< std::vector< double > >& name)
{
	// return a pair structure of the form (nrows, ncols)
	// telling you the size of the array
	// R. Sheehan 22 - 6 - 2018

	try {
		if ((int)(name.size()) > 0 && (int)(name[0].size()) > 0) {
			return std::pair<int, int>((int)(name.size()), (int)(name[0].size()));
		}
		else {
			return std::pair<int, int>(-1, -1);
			std::string reason = "Error: std::pair<int, int> lin_alg::array_2D_size()\n";
			reason += "Input array dimensions not assigned correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		return std::pair<int, int>(-1, -1);
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
			std::string reason = "Error: std::vector< std::vector< double > > lin_alg::array_2D\n";
			if (nrows <= 1) reason += "nrows = " + toStringInt(nrows) + " too small\n";
			if (ncols <= 1) reason += "ncols = " + toStringInt(ncols) + " too small\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::vector< std::vector< double > > temp;
		return temp;
		std::cerr << e.what();
	}
}

// Linear Algebra
void gaussj(std::vector< std::vector< double > >& a, int& n, std::vector< std::vector< double > >& b, int& m)
{
	// Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
	// is the input matrix.b[1..n][1..m] is input containing the m right - hand side vectors. On
	// output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
	// vectors

	try {
		std::pair<int, int> sizea, sizeb;

		sizea = array_2D_size(a); sizeb = array_2D_size(b);

		bool c1 = sizea.first == sizea.second ? true : false; // test for squareness of a
		bool c2 = sizea.second == sizeb.first ? true : false; // columns of A equal to rows of B? 

		if (c1 && c2) {
			int i, icol, irow, j, k, l, ll;
			double big, dum, pivinv;

			std::vector<int> indxc(n, 0);
			std::vector<int> indxr(n, 0);
			std::vector<int> ipiv(n, 0);

			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (fabs(a[j][k]) >= big) {
									big = fabs(a[j][k]);
									irow = j;
									icol = k;
								}
							}
							else if (ipiv[k] > 1) {
								std::string reason = "Error: void lin_alg::gaussj()\nSingular Matrix-1\n";
								throw std::runtime_error(reason);
							}
						}
				++(ipiv[icol]);

				if (irow != icol) {
					for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l]);
					for (l = 0; l < m; l++) SWAP(b[irow][l], b[icol][l]);
				}

				indxr[i] = irow;

				indxc[i] = icol;

				if (a[icol][icol] == 0.0) {
					std::string reason = "Error: void lin_alg::gaussj()\nSingular Matrix-2\n";
					throw std::runtime_error(reason);
				}

				pivinv = 1.0 / a[icol][icol];

				a[icol][icol] = 1.0;

				for (l = 0; l < n; l++) a[icol][l] *= pivinv;

				for (l = 0; l < m; l++) b[icol][l] *= pivinv;

				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
						for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
					}
			}

			// End of Elimination procedure, unscramble row permuations
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						SWAP(a[k][indxr[l]], a[k][indxc[l]]);
			}

			ipiv.clear(); indxr.clear(); indxc.clear();
		}
		else {
			std::string reason = "Error: void lin_alg::gaussj()\n";
			reason += "Input array dimensions are not correct\n";
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

// Fitting Models

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

	try {
		double t1 = x - a[1]; // ( x - x_{centre} )
		double t2 = (a[0] * a[2]); // A (G/2)
		double denom = (DSQR(t1) + DSQR(a[2])); // ( x - x_{centre} )^{2} + (G/2)^{2}
		*y = t2 / denom; // A (G/2) / [ ( x - x_{centre} )^{2} + (G/2)^{2} ]
		dyda[0] = (*y) / a[0]; // \partial L / \partial A
		dyda[1] = (2.0 * t1 * DSQR(*y)) / t1; // \partial L / \partial x_{centre}
		dyda[2] = (*y) * ((1.0 / a[2]) - ((2.0 * (*y)) / a[0])); // \partial L / \partial G
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void Gaussian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Gaussian function to be fitted
	// a stores Gaussian parameters a = { amplitude, mean, standard deviation} = {A, b, c}
	// a[0] = A, a[1] = b, a[2] = c
	// c is related to the Gaussian half-width at half-maximum (i.e. linewidth) HWHM = sqrt{ 2 log(2) } c ~ 1.17741 c
	// Gaussian value is given by *y
	// dyda is array that stores value of derivative of Lorentzian function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 27 - 10 - 2021
	
	try {
		double t1 = x - a[1]; // ( x - b )
		double t1sqr = DSQR(t1); // ( x - b )^{2}
		double csqr = DSQR(a[2]); // c^{2}
		double arg = (-1.0 * t1sqr) / (2.0 * csqr); // -( x - b )^{2} / 2 c^{2}
		*y = a[0] * exp(arg); // A exp( -( x - b )^{2} / 2 c^{2} )
		dyda[0] = (*y) / a[0]; // \partial G / \partial A
		dyda[1] = (t1 / csqr) * (*y); // \partial G / \partial b
		dyda[2] = (t1sqr / (a[2] * csqr)) * (*y); // \partial G / \partial c
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

// Probability Functions

double gammln(double xx)
{
	// Gamma Function
	// Return the value of ln[gamma(xx)] for xx>0

	try {
		if (xx > 0.0) {
			double x, tmp, ser;
			static double cof[6] = { 76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5 };
			int j;

			x = xx - 1.0;
			tmp = x + 5.5;
			tmp -= (x + 0.5) * log(tmp);
			ser = 1.0;
			for (j = 0; j <= 5; j++) {
				x += 1.0;
				ser += (cof[j] / x);
			}
			return -tmp + log(2.50662827465 * ser);
		}
		else {
			std::string reason;
			reason = "Error: double gammln(double xx)\n";
			reason += "xx input with value = " + toStringDouble(xx, 2) + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		return 0.0; 
		std::cerr << e.what();
	}
}

double gammq(double a, double x)
{
	// Incomplete Gamma Function
	// Computes the incomplete gamma function Q(a,x)=1-P(a,x) from the functions gser and gcf

	try {
		if (x > 0.0 && a > 0.0) {
			double gamser, gammcf, gln;
			if (x < (a + 1.0)) {
				gser(&gamser, a, x, &gln);
				return 1.0 - gamser;
			}
			else {
				gcf(&gammcf, a, x, &gln);
				return gammcf;
			}
		}
		else {
			std::string reason;
			reason = "Error: double gammq(double a,double x)\n";
			if (x < 0.0) reason += "x input with value = " + toStringDouble(x, 2) + "\n";
			if (a <= 0.0) reason += "a input with value = " + toStringDouble(a, 2) + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		return 0.0; 
		std::cerr << e.what();
	}
}

void gser(double* gamser, double a, double x, double* gln)
{
	// This function returns the incomplete gamma function P(a,x), calculated by its series representation
	// Also returns ln[gamma(a)] as gln

	try {
		if (x > 0.0 && a > 0.0) {

			int n;
			double sum, del, ap;

			static const int ITMAX = (150);
			static const double EPS = (3.0e-7);

			*gln = gammln(a);
			if (x <= 0.0) {
				if (x < 0.0)
					std::cerr << "x less than 0 in routine GSER" << "\n";
				*gamser = 0.0;
				return;
			}
			else {
				ap = a;
				del = sum = 1.0 / a;
				for (n = 1; n <= ITMAX; n++) {
					ap += 1.0;
					del *= x / ap;
					sum += del;
					if (fabs(del) < fabs(sum) * EPS) {
						*gamser = sum * exp(-x + a * log(x) - (*gln));
						return;
					}
				}
				std::cerr << "a too large, ITMAX too small in routine GSER" << "\n";
				return;
			}
		}
		else {
			std::string reason;
			reason = "Error: void gser(double *gamser,double a,double x,double *gln)\n";
			if (x < 0.0) reason += "x input with value = " + toStringDouble(x, 2) + "\n";
			if (a <= 0.0) reason += "a input with value = " + toStringDouble(a, 2) + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void gcf(double* gammcf, double a, double x, double* gln)
{
	// This function returns the incomplete gamma function Q(a,x), calculated by its continued fraction representation
	// Also returns ln[gamma(a)] as gln

	try {
		if (x > 0.0 && a > 0.0) {
			int n;
			double gold = 0.0, g, fac = 1.0, b1 = 1.0;
			double b0 = 0.0, anf, ana, an, a1, a0 = 1.0;

			static const int ITMAX = (100);
			static const double EPS = (3.0e-7);

			*gln = gammln(a);
			a1 = x;
			for (n = 1; n <= ITMAX; n++) {
				an = static_cast<double>(n);
				ana = an - a;
				a0 = (a1 + a0 * ana) * fac;
				b0 = (b1 + b0 * ana) * fac;
				anf = an * fac;
				a1 = x * a0 + anf * a1;
				b1 = x * b0 + anf * b1;
				if (a1) {
					fac = 1.0 / a1;
					g = b1 * fac;
					if (fabs((g - gold) / g) < EPS) {
						*gammcf = exp(-x + a * log(x) - (*gln)) * g;
						return;
					}
					gold = g;
				}
			}
			std::cerr << "a too large, ITMAX too small in routine GCF" << "\n";
		}
		else {
			std::string reason;
			reason = "Error: void gcf(double *gammcf, double a, double x, double *gln)\n";
			if (x < 0.0) reason += "x input with value = " + toStringDouble(x, 2) + "\n";
			if (a <= 0.0) reason += "a input with value = " + toStringDouble(a, 2) + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

// Fitting Functions
void covsrt(std::vector< std::vector< double > >& covar, int& ma, std::vector<int>& ia, int& mfit)
{
	// Expand in storage the covariance matrix covar, so as to take into account parameters that are
	// being held fixed. (For the latter, return zero covariances.)

	try {

		if (array_2D_square(covar)) {
			int i, j, k;

			for (i = mfit; i < ma; i++)
				for (j = 0; j < i + 1; j++) covar[i][j] = covar[j][i] = 0.0;

			k = mfit - 1;
			for (j = ma - 1; j >= 0; j--) {
				if (ia[j]) {
					for (i = 0; i < ma; i++) SWAP(covar[i][k], covar[i][j]);
					for (i = 0; i < ma; i++) SWAP(covar[k][i], covar[j][i]);
					k--;
				}
			}

		}
		else {
			std::string reason = "Error: void fit::covsrt()\nInput array is not square\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what(); 
	}
}

void non_lin_fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& covar, std::vector<std::vector<double>>& alpha, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), int& itmax, double& toler, bool loud)
{
	// Levenberg-Marquardt method, attempting to reduce the value chi^{2} of a fit between a set of data points x[0..na-1], y[na-1]
	// with individual standard deviations sig[0..na-1] and a nonlinear function dependent on ma coefficients a[0..ma-1]
	// R. Sheehan 13 - 7 - 2018
	// 
	// Fixed the covariance matrix error and the error which meant that single-parameter fits could not be performed. 
	// Feeling pleased
	// R. Sheehan 22 - 11 - 2019

	// sig cannot contain any zeroes, all elements must be strictly non-zero!
	// R. Sheehan 19 - 10 - 2021

	// What can you do if you do not know the measurement error \sigma? 
	// Method will not work without some estimate of sigma
	// R. Sheehan 26 - 10 - 2021

	try {
		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;
		bool c7 = (int)(ia.size()) == ma ? true : false;
		std::pair<int, int> sze = array_2D_size(alpha);
		bool c9 = sze.first == ma ? true : false;
		bool c10 = sze.second == ma ? true : false;
		std::pair<int, int> szea = array_2D_size(covar);
		bool c9a = szea.first == ma ? true : false;
		bool c10a = szea.second == ma ? true : false;
		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c9 && c10 && c9a && c10a ? true : false;

		if (c11) {

			int itnum = 0;
			double ochisq, alamda;

			ochisq = *chisq;

			while (itnum < itmax) {

				if (itnum == 0) alamda = -1.0;

				// call mrqmin routine until convergence is achieved
				mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, chisq, funcs, &alamda);

				// check for convergence
				if (fabs(*chisq - ochisq) < toler) {

					if (loud) std::cout << "Iterations converged to within tolerance after " << itnum << " steps\n";

					break;
				}

				ochisq = *chisq;

				itnum++;
			}

			// Final call to mrqmin with alamda = 0, so that covar[1..ma][1..ma] returns the 
			// covariance matrix, and alpha the curvature matrix
			alamda = 0.0;
			mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, chisq, funcs, &alamda);

			if (loud) {
				std::cout << "The computed fit parameters are:\n";
				int count = 0;
				for (int i = 0; i < ma; i++) {
					if (ia[i]) {
						std::cout << "a[" << i << "] = " << a[i] << " +/- " << sqrt(covar[count][count]) << "\n";
						count++;
					}
					else {
						std::cout << "a[" << i << "] = " << a[i] << "\n";
					}
				}
				double nu = ndata - count; // nu must be computed on basis of number of fitted parameters
				double q = gammq(0.5 * nu, 0.5 * (*chisq)); // goodness of fit
				std::cout << "\nThe chi-sq value for the fit is " << *chisq << "\n";
				std::cout << "nu for the fit is " << nu << "\n";
				std::cout << "chi-sq / nu = " << *chisq / nu << "\n";
				std::cout << "goodness of fit = " << q << "\n\n";
			}
		}
		else {
			std::string reason = "Error: fit::mrqmin()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + toStringInt(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + toStringInt(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + toStringInt(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + toStringInt(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + toStringInt(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + toStringInt(a.size()) + "\n";
			if (!c7) reason += "ia does not have correct size ia.size() = " + toStringInt(ia.size()) + "\n";
			if (!c9 || !c10) reason += "alpha does not have correct size alpha.size() = ( " + toStringInt(sze.first) + " , " + toStringInt(sze.second) + " )\n";
			if (!c9a || !c10a) reason += "covar does not have correct size covar.size() = ( " + toStringInt(szea.first) + " , " + toStringInt(szea.second) + " )\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void mrqmin(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& covar, std::vector<std::vector<double>>& alpha, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), double* alamda)
{
	// Levenberg-Marquardt method, attempting to reduce the value chi^{2} of a fit between a set of data points x[0..na-1], y[na-1]
	// with individual standard deviations sig[0..na-1] and a nonlinear function dependent on ma coefficients a[0..ma-1]

	// The input array ia[0..ma-1] indicates by nonzero entries those components of a that should be fitted for, and by zero
	// entries those components that should be held fixed at their input values. 

	// The program returns current best - fit values for the parameters in a[0..ma-1] and chi^{2}.  
	// The arrays covar[0..ma-1][0..ma-1] and alpha[0..ma-1][0..ma-1] are used as working space during most iterations. 

	// Routine funcs evaluates the fitting function and its derivatives dyda[0..ma-1] with respect to the fitting parameters a at x. 

	// On the first call provide an initial guess for the parameters a, and set alamda<0 for initialization (which then sets 
	// alamda = .001). If a step succeeds chisq becomes smaller and alamda decreases by a factor of 10. If a step fails alamda 
	// grows by a factor of 10

	// You must call this routine repeatedly until convergence is achieved. Then, make one final 
	// call with alamda = 0, so that covar[1..ma][1..ma] returns the covariance matrix, and alpha the curvature matrix. 
	// (Parameters held fixed will return zero covariances.)

	// R. Sheehan 13 - 7 - 2018

	try {

		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;
		bool c7 = (int)(ia.size()) == ma ? true : false;
		std::pair<int, int> sze = array_2D_size(alpha);
		bool c9 = sze.first == ma ? true : false;
		bool c10 = sze.second == ma ? true : false;
		std::pair<int, int> szea = array_2D_size(covar);
		bool c9a = szea.first == ma ? true : false;
		bool c10a = szea.second == ma ? true : false;
		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c9 && c10 && c9a && c10a ? true : false;

		if (c11) {

			int j, k, l, ncols = 1;
			static int mfit;
			static double ochisq;

			std::vector<double> da;
			std::vector<double> beta;
			std::vector<double> atry;

			std::vector<std::vector<double>> oneda;
			std::vector<std::vector<double>> temp;

			atry = std::vector<double>(ma, 0.0);
			beta = std::vector<double>(ma, 0.0);
			da = std::vector<double>(ma, 0.0);

			for (mfit = 0, j = 0; j < ma; j++) if (ia[j]) mfit++;

			oneda = array_2D(mfit, ncols);
			temp = array_2D(mfit, mfit);

			*alamda = 0.001;

			mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, funcs);

			ochisq = (*chisq);

			for (j = 0; j < ma; j++) atry[j] = a[j];

			// Alter linearised fit matrix along its diagonal elements
			for (j = 0; j < mfit; j++) {
				for (k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
				covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
				for (k = 0; k < mfit; k++) temp[j][k] = covar[j][k];
				oneda[j][0] = beta[j];
			}

			// solve the system of equations covar . x = oneda, store solution in oneda
			// do calculation for covar using array temp
			gaussj(temp, mfit, oneda, ncols);

			// update the vector da
			for (j = 0; j < mfit; j++) {
				for (k = 0; k < mfit; k++) covar[j][k] = temp[j][k]; // retrieve computed values for covar from temp
				da[j] = oneda[j][0];
			}

			// if solution is converged evaluate covariance matrix
			if (*alamda == 0.0) {
				covsrt(covar, ma, ia, mfit);
				covsrt(alpha, ma, ia, mfit);

				oneda.clear(); temp.clear(); da.clear(); beta.clear(); atry.clear();

				return;
			}

			// Examine whether or not the trial is successful
			for (j = 0, l = 0; l < ma; l++)	if (ia[l]) atry[l] = a[l] + da[j++];

			// mrqcof overwrites covar, so work with temp at this point instead of covar
			// that way the correct values for covar are stored in covar
			mrqcof(x, y, sig, ndata, atry, ia, ma, temp, da, chisq, funcs);

			// If success, accept the new solution
			if (*chisq < ochisq) {

				*alamda *= 0.1;

				ochisq = (*chisq);

				for (j = 0; j < mfit; j++) {
					for (k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
					beta[j] = da[j];
				}

				for (l = 0; l < ma; l++) a[l] = atry[l];
			}
			else { // solution is not accurate repeat process with updated alamda value
				*alamda *= 10.0;
				*chisq = ochisq;
			}

			oneda.clear(); temp.clear(); da.clear(); beta.clear(); atry.clear();
		}
		else {
			std::string reason = "Error: fit::mrqmin()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + toStringInt(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + toStringInt(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + toStringInt(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + toStringInt(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + toStringInt(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + toStringInt(a.size()) + "\n";
			if (!c7) reason += "ia does not have correct size ia.size() = " + toStringInt(ia.size()) + "\n";
			if (!c9 || !c10) reason += "alpha does not have correct size alpha.size() = ( " + toStringInt(sze.first) + " , " + toStringInt(sze.second) + " )\n";
			if (!c9a || !c10a) reason += "covar does not have correct size covar.size() = ( " + toStringInt(szea.first) + " , " + toStringInt(szea.second) + " )\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void mrqcof(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, std::vector<int>& ia, int& ma, std::vector<std::vector<double>>& alpha, std::vector<double>& beta, double* chisq, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&))
{
	// Used by mrqmin to evaluate the linearised alpha and beta arrays and to compute chi^{2} of a fit between a set of data points 
	// x[0..na-1], y[na-1] with individual standard deviations sig[0..na-1] and a nonlinear function dependent on ma coefficients 
	// a[0..ma-1]

	// The input array ia[0..ma-1] indicates by nonzero entries those components of a that should be fitted for, and by zero
	// entries those components that should be held fixed at their input values. 

	// This functions computes the system of linear equations that must be solved for at each step
	// alpha . delta = beta

	// Routine funcs evaluates the fitting function and its derivatives dyda[0..ma-1] with respect to the fitting parameters a at x. 

	// R. Sheehan 13 - 7 - 2018

	try {

		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;
		bool c7 = (int)(ia.size()) == ma ? true : false;
		bool c8 = (int)(beta.size()) == ma ? true : false;
		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 ? true : false;

		if (c11) {

			int i, j, k, l, m, mfit = 0;
			double ymod, wt, sig2i, dy;

			std::vector<double> dyda(ma, 0.0);

			for (j = 0; j < ma; j++) if (ia[j]) mfit++;

			// initialise symmetric alpha, beta	
			for (j = 0; j < mfit; j++) {
				for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
				beta[j] = 0.0;
			}

			*chisq = 0.0;

			// sum loop over all data
			for (i = 0; i < ndata; i++) {

				(*funcs)(x[i], a, &ymod, dyda, ma);

				sig2i = 1.0 / (DSQR(sig[i]));

				dy = y[i] - ymod;

				for (j = 0, l = 0; l < ma; l++) {
					if (ia[l]) {
						wt = dyda[l] * sig2i;
						for (k = 0, m = 0; m <= l; m++)
							if (ia[m]) alpha[j][k++] += wt * dyda[m];
						beta[j++] += dy * wt;
					}
				}

				*chisq += DSQR(dy) * sig2i; // find chi^{2}
				//*chisq += DSQR(dy); // find chi^{2}
			}

			// fill the symmetric side
			for (j = 1; j < mfit; j++)
				for (k = 0; k < j; k++) alpha[k][j] = alpha[j][k];

			/*array_2d_print(alpha);
			std::cout << "\n";
			print_vec(beta);
			std::cout << "\n";*/

			dyda.clear();
		}
		else {
			std::string reason = "Error: fit::mrqcof()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + toStringInt(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + toStringInt(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + toStringInt(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + toStringInt(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + toStringInt(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + toStringInt(a.size()) + "\n";
			if (!c7) reason += "ia does not have correct size ia.size() = " + toStringInt(ia.size()) + "\n";
			if (!c8) reason += "beta does not have correct size beta.size() = " + toStringInt(beta.size()) + "\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void goodness_of_fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, int& ma, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), double* chisq, double* nu, double* rsqr, double* gof, bool loud)
{
	// Given a set of data points x[0..na - 1], y[na - 1] with individual standard deviations sig[0..na-1] and a nonlinear function dependent on ma coefficients a[0..ma-1]
	// Compute the value of chi^{2} for the fit, degrees of freedom nu = ndata - no. fitted parameters, the goodness of fit probability, R^{2} coefficient

	// The following contains some useful info	
	// https://en.wikipedia.org/wiki/All_models_are_wrong
	// https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
	// https://en.wikipedia.org/wiki/Coefficient_of_determination
	// https://en.wikipedia.org/wiki/Goodness_of_fit
	// https://en.wikipedia.org/wiki/Kitchen_sink_regression
	// https://en.wikipedia.org/wiki/Data_dredging

	// R. Sheehan 1 - 11 - 2021

	// To characterise a fit as good
	// nu = no. data points - no. fit parameters = degrees of freedom
	// 
	// In general you want chi^{2} / nu ~ 1
	// chi^{2} / nu is the reduced chi^{2} statistic chi^{2} per degree of freedom
	// As a rule of thumb, when measurement variance is known a priori chi^{2} / nu >> 1 implies a poor fit, chi^{2} / nu > 1 indicates fit has not fully
	// captured the data, or that the measurement variance has been underestimated. chi^{2} / nu ~ 1 indicates that the match between measurement and
	// estimates is in line with measurement variance, aka good fit. chi^{2} / nu < 1 indicates the model is over-fitting the data, either the model is 
	// improperly fitting the noise or measurement variance has been overestimated. When the variance of the measurement error is only partially known, 
	// the reduced chi-squared may serve as a correction estimated a posteriori. 
	// 
	// In general you want R^{2} ~ 1
	// The coefficient of determination or R^{2}, is the proportion of the variation in the dependent variable that is predictable from the independent variable(s). 
	// It provides a measure of how well observed outcomes are replicated by the model, based on the proportion of total variation of outcomes explained by the model. 
	// R^{2} ~ 1 indicates a good fit, R^{2} ~ 0 idicates a bad fit, R^{2} outside the range [0, 1] also indicates a bad fit
	// 
	// q = gammq(nu/2, chi^{2}/2) > 0.1 // goodness of fit probability
	// "q measures the probability that a value of chi^{2} as poor as the one computed for the fit should occur by chance. 
	// q > 0.1 implies the gof is believable, q > 0.001  implies the fit may be acceptable if errors are nonnormal or have been underestimated.
	// q < 0.001 implies the fitting procedure may be called into question. " - NRinC, sect 15.2. 

	try {
		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;

		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 ? true : false;

		if (c11) {

			// declare parameters
			double ei = 0.0, ybar = 0.0, yval = 0.0, ssres = 0.0, sstot = 0.0;
			std::vector<double> dyda(ma, 0.0);

			// Compute the mean observed value
			for (int i = 0; i < ndata; i++) ybar += y[i];
			ybar /= ndata;

			// compute the values of the sums needed to compute the statistics
			for (int i = 0; i < ndata; i++) {
				(*funcs)(x[i], a, &yval, dyda, ma);
				ei = y[i] - yval; // difference between observed value and model value
				*chisq += DSQR(ei / sig[i]); // compute the value of the chi^{2} coefficient
				ssres += DSQR(ei); // compute the residual sum of squares
				sstot += DSQR(y[i] - ybar); // compute the total sum of squares
			}
			*rsqr = 1.0 - (ssres / sstot);
			*gof = gammq(0.5 * (*nu), 0.5 * (*chisq)); // goodness of fit probability

			if (loud) {
				std::cout << "\nGoodness-of-fit statistics\n";
				std::cout << "The chi-sq value for the fit is " << *chisq << "\n";
				std::cout << "nu for the fit is " << *nu << "\n";
				std::cout << "chi-sq / nu = " << *chisq / *nu << "\n";
				std::cout << "goodness of fit = " << *gof << "\n";
				std::cout << "coefficient of determination = " << *rsqr << "\n\n";
			}
		}
		else {
			std::string reason = "Error: fit::residuals()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + toStringInt(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + toStringInt(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + toStringInt(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + toStringInt(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + toStringInt(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + toStringInt(a.size()) + "\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void residuals(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, int& ma, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), std::vector<std::vector<double>>& data)
{
	// Compute the residuals of a fitted function
	// i.e. compute the difference between the original data and the fitted function computed using non_lin_fit
	// store the results in data for later graphical analysis
	// output array has the form {x, y, sig, f(x), res}
	// R. Sheehan 21 - 10 - 2021

	try {
		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;

		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 ? true : false;

		if (c11) {
			int ncols = 5;

			data = array_2D(ncols, ndata);

			// store the data from the calculation in an array
			std::copy(x.begin(), x.end(), data[0].begin());
			std::copy(y.begin(), y.end(), data[1].begin());
			std::copy(sig.begin(), sig.end(), data[2].begin());

			double fit_val = 0.0;
			std::vector<double> dyda(ma, 0.0);
			for (int i = 0; i < ndata; i++) {
				funcs(x[i], a, &fit_val, dyda, ma);
				data[3][i] = fit_val;
				data[4][i] = fit_val - y[i];
			}
		}
		else {
			std::string reason = "Error: fit::residuals()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + toStringInt(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + toStringInt(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + toStringInt(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + toStringInt(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + toStringInt(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + toStringInt(a.size()) + "\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

// Actual functions that will be called by external user of DLL

// Fit a Lorentzian to a data set
void Lorentz_Fit(int n_data, double freq_data[], double spctrm_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[])
{
	// Use non-lin-fit to fit the Lorentz model to a set of RSA spectrum data
	// n_data is no. of data points in measurement
	// freq_data[] is an array holding the frequency data, assumed to be in units of MHz
	// spctrm_data[] is an array holding the measured ESA spectral data, assumed to be in units of dBm
	// fit_data[] is an array that will store the computed values of the fitted model on output, assumed to be in units of dBm
	// n_pars is the no. of fitting parameters, 3 in the case of the Lorentz fit
	// a[] = {A, f_{0}, HWHM} is an array holding initial estimates of the fit parameters, this will be overwritten 
	// with the fitted values on output
	// n_stats is the number of goodness of fit statistics that are computed
	// gof_stats[] = {chi^{2} value for fit, chi^{2} / nu, R^{2} coefficient, gof probability } is an array that will store the computed goodness of fit stats on output
	// R. Sheehan 4 - 11 - 2021

	// Declare various parameters
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Create std::vector for computing the fits
	std::vector<double> x(n_data, 0.0); 
	std::vector<double> y(n_data, 0.0); 
	std::vector<double> sig(n_data, 0.0); 

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = array_2D(n_pars, n_pars);
	std::vector<std::vector<double>> alpha = array_2D(n_pars, n_pars);

	std::vector<std::vector<double>> data;

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(n_pars, 0.0);
	std::vector<int> ia(n_pars, 1); // tell the algorithm that you want to locate all parameters 

	// store the data in the vector containers
	double spread = 0.05; 
	double scale_fac = 1.0e+6; 
	for (int i = 0; i < n_data; i++) {
		x[i] = freq_data[i]; 
		y[i] = scale_fac * convert_dBm_to_mW( spctrm_data[i] ); // convert from dBm scale to mW scale
		sig[i] = spread * y[i];
	}

	// store the initial guesses
	ia[1] = 0; // no sense in trying to fit to the centre frequency since you know its value already
	for (int i = 0; i < n_pars; i++) {
		a_guess[i] = a_pars[i]; 
	}

	// perform the fit process
	bool loud = false; 

	non_lin_fit(x, y, sig, n_data, a_guess, ia, n_pars, covar, alpha, &chisq, Lorentzian, ITMAX, TOL, loud);

	// compute the residuals
	residuals(x, y, sig, n_data, a_guess, n_pars, Lorentzian, data);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(n_data - n_pars), gof = 0.0;
	goodness_of_fit(x, y, data[4], n_data, a_guess, n_pars, Lorentzian, &chisqr, &dof, &rsqr, &gof, loud);

	// store the computed model values
	for (int i = 0; i < n_data; i++) {
		fit_data[i] = convert_mW_to_dBm( data[3][i] / scale_fac ); // convert from mW to dBm scale
	}

	// store the computed fit parameters
	for (int i = 0; i < n_pars; i++) {
		a_pars[i] = a_guess[i];
	}

	// store the computed fit statistics
	gof_stats[0] = chisqr; gof_stats[1] = chisqr / dof; gof_stats[2] = rsqr; gof_stats[3] = gof;

	// release memory
	x.clear(); y.clear(); sig.clear(); data.clear(); 
	a_guess.clear(); ia.clear(); 
	covar.clear(); alpha.clear(); 
}

void Gauss_Fit(int n_data, double freq_data[], double spctrm_data[], double fit_data[], int n_pars, double a_pars[], int n_stats, double gof_stats[])
{
	// Use non-lin-fit to fit the Lorentz model to a set of RSA spectrum data
	// n_data is no. of data points in measurement
	// freq_data[] is an array holding the frequency data, assumed to be in units of MHz
	// spctrm_data[] is an array holding the measured ESA spectral data, assumed to be in units of dBm
	// fit_data[] is an array that will store the computed values of the fitted model on output, assumed to be in units of dBm
	// n_pars is the no. of fitting parameters, 3 in the case of the Gaussian fit
	// a[] = {A, f_{0}, HWHM} is an array holding initial estimates of the fit parameters, this will be overwritten 
	// with the fitted values on output
	// n_stats is the number of goodness of fit statistics that are computed
	// gof_stats[] = {chi^{2} value for fit, chi^{2} / nu, R^{2} coefficient, gof probability } is an array that will store the computed goodness of fit stats on output
	// R. Sheehan 19 - 11 - 2021

	// Declare various parameters
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Create std::vector for computing the fits
	std::vector<double> x(n_data, 0.0);
	std::vector<double> y(n_data, 0.0);
	std::vector<double> sig(n_data, 0.0);

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = array_2D(n_pars, n_pars);
	std::vector<std::vector<double>> alpha = array_2D(n_pars, n_pars);

	std::vector<std::vector<double>> data;

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(n_pars, 0.0);
	std::vector<int> ia(n_pars, 1); // tell the algorithm that you want to locate all parameters 

	// store the data in the vector containers
	double spread = 0.05;
	double scale_fac = 1.0e+6;
	for (int i = 0; i < n_data; i++) {
		x[i] = freq_data[i];
		y[i] = scale_fac * convert_dBm_to_mW(spctrm_data[i]); // convert from dBm scale to mW scale
		sig[i] = spread * y[i];
	}

	// store the initial guesses
	ia[1] = 0; // no sense in trying to fit to the centre frequency since you know its value already
	for (int i = 0; i < n_pars; i++) {
		a_guess[i] = a_pars[i];
	}

	// perform the fit process
	bool loud = false;

	non_lin_fit(x, y, sig, n_data, a_guess, ia, n_pars, covar, alpha, &chisq, Gaussian, ITMAX, TOL, loud);

	// compute the residuals
	residuals(x, y, sig, n_data, a_guess, n_pars, Gaussian, data);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(n_data - n_pars), gof = 0.0;
	goodness_of_fit(x, y, data[4], n_data, a_guess, n_pars, Gaussian, &chisqr, &dof, &rsqr, &gof, loud);

	// store the computed model values
	for (int i = 0; i < n_data; i++) {
		fit_data[i] = convert_mW_to_dBm(data[3][i] / scale_fac); // convert from mW to dBm scale
	}

	// store the computed fit parameters
	for (int i = 0; i < n_pars; i++) {
		a_pars[i] = a_guess[i];
	}

	// store the computed fit statistics
	gof_stats[0] = chisqr; gof_stats[1] = chisqr / dof; gof_stats[2] = rsqr; gof_stats[3] = gof;

	// release memory
	x.clear(); y.clear(); sig.clear(); data.clear();
	a_guess.clear(); ia.clear();
	covar.clear(); alpha.clear();
}

void Testing()
{
	std::cout << "Hello World\n"; 
	std::cout << "I am being run from the DLL\n\n";
}