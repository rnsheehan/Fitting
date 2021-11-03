// Non_Lin_Fit.h - Contains declarations of non-linear fitting functions
// This header file declares some functions perform various non-linear fits
#pragma once

#include <vector>

// preprocessor statements
// The new project template for a DLL project adds PROJECTNAME_EXPORTS to the defined preprocessor macros. 
// In this example, Visual Studio defines MATHLIBRARY_EXPORTS when your MathLibrary DLL project is built.
// When the MATHLIBRARY_EXPORTS macro is defined, the MATHLIBRARY_API macro sets the __declspec(dllexport) modifier on the function declarations. 
// This modifier tells the compiler and linker to export a function or variable from the DLL for use by other applications.
// When MATHLIBRARY_EXPORTS is undefined, for example, when the header file is included by a client application, MATHLIBRARY_API applies the __declspec(dllimport) modifier to the declarations.
// This modifier optimizes the import of the function or variable in an application.

#ifdef NON_LIN_FIT_EXPORTS
#define NON_LIN_FIT_API __declspec(dllexport)
#else
#define NON_LIN_FIT_API __declspec(dllimport)
#endif

// Function Declarations

extern "C" NON_LIN_FIT_API double DSQR(double a);

extern "C" NON_LIN_FIT_API void Lorentzian(double x, std::vector<double>&a, double* y, std::vector<double>&dyda, int& na);

extern "C" NON_LIN_FIT_API void Gaussian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na);

