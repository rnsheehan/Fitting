#ifndef ATTACH_H
#include "Attach.h"
#endif

// Implementation of a namespace that contains functions used to compute probabilities in Statistics
// R. Sheehan 5 - 9 - 2017

double probability::gammln(double xx)
{
	// Return the value of ln[gamma(xx)] for xx>0

	try{
		if(xx > 0.0){
			double x,tmp,ser;
			static double cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
			int j;

			x=xx-1.0;
			tmp=x+5.5;
			tmp-=(x+0.5)*log(tmp);
			ser=1.0;
			for (j=0;j<=5;j++){
				x+=1.0;
				ser+=(cof[j]/x);
			}
			return -tmp+log(2.50662827465*ser);		
		}
		else{
			std::string reason; 
			reason = "Error: double probability::gammln(double xx)\n"; 
			reason += "xx input with value = " + template_funcs::toString(xx, 2) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double probability::gammp(double a,double x)
{
	// Computes the incomplete gamma function P(a,x) from the functions gser and gcf

	try{
		if(x > 0.0 && a > 0.0){
			double gamser,gammcf,gln;

			if(x<(a+1.0)){
				gser(&gamser,a,x,&gln);
				return gamser;
			} 
			else{
				gcf(&gammcf,a,x,&gln);
				return 1.0-gammcf;
			}
		}
		else{
			std::string reason; 
			reason = "Error: double probability::gammp(double a,double x)\n"; 
			if(x < 0.0) reason += "x input with value = " + template_funcs::toString(x, 2) + "\n"; 
			if(a <= 0.0) reason += "a input with value = " + template_funcs::toString(a, 2) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}	
}

double probability::gammq(double a,double x)
{
	// Computes the incomplete gamma function Q(a,x)=1-P(a,x) from the functions gser and gcf

	try{
		if(x > 0.0 && a > 0.0){
			double gamser,gammcf,gln;
			if(x<(a+1.0)){
				gser(&gamser,a,x,&gln);
				return 1.0-gamser;
			} 
			else{
				gcf(&gammcf,a,x,&gln);
				return gammcf;
			}
		}
		else{
			std::string reason; 
			reason = "Error: double probability::gammq(double a,double x)\n"; 
			if(x < 0.0) reason += "x input with value = " + template_funcs::toString(x, 2) + "\n"; 
			if(a <= 0.0) reason += "a input with value = " + template_funcs::toString(a, 2) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void probability::gser(double *gamser,double a,double x,double *gln)
{
	// This function returns the incomplete gamma function P(a,x), calculated by its series representation
	// Also returns ln[gamma(a)] as gln

	try{
		if(x > 0.0 && a > 0.0){

			int n;
			double sum,del,ap;

			static const int ITMAX=(150);
			static const double EPS=(3.0e-7);

			*gln=gammln(a);
			if(x<=0.0){
				if(x<0.0)
					std::cerr<<"x less than 0 in routine GSER"<<"\n";
				*gamser=0.0;
				return;
			}else{
				ap=a;
				del=sum=1.0/a;
				for(n=1;n<=ITMAX;n++){
					ap+=1.0;
					del*=x/ap;
					sum+=del;
					if(fabs(del)<fabs(sum)*EPS){
						*gamser=sum*exp(-x+a*log(x)-(*gln));
						return;
					}
				}
				std::cerr<<"a too large, ITMAX too small in routine GSER"<<"\n";
				return;
			}			
		}
		else{
			std::string reason; 
			reason = "Error: void probability::gser(double *gamser,double a,double x,double *gln)\n"; 
			if(x < 0.0) reason += "x input with value = " + template_funcs::toString(x, 2) + "\n"; 
			if(a <= 0.0) reason += "a input with value = " + template_funcs::toString(a, 2) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void probability::gcf(double *gammcf, double a, double x, double *gln)
{
	// This function returns the incomplete gamma function Q(a,x), calculated by its continued fraction representation
	// Also returns ln[gamma(a)] as gln

	try{
		if(x > 0.0 && a > 0.0){
			int n;
			double gold=0.0,g,fac=1.0,b1=1.0;
			double b0=0.0,anf,ana,an,a1,a0=1.0;

			static const int ITMAX=(100);
			static const double EPS=(3.0e-7);

			*gln=gammln(a);
			a1=x;
			for(n=1;n<=ITMAX;n++){
				an=static_cast<double>(n);
				ana=an-a;
				a0=(a1+a0*ana)*fac;
				b0=(b1+b0*ana)*fac;
				anf=an*fac;
				a1=x*a0+anf*a1;
				b1=x*b0+anf*b1;
				if(a1){
					fac=1.0/a1;
					g=b1*fac;
					if(fabs((g-gold)/g)<EPS){
						*gammcf=exp(-x+a*log(x)-(*gln))*g;
						return;
					}
					gold=g;
				}
			}
			std::cerr<<"a too large, ITMAX too small in routine GCF"<<"\n";
		}
		else{
			std::string reason; 
			reason = "Error: void probability::gcf(double *gammcf, double a, double x, double *gln)\n"; 
			if(x < 0.0) reason += "x input with value = " + template_funcs::toString(x, 2) + "\n"; 
			if(a <= 0.0) reason += "a input with value = " + template_funcs::toString(a, 2) + "\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double probability::erff(double x)
{
	// Return the error function erf(x)

	return ( x < 0.0 ? -gammp(0.5, template_funcs::DSQR(x) ) : gammp(0.5, template_funcs::DSQR(x) ) ) ; 
}

double probability::erffc(double x)
{
	// Return the complementary error function erfc(x)

	return ( x < 0.0 ? 1.0 + gammp(0.5, template_funcs::DSQR(x) ) : gammq(0.5, template_funcs::DSQR(x) ) ) ; 
}

double probability::erfcc(double x)
{
	// Returns the complementary error function erfc(x) with fractional error everywhere less than 1.2e-7

	double t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}

double probability::factorial(int n)
{
	// return the factorial of n by recursion
	// It is possible to represent 170! as a floating point number
	// Attempting to compute 171! will cause numerical overflow
	// Factorials up to 22! are exact, 23! and higher are numerically approximate
	// See NRinC, sect 6.1. 
	// R. Sheehan 15 - 4 - 2014

	if(n >-1 && n < 170){

		if(n == 0 || n == 1){
			return 1.0; 
		}
		else if(n == 2){
			return 2.0; 
		}
		else if(n == 3){
			return 6.0; 
		}
		else if(n == 4){
			return 24.0; 
		}
		else if(n == 5){
			return 120.0; 
		}
		else if(n == 6){
			return 720.0; 
		}
		else{
			// recursively compute n!
			return ( static_cast<double>(n)*factorial(n-1) ); 
		}
	}
	else{
		std::cerr<<"Cannot compute "<<n<<"! by this method\nc.f. gammln routine for factorials of large numbers"; 
		return 0.0; 
	}
}

double probability::betacf(double a, double b, double x)
{
	// Evaluation of continued fraction expansion for incomplete beta function
	// NRinC equation 6.4.5, 6.4.6

	try{
		int m,m2;
		int MAXIT = 100; 
		double EPS1 = 3.0E-7; 
		static const double FPMIN=(1.0e-30);
		double aa,c,d,del,h,qab,qam,qap;

		qab=a+b;
		qap=a+1.0;
		qam=a-1.0;
		c=1.0;
		d=1.0-qab*x/qap;
		if (fabs(d) < FPMIN) d=FPMIN;
		d=1.0/d;
		h=d;
		for (m=1;m<=MAXIT;m++) {
			m2=2*m;
			aa=m*(b-m)*x/((qam+m2)*(a+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			h *= d*c;
			aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			del=d*c;
			h *= del;
			if (fabs(del-1.0) < EPS1) break;
		}
		if (m > MAXIT){
			std::string reason; 
			reason = "Error: probability::betacf()\na or b too big, or MAXIT too small in betacf\n"; 
			throw std::runtime_error(reason); 
		}
		return h;
	}
	catch(std::runtime_error &e){
		std::cerr<<e.what(); 
		return 0.0; 
	}
}

double probability::betai(double a, double b, double x)
{
	// Return the value of the incomplete Beta function I_{x}(a, b)
	// NRinC equation 6.4.1
	// this is used to compute probabilities related to Student's distribution

	try{
		if(x < 0.0 || x > 1.0){
			std::string reason; 
			reason = "Error: probability::betai()\n"; 
			reason += "Bad x value\n"; 
			throw std::invalid_argument(reason); 
		}
		else{
			double bt; 
			if (x == 0.0 || x == 1.0){
				bt=0.0;
			}
			else{
				bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
			}

			if (x < (a+1.0)/(a+b+2.0)){
				return bt*betacf(a,b,x)/a;
			}
			else{
				return 1.0-bt*betacf(b,a,1.0-x)/b;		
			}
		}
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what(); 
		return 0.0; 
	}
}

double probability::probks(double alam)
{
	// Kolmorgorov-Smirnov Probability Function

	try{

		if(alam < 0.0){
			std::string reason; 
			reason = "Error: double probability::probks(double alam)\n"; 
			reason += "alam = " + template_funcs::toString(alam) + " < 0\n"; 
			throw std::invalid_argument(reason); 
		}
		else{
			int j;
			static const double EPS1 = 0.001; 
			static const double EPS2 = 1.0e-8; 
			double a2,fac=2.0,sum=0.0,term,termbf=0.0;

			a2 = -2.0*alam*alam;
			for (j=1;j<=100;j++) {
				term=fac*exp(a2*j*j);
				sum += term;
				if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
				fac = -fac;
				termbf=fabs(term);
			}
			return 1.0;
		}	
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what(); 
		return 0.0; 
	}
}