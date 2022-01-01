#ifndef _M_ESTIMATORS_H_
#define _M_ESTIMATORS_H_

#include <cmath>

inline double get_l1_rho(double r){
	return fabs(r);
}

inline double get_l1_upsilon(double r){
	return (r > 0) ? 1 : ((r < 0) ? -1 : 0);
}

inline double get_l1_w(double r){
	return 1/fabs(r);
}

inline double get_l2_rho(double r){
	return (r*r)/2;
}

inline double get_l2_upsilon(double r){
	return r;
}

inline double get_l2_w(double r){
	return 1;
}

inline double get_l1l2_rho(double r){
	return 2*(sqrt(1 + (r*r)/2.0) - 1);
}

inline double get_l1l2_upsilon(double r){
	return r/(sqrt(1+(r*r)/2.0));
}

inline double get_l1l2_w(double r){
	return 1/(sqrt(1+(r*r)/2.0));
}

inline double get_lp_rho(double r, double nu){
	return (pow(fabs(r),nu))/nu;
}

inline double get_lp_upsilon(double r, double nu){
	return ((r > 0) ? 1 : ((r < 0) ? -1 : 0)) * pow(fabs(r),nu-1);
}

inline double get_lp_w(double r, double nu){
	return pow(fabs(r),nu-2);
}

inline double get_fair_rho(double r, double c){
	return c*c*(fabs(r)/c - log(1 + fabs(r)/c));
}

inline double get_fair_upsilon(double r, double c){
	return r/(1+fabs(r)/c);
}

inline double get_fair_w(double r, double c){
	return 1/(1+fabs(r)/c);
}

inline double get_huber_rho(double r, double k){
	if(fabs(r) < k){
		return (r*r)/2;
	}else{
		return k*(fabs(r)-k/2);
	}
}

inline double get_huber_upsilon(double r, double k){
	if(fabs(r) < k){
		return r;
	}else{
		return k*((r > 0) ? 1 : ((r < 0) ? -1 : 0));
	}
}

inline double get_huber_w(double r, double k){
	if(fabs(r) < k){
		return 1;
	}else{
		return k/fabs(r);
	}
}

inline double get_cauchy_rho(double r, double c){
	return (c*c)/2 * log(1 + (r/c)*(r/c));
}

inline double get_cauchy_upsilon(double r, double c){
	return r/(1 + (r/c)*(r/c));
}

inline double get_cauchy_w(double r, double c){
	return 1/(1 + (r/c)*(r/c));
}

inline double get_geman_mcclure_rho(double r){
	return ((r*r)/2)/(1 +r*r);
}

inline double get_geman_mcclure_upsilon(double r){
	return r/( (1+r*r)*(1+r*r) );
}

inline double get_geman_mcclure_w(double r){
	return 1/( (1+r*r)*(1+r*r) );
}

inline double get_welsch_rho(double r, double c){
	return (c*c)/2 * (1 - exp(-(r/c)*(r/c)));
}

inline double get_welsch_upsilon(double r, double c){
	return r*exp(-(r/c)*(r/c));
}

inline double get_welsch_w(double r, double c){
	return exp(-(r/c)*(r/c));
}

inline double get_tukey_rho(double r, double c){
	if(fabs(r) < c){
		double temp = 1-(r/c)*(r/c);
		return (c*c)/6 * (1 - temp*temp*temp);
	}else{
		return (c*c)/6;
	}
}

inline double get_tukey_upsilon(double r, double c){
	if(fabs(r) < c){
		double temp = 1-(r/c)*(r/c);
		return r*temp*temp;
	}else{
		return 0;
	}
}

inline double get_tukey_w(double r, double c){
	if(fabs(r) < c){
		double temp = 1-(r/c)*(r/c);
		return temp*temp;
	}else{
		return 0;
	}
}

#endif
