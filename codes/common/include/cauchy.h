#ifndef _CAUCHY_H_
#define _CAUCHY_H_

inline double cauchy(double delta, double b){
	return 1.0 / (M_PI * b *( 1.0 + ((delta)/b) * ((delta)/b) ) );
}

#endif
