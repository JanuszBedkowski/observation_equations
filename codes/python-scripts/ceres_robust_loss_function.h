inline double trivial_loss(double s)
{return s;
}
inline double trivial_loss_first_derivative(double s)
{return 1;
}
inline double trivial_loss_second_derivative(double s)
{return 0;
}
inline double huber_loss_less_1(double s)
{return s;
}
inline double huber_loss_less_1_first_derivative(double s)
{return 1;
}
inline double huber_loss_less_1_second_derivative(double s)
{return 0;
}
inline double huber_loss_more_1(double s)
{return 2*sqrt(s) - 1;
}
inline double huber_loss_more_1_first_derivative(double s)
{return pow(s, -1.0/2.0);
}
inline double huber_loss_more_1_second_derivative(double s)
{return -1.0/2.0/pow(s, 3.0/2.0);
}
inline double soft_lone_loss(double s)
{return 2*sqrt(s + 1) - 2;
}
inline double soft_lone_loss_first_derivative(double s)
{return pow(s + 1, -1.0/2.0);
}
inline double soft_lone_loss_second_derivative(double s)
{return -1.0/2.0/pow(s + 1, 3.0/2.0);
}
inline double cauchy_loss(double s)
{return log(s + 1);
}
inline double cauchy_loss_first_derivative(double s)
{return 1.0/(s + 1);
}
inline double cauchy_loss_second_derivative(double s)
{return -1/pow(s + 1, 2);
}
inline double arctan_loss(double s)
{return atan(s);
}
inline double arctan_loss_first_derivative(double s)
{return 1.0/(pow(s, 2) + 1);
}
inline double arctan_loss_second_derivative(double s)
{return -2*s/pow(pow(s, 2) + 1, 2);
}
inline double tolerant_loss(double s, double a, double b)
{return -b*log(1 + exp(-a/b)) + b*log(exp((-a + s)/b) + 1);
}
inline double tolerant_loss_first_derivative(double s, double a, double b)
{return exp((-a + s)/b)/(exp((-a + s)/b) + 1);
}
inline double tolerant_loss_second_derivative(double s, double a, double b)
{return exp((-a + s)/b)/(b*(exp((-a + s)/b) + 1)) - exp(2*(-a + s)/b)/(b*pow(exp((-a + s)/b) + 1, 2));
}
