#include <fstream>

void ciddor(double lambda, double RH, double p, double t, double CO2){

    // Saturatiokn Vapor Pressure //

    cout << setprecision(10) ;

    // For saturation vapor pressure over water
    const double K1  = 1.16705214528E+03;
    const double K2  = -7.24213167032E+05;
    const double K3  = -1.70738469401E+01;
    const double K4  = 1.20208247025E+04;
    const double K5  = -3.23255503223E+06;
    const double K6  = 1.49151086135E+01;
    const double K7  = -4.82326573616E+03;
    const double K8  = 4.05113405421E+05;
    const double K9  = -2.38555575678E-01;
    const double K10 = 6.50175348448E+02;
    const double exp = 2.71828182846;

    double T = t + 273.15 ;
    double Omega = T + K9/(T-K10);
    double A = pow(Omega, 2) + K1 * Omega + K2;
    double B = K3 * pow(Omega, 2) + K4 * Omega + K5;
    double C = K6 * pow(Omega, 2) + K7 * Omega + K8;
    double X = -B + sqrt(pow(B,2) - 4*A*C) ;
    double psv = pow(10,6) * pow((2*C/X),4);

    /// Humidity for Edlen Equation ///
    // doube pv = (RH/100)*psv ;

    ///// For saturation vapor pressure over ice

    // const A1 =  -13.928169;
    // const A2 =  34.7078238;

    // double Theta = T/273.16;
    // double Y = A1 * (1 - Theta^(-1.5)) + A2 * (1 - Theta^(-1.25));
    // double psv = 611.657 * exp^Y;

    // Humidity for Ciddor Equation //

    // First calculate the enhancement factor f using

    const double alpha = 1.00062;
    const double beta  = 3.14 * pow(10,-8);
    const double gamma = 5.60 * pow(10,-7);
    double f = alpha + beta * p + gamma*pow(t, 2);

    // If you are given dew point td (or frost point) in degrees Celsius, and air pressure p in Pascals, calculate mole fraction xv using

    // xv = f * psv * / p;

    // If you are given relative humidity RH in percent (ranging numerically from 0 to 100), air pressure p (Pascals), and air temperature t, calculate the mole fraction xv using

    double xv = (RH/100) * f * psv / p;

    // If you are given partial pressure pv (Pascals), air pressure p (Pascals), and air temperature t (Celsius), calculate mole fraction xv using

    // double xv = f * pv / p;



    ///// Ciddor Calculation of Index of Refraction /////

    // Convert all temperatures to Celsius.
    // Convert all pressures to Pascal.
    // Calculate the mole fraction xv as described previously.

    const double w0 = 295.235;
    const double w1 = 2.6422;
    const double w2 = -0.0323;
    const double w3 = 0.004028;

    const double k0 = 238.0185;
    const double k1 = 5792105;
    const double k2 = 57.362;
    const double k3 = 167917;

    const double a0 = 1.58123 * pow(10,-6);
    const double a1 = -2.9331 * pow(10,-8);
    const double a2 = 1.1043 * pow(10,-10);

    const double b0 = 5.707 * pow(10,-6);
    const double b1 = -2.051 * pow(10,-8);

    const double c0 = 1.9898 * pow(10,-4);
    const double c1 = -2.376 * pow(10,-6);

    const double d =1.83 * pow(10,-11);
    const double e = -0.765 * pow(10,-8);

    const double pR1 = 101325;
    const double TR1 = 288.15;
    const double Za = 0.9995922115;

    const double rhovs = 0.00985938;
    const double R = 8.314472;
    const double Mv = 0.018015;
    const double S = 1 / pow(lambda,2);

    const double ras = pow(10,-8) * (( k1 / (k0 - S)) + (k3 / (k2 - S)));
    const double rvs = 1.022 * pow(10,-8) * (w0 + w1 * S + w2 * pow(S,2) + w3 * pow(S,3));

    const double Ma = 0.0289635 + 1.2011 * pow(10,-8) * ( CO2 - 400 );
    const double raxs = ras * (1 + 5.34 * pow(10,-7) * (CO2 - 450));

    const double Zm = 1 - (p / T) * ( a0 + a1 * t + a2 * pow(t,2) + (b0 + b1 * t) * xv + (c0 + c1 * t) * pow(xv,2) ) + pow((p / T),2) * (d + e * pow(xv,2));

    const double rhoaxs = pR1 * Ma / (Za * R * TR1);
    const double rhov = xv * p * Mv / (Zm * R * T);
    const double rhoa = (1 - xv) * p * Ma / (Zm * R * T);

    const double n = 1 + (rhoa / rhoaxs) * raxs + (rhov / rhovs) * rvs;

    cout << "Refractive Index : " << n << endl;


}
