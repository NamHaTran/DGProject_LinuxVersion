#include <cmath>
#include <vector>
#include <iostream>

double newtonRaphsonForPolynomial(std::vector<double> &power, std::vector<double> &coefs, double initialValue)
{
    int polySize(static_cast<int>(power.size()));
    std::vector<double> power_deriv(polySize, 0.0),
            coefs_deriv(polySize, 0.0);
    double error(1), convergence_cri(1e-6), output(0.0), fValue(0.0), derivfValue(0.0);

    //find 1st derivative of input polynomial
    for (int polyOrder = 0; polyOrder < polySize; polyOrder++) {
        power_deriv[polyOrder]=power[polyOrder] - 1.0;
        coefs_deriv[polyOrder]=power[polyOrder]*coefs[polyOrder];
    }

    //solve equation
    while (error > convergence_cri) {
        fValue = 0.0;
        derivfValue = 0.0;
        for (int polyOrder = 0; polyOrder < polySize; polyOrder++) {
            fValue+=pow(initialValue,power[polyOrder])*coefs[polyOrder];
            derivfValue+=pow(initialValue,power_deriv[polyOrder])*coefs_deriv[polyOrder];
        }
        output = initialValue - fValue/derivfValue;
        error = fabs(output - initialValue)/initialValue;
        initialValue = output;
    }
    return output;
}

int main()
{
    std::vector<double> power{4.0, 3.0, 2.0, 1.0, 0.0},
    coefs{20, -60, 130.0, 122.0, -197.0};
    double iniVal(-1.2), root(0.0);
    root = newtonRaphsonForPolynomial(power,coefs,iniVal);
    std::cout<<"Root is "<<root<<std::endl;
    return 0;
}
