

#include "hyperlog.h"
#include "zeroin.h"
#include <memory.h>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <limits>

const double Hyperlog::DEFAULT_DECADES = 4.5;

const double Hyperlog::LN_10 = log(10.);
const double Hyperlog::EPSILON = std::numeric_limits<double>::epsilon();
const double Hyperlog::NaN = std::numeric_limits<double>::quiet_NaN();

const int Hyperlog::TAYLOR_LENGTH = 16;

Hyperlog::Exception::Exception()
{
    buffer = 0;
}

Hyperlog::Exception::Exception(const Hyperlog::Exception & e)
{
    buffer = strdup(e.buffer);
}

Hyperlog::Exception::Exception (const char * const message)
{
    buffer = strdup(message);
}

Hyperlog::Exception::~Exception ()
{
    delete buffer;
}

const char * Hyperlog::Exception::message () const
{
    return buffer;
}

Hyperlog::IllegalArgument::IllegalArgument (double value)
{
    buffer = new char[128];
    snprintf(buffer, sizeof(buffer), "Illegal argument value %.17g", value);
}

Hyperlog::IllegalArgument::IllegalArgument (int value)
{
    buffer = new char[128];
    snprintf(buffer, sizeof(buffer), "Illegal argument value %d", value);
}

Hyperlog::IllegalParameter::IllegalParameter (const char * const message) : Exception(message)
{
}

Hyperlog::DidNotConverge::DidNotConverge(const char * const message) : Exception(message)
{
}

void Hyperlog::initialize (double T, double W, double M, double A, int bins)
{
    // allocate the parameter structure
    p = new hyperlog_params;
    p->taylor = 0;
    if (T <= 0)
        throw "IllegalParameter: T is not positive";
    if (W < 0)
        throw "IllegalParameter: W is negative";
    if (W <= 0)
        throw "IllegalParameter: W is not positive";
    if (M <= 0)
        throw "IllegalParameter: M is not positive";
    if (2 * W > M)
        throw "IllegalParameter: W is too large";
    if (-A > W || A + W > M - W)
       throw "IllegalParameter: A is too large";

    // if we're going to bin the data make sure that
    // zero is on a bin boundary by adjusting A
    if (bins > 0)
    {
        double zero = (W + A) / (M + A);
        zero = floor(zero * bins + .5) / bins;
        A = (M * zero - W) / (1 - zero);
    }

    // standard parameters
    p->T = T;
    p->M = M;
    p->W = W;
    p->A = A;

    // actual parameters
		// choose the data zero location and the width of the linearization region
		// to match the corresponding logicle scale
    p->w = W / (M + A);
    p->x2 = A / (M + A);
    p->x1 = p->x2 + p->w;
    p->x0 = p->x2 + 2 * p->w;
    p->b = (M + A) * LN_10;
    double e2bx0 = exp(p->b * p->x0);
    double c_a = e2bx0 / p->w;
    double f_a = exp(p->b * p->x1) + c_a * p->x1;

    p->a = T / ((exp(p->b) + c_a) - f_a);
    p->c = c_a * p->a;
    p->f = f_a * p->a;

    // use Taylor series near x1, i.e., data zero to
    // avoid round off problems of formal definition
    p->xTaylor = p->x1 + p->w / 4;

    // compute coefficients of the Taylor series
    double coef = p->a * exp(p->b * p->x1);
    // 16 is enough for full precision of typical scales
    p->taylor = new double[TAYLOR_LENGTH];
    for (int i = 0; i < TAYLOR_LENGTH; ++i)
    {
        coef *= p->b / (i + 1);
        (p->taylor)[i] = coef;
    }
    p->taylor[0] += p->c; // hyperlog condition
    p->inverse_x0 = inverse(p->x0);
}

Hyperlog::Hyperlog (double T, double W, double M, double A)
{
    initialize(T, W, M, A, 0);
}

Hyperlog::Hyperlog (double T, double W, double M, double A, int bins)
{
    initialize(T, W, M, A, bins);
}

Hyperlog::Hyperlog (const Hyperlog & hyperlog)
{
    p = new hyperlog_params;
    memcpy(p, hyperlog.p, sizeof(hyperlog_params) );
    p->taylor = new double[TAYLOR_LENGTH];
    memcpy(p->taylor, hyperlog.p->taylor, TAYLOR_LENGTH * sizeof(double));
}

Hyperlog::~Hyperlog ()
{
    delete[] p->taylor;
    delete p;
}

double Hyperlog::slope (double scale) const
{
    // reflect negative scale regions
    if (scale < p->x1)
        scale = 2 * p->x1 - scale;

    // compute the slope of the biexponential
    return p->a * p->b * exp(p->b * scale) + p->c;
}

double Hyperlog::taylorSeries (double scale) const
{
    // Taylor series is around x1
    double x = scale - p->x1;
    double sum = p->taylor[TAYLOR_LENGTH - 1] * x;
		for (int i = TAYLOR_LENGTH - 2; i >= 0; --i)
			  sum = (sum + p->taylor[i]) * x;
		return sum;
}

double Hyperlog::scale (double value) const
{
    // handle true zero separately
    if (value == 0)
        return p->x1;

    // reflect negative values
    bool negative = value < 0;
    if (negative)
        value = -value;

    // initial guess at solution
    double x;
    if (value < p->inverse_x0)
        // use linear approximation in the quasi linear region
        x = p->x1 + value * p->w / p->inverse_x0;
    else
        // otherwise use ordinary logarithm
        x = log(value / p->a) / p->b;

    // try for double precision unless in extended range
    double tolerance = 3 * EPSILON;
    if (x > 1)
        tolerance = 3 * x * EPSILON;

    for (int i = 0; i < 10; ++i)
    {
        // compute the function and its first two derivatives
        double ae2bx = p->a * exp(p->b * x);
        //double ce2mdx = p->c / exp(p->d * x);
        double y;
        if (x < p->xTaylor)
            // near zero use the Taylor series
            y = taylorSeries(x) - value;
        else
            // this formulation has better roundoff behavior
            y = (ae2bx + p->c * x) - (p->f + value);

        double abe2bx = p->b * ae2bx;
        double dy = abe2bx + p->c;
        double ddy = p->b * abe2bx;

        // this is Halley's method with cubic convergence
        double delta = y / (dy * (1 - y * ddy / (2 * dy * dy)));
        x -= delta;

        // if we've reached the desired precision we're done
        if (std::abs(delta) < tolerance) {
            // handle negative arguments
            if (negative)
                return 2 * p->x1 - x;
            else
                return x;
        }
    }

    throw "DidNotConverge: scale() didn't converge";
}

double Hyperlog::inverse (double scale) const
{
    // reflect negative scale regions
    bool negative = scale < p->x1;
    if (negative)
        scale = 2 * p->x1 - scale;

    double inverse;
    if (scale < p->xTaylor)
        // near x1, i.e., data zero use the series expansion
        inverse = taylorSeries(scale);
    else
        // this formulation has better roundoff behavior
        inverse = (p->a * exp(p->b * scale) + p->c * scale) - p->f;

    // handle scale for negative values
    if (negative)
        return -inverse;
    else
        return inverse;
}

