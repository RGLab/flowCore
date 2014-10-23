#include "logicle.h"
#include "zeroin.h"
#include <memory.h>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <limits>

const double Logicle::DEFAULT_DECADES = 4.5;

const double Logicle::LN_10 = log(10.);
const double Logicle::EPSILON = std::numeric_limits<double>::epsilon();
const double Logicle::NaN = std::numeric_limits<double>::quiet_NaN();

const int Logicle::TAYLOR_LENGTH = 16;

Logicle::Exception::Exception()
{
	buffer = 0;
}

Logicle::Exception::Exception(const Logicle::Exception & e)
{
	buffer = strdup(e.buffer);
}

Logicle::Exception::Exception (const char * const message)
{
	buffer = strdup(message);
}

Logicle::Exception::~Exception ()
{
	delete buffer;
}

const char * Logicle::Exception::message () const
{
	return buffer;
}

Logicle::IllegalArgument::IllegalArgument (double value)
{
	buffer = new char[128];
	sprintf(buffer, "Illegal argument value %.17g", value);
}

Logicle::IllegalArgument::IllegalArgument (int value)
{
	buffer = new char[128];
	sprintf(buffer, "Illegal argument value %d", value);
}

Logicle::IllegalParameter::IllegalParameter (const char * const message) : Exception(message)
{	}

Logicle::DidNotConverge::DidNotConverge(const char * const message) : Exception(message)
{	}

void Logicle::initialize (double T, double W, double M, double A, int bins)
{
	// allocate the parameter structure
	p = new logicle_params;
	p->taylor = 0;
  	if (T <= 0)
		throw "IllegalParameter: T is not positive";
		//throw IllegalParameter("T is not positive");
	if (W < 0)
        throw "IllegalParameter: W is not positive";
        //throw IllegalParameter("W is not positive");
	if (M <= 0)
		throw "IllegalParameter: M is not positive";
        //throw IllegalParameter("M is not positive");
	if (2 * W > M)
		throw "IllegalParameter: W is too large";
        //throw IllegalParameter("W is too large");
	if (-A > W || A + W > M - W)
        throw "IllegalParameter: A is too large";
        //throw IllegalParameter("A is too large");

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
	// formulas from biexponential paper
	p->w = W / (M + A);
	p->x2 = A / (M + A);
	p->x1 = p->x2 + p->w;
	p->x0 = p->x2 + 2 * p->w;
	p->b = (M + A) * LN_10;
	p->d = solve(p->b, p->w);
	double c_a = exp(p->x0 * (p->b + p->d));
	double mf_a = exp(p->b * p->x1) - c_a / exp(p->d * p->x1);
	p->a = T / ((exp(p->b) - mf_a) - c_a / exp(p->d));
	p->c = c_a * p->a;
	p->f = -mf_a * p->a;

	// use Taylor series near x1, i.e., data zero to
	// avoid round off problems of formal definition
	p->xTaylor = p->x1 + p->w / 4;
	// compute coefficients of the Taylor series
	double posCoef = p->a * exp(p->b * p->x1);
	double negCoef = -p->c / exp(p->d * p->x1);
	// 16 is enough for full precision of typical scales
	p->taylor = new double[TAYLOR_LENGTH];
	for (int i = 0; i < TAYLOR_LENGTH; ++i)
	{
		posCoef *= p->b / (i + 1);
		negCoef *= -p->d / (i + 1);
		(p->taylor)[i] = posCoef + negCoef;
	}
	p->taylor[1] = 0; // exact result of Logicle condition
}

Logicle::Logicle (double T, double W, double M, double A)
{
	initialize(T, W, M, A, 0);
}

Logicle::Logicle (double T, double W, double M, double A, int bins)
{
	initialize(T, W, M, A, bins);
}

Logicle::Logicle (const Logicle & logicle)
{
	p = new logicle_params;
	memcpy(p, logicle.p, sizeof(logicle_params) );
	p->taylor = new double[TAYLOR_LENGTH];
	memcpy(p->taylor, logicle.p->taylor, TAYLOR_LENGTH * sizeof(double));
}

Logicle::~Logicle ()
{
	delete[] p->taylor;
	delete p;
}

// f(w,b) = 2 * (ln(d) - ln(b)) + w * (b + d)
double logicle_fn(double x,void*info) {
	struct sfun_info *p = (struct sfun_info *)info;
	double B = 2 * (log(x) - log(p->b)) + p->w * (p->b + x);
    return (B);
}

/*
 * root finder routines are copied from stats/src/zeroin.c
 */
double Logicle::R_zeroin(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit)				/* Max # of iterations */
{
    double fa = (*f)(ax, info);
    double fb = (*f)(bx, info);
    return R_zeroin2(ax, bx, fa, fb, f, info, Tol, Maxit);
}

/* R_zeroin2() is faster for "expensive" f(), in those typical cases where
 *             f(ax) and f(bx) are available anyway : */

double Logicle::R_zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit)				/* Max # of iterations */
{
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return a;
    }
    if(fb ==  0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return b;
    }

    while(maxit--)		/* Main iteration loop	*/
    {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	{				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	}
	tol_act = 2*EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	{
	    *Maxit -= maxit;
	    *Tol = fabs(c-b);
	    return b;			/* Acceptable approx. is found	*/
	}

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	    register double t1,cb,t2;
	    cb = c-b;
	    if( a==c ) {		/* If we have only two distinct	*/
					/* points linear interpolation	*/
		t1 = fb/fa;		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	    }
	    else {			/* Quadric inverse interpolation*/

		q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	    if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	    else			/* and assign possible minus to	*/
		p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	    if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	    else
		new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;	fb = (*f)(b, info);	/* Do step to a new approxim. */
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	    /* Adjust c for it to have a sign opposite to that of b */
	    c = a;  fc = fa;
	}

    }
    /* failed! */
    *Tol = fabs(c-b);
    *Maxit = -1;
    return b;
}

/*
 * use R built-in root finder API :R_zeroin
 */
double Logicle::solve (double b, double w)
{
	// w == 0 means its really arcsinh
	if (w == 0)
		return b;

	// precision is the same as that of b
	double tolerance = 2 * b * EPSILON;
	struct sfun_info params;
	params.b=b;
	params.w=w;

	// bracket the root
	double d_lo = 0;
	double d_hi = b;


	int MaxIt = 20;
	double d ;
	d= R_zeroin(d_lo,d_hi,logicle_fn,(void*)&params,&tolerance,&MaxIt);
	return d;
}

double Logicle::slope (double scale) const
{
	// reflect negative scale regions
	if (scale < p->x1)
		scale = 2 * p->x1 - scale;

	// compute the slope of the biexponential
	return p->a * p->b * exp(p->b * scale) + p->c * p->d / exp(p->d * scale);
}

double Logicle::seriesBiexponential (double scale) const
{
	// Taylor series is around x1
	double x = scale - p->x1;
	// note that taylor[1] should be identically zero according
	// to the Logicle condition so skip it here
	double sum = p->taylor[TAYLOR_LENGTH - 1] * x;
	for (int i = TAYLOR_LENGTH - 2; i >= 2; --i)
		sum = (sum + p->taylor[i]) * x;
	return (sum * x + p->taylor[0]) * x;
}

double Logicle::scale (double value) const
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
	if (value < p->f)
		// use linear approximation in the quasi linear region
		x = p->x1 + value / p->taylor[0];
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
		double ce2mdx = p->c / exp(p->d * x);
		double y;
		if (x < p->xTaylor)
			// near zero use the Taylor series
			y = seriesBiexponential(x) - value;
		else
			// this formulation has better roundoff behavior
			y = (ae2bx + p->f) - (ce2mdx + value);
		double abe2bx = p->b * ae2bx;
		double cde2mdx = p->d * ce2mdx;
		double dy = abe2bx + cde2mdx;
		double ddy = p->b * abe2bx - p->d * cde2mdx;

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
	//throw DidNotConverge("scale() didn't converge");
}

double Logicle::inverse (double scale) const
{
	// reflect negative scale regions
	bool negative = scale < p->x1;
	if (negative)
		scale = 2 * p->x1 - scale;

	// compute the biexponential
	double inverse;
	if (scale < p->xTaylor)
		// near x1, i.e., data zero use the series expansion
		inverse = seriesBiexponential(scale);
	else
		// this formulation has better roundoff behavior
		inverse = (p->a * exp(p->b * scale) + p->f) - p->c / exp(p->d * scale);

	// handle scale for negative values
	if (negative)
		return -inverse;
	else
		return inverse;
}

double Logicle::dynamicRange () const
{
	return slope(1) / slope(p->x1);
}

int PullInMyLibrary () { return 0; }
