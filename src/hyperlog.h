// Copyright (c) 2021 Ozette Technologies
// Copyright (c) 2021 Fred Hutch Cancer Center
// Hyperlog transformation added by Josef Spidlen.
// This hyperlog implementation is based on Java reference
// implementation that is part of the full Gating-ML 2.0
// specification. The Java reference implementation has
// been provided by Wayne Moore, see hyperlog.notice.html
// for details. Josef Spidlen ported it to C/CPP and
// integrated it with R/flowCore.

#ifdef __cplusplus

extern "C" {

#endif

    struct hyperlog_params
    {
        double T, W, M, A;

    double a, b, c, f;
        double w, x0, x1, x2;
    double inverse_x0;

        double xTaylor;
        double *taylor;

        double *lookup;
        int bins;
    };

    struct sfun_info{
        double b,w;
    };

    double hyperlog_fn(double x,void*info);
    const char * hyperlog_error ();
    const struct hyperlog_params * hyperlog_initialize (double T, double W, double M, double A, int bins);
    void hyperlog_destroy (const struct hyperlog_params * params);
    double hyperlog_scale (const struct hyperlog_params * hyperlog, double value);
    int hyperlog_int_scale (const struct hyperlog_params * hyperlog, double value);
    double hyperlog_inverse (const struct hyperlog_params * hyperlog, double scale);

#ifdef __cplusplus

}

class Hyperlog
{
public:
    static const double DEFAULT_DECADES;

    class Exception
    {
    public:
        Exception (const Exception & e);
        virtual ~Exception ();
        const char * message () const;

    protected:
        char * buffer;
        Exception ();
        Exception (const char * const message);

    private:
        Exception & operator= (const Exception & e);
        friend class Hyperlog;
    };

    class IllegalArgument : public Exception
    {
    private:
        IllegalArgument (double value);
        IllegalArgument (int value);
        friend class Hyperlog;
    };

    class IllegalParameter : public Exception
    {
    private:
        IllegalParameter (const char * const message);
        friend class Hyperlog;
    };

    class DidNotConverge : public Exception
    {
    private:
        DidNotConverge (const char * const message);
        friend class Hyperlog;
    };

    Hyperlog (double T, double W, double M = DEFAULT_DECADES, double A = 0);
    Hyperlog (const Hyperlog & hyperlog);

    virtual ~Hyperlog ();

    inline double T() const { return p->T; };
    inline double W() const { return p->W; };
    inline double M() const { return p->M; };
    inline double A() const { return p->A; };

    inline double a() const { return p->a; };
    inline double b() const { return p->b; };
    inline double c() const { return p->c; };
    inline double f() const { return p->f; };

    inline double w() const { return p->w; };
    inline double x0() const { return p->x0; };
    inline double x1() const { return p->x1; };
    inline double x2() const { return p->x2; };
    inline double inverse_x0() const { return p->inverse_x0; };

    virtual double scale (double value) const;
    virtual double inverse (double scale) const;

protected:
    static const double LN_10;
    static const double EPSILON;
    static const double NaN;
    static const int TAYLOR_LENGTH;

    hyperlog_params * p;

    Hyperlog (double T, double W, double M, double A, int bins);

    double slope (double scale) const;
    double taylorSeries (double scale) const;

private:
    Hyperlog & operator= (const Hyperlog & hyperlog);
    void initialize (double T, double W, double M, double A, int bins);

};

int PullInMyLibrary ();

#endif

