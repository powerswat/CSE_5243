
#include <all.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>
#include <assert.h>
#include <float.h>

// $Id: bayeslsh.c,v 1.2 2011-12-14 22:32:35 venu Exp $

static inline float convertJaccardForHashing(float sim)
{
	return (1.0+sim)/2;
}

static inline float convertHashingForJaccard(float sim)
{
	return (2*sim-1);
}

static inline float convertCosineForHashing(float sim)
{
	return 1.0-acos(sim)/PI;
}

static inline float convertHashingForCosine(float x)
{
	return cos(PI*(1.0-x));
}

static inline uint32_t popcount ( uint32_t v )
{
	// http://graphics.stanford.edu/~seander/bithacks.html
	v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
	uint32_t c = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count

	return c;
}


void my_gsl_error_handler(const char * reason, const char * file,
int line, int gsl_errno)
{
	fprintf(stderr, "Error: %s, in file:%s, line:%d\n", 
			reason, file, line);
}

static inline double gsl_lnBetaFunction(double a, double b)
{
	gsl_sf_result r;
	int eno = gsl_sf_lnbeta_e(a, b, &r);	
	if ( !eno || eno == GSL_EUNDRFLW )
	{
		return r.val;
	}
	else
	{
		fprintf(stderr, "Error getting beta for a %f, b %f\n", 
						a, b);
		abort();
	}
}

static inline double gsl_betaCdf(double a, double b, double x)
{
	gsl_sf_result r;
	int eno = gsl_sf_beta_inc_e(a, b, x, &r);	
	if ( !eno || eno == GSL_EUNDRFLW )
	{
		return r.val;
	}
	else
	{
			fprintf(stderr, "Error getting betaCdf for a %f, b %f, x %f\n", 
						a, b, x);
			abort();
	}
}

inline float posteriorCdf_unifHalfTo1(void *vprior, float s,
float n, float x)
{
	if ( x >= 1.0 )
		return 1.0;
	if ( x <= 0.5 )
		return 0;
	float b1 = 1.0;
	float bHalf = gsl_betaCdf(s+1, n-s+1, 0.5);
	float bx = gsl_betaCdf(s+1, n-s+1, x);
	float den = b1 - bHalf;
	if ( den < 1.0e-15 )
		return exp( log(bx-bHalf) - log(den) );
	else
		return (bx-bHalf)/den;
}


inline float posteriorCdf_unif(void *vprior, float s, float n,
float x)
{
	Unif prior = *((Unif*) vprior);

	if ( x <= prior.lb )
		return 0;
	if ( x >= prior.ub )
		return 1;

	return gsl_betaCdf(s+1.0, n-s+1.0, x);
}

inline float posteriorMode_unif_dummyCache(void *vprior, int s, int
n, float **pmCache, const int MIN_MATCHES_STEP, const int
*minMatches)
{
	Unif prior = *((Unif*)vprior);
//	float m = 1.0*s/n;
	if ( prior.lb <= 0.0001 )
		s -= 6;
	else 
		s -= 12;

	float m = 1.0*(s)/n;

	if ( m <= prior.lb )
		return prior.lb;
	
	return m;
}

inline float posteriorMode_unif(void *vprior, int s, int n )
{
	Unif prior = *((Unif*)vprior);
	float m = 1.0*s/n;

	if ( m <= prior.lb )
		return prior.lb;
	
	return m;
}

inline float posteriorCdf_beta(void *vprior, float s,
float n, float x)
{
	Beta prior = *((Beta*) vprior);

	if ( x <= 0 )
		return 0;
	if ( x >= 1 )
		return 1;

	return gsl_betaCdf(s+prior.alpha, n-s+prior.beta, x);
}

inline float posteriorMode_beta(void *vprior, int s, int
n)
{
	Beta prior = *((Beta*) vprior);
	s -= 5;
	float mode = ((prior.alpha+s-1)/(prior.alpha+n+prior.beta-2));
	return (mode > 0) ? mode : 0;
}

inline float posteriorMode_beta_dummyCache(void *vprior, int s, int
n, float **pmc, const int MMS, const int *minMatches)
{
	return posteriorMode_beta(vprior, s, n);
}


float posteriorCdf_pwUnif0to1(void *vprior, float nSuccesses, float nTrials,
float x)
{
	PwUnif prior = *((PwUnif*) vprior);

	if ( x <= prior.lb )
		return 0;
	if ( x >= prior.ub )
		return 1;
	float s = nSuccesses;
	float n = nTrials;

	float bx = gsl_betaCdf(s+1.0, n-s+1.0, x);
	float bChangePt = gsl_betaCdf(s+1.0, n-s+1.0,
							prior.changePt);
	float den = prior.p1 * (bChangePt) + 
					prior.p2*(1.0-bChangePt);
	float p;
	if ( x < prior.changePt )
	{
		p = prior.p1 * bx;
	}
	else
	{
		p = prior.p1 * bChangePt + prior.p2 * (bx-bChangePt) ;	
	}
	if ( den < 1.0e-15 )
	{
		p = (float) exp(log(p)-log(den));
	//	p = 1;
	}
	else
	{
		p /= den;
	}

	return p;
}

double sign_dx2_pdf_linear(const double a, const double b, const
double s, const double n, const double x)
{
	// http://www.wolframalpha.com/input/?i=d^2%2Fdx^2+%28x^s%29*%28%281-x%29^%28n-s%29%29*%28a*x%2Bb%29

	const double ax_mult = -2*x*n*s + n*x*(n*x+x-2) + s*(s+1);
	const double b_mult = s*(-2*n*x + 2*x -1) + (n-1)*n*x*x +
	s*s;
	const double ret = a*x*ax_mult + b*b_mult;
	return ret;
}

double posteriorMode_linear(const double a, const double b, const
double s, const double n, const double lb, const double ub,
double * const mode_pdf)
{
	// http://www.wolframalpha.com/input/?i=d%2Fdx+%28x^s%29*%28%281-x%29^%28n-s%29%29*%28a*x%2Bb%29
	const double a_eqn = -a*(n+1);
	const double b_eqn = a*(s+1)-n*b;
	const double c_eqn = b*s;
	const double sq = b_eqn*b_eqn - 4*a_eqn*c_eqn;

	if ( sq > 0 )
	{
		const double sqt = sqrt(sq);
		const double x1 = (-1*b_eqn + sqt)/(2*a_eqn);
		const double x2 = (-1*b_eqn - sqt)/(2*a_eqn);
		// check that double derivative is negative
		if ( x1 > lb && x1 < ub && (sign_dx2_pdf_linear(a, b, s,
		n, x1) < 0) )
		{
			*mode_pdf = s*log(x1)+(n-s)*log(1-x1)+log(a*x1+b);
			return x1;
		}
		if ( x2 > lb && x2 < ub && (sign_dx2_pdf_linear(a, b, s,
		n, x2) < 0) )
		{
			*mode_pdf = s*log(x2)+(n-s)*log(1-x2)+log(a*x2+b);
			return x2;
		}
	}

	double ub_pdf = -DBL_MAX, lb_pdf = -DBL_MAX;
	if ( ub > 0 && ub < 1 )
		ub_pdf = s*log(ub)+(n-s)*log(1-ub)+log(a*ub+b);
	if ( lb > 0 && lb < 1 )
		lb_pdf = s*log(lb)+(n-s)*log(1-lb)+log(a*lb+b);

	if ( ub_pdf > lb_pdf )
	{
		*mode_pdf = ub_pdf;
		return ub;
	}
	
	*mode_pdf = lb_pdf;
	return lb;
}

float posteriorMode_PwLinear(void *vprior, float s, float n)
{
	PwLinear prior = *((PwLinear*) vprior);

	double bestMode = -1;
	double bestModeP = -DBL_MAX;

	double lb, ub, mode, modeP;
	for ( int i=0; i<prior.k; i++ )
	{
		if ( i == prior.k-1 )
			ub = prior.ub;
		else
			ub = prior.changePts[i];

		if ( i > 0 )
			lb = prior.changePts[i-1];
		else
			lb = prior.lb;

		mode = posteriorMode_linear(prior.as[i], prior.bs[i], 
					s, n, lb, ub, &modeP);
		if ( modeP > bestModeP )
		{
			bestMode = mode;
			bestModeP = modeP;
		}
	}


	return bestMode;
}

float posteriorMode_PwLinear_cache(void *vprior, int s, int n, float **posteriorModeCache, const int MIN_MATCHES_STEP, const int *minMatches)
{
	int step = (n/MIN_MATCHES_STEP) - 1;
	float c = posteriorModeCache[step][s-minMatches[step]];

	if ( c > 0 )
		return c;

	c = posteriorMode_PwLinear(vprior, s, n);
	posteriorModeCache[step][s-minMatches[step]] = c;

	return c;
}

double posteriorCdf_linear(const double a, const double b, 
	const double s, const double n, const double x, 
	const double beta_a_lower, const double beta_b_lower, const
	double lnBetaSPlus2, const double lnBetaSPlus1, 
	double * const beta_a_upper, double * const beta_b_upper)
{
	*beta_a_upper = exp( log(gsl_betaCdf(s+2.0, n-s+1.0, x)) + 
						lnBetaSPlus2);
	*beta_b_upper = exp( log(gsl_betaCdf(s+1.0,n-s+1.0,x)) +
						lnBetaSPlus1);

	double ret = a*(*beta_a_upper-beta_a_lower) + 
				b*(*beta_b_upper-beta_b_lower);

	assert(ret > -0.0000001);
	return ret;
}

float posteriorCdf_PwLinear(void *vprior, float s, float n, float
x)
{
	PwLinear prior = *((PwLinear*) vprior);

	assert ( x >= prior.lb && x <= prior.ub);

	double num=0, den=0;
//	const double betaSPlus2 = gsl_betaFunction(s+2.0, n-s+1.0);
//	const double lnBetaSPlus2 = gsl_lnBetaFunction(s+2.0, n-s+1.0);
	const double lnBetaSPlus2 = gsl_lnBetaFunction(s+2.0,
	n-s+1.0) - gsl_lnBetaFunction(s+1.0, n-s+1.0);
//	const double betaSPlus1 = gsl_betaFunction(s+1.0, n-s+1.0);
//	const double lnBetaSPlus1 = gsl_lnBetaFunction(s+1.0, n-s+1.0);
	const double lnBetaSPlus1 = 0;

	// dividing both numerator and denominator by
	// lnBeta(s+1,n-s+1)

	double beta_a_lower = exp(log(gsl_betaCdf(s+2.0, n-s+1.0,
		prior.lb)) + lnBetaSPlus2);
//	double beta_a_lower = gsl_betaCdf(s+2.0, n-s+1.0,
//	prior.lb) * betaSPlus2;
	double beta_b_lower = exp(log(gsl_betaCdf(s+1.0, n-s+1.0,
		prior.lb)) + lnBetaSPlus1);
//	double beta_b_lower = gsl_betaCdf(s+1.0, n-s+1.0, prior.lb) *
//	betaSPlus1;
	double beta_a_upper, beta_b_upper, 
		x_beta_a_lower = beta_a_lower, x_beta_b_lower = beta_b_lower;

	int xChangePt = -1;
	for ( int i=0; i<prior.k-1; i++ )
	{
		double p = posteriorCdf_linear(prior.as[i], prior.bs[i],
		s, n, prior.changePts[i], beta_a_lower, beta_b_lower,
		lnBetaSPlus2, lnBetaSPlus1, &beta_a_upper, &beta_b_upper);
		if ( prior.changePts[i] < x )
		{
			num += p;
			x_beta_a_lower = beta_a_upper;
			x_beta_b_lower = beta_b_upper;
		}
		else if ( xChangePt == -1 )
		{
			xChangePt = i;
		}

		den += p;
		beta_a_lower = beta_a_upper;
		beta_b_lower = beta_b_upper;
	}
	if ( xChangePt == -1 )
		xChangePt = prior.k-1;
	
	num += posteriorCdf_linear(prior.as[xChangePt],
	prior.bs[xChangePt], s, n, x, x_beta_a_lower, x_beta_b_lower,
	lnBetaSPlus2, lnBetaSPlus1, &beta_a_upper, &beta_b_upper);
	den += posteriorCdf_linear(prior.as[prior.k-1],
	prior.bs[prior.k-1], s, n, prior.ub, beta_a_lower,
	beta_b_lower, lnBetaSPlus2, lnBetaSPlus1, &beta_a_upper, &beta_b_upper);

	float r = (float) exp( log(num) - log(den));
	assert ( r < 1.00000001 );
	return r;
}

float posteriorCdf_KPwUnif(void *vprior, float s, float n, float
x)
{
	KPwUnif prior = *((KPwUnif*) vprior);

	if ( x <= prior.lb )
		return 0;
	if ( x >= prior.ub )
		return 1;

	double num=0; // we will accumulate probability in this
	double bPreviousChangePt = gsl_betaCdf(s+1.0, n-s+1.0,
	prior.lb);
	double bxChangePt = bPreviousChangePt;
	double bCurrentChangePt;
	double den=0;
	int j = -1;
	for ( int i=0; i < prior.k-1; i++ )
	{
		bCurrentChangePt = gsl_betaCdf(s+1.0, n-s+1.0,
							prior.changePts[i]);
		double p = prior.ps[i] * (bCurrentChangePt -
				bPreviousChangePt);
		if ( prior.changePts[i] < x )
		{
			num += p;
			bxChangePt = bCurrentChangePt; // the last of these
			// will be stored in bxChangePt
	//		bxChangePt = bPreviousChangePt;
		}
		else if ( j == -1 )
			j = i;

		den += p;
		bPreviousChangePt = bCurrentChangePt;
	}
	if ( j == -1 )
		j = prior.k-1;

	num += prior.ps[j] * (gsl_betaCdf(s+1.0,n-s+1.0,x) -
	bxChangePt);
	den += prior.ps[prior.k-1] * (1.0 - bPreviousChangePt);

	float r = (float) exp( log(num) - log(den));
	return r;
}

float posteriorMode_KPwUnif(void *vprior, int s, int n)
{
	KPwUnif prior = *((KPwUnif*) vprior);
	// first locate interval containing x
	float x = s*1.0/n;
	int xInterval = 0;
	for (; prior.changePts[xInterval] < x; xInterval++ );
	float bestYet = x;
	float bestYetP = prior.logps[xInterval] + s*log(x) +
						(n-s)*log(1-x);

	// strictly speaking, should examine all changepts, but we
	// will only examine changept to the left and to the right of
	// x.  Assuming our intervals are reasonably spaced out, this
	// should be fine.

	// check left intervals
/*	int numLeftIntervalsToCheck=10;
	for( int leftInterval=xInterval-1; leftInterval>=0 &&
	xInterval-leftInterval
*/
	if ( xInterval > 0 && prior.ps[xInterval-1] > prior.ps[xInterval])
	{
		float p = prior.logps[xInterval-1] +
		s*prior.logChangePts[xInterval-1] +
		(n-s)*prior.log1minusChangePts[xInterval-1];
		if ( p > bestYetP )
		{
			bestYetP = p;
			bestYet = prior.changePts[xInterval-1];
		}
	}

	// check right interval
	if ( xInterval < prior.k && prior.ps[xInterval+1] >
	prior.ps[xInterval] )
	{
		float p = prior.logps[xInterval+1] +
		s*prior.logChangePts[xInterval+1] +
		(n-s)*prior.log1minusChangePts[xInterval+1];
		if ( p > bestYetP )
		{
			bestYetP = p;
			bestYet = prior.changePts[xInterval];
		}
	}

	return bestYet;
}

float posteriorMode_KPwUnif_cache(void *vprior, int s, int n, float **posteriorModeCache, const int MIN_MATCHES_STEP, const int *minMatches)
{
	int step = (n/MIN_MATCHES_STEP) - 1;
	float c = posteriorModeCache[step][s-minMatches[step]];

	if ( c > 0 )
		return c;

	c = posteriorMode_KPwUnif(vprior, s, n);
	posteriorModeCache[step][s-minMatches[step]] = c;

	return c;
}

int fitToLinear(const float * const samples,  
const double probPerSample, const double lb, const double ub, const
double valuesToTry, float * const ret_a, float * const ret_b)
{
	int numSamples = 0;
	double totalProb = 0;
	float maxSampleValue = ub;

	assert ( samples[0] < maxSampleValue );

	for ( numSamples=0; samples[numSamples] < maxSampleValue;
	numSamples++ )
		totalProb += probPerSample;

	assert( totalProb > 0 );

	const double a_lb = -2*totalProb/((ub-lb)*(ub-lb)) + 0.000001;
	const double a_ub = -1 * a_lb;
	const double b_conv_add = totalProb/(ub-lb);
	const double b_conv_mult = -1*(ub+lb)/2;

	double bestA = -DBL_MAX, bestB = -DBL_MAX;
	double bestA_ll = -DBL_MAX;

	const double inc = (a_ub - a_lb)/valuesToTry;

	for ( double a=a_lb; a<a_ub; a += inc)
	{
		double b = a*b_conv_mult + b_conv_add;

		double ll=0;
		int i;
		for ( i=0; samples[i] < maxSampleValue; i++ )
		{
			ll += log(a*samples[i] + b);
		}

		if ( ll > bestA_ll )
		{
			bestA = a;
			bestA_ll = ll;
			bestB = b;
		}
	}

	*ret_a = bestA;
	*ret_b = bestB;

	return numSamples;
}

PwLinear::PwLinear(float *samples, const int
numSamples, const int kprime, const float lba, const float uba,
const float threshold)
{
	k = kprime;
	lb = lba;
	ub = uba;

	int numFloats = 2*(k) + k-1;
	base =(float*) calloc(numFloats, sizeof(float));
	as = base;
	bs = as + k;
	changePts = bs + k;

	const double a_fitting_inc = 0.25; // at what spacing we'll
	// possible values of a, in the max. likelihood estimation

	qsort(samples, numSamples, sizeof(float), float_compare);

	for ( int i=0; i<k-1; i++ )
	{
		changePts[i] = lb + (i+1)*(ub-lb)/k;
//		changePts[i] = samples[(i+1)*numSamples/k]+0.00001;
	assert(changePts[i]>lb && changePts[i] <ub);
	}
	
	float incPerSample = 1.0/numSamples;

	int numSamplesSoFar = 0;
	for ( int interval=0; interval < k; interval++ )
	{
		float ub = uba + 0.000001;
		if ( interval < k-1)
			ub = changePts[interval]+0.000001;
		float lb = lba;
		if ( interval > 0 )
			lb = changePts[interval-1];
		numSamplesSoFar += 
			fitToLinear(samples+numSamplesSoFar, incPerSample, lb, 
				ub, 1000, &as[interval], &bs[interval]);
//		printf("Interval %d, lb: %f, ub: %f, a: %f, b: %f\n",
//		interval, lb, ub, as[interval], bs[interval]);
	}

}


KPwUnif::KPwUnif(float *sims, const int numSamples,
const int kprime, const float lba, const float uba, const float
firstChangePt)
{
	k = kprime;
	int numFloats = 2*(k) + 3*(k-1);
	base =(float*) calloc(numFloats, sizeof(float));
	ps = base;
	logps = ps + (k);
	changePts = logps + (k);
	logChangePts = changePts + k-1;
	log1minusChangePts = logChangePts +  k-1;

	lb = lba;
	ub = uba;

	qsort(sims, numSamples, sizeof(float), float_compare);

	float incPerSample = 1.0/numSamples;
	for ( int i=0; i<k-1; i++ )
	{
//		changePts[i] = lb + (i+1)*(ub-lb)/(k);
		changePts[i] = sims[(i+1)*numSamples/k]+0.00001;
	}

/*	changePts[0] = firstChangePt;
	for ( int i=1; i<k-1; i++ )
	{
		changePts[i] = firstChangePt +
		i*(ub-firstChangePt)/(k-1);
	}
*/
	for ( int i=0; i<k-1; i++ )
	{
		logChangePts[i] = log(changePts[i]);
		log1minusChangePts[i] = log(1.0-changePts[i]);
	}

	int iSample=0, iChangePt=0;
	for ( ; iChangePt < k-1; iChangePt++ )
	{
		for ( ; sims[iSample] < changePts[iChangePt] && iSample <
				numSamples; iSample++ )
			ps[iChangePt] += incPerSample;
	}
	for (; iSample<numSamples; iSample++)
		ps[iChangePt] += incPerSample;
	
	for ( int i=0; i<k; i++ )
		logps[i] = log(ps[i]);
}

void KPwUnif::print()
{
	for ( int i=0; i<k-1; i++ )
	{
		printf("changePt: %.3f, prob: %.4f\n", changePts[i],
		ps[i]);
	}
	printf("changePt: 1.0, prob: %.4f\n", ps[k-1]);
}
/*
KPwUnif::~KPwUnif()
{
	free(base);
	base = ps = logps = changePts = logChangePts =
	log1minusChangePts = NULL;
}
*/
//static inline float posteriorCdf_pwUnif(void *vprior, float nSuccesses, float nTrials,
float posteriorCdf_pwUnif(void *vprior, float nSuccesses, float nTrials,
float x)
{
	PwUnif prior = *((PwUnif*) vprior);

	if ( x <= prior.lb )
		return 0;
	if ( x >= prior.ub )
		return 1;
	float s = nSuccesses;
	float n = nTrials;

	float bx = gsl_betaCdf(s+1.0, n-s+1.0, x);
	float bChangePt = gsl_betaCdf(s+1.0, n-s+1.0,
							prior.changePt);
	float bHalf = gsl_betaCdf(s+1.0, n-s+1.0, 0.5);
	float den = prior.p1 * (bChangePt-bHalf) + 
					prior.p2*(1.0-bChangePt);
	float p;
	if ( x < prior.changePt )
	{
		p = prior.p1 * (bx - bHalf);
	}
	else
	{
		p = prior.p1 * (bChangePt-bHalf) + prior.p2 *
				(bx - bChangePt);	
	}
	if ( den < 1.0e-15 )
	{
		p = (float) exp(log(p)-log(den));
	//	p = 1;
	}
	else
	{
		p /= den;
	}

	return p;
}

float posteriorMean(PwUnif prior, float nSuccesses, float nTrials)
{
	float s = nSuccesses;
	float n = nTrials;
	float num_bChangePt = gsl_betaCdf(s+2.0, n-s+1.0, prior.changePt);
	float num_bHalf = gsl_betaCdf(s+2.0, n-s+1.0, 0.5);
	
	float den_bChangePt = gsl_betaCdf(s+1.0, n-s+1.0, prior.changePt);
	float den_bHalf = gsl_betaCdf(s+1.0, n-s+1.0, 0.5);

	float numerator = prior.p1*(num_bChangePt - num_bHalf) +
						prior.p2*(1.0 - num_bChangePt);
	float denominator = prior.p1*(den_bChangePt - den_bHalf) +
						prior.p2*(1.0 - den_bChangePt);

	if ( denominator < 1.0e-25 )
	{
		return exp ( log(numerator) - log(denominator) );
	//	return 1;
	}
	float p = numerator / denominator;
	return p;

}

float posteriorMode_pwUnif(void *vprior, int is, int in)
{
	PwUnif prior = *((PwUnif*)vprior);
	float s = (float) is;
	float n = (float) in;
	// s is no. of successes
	// n is no. of trials
	float m = s/n;
	if ( m <= prior.lb )
		return prior.lb;

	float p_m;
	if ( is == in)
		p_m = 0;
	else
		p_m = s*log(m) + (n-s)*log(1.0-m);

	float p_c = s*prior.logChangePt +
				(n-s)*prior.log1minusChangePt;
	if ( m < prior.changePt )
	{
		if ( p_m + prior.logp1 > p_c + prior.logp2 )
			return m;
		else
			return prior.changePt;
	}
	else
	{
		if ( p_m + prior.logp2 > p_c + prior.logp1 )
			return m;
		else
			return prior.changePt;
	}
}

inline float posteriorMode_pwUnif_dummyCache(void *vprior, int matches, int
nTrials, float **posteriorModeCache, const int
MIN_MATCHES_STEP, const int *minMatches)
{
	return posteriorMode_pwUnif(vprior, matches, nTrials);
}

inline float posteriorMode_pwUnifcache(void *vprior, int matches, int
nTrials, float **posteriorModeCache, const int
MIN_MATCHES_STEP, const int *minMatches)
{
	int step = (nTrials/MIN_MATCHES_STEP) - 1;
	float c = posteriorModeCache[step][matches-minMatches[step]];
	if ( c > 0 )
		return c;

	c = posteriorMode_pwUnif(vprior, matches, nTrials);
	posteriorModeCache[step][matches-minMatches[step]] = c;

	return c;

}

int worstCaseHashesChebyshev(float threshold, float delta, float
epsilon, int factor)
{
	if ( threshold < 0.5 )
		threshold = 0.5;
	int n = factor * (int) round( (threshold*(1-threshold)) 
						/(delta*epsilon*epsilon*factor) );
	return n;
}

float confidenceOfConcentration_jacBits(float mode, const int nMatches,
const int nTrials, float (*posteriorCdf)(void*, float, float,
float), void *prior, float cosDelta)
{
	float x = mode;
	float y = convertHashingForJaccard(x);
	float z = y + cosDelta;
	float rl, prl=1.0, ll, pll=0.0;
	if ( z < 1.0 )
	{
		rl = convertJaccardForHashing( z ); 
		prl = posteriorCdf(prior, nMatches, nTrials, rl);
	}

	z = y - cosDelta;
	if ( z > 0 )
	{
		ll = convertJaccardForHashing( z );
		pll = posteriorCdf(prior, nMatches, nTrials, ll);
	}

	return (prl - pll);
}

float confidenceOfConcentration_cos(float mode, const int nMatches,
const int nTrials, float (*posteriorCdf)(void*, float, float,
float), void *prior, float cosDelta)
{
	float x = mode;
	float y = convertHashingForCosine(x);
	float z = y + cosDelta;
	float rl, prl=1.0, ll, pll=0.0;
	if ( z < 1.0 )
	{
		rl = convertCosineForHashing( z ); 
		prl = posteriorCdf(prior, nMatches, nTrials, rl);
	}

	z = y - cosDelta;
	if ( z > 0 )
	{
		ll = convertCosineForHashing( z );
		pll = posteriorCdf(prior, nMatches, nTrials, ll);
	}

	return (prl - pll);
}

int worstCaseHashesToPrune_jacBits(void *prior, BssOptions opts, float
(*posteriorCdf)(void*, float, float, float))
{
	float pct = 1.0;
	float nt = pct * opts.angleThreshold;
	int nHashes = opts.MIN_MIN_HASHES;
	float p;
	for ( ; nHashes < BayesLSH::BITS_MAX_MAX_HASHES; nHashes += opts.HASHES_STEP)
	{
		int s = (int) round(nt*nHashes); 
		if ( s < nHashes/2 )
			s = nHashes/2;
		p = confidenceOfConcentration_jacBits(nt, s, nHashes,
		posteriorCdf, prior, opts.cosDelta);
		if ( p > 1.0 - opts.eps2 )
			break;
	}

	if ( opts.HASHES_STEP < 256 )
	// add a number that's a multiple of HASHES_STEP and at least
	// 256
		nHashes += opts.HASHES_STEP * (int) ceil(256.0/opts.HASHES_STEP);
	else
		nHashes += opts.HASHES_STEP;

	#ifndef NDEBUG
	printf("Worst-case hashes: %d\n", nHashes);
	printf("Worst-case confidence: %f\n", p);
	#endif
//	exit(0);

	return nHashes;

}

int worstCaseHashesToPrune_cos(void *prior, BssOptions opts, float
(*posteriorCdf)(void*, float, float, float))
{
	float pct = 1.0;
	float nt = pct * opts.angleThreshold;
	int nHashes = opts.MIN_MIN_HASHES;
	float p;
	for ( ; nHashes < BayesLSH::BITS_MAX_MAX_HASHES; nHashes += opts.HASHES_STEP)
	{
/*		int s = (int) round(nt*nHashes); 
		if ( s < nHashes/2 )
			s = nHashes/2;
		p = confidenceOfConcentration_cos(nt, s, nHashes,
		posteriorCdf, prior, opts.cosDelta);
		if ( p > 1.0 - opts.eps2 )
			break;
*/	
		float true1 = opts.angleThreshold;
		float true2 = opts.angleThreshold+0.05;
		float true3 = opts.angleThreshold-0.05;
		int s1 = (int) round(true1*nHashes);
		int s2 = (int) round(true2*nHashes);
		int s3 = (int) round(true3*nHashes);
		if ( s2 > nHashes )
			s2 = nHashes;
		if ( s3 < nHashes/2 )
			s3 = nHashes/2;
		float p1, p2, p3;
		p1 = confidenceOfConcentration_cos(true1,s1,nHashes,
		posteriorCdf, prior, opts.cosDelta);
		p2 = confidenceOfConcentration_cos(true2,s2,nHashes,
		posteriorCdf, prior, opts.cosDelta);
		p3 = confidenceOfConcentration_cos(true3,s3,nHashes,
		posteriorCdf, prior, opts.cosDelta);
		p=p1;
		if ( p1 > 1.0 - opts.eps2 )
			break;
//		if ( p1 > 1.0 - opts.eps2 && p2 > 1.0 - opts.eps2 && p3 >
//		1.0 - opts.eps2 )
//			break;
			
	}

	if ( opts.HASHES_STEP < 256 )
	// add a number that's a multiple of HASHES_STEP and at least
	// 256
		nHashes += opts.HASHES_STEP * (int) ceil(256.0/opts.HASHES_STEP);
	else
		nHashes += opts.HASHES_STEP;

	#ifndef NDEBUG
	printf("Worst-case hashes: %d\n", nHashes);
	printf("Worst-case confidence: %f\n", p);
	#endif

	return nHashes;

}

float confidenceOfConcentration_jac(float mode, const int nMatches,
const int nTrials, float (*posteriorCdf)(void*, float, float,
float), void *prior, float delta)
{
	float y = mode;
	float rl = y + delta;
	float prl=1.0, ll, pll=0;
	if ( rl < 1.0 )
		prl = posteriorCdf(prior, nMatches, nTrials, rl);

	ll = y - delta;
	if ( ll > 0 )
		pll = posteriorCdf(prior, nMatches, nTrials, ll);

	return (prl - pll);
}

int worstCaseHashesToPrune_jac(void *prior, BssOptions opts, 
	float (*posteriorCdf)(void*, float, float, float))
{
	int nHashes = opts.MIN_MIN_HASHES;

	// maximum uncertainty occurs at 0.5,
	float p = 0;
	float nt = 0.5;
	if ( opts.cosThreshold > nt )
		nt = opts.cosThreshold;
	for ( ; nHashes < BayesLSH::INTS_MAX_MAX_HASHES; nHashes +=
				opts.HASHES_STEP)
	{
		// assume no. of successes is proportional to
		// similarity-threhold.
		int s = (int) round(nt*nHashes); 
		p = confidenceOfConcentration_jac(nt, s, nHashes,
		posteriorCdf, prior, opts.cosDelta);
		if ( p > 1.0 - opts.eps2 )
			break;
	}

	if ( opts.HASHES_STEP < 32 )
	// add a number that's a multiple of HASHES_STEP and at least
	// 32
		nHashes += opts.HASHES_STEP * (int) ceil(32.0/opts.HASHES_STEP);
	else
		nHashes += opts.HASHES_STEP;

	#ifndef NDEBUG
	printf("Worst-case hashes: %d\n", nHashes);
	printf("Worst-case confidence: %f\n", p);
	#endif

	return nHashes;
}


void allocPosteriorModeCache(long concentratedEnough_size, const int
*minMatches, const int maxSketchSize, const int MIN_MATCHES_STEP, float
**posteriorModeCache_base, float ***posteriorModeCache)
{
	(*posteriorModeCache_base) = (float*)
	calloc (concentratedEnough_size, sizeof(float));

	int t = maxSketchSize/MIN_MATCHES_STEP;
	(*posteriorModeCache) = (float**) malloc( sizeof(float*) * t);

	(*posteriorModeCache)[0] = *posteriorModeCache_base;
	for ( int i=1; i<t; i++ )
	{
		(*posteriorModeCache)[i] = (*posteriorModeCache)[i-1] + 
						i*MIN_MATCHES_STEP - minMatches[i-1] + 1;
	}
}

long allocConcentratedEnoughs(const int *minMatches, const int
MIN_MATCHES_STEP, const int maxSketchSize, int8_t**
concentratedEnough_base, int8_t*** concentratedEnough )
{
	long concentratedEnough_size = 0;
	for ( int i=MIN_MATCHES_STEP; i<=maxSketchSize; i+=MIN_MATCHES_STEP )
	{
		concentratedEnough_size += (i+1) - minMatches[i/MIN_MATCHES_STEP - 1];
	}

	(*concentratedEnough_base) = (int8_t*)
		calloc( concentratedEnough_size, sizeof(int8_t) );
	
	int t = maxSketchSize/MIN_MATCHES_STEP;
	(*concentratedEnough) = (int8_t**) malloc( sizeof(int8_t*) * t);

	(*concentratedEnough)[0] = *concentratedEnough_base;
	for ( int i=1; i<t; i++ )
	{
		(*concentratedEnough)[i] = (*concentratedEnough)[i-1] +
						i*MIN_MATCHES_STEP - minMatches[i-1] + 1;
	}
	
	return concentratedEnough_size;
}

int computeMinMatches(int nTrials, float (*posteriorCdf)(void*,
float, float, float), float threshold, void *prior, float eps1 )
{
	float pct = 0.7;
	int guess = nTrials/2;
	float p;
	int start = 0, end = nTrials;

	int retValue = -1;
	do
	{
		p = 1.0 - posteriorCdf(prior, guess, nTrials,
			threshold);
		if ( p > eps1 )
		{
			if ( start == end )
			{
				retValue = start-1;
				break;
			}
			end = guess - 1;
			if ( end < start )
			{
				retValue = guess -1;
				break;
			}
		}
		else if ( p < eps1 )
		{
			if ( start == end)
			{
				retValue = start;
				break;
			}
			start = guess + 1;
			if ( start > end )
			{
				retValue = guess;
				break;
			}
		}
		else
		{
			retValue = guess;
			break;
		}
		guess = (start+end)/2;
		if ( end < 0 || start > nTrials )
			errexit("Yikes! mistake in computeMinMatches start %d end %d\n", start, end);
	} while ( 1 );
	
	if ( retValue < 0 || retValue > nTrials )
	{
		if ( retValue < 0 )
			retValue = 0;
		else
			retValue = nTrials;
	}

	#ifndef NDEBUG
	printf("minMatches for %d is %d, cdf: %f", nTrials,
	retValue, posteriorCdf(prior, retValue, nTrials, threshold));
	if ( retValue < nTrials)
		printf(", cdf of %d: %f\n", retValue+1, posteriorCdf(prior,
		retValue+1, nTrials, threshold));
	else 
		printf("\n");
	#endif

	return retValue;

	return -1;
}

void computeAllMinMatches(const int maxSketchSize, const int
MIN_MATCHES_STEP, int **minMatches, long **numPruned_stages, long
**numEarlyConcentrated_stages, float (*posteriorCdf)(void*,
float, float, float), float threshold, void *prior, float eps1)
{
	int	numMinMatches = maxSketchSize/MIN_MATCHES_STEP;
	(*minMatches) = (int*)malloc(sizeof(int)*numMinMatches);
	(*numPruned_stages) = (long*) calloc( numMinMatches, sizeof(long));
	(*numEarlyConcentrated_stages) = (long*) calloc( numMinMatches, sizeof(long));

	for ( int i=0; i<numMinMatches; i++ )
	{
		int nTrials = (i+1)*MIN_MATCHES_STEP;
		(*minMatches)[i] = computeMinMatches(nTrials, posteriorCdf,
		threshold, prior, eps1);
	}
	fflush(stdout);
}

IntVector::~IntVector()
{
	if ( a != NULL )
		free(a);
	a = NULL;
}

MinwisePerms::MinwisePerms(int nh, int hr)
{
	numHashes = nh;
	nHashes = nh;
	randoms_base = (int*)malloc(sizeof(int)*numHashes*2);
	a = randoms_base;
	b = randoms_base + numHashes;

	if ( hr )
		InitRandom(time(NULL));
	else
		InitRandom(-1);

	generateRandoms(2*numHashes, randoms_base);
}

MinwisePerms::~MinwisePerms()
{
	if ( randoms_base != NULL )
		free(randoms_base);
	randoms_base = NULL;
	a = NULL;
	b = NULL;
}

MinwiseIntSketches::MinwiseIntSketches(int np, int ms, int ah,
int as)
{
	nPoints = np;
	nVectors = np;
	appendHashes = ah;
	appendAllocSize = as;
	minSketchSize = ms;

	long s = minSketchSize * ((long) nPoints);
	minSketches = (int*)calloc(s, sizeof(int));
	int max = (minSketchSize > appendHashes) ? minSketchSize :
					appendHashes;							
	mins = (int*) calloc(max, sizeof(int));

	extraSketches = new IntVector[nPoints];

}

MinwiseIntSketches::~MinwiseIntSketches()
{
	free(minSketches);
	free(mins);
	delete[] extraSketches;
}

inline int MinwisePerms::hash(int key, int nHash) const 
{
	return (int) ((a[nHash]*(long)key + b[nHash])%largePrime);
}

void sketch(const int* set, int length, int
startHash, int numHashes, const MinwisePerms *mp, int* mins, int* out)
{
	int t;
	for ( int i=0; i<length; i++ )
	{
		int feature = set[i];
		for ( int k=startHash, l=0; k<startHash+numHashes; k++,
		l++ )
		{
			t = mp->hash(feature, k);
			if ( t < mins[l] || i==0 )
			{
				out[l] = feature;
				mins[l] = t;
			}
		}
	}
}


inline int CosineSketches::hashIntegerToBit(int key, int hashId)
{
	int l = (integerToBitPerms->hash(key, hashId));
	return l%2;
}

void CosineSketches::buildMinSketches(const RectMatrix* m, 
MinwisePerms *mp)
{
	long offset = 0;
	uint8_t* dest = minSketches;
	for ( int i=0; i<m->nVectors; i++ )
	{
		sketch(m->adjncy+m->xadj[i], m->xadj[i+1]-m->xadj[i], 0,
		minSketchSize, mp, mins, minHashes);

		for ( int k=0, j=0; k<minSketchSize/8; k++ )
		{
			for ( int l=0; l<8; l++, j++ )
			{
				dest[k] <<= 1;
				if ( hashIntegerToBit(minHashes[j], j) > 0 )
					dest[k]++;
			}
		}
	
		dest += minSketchSize/8;
	}

}

void CosineSketches::appendToSketch_mp(const RectMatrix *m,
const int vId, void *mp, const int numAppendHashes)
{
	appendToSketch(m, vId, (MinwisePerms*)mp, numAppendHashes);
}

void CosineSketches::appendToSketch(const RectMatrix* m, const
int vId, MinwisePerms *mp, const int numAppendHashes)
{
	int i = vId;	
	int curr_extraLength = extraSketches[i].length;
	int curr_totalLength = minSketchSize + curr_extraLength;
	int future_totalLength = curr_totalLength +	numAppendHashes;

	if ( extraSketches[i].allocSize + minSketchSize < future_totalLength )
	{
		extraSketches[i].allocSize += numAppendHashes;
		extraSketches[i].a = (uint8_t*) realloc( extraSketches[i].a, 
								extraSketches[i].allocSize/8);
	}

	uint8_t* dest = extraSketches[i].a + curr_extraLength/8;

	sketch(m->adjncy+m->xadj[i], m->xadj[i+1]-m->xadj[i],
	curr_totalLength, numAppendHashes, mp, mins, minHashes);

	for ( int k=0, j=0; k<numAppendHashes/8; k++ )
	{
		for ( int l=0; l<8; l++, j++ )
		{
			dest[k] <<= 1;
			if ( hashIntegerToBit(minHashes[j],
			j+curr_totalLength) > 0 )
				dest[k]++;
		}
	}

	extraSketches[i].length += numAppendHashes;
}

void MinwiseIntSketches::buildMinSketches(const RectMatrix* m,
const MinwisePerms *mp)
{
	assert ( minSketchSize <= mp->nHashes );
/*	if ( minSketchSize > mp->nHashes )
	{
		printf("Yikes! minSketchSize: %d > mp->nHashes: %d\n",
		minSketchSize, mp->nHashes);
		exit(1);
	}
*/
	long offset = 0; // offset into minSketches
	int t;
	int *offset_sketches = minSketches;
	for ( int i=0; i<m->nVectors; i++ )
	{
		sketch(m->adjncy+m->xadj[i], m->xadj[i+1]-m->xadj[i], 0,
		minSketchSize, mp, mins, offset_sketches);

		offset_sketches += minSketchSize;
	}
}

void MinwiseIntSketches::appendToSketch(const RectMatrix* m,
const int vId, const MinwisePerms *mp, const int numAppendHashes)
{
	int i=vId;
	int curr_extraLength = extraSketches[i].length;
	int curr_totalLength = minSketchSize + curr_extraLength;
	int future_totalLength = curr_totalLength +	numAppendHashes;
	if ( extraSketches[i].allocSize + minSketchSize < future_totalLength )
	{
		extraSketches[i].allocSize += numAppendHashes;
		extraSketches[i].a = (int*) realloc( extraSketches[i].a, 
							extraSketches[i].allocSize*sizeof(int));
	}
	int* dest = extraSketches[i].a + curr_extraLength;

	sketch(m->adjncy+m->xadj[i], m->xadj[i+1]-m->xadj[i],
	curr_totalLength, numAppendHashes, mp, mins, dest);
	extraSketches[i].length += numAppendHashes;

}

void MinwiseIntSketches::appendToSketch(const RectMatrix* m,
const int vId, const MinwisePerms *mp)
{
	int i=vId;
	int curr_extraLength = extraSketches[i].length;
	int curr_totalLength = minSketchSize + curr_extraLength;
	int future_totalLength = curr_totalLength +	appendHashes;
	if ( extraSketches[i].allocSize + minSketchSize < future_totalLength )
	{
		extraSketches[i].allocSize += appendAllocSize;
		extraSketches[i].a = (int*) realloc( extraSketches[i].a, 
							extraSketches[i].allocSize*sizeof(int));
	}
	int* dest = extraSketches[i].a + curr_extraLength;

	sketch(m->adjncy+m->xadj[i], m->xadj[i+1]-m->xadj[i],
	curr_totalLength, appendHashes, mp, mins, dest);
	extraSketches[i].length += appendHashes;

}

RandomIntGaussians::RandomIntGaussians(int nd, int nh, 
int maxMemoryInMB, int hr)
{
	if ( nh % 8 != 0 )
	{
		errexit("Yikes! nHashes %d not a multiple of 8\n", nh);
	}

	// Going to need 2*nd bytes per hash. Ensure hashes is a
	// multiple of 32.
	int maxHashes = 32 * (int)
	round((1024*1024.0)/(nd) * (maxMemoryInMB*1.0)/(2*32));

	if ( maxHashes > nh )
		maxHashes = nh;

	#ifndef NDEBUG
	if ( maxHashes < nh )
	{
		printf("Setting maxHashes to %d; can't store more with %d MB.\n",
				maxHashes, maxMemoryInMB);
	}
	#endif

	initRandom = hr;

	setup(nd, maxHashes, nh, 16);
}

void RandomIntGaussians::fill(int startHash, int numHashes)
{
	currentStartHash = startHash;
	currentNumHashes = numHashes;
	long offset = 0;
	for ( int i=0; i<nDimensions; i++ )
	{
		getRandomIntGaussians(i, startHash, numHashes, intGaussians+offset);
		offset += numHashes;
	}
	
}

void RandomIntGaussians::setup(int nd, int mh, int nh, int rge)
{
	assert( nh % 8 == 0 );

	nDimensions = nd;
	nHashes = nh;
	maxHashes = mh;

	range = rge;
	leftLimit = (0-range)/2.0;
	rightLimit = leftLimit+range;
	floatToIntFactor = (65535)*1.0/range;
	intToFloatFactor = range*1.0/(65535);

	long s = sizeof(uint16_t) * maxHashes * nDimensions;
	intGaussians = (uint16_t*)malloc(s);

	r = gsl_rng_alloc(gsl_rng_taus2);
	if ( initRandom )
		gsl_rng_set(r, time(NULL));
	else
		gsl_rng_set(r, 0);

	tempMurmurHash= (void*)malloc(16); // 16 bytes = 128 bits

	int cacheSize = 65536;
	intToFloatCache = (float*)malloc(sizeof(float)*cacheSize);
	intToFloatCache[0] = leftLimit;
	for ( int i=1; i<cacheSize; i++ )
	{
		intToFloatCache[i] = intToFloatCache[i-1] +
		intToFloatFactor;
	}

}


RandomIntGaussians::RandomIntGaussians(int nd, int mh, int hr) 
{
	int rge = 16;
	initRandom = hr;
	setup(nd, mh, mh, rge);
}


void RandomIntGaussians::getRandomIntGaussians(int feature, int
startHash, int numHashes, uint16_t* result)
{
	float rf;
	for ( int i = startHash; i<startHash+numHashes;	i+=numGaussiansPerSeed )
	{
		for ( int j=i; j < i+numGaussiansPerSeed &&
						j < startHash+numHashes; j++ )
		{
			rf = (float) gsl_ran_gaussian_ziggurat(r, 1.0);
			if ( rf < leftLimit )
			{
				result[j-startHash] = (uint16_t) 0;
				continue;
			}
			if ( rf > rightLimit )
			{
				result[j-startHash] = (uint16_t) 65535;
				continue;
			}
			result[j-startHash] = floatToInt(rf);
		}

	}
}

inline float RandomIntGaussians::intToFloat(uint16_t a)
{
	return (float) (a*intToFloatFactor + leftLimit);
}

inline uint16_t RandomIntGaussians::floatToInt(float a)
{
	return (uint16_t) round((a-leftLimit)*floatToIntFactor);
}

RandomIntGaussians::~RandomIntGaussians()
{
	free(intGaussians);
	free(tempMurmurHash);
	gsl_rng_free(r);
	free(intToFloatCache);
}

void CosineSketches::buildBatch(const RectMatrix* m, 
	RandomIntGaussians* rig, uint8_t* tempSketches, int startHash,
	int numHashes)
{
	rig->fill(startHash, numHashes);

	long offset = startHash/8;
	for ( int i=0; i<m->nVectors; i++ )
	{
		// sums must be completely zero at this point.
		for ( int j=m->xadj[i]; j<m->xadj[i+1]; j++ )
		{
			int feature = m->adjncy[j];
			long f_offset = ((long)feature)*rig->currentNumHashes;
			if ( m->adjwgt != NULL )
			{
				wgttype wt = m->adjwgt[j];
				for ( long k=f_offset,l=0; k<f_offset+numHashes; k++,l++ )
				{
					sums[l] += wt *	rig->intToFloatCache[rig->intGaussians[k]];
				}
			}
			else
			{
				for ( long k=f_offset,l=0; k<f_offset+numHashes; k++,l++ )
				{
					sums[l] += rig->intToFloatCache[rig->intGaussians[k]];
				}
				
			}
		}

		for ( int k=0, n=0; k<numHashes/8; k++ )
		{
			for ( int l=0; l<8; l++, n++ )
			{
				tempSketches[k] <<= 1; 
				if ( sums[n] > 0 )
					tempSketches[k]++;
				sums[n] = 0; // clear out sums
			}
		}

		memcpy( (void*)(minSketches+offset), (void*)tempSketches, 
				numHashes*sizeof(uint8_t)/8);

		offset += minSketchSize/8;
	}

}

void CosineSketches::buildMinSketchesInBatches(const RectMatrix*
m, RandomIntGaussians *rig)
{
	int batchSize = rig->maxHashes;
	uint8_t* tempSketches =	(uint8_t*) calloc(batchSize/8, sizeof(uint8_t));

	if ( appendHashes > rig->maxHashes )
	{
		appendHashes = rig->maxHashes;
		appendAllocSize = appendHashes;
		#ifndef NDEBUG
		printf("setting appendHashes to %d\n", appendHashes);
		#endif
	}

	int numBatches = (int) ceil(minSketchSize*1.0/batchSize);
	#ifndef NDEBUG
	printf("Going to build minSketches in %d batches\n",
	numBatches);
	fflush(stdout);
	#endif
	int numHashes = batchSize;
	for( int done=0; done<minSketchSize; )
	{
		if ( done + numHashes > minSketchSize )
			numHashes = minSketchSize - done;
		buildBatch(m, rig, tempSketches, done, numHashes);
		done += numHashes;
	}

}

void CosineSketches::setupMinwiseBitSketches(int maxSketchSize,
int hr)
{
	int max = (minSketchSize > appendHashes) ? minSketchSize :
					appendHashes;							

	minHashes = (int*)calloc(max, sizeof(int));
	mins = (int*)calloc(max, sizeof(int));

	integerToBitPerms = new MinwisePerms(maxSketchSize, hr);
}

CosineSketches::CosineSketches(int np, int ms, int ah, int as)
{
	assert( ms % 8 == 0 );

	nPoints = np;
	appendHashes = ah;
	appendAllocSize = as;
	minSketchSize = ms;
	gaussians = 1;
	long s = (minSketchSize/8)*nPoints;
	minSketches = (uint8_t*)calloc(s,sizeof(uint8_t));
	int max = (minSketchSize > appendHashes) ? minSketchSize :
					appendHashes;							
	sums = (wgttype*) calloc(max, sizeof(wgttype));

	if ( !gaussians )
		tempFeatureBits = (uint8_t*) calloc(appendHashes,
					sizeof(uint8_t));
	else
	{
		tempRandomGaussians = (float*) calloc(appendHashes,
		sizeof(float));
		tempRandomIntGaussians = (uint16_t*) calloc(appendHashes,
		sizeof(uint16_t));
	}

		
	tempMurmurHash = (void*)malloc(16); // 16 bytes = 128 bits

	extraSketches = new BoolVector[np];
}

CosineSketches::~CosineSketches()
{
	free(minSketches);
	free(sums);
	if ( !gaussians) 
		free(tempFeatureBits);
	else
	{
		free(tempRandomGaussians);
		free(tempRandomIntGaussians);
	}
	free(tempMurmurHash);
	delete[] extraSketches;
}

BoolVector::~BoolVector()
{
	if ( a != NULL )
		free(a);
	a = NULL;
}

void CosineSketches::printSketch(int vId, int startBytes, int
numBytes)
{
	printf("Sketch for %d, bytes %d to %d\n", vId, startBytes,
	numBytes);
	int i=startBytes;
	for ( i=startBytes; i<minSketchSize/8 &&
	i<startBytes+numBytes; i++ )
	{
		printf("%d,", minSketches[vId*minSketchSize/8 + i]);
	}

	if ( i == startBytes+numBytes )
	{
		printf("\n");
		return;
	}
	
	numBytes -= (minSketchSize/8 - startBytes);
	for ( i=0; i<extraSketches[vId].length/8 &&
				i<numBytes; i++ )
	{
		printf("%d, ", extraSketches[vId].a[i]);
	}
	printf("\n\n");
	fflush(stdout);

}

void CosineSketches::appendToSketch_rig(const RectMatrix *m,
const int vId)
{
	appendToSketch(m, vId, rig);
}

void CosineSketches::appendToSketch_rig2(const RectMatrix *m,
const int vId, void *rg, const int numAppendHashes)
{
	appendToSketch(m, vId, (RandomIntGaussians*)rg,
	numAppendHashes);
}

void CosineSketches::appendToSketch(const RectMatrix *m, const
int vId, RandomIntGaussians *rg, const int numAppendHashes)
{
	const uint8_t msbmask = 128;

	int i = vId;	
	int curr_extraLength = extraSketches[i].length;
	int curr_totalLength = minSketchSize + curr_extraLength;
	int future_totalLength = curr_totalLength +	numAppendHashes;

	assert( ! (curr_totalLength < rg->currentStartHash ||
	future_totalLength >
	rg->currentStartHash+rg->currentNumHashes) );
/*	if ( 
	)
	{
		printf("Yikes! Append hashes not aligned\n");
		printf("vId:%d, curr_totalLength:%d, numAppendHashes:%d\n", 
				vId, curr_totalLength, numAppendHashes);
		printf("rg.currentStartHash:%d,	rg.currentNumHashes:%d\n", 
				rg->currentStartHash, rg->currentNumHashes);
		exit(-1);
	}
*/
	if ( extraSketches[i].allocSize + minSketchSize < future_totalLength )
	{
		extraSketches[i].allocSize += numAppendHashes;
		extraSketches[i].a = (uint8_t*) realloc( extraSketches[i].a, 
								extraSketches[i].allocSize/8);
	}

	uint8_t* dest = extraSketches[i].a + curr_extraLength/8;

	// sums should be zero at this point
	for ( int j=m->xadj[vId]; j<m->xadj[vId+1]; j++ )
	{
		int feature = m->adjncy[j];

		uint16_t* intGaussians;
		intGaussians = rg->intGaussians + (long)
		( ((long)feature)*rg->currentNumHashes
			+ (curr_totalLength-rg->currentStartHash));

		if ( m->adjwgt != NULL )
		{
			wgttype wt = m->adjwgt[j];
			for ( int k=0; k < numAppendHashes; k++ )
			{
				sums[k] += wt *	rg->intToFloatCache[intGaussians[k]];
			}
		}
		else
		{
			for ( int k=0; k < numAppendHashes; k++ )
			{
				sums[k] += rg->intToFloatCache[intGaussians[k]];
			}
		}
	}

	for ( int k=0, offset=0; k < numAppendHashes/8; k++ )
	{
		for ( int l=0; l<8; l++, offset++ )
		{
			dest[k] <<= 1; 
			if ( sums[offset] >= 0 )
				dest[k]++;
			sums[offset] = 0; // clear out sums
		}
	}

	extraSketches[i].length += numAppendHashes;
//	memset((void*)sums, 0, sizeof(wgttype)*numAppendHashes);
}

void CosineSketches::appendToSketch(const RectMatrix *m,const int
vId, RandomIntGaussians *rg)
{
	appendToSketch(m, vId, rg, appendHashes);
}

void BayesLSH::setupUnifPrior()
{
	uprior.ub = 1;
	uprior.lb = 0;
	posteriorMode = &posteriorMode_unif_dummyCache;
	posteriorCdf = &posteriorCdf_unif;
	if ( !opts.jaccard || opts.useMinwiseBitSketches )
	{
		uprior.lb = 0.5;
		posteriorCdf = &posteriorCdf_unifHalfTo1;
	}
	prior = (void*) &uprior;
}

void BayesLSH::setupBetaPrior(Beta *bPrior)
{
	prior = (void*)bPrior;
	posteriorMode = &posteriorMode_beta_dummyCache;
	posteriorCdf = &posteriorCdf_beta;
}

void BayesLSH::setupPwUnifPrior(PwUnif *pup)
{
	prior = (void*)pup;
	posteriorMode = &posteriorMode_pwUnifcache;
	posteriorCdf = &posteriorCdf_pwUnif;
	if ( opts.jaccard && !opts.useMinwiseBitSketches)
		posteriorCdf = &posteriorCdf_pwUnif0to1;
}

void BayesLSH::setupKPwUnifPrior(KPwUnif *kp)
{
	prior = (void*)kp;
	posteriorMode = &posteriorMode_KPwUnif_cache;
	posteriorCdf = &posteriorCdf_KPwUnif;
}

void BayesLSH::setupPwLinearPrior(PwLinear *pl)
{
	prior = (void*)pl;
	posteriorMode = &posteriorMode_PwLinear_cache;
	posteriorCdf = &posteriorCdf_PwLinear;
}

void BayesLSH::setupFunctionPointers()
{
	if ( opts.jaccard && !opts.useMinwiseBitSketches )
	{
		worstCaseHashesToPrune = &worstCaseHashesToPrune_jac;
		confidenceOfConcentration =
		&confidenceOfConcentration_jac;
		processPair = &BayesLSH::processPair_jac;
	}
	else
	{
		if ( opts.useMinwiseBitSketches )
		{
			worstCaseHashesToPrune =
			&worstCaseHashesToPrune_jacBits;
			confidenceOfConcentration =
			&confidenceOfConcentration_jacBits;
		}
		else
		{
			worstCaseHashesToPrune = &worstCaseHashesToPrune_cos;
			confidenceOfConcentration =
			&confidenceOfConcentration_cos;
		}
		processPair = &BayesLSH::processPair_cos;
	}
}

/* need to call setupPrior and setMinMinSketchSize on your own,
 * before this */
void BayesLSH::setup()
{
	setupFunctionPointers();
	setupParameters();
	setupCaches();
	if ( !doneGeneratingMinSketches )
		setupSketches();
	numNotConcentrated = 0;
}

int BayesLSH::getMaxMaxSketchSize()
{
	long nnz = fvs->xadj[fvs->nVectors];
	int datasetMemoryInMB = sizeof(int) * (int) round( 
		(nnz/(1024*1024.0)) );
	if ( fvs->adjwgt != NULL )
	{
		datasetMemoryInMB += sizeof(float) * (int) round( 
		(nnz/(1024*1024.0)) );
	}
	int maxSketchMemoryInMB = opts.availableMemoryInMB -
	opts.availableRandomGaussiansMemoryInMB - datasetMemoryInMB;

	#ifndef NDEBUG
	printf("datasetMemoryInMB: %d, ", datasetMemoryInMB);
	printf("available maxSketchMemoryInMB: %d\n", maxSketchMemoryInMB);
	#endif

	int maxMaxSketchSize = 32 * (int) round(
	(maxSketchMemoryInMB*1024.0*8)/(fvs->nVectors) * (1024/32) );

	return maxMaxSketchSize;
}

void BayesLSH::setMinMinSketchSize(int mms)
{
	minMinSketchSize = opts.MIN_MATCHES_STEP * (int)
	ceil(mms*1.0/opts.MIN_MATCHES_STEP);
}

int BayesLSH::getDefaultMinSketchSize()
{
	return (MIN_MIN_HASHES < minMinSketchSize) ?
	minMinSketchSize: MIN_MIN_HASHES;
	// return minSketchSize/2;
}

void BayesLSH::setupParameters()
{
	minMatchesThreshold = opts.cosThreshold;
	if ( opts.jaccard && opts.useMinwiseBitSketches )
		minMatchesThreshold =
		convertJaccardForHashing(opts.cosThreshold);
	else if ( !opts.jaccard )
		minMatchesThreshold = convertCosineForHashing(opts.cosThreshold);

	MIN_MIN_HASHES = opts.MIN_MIN_HASHES;
	HASHES_STEP = opts.HASHES_STEP;

	MIN_MATCHES_STEP = opts.MIN_MATCHES_STEP;

	if ( opts.pruneOnly )
	{
		minSketchSize = (MIN_MIN_HASHES < minMinSketchSize) ?
		minMinSketchSize : MIN_MIN_HASHES;
		maxSketchSize = minSketchSize;
		opts.onDemandSketches = 0;
	}
	else
	{
		maxSketchSize = worstCaseHashesToPrune(prior, opts, posteriorCdf);
	}

	assert ( maxSketchSize >= minMinSketchSize );
/*	if ( maxSketchSize < minMinSketchSize )
	{
		printf("Yikes! worstCaseHashesToPrune %d < minMinSketchSize %d\n"
			, maxSketchSize, minMinSketchSize);	
		exit(-1);
	}
*/
	if ( opts.jaccard )
	{
		maxHashesAtOnce = maxSketchSize;
		if ( opts.onDemandSketches )
		{
			minSketchSize = getDefaultMinSketchSize();
			APPEND_HASHES = opts.APPEND_HASHES;
		}
		else
		{
			minSketchSize = maxSketchSize;
			APPEND_HASHES = 1;
		}
	}
	else
	{
		int maxMaxSketchSize = getMaxMaxSketchSize();

		if ( maxSketchSize > maxMaxSketchSize )
			maxSketchSize = maxMaxSketchSize;

		maxHashesAtOnce = 32 * (int)
			round((1024*1024.0)/(fvs->nDimensions) *
			(opts.availableRandomGaussiansMemoryInMB*1.0)/(2*32));

		assert( maxHashesAtOnce >= HASHES_STEP );
/*		if ( maxHashesAtOnce < HASHES_STEP )
		{
			printf("Yikes! maxHashesAtOnce %d < HASHES_STEP	%d\n", maxHashesAtOnce, HASHES_STEP);
			printf("Either reduce HASHES_STEP or increase");
			printf(" opts.availableRandomGaussiansMemoryInMB\n");
			exit(-1);
		}
*/
		if ( maxHashesAtOnce > maxSketchSize )
			maxHashesAtOnce = maxSketchSize;

		if ( opts.onDemandSketches )
		{
			if ( maxHashesAtOnce < maxSketchSize )
			{
				minSketchSize = maxSketchSize - maxHashesAtOnce;
				if ( minSketchSize < getDefaultMinSketchSize() )
					minSketchSize = getDefaultMinSketchSize();
				APPEND_HASHES = opts.APPEND_HASHES;
				if ( opts.APPEND_HASHES > maxHashesAtOnce )
				{
					APPEND_HASHES = HASHES_STEP * (int)
					floor(maxHashesAtOnce*1.0/HASHES_STEP);
				}
			}
			else
			{
				minSketchSize = getDefaultMinSketchSize();
				APPEND_HASHES = opts.APPEND_HASHES;
			}
		}
		else
		{
			minSketchSize = maxHashesAtOnce;
			APPEND_HASHES = 8;
		}
	}
	
	#ifndef NDEBUG
	printf("minSketchSize: %d\n", minSketchSize);
	printf("maxSketchSize: %d\n", maxSketchSize);
	printf("maxHashesAtOnce: %d\n", maxHashesAtOnce);
	printf("MIN_MIN_HASHES: %d\n", MIN_MIN_HASHES);
	printf("HASHES_STEP: %d\n", HASHES_STEP);
	printf("APPEND_HASHES: %d\n", APPEND_HASHES);
	printf("MIN_MATCHES_STEP: %d\n", MIN_MATCHES_STEP);
	fflush(stdout);
	#endif
}

void BayesLSH::setupCaches()
{
	computeAllMinMatches(maxSketchSize, MIN_MATCHES_STEP,
	&minMatches, &numPruned_stages, &numEarlyConcentrated_stages,
	posteriorCdf, minMatchesThreshold, prior, opts.eps1);

	long concentratedEnough_size = allocConcentratedEnoughs(minMatches, 
		MIN_MATCHES_STEP, maxSketchSize, &concentratedEnough_base,
		&concentratedEnough);

	if ( !opts.uniformPrior )
		allocPosteriorModeCache(concentratedEnough_size,
		minMatches, maxSketchSize, MIN_MATCHES_STEP,
		&posteriorModeCache_base, &posteriorModeCache);

}

void BayesLSH::jaccardLSH_setupSketchesFirst()
{
	timer skTimer;
	cleartimer(skTimer);
	starttimer(skTimer);

	#ifndef NDEBUG
	printf("Going to generate sketches beforehand for jaccard beta\n");
	fflush(stdout);
	#endif

	mp = new MinwisePerms(BayesLSH::INTS_MAX_MAX_HASHES,
	opts.hash_initRandomSeed);
	minSketchSize = getDefaultMinSketchSize();
	mis = new MinwiseIntSketches(fvs->nVectors,
	minSketchSize, APPEND_HASHES, APPEND_HASHES);
	mis->buildMinSketches(fvs, mp);

	stoptimer(skTimer);
	#ifndef NDEBUG
	printf("Time to generate %d hash functions and %d*%d sketches: %.3f\n",
			BayesLSH::INTS_MAX_MAX_HASHES, fvs->nVectors, minSketchSize, gettimer(skTimer));
	fflush(stdout);
	#endif

	doneGeneratingMinSketches = 1;
}

void BayesLSH::setupSketches()
{
	timer skTimer;
	cleartimer(skTimer);
	starttimer(skTimer);
	int numHashesGenerated = 0;

	#ifndef NDEBUG
	printf("Going to generate sketches\n");
	fflush(stdout);
	#endif 

	if ( opts.jaccard )
	{
		mp = new MinwisePerms(maxSketchSize,
		opts.hash_initRandomSeed);
		if ( opts.useMinwiseBitSketches )
		{
			cs = new CosineSketches(fvs->nVectors, minSketchSize,
			APPEND_HASHES, APPEND_HASHES);
			cs->setupMinwiseBitSketches(maxSketchSize,
			opts.hash_initRandomSeed);
			cs->buildMinSketches(fvs, mp);
		}
		else
		{
			mis = new MinwiseIntSketches(fvs->nVectors,
			minSketchSize, APPEND_HASHES, APPEND_HASHES);
			mis->buildMinSketches(fvs, mp);
		}
		numHashesGenerated = maxSketchSize;
	}
	else
	{
		rig = new RandomIntGaussians(fvs->nDimensions,
		maxHashesAtOnce, opts.hash_initRandomSeed);
		cs = new CosineSketches(fvs->nVectors, minSketchSize,
		APPEND_HASHES, APPEND_HASHES);
		cs->rig = rig;
		cs->buildMinSketchesInBatches(fvs, rig);
		numHashesGenerated = minSketchSize;

		if ( opts.onDemandSketches )
		{
			int remainingHashes = maxSketchSize - minSketchSize;
			if ( maxHashesAtOnce < remainingHashes)
			{
				#ifndef NDEBUG
				printf("maxHashesAtOnce %d + minSketchSize %d =",
				maxHashesAtOnce, minSketchSize);
				printf(" %d < maxSketchSize %d\n",
				maxHashesAtOnce+minSketchSize, maxSketchSize); 
				#endif
				rig->fill(minSketchSize, maxHashesAtOnce);
				numHashesGenerated += maxHashesAtOnce;
			}
			else
			{
				rig->fill(minSketchSize,remainingHashes);
				numHashesGenerated += remainingHashes;
			}
		}
	}

	if ( !opts.onDemandSketches )
	{
		if ( opts.jaccard )
			delete mp;
		else
			delete rig;
	}

	stoptimer(skTimer);
	#ifndef NDEBUG
	printf("Time to generate %d hash functions and %d*%d sketches: %.3f\n",
			numHashesGenerated, fvs->nVectors, minSketchSize, gettimer(skTimer));
	fflush(stdout);
	#endif

}

inline int8_t BayesLSH::isConcentratedEnough(const int nMatches, const int
step)
{
	int a = step;
	if ( concentratedEnough[a][nMatches-minMatches[a]] != 0 )
		return concentratedEnough[a][nMatches-minMatches[a]];
	
	int nTrials = (step+1)*(MIN_MATCHES_STEP);
	float mode = posteriorMode(prior, nMatches, nTrials,
			posteriorModeCache, MIN_MATCHES_STEP, minMatches);
	float p = confidenceOfConcentration(mode, nMatches, nTrials,
			posteriorCdf, prior, opts.cosDelta);

	if ( p > 1.0 - opts.eps2 )
		concentratedEnough[a][nMatches-minMatches[a]] = (int8_t) 1;
	else
		concentratedEnough[a][nMatches-minMatches[a]] = (int8_t) -1;

	return concentratedEnough[a][nMatches-minMatches[a]];
	
}

float BayesLSH::processPair_cos(const int vId, const int uId)
{
	if ( fvs->xadj[vId+1] == fvs->xadj[vId] || fvs->xadj[uId+1]
	== fvs->xadj[uId] )
		return -1;

	int matches = 0;
	int pruned = 0, concentrated = 0;
	int nTrials = 0;
	int step = 0;
	uint32_t *uw = (uint32_t*) ( (uint8_t*)cs->minSketches + 
					((long) uId)*(cs->minSketchSize/8) );
	uint32_t *vw =  (uint32_t*) ( (uint8_t*)cs->minSketches + 
					((long) vId)*(cs->minSketchSize/8) );

	float sim_mode = -1;

	for ( int k=0; k<MIN_MIN_HASHES; k+=MIN_MATCHES_STEP, step++)
	{
		for ( int l=0; l<MIN_MATCHES_STEP/32; l++ )
			matches += 32 - (int) popcount(*(vw++) ^ *(uw++));
		if ( matches < minMatches[step] )
		{
			numPruned_stages[step]++;
			return -1;
		}	
	}
	
	nTrials += MIN_MIN_HASHES;
	if ( ! opts.pruneOnly )
	{
		step = nTrials/MIN_MATCHES_STEP-1;
		if ( isConcentratedEnough(matches, step) > 0 )
		{
			concentrated = 1;
			numEarlyConcentrated_stages[step]++;
		}

		while( !concentrated && nTrials < cs->minSketchSize )
		{
			for ( int l = 0; l < HASHES_STEP/32 && nTrials <
			cs->minSketchSize; l++, nTrials+=32 )
			{
				matches += 32 - (int) popcount(*(vw++) ^ *(uw++));
			}

			step = nTrials/MIN_MATCHES_STEP-1;

			if ( matches < minMatches[step] )
			{
				numPruned_stages[step]++;
				return -1;
			}

			if ( isConcentratedEnough(matches, step) > 0 )
			{
				concentrated = 1;
				numEarlyConcentrated_stages[step]++;
				break;
			}

		}
		
		if ( !concentrated && nTrials < maxSketchSize )
		{
			int curr_HASHES_STEP = HASHES_STEP;

			vw = (uint32_t*) (cs->extraSketches[vId]).a;
			uw = (uint32_t*) (cs->extraSketches[uId]).a;

			int extraTrials = 0;
					
			while ( !concentrated && nTrials < maxSketchSize )
			{
				if ( nTrials + curr_HASHES_STEP > maxSketchSize )
				{
					curr_HASHES_STEP = maxSketchSize - nTrials;
				}

				if ( (cs->extraSketches[vId]).length <
				extraTrials + curr_HASHES_STEP)
				{
					int curr_APPEND_HASHES = APPEND_HASHES;

					if ( (cs->extraSketches[vId]).length +
					minSketchSize +	curr_APPEND_HASHES > maxSketchSize)
					{
						curr_APPEND_HASHES = maxSketchSize -
						(cs->extraSketches[vId]).length -
						minSketchSize;
					}

					if ( !opts.jaccard )
						cs->appendToSketch(fvs, vId, cs->rig, curr_APPEND_HASHES);
					else
						cs->appendToSketch(fvs, vId, mp,
						curr_APPEND_HASHES);

					vw = ((uint32_t*) (cs->extraSketches[vId]).a) 
							+ extraTrials/32;
				}

				if ( (cs->extraSketches[uId]).length <
				extraTrials + curr_HASHES_STEP )
				{
					int curr_APPEND_HASHES = APPEND_HASHES;

					if ( (cs->extraSketches[uId]).length +
					minSketchSize +	curr_APPEND_HASHES > maxSketchSize)
					{
						curr_APPEND_HASHES = maxSketchSize -
						(cs->extraSketches[uId]).length -
						minSketchSize;
					}

					if ( !opts.jaccard )
						cs->appendToSketch(fvs, uId, cs->rig, curr_APPEND_HASHES);
					else
						cs->appendToSketch(fvs, uId, mp,
						curr_APPEND_HASHES);

					uw = ((uint32_t*) (cs->extraSketches[uId]).a) 
								+ extraTrials/32;
				}
				
				for ( int l = 0; l < curr_HASHES_STEP/32; l++ )
					matches += 32 - (int)popcount(*(vw++) ^ *(uw++));

				nTrials += curr_HASHES_STEP;
				extraTrials += curr_HASHES_STEP;
				step = nTrials/MIN_MATCHES_STEP-1;

				if ( matches < minMatches[step] )
				{
					numPruned_stages[step]++;
					return -1;
				}

				if ( isConcentratedEnough(matches, step) > 0 )
				{
					concentrated = 1;
					numEarlyConcentrated_stages[step]++;
					break;
				}
			}
		}
	}

	if ( !concentrated || opts.pruneOnly )
	{
		numNotConcentrated++;
		if ( !opts.binaryCosine)
		{
			sim_mode = dotProduct(fvs, vId, uId);
		}
		else
		{
			sim_mode = cosine_binary(fvs, vId, uId);
		}
		if ( sim_mode < opts.cosThreshold )
			return -1;
	}
	else
	{
		sim_mode = posteriorMode(prior, matches, nTrials,
		posteriorModeCache, MIN_MATCHES_STEP, minMatches);
		if ( !opts.jaccard )
			sim_mode = convertHashingForCosine(sim_mode);
		else
			sim_mode = convertHashingForJaccard(sim_mode);
	}

	return sim_mode;
}

float BayesLSH::processPair_jac(const int vId, const int uId)
{
	if ( fvs->xadj[vId+1] == fvs->xadj[vId] || fvs->xadj[uId+1]
	== fvs->xadj[uId] )
		return -1;

	int matches = 0;
	int pruned = 0, concentrated = 0;
	int nTrials = 0;
	int step = 0;
	int* uw = mis->minSketches + ((long) uId)*(mis->minSketchSize);
	int* vw = mis->minSketches + ((long) vId)*(mis->minSketchSize);
	float sim_mode = -1;


	for ( int k=0; k<MIN_MIN_HASHES; k+=MIN_MATCHES_STEP, step++)
	{
		for ( int l=0; l<MIN_MATCHES_STEP; l++ )
		{
			matches += (*(vw++) == *(uw++));
		}
		if ( matches < minMatches[step] )
		{
			numPruned_stages[step]++;
			return -1;
		}	
	}

	nTrials += (MIN_MIN_HASHES);
	if ( ! opts.pruneOnly )
	{
		step = nTrials/MIN_MATCHES_STEP-1;
		if ( isConcentratedEnough(matches, step) > 0 )
		{
			concentrated = 1;
			numEarlyConcentrated_stages[step]++;
		}

		for ( ; nTrials < (mis)->minSketchSize;  )
		{
			for ( int l = 0; l < (HASHES_STEP) && nTrials <
			(mis)->minSketchSize; l++, nTrials++ )
			{
				matches += (*(vw++) == *(uw++));
			}

			step = nTrials/(MIN_MATCHES_STEP)-1;

			if ( matches < (minMatches)[step] )
			{
				(numPruned_stages)[step]++;
				return -1;
			}

			if ( isConcentratedEnough(matches, step) > 0 )
			{
				concentrated = 1;
				(numEarlyConcentrated_stages)[step]++;
				break;
			}
		}

		if ( !concentrated && nTrials < maxSketchSize )
		{
			int curr_HASHES_STEP = HASHES_STEP;
			vw = (mis->extraSketches[vId]).a;
			uw = (mis->extraSketches[uId]).a;
			int extraTrials = 0;
					
			while ( !concentrated && nTrials < maxSketchSize )
			{
				if ( nTrials + curr_HASHES_STEP > maxSketchSize )
				{
					curr_HASHES_STEP = maxSketchSize - nTrials;
				}

				if ( (mis->extraSketches[vId]).length <
				extraTrials + HASHES_STEP)
				{
					int curr_APPEND_HASHES = APPEND_HASHES;

					if ( (mis->extraSketches[vId]).length +
					minSketchSize +	curr_APPEND_HASHES > maxSketchSize)
					{
						curr_APPEND_HASHES = maxSketchSize -
						(mis->extraSketches[vId]).length -
						minSketchSize;
					}

					mis->appendToSketch(fvs, vId, mp,
					curr_APPEND_HASHES);
					vw = mis->extraSketches[vId].a + extraTrials;
				}
				if ( (mis->extraSketches[uId]).length <
				extraTrials + HASHES_STEP)
				{
					int curr_APPEND_HASHES = APPEND_HASHES;

					if ( (mis->extraSketches[uId]).length +
					minSketchSize +	curr_APPEND_HASHES > maxSketchSize)
					{
						curr_APPEND_HASHES = maxSketchSize -
						(mis->extraSketches[uId]).length -
						minSketchSize;
					}

					mis->appendToSketch(fvs, uId, mp,
					curr_APPEND_HASHES);
					uw = mis->extraSketches[uId].a + extraTrials;
				}

				for ( int l=0; l < curr_HASHES_STEP; l++ )
				{
					matches += (*(vw++) == *(uw++));
				}

				nTrials += curr_HASHES_STEP;
				extraTrials += curr_HASHES_STEP;
				step = nTrials/MIN_MATCHES_STEP-1;

				if ( matches < minMatches[step] )
				{
					numPruned_stages[step]++;
					return -1;
				}

				if ( isConcentratedEnough(matches, step) > 0 )
				{
					concentrated = 1;
					numEarlyConcentrated_stages[step]++;
					break;
				}
			}
		}
	}

	if ( !concentrated )
	{
		numNotConcentrated++;
		sim_mode = jaccard(fvs, vId, uId);
		if ( sim_mode < opts.cosThreshold )
			return -1;
	}
	else
	{
		sim_mode = posteriorMode(prior, matches, nTrials,
		posteriorModeCache, MIN_MATCHES_STEP, minMatches);
	}

	return sim_mode;
}

void BayesLSH::printPruneStatistics()
{
	printf("Num. not concentrated: %ld\n", numNotConcentrated);
	long totalEarlyPruned = 0;
	int i;
	for ( i=MIN_MATCHES_STEP; i<=MIN_MIN_HASHES; i+=MIN_MATCHES_STEP )
	{
		printf("Num. early pruned at %d hashes: %ld\n",
		i,
		numPruned_stages[i/MIN_MATCHES_STEP-1]);
		totalEarlyPruned += numPruned_stages[i/MIN_MATCHES_STEP-1];
	}
	i -= MIN_MATCHES_STEP;

	printf("Num concentrated at %d trials: %ld\n", i,
				numEarlyConcentrated_stages[i/MIN_MATCHES_STEP-1]);

	long totalLatePruned = 0;
	while ( i < minSketchSize )
	{
		if ( i+HASHES_STEP > minSketchSize )
			i = minSketchSize;
		else
			i += HASHES_STEP;

		if ( i % 1 == 0 )
		{
			printf("Num pruned at %d trials: %ld\n", 
				i, numPruned_stages[i/MIN_MATCHES_STEP-1]);
			printf("Num concentrated at %d trials: %ld\n", 
				i,
				numEarlyConcentrated_stages[i/MIN_MATCHES_STEP-1]);
			totalLatePruned +=
			numPruned_stages[i/MIN_MATCHES_STEP-1];
		}
	}

	while ( i < maxSketchSize )
	{
		if ( i+HASHES_STEP > maxSketchSize )
			i = maxSketchSize;
		else
			i += HASHES_STEP;

		if ( i % 1 == 0 )
		{
			printf("Num pruned at %d trials: %ld\n", 
				i, numPruned_stages[i/MIN_MATCHES_STEP-1]);
			printf("Num concentrated at %d trials: %ld\n", 
				i,
				numEarlyConcentrated_stages[i/MIN_MATCHES_STEP-1]);
			totalLatePruned +=
			numPruned_stages[i/MIN_MATCHES_STEP-1];
		}
	}

	printf("Total early pruned: %ld\n", totalEarlyPruned);
	printf("Total late pruned: %ld\n", totalLatePruned);
	fflush(stdout);

	if ( opts.onDemandSketches )
	{
		long numExtraSketches = 0;
		if ( !opts.jaccard || opts.useMinwiseBitSketches )
		{
			for ( int i=0; i<fvs->nVectors; i++ )
				numExtraSketches += cs->extraSketches[i].length;
		}
		else
		{
			for ( int i=0; i<fvs->nVectors; i++ )
				numExtraSketches += mis->extraSketches[i].length;
		}
		printf("Avg. hashes per point: %.2f\n", minSketchSize +
		(numExtraSketches*1.0/fvs->nVectors));
		fflush(stdout);
	}

}

BayesLSH::~BayesLSH()
{
	if ( opts.onDemandSketches )
	{
		if ( opts.jaccard )
			delete mp;
		else
			delete rig;
	}

	if ( opts.jaccard && !opts.useMinwiseBitSketches )
	{
		delete mis;
	}
	else
	{
		delete cs;
	}

}

