/*
 * $Id: struct.h,v 1.1.1.1 2011-12-14 18:42:37 venu Exp $
 */

/* Undefine the following #define in order to use short int as the idxtype */
#define IDXTYPE_INT

#include <stdint.h>

const int largePrime = 2147483647;

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
typedef int idxtype;
#else
typedef short idxtype;
#endif

#define WGTTYPE_FLOAT
#ifdef WGTTYPE_FLOAT
typedef float wgttype;
#else
typedef double wgttype;
#endif

#define MAXIDX	(1<<8*sizeof(idxtype)-2)

/*************************************************************************
* The following data type implements a timer
**************************************************************************/
typedef double timer;


/****************
* This data structures holds a (square) sparse matrix
*/
struct matrixdef{
	int nvtxs, nnz; /* nnz stands for number of non-zero
	entries*/
	idxtype* xadj; /* xadj[i+1]-xadj[i] gives the number of
	non-zero entries in column i */
	idxtype* adjncy;
	wgttype* adjwgt; /* array that stores weights of the
	adjacency lists */
	wgttype* adjwgtsum; /* sum of adjacency weights of each node,
	or the sum of a column, basically. assigned only during the
	expandAndInflate phase */
	wgttype* maxwgt; /* max wgt of each column, assigned only
	during the expandAndInflate phase */
	idxtype* attractors; /* the row with the maximum weight in a
	column, and which has weight > 0.50 */

	idxtype* rmap; /* This is for mis_coarsen */
	long currentSize;
	int sizeIncrement;
};

typedef struct matrixdef Matrix;

struct rectmatrixdef{
	int nVectors;
	int nDimensions;
	idxtype *xadj, *adjncy;
	wgttype *adjwgt;

	int currentSize;
	int sizeIncrement;
};

typedef struct rectmatrixdef RectMatrix;

struct AllPairsOptsDef
{
	float threshold;
};

typedef struct AllPairsOptsDef AllPairsOpts;

class AllPairsSimSearch{
	
	long addPairsForVertex(int, Matrix*);
	Matrix* simPairs_reordered();
	void initSimPairs();
	void finalizeSimPairs();
	Matrix* binarySimPairs();
	long binaryAddPairsForVertex(int, Matrix*);
	
	RectMatrix* reordered;
	idxtype* inDegrees, *outDegrees;
	wgttype* maxInWeights;
	wgttype* maxOutWeights;
	idxtype* descMaxOutWeightOrder, *ascOutDegreeOrder;

	int nVectors, nDimensions;
	long numTotalCandidates2, numTotalCandidates1;

	int binary;
	int binaryCosine;

	idxtype *xadj, *adjncy;
	wgttype *adjwgt;
	idxtype *invIndexAdjncy, *invIndex_starts, *invIndex_ends,
	*xadj_ends, *candidates;
	wgttype* invIndexWgts, *temp;
	float threshold;

	public:

	RectMatrix *featureVectors;
	AllPairsOpts opts;

	AllPairsSimSearch(RectMatrix *fvs, AllPairsOpts opts) :
	featureVectors(fvs), opts(opts){}
	AllPairsSimSearch(RectMatrix *fvs, AllPairsOpts opts, int
	b) :
	featureVectors(fvs), opts(opts), binary(b) 
	{}

	Matrix* run();
};

// forward declaration of struct
#include <gsl/gsl_rng.h>
class RandomIntGaussians
{
	public:

	gsl_rng *r;
	void *tempMurmurHash;
	uint16_t* intGaussians;
	int range; 
	// if we only want gaussian samples between	-8 and 
	// +8, set range to 16.
	float leftLimit, rightLimit;
	float floatToIntFactor;
	float intToFloatFactor;
	float *intToFloatCache;

	int initRandom;

	int nDimensions;
	int nHashes;

	int maxHashes;
	int currentStartHash, currentNumHashes;

	void setup(int, int, int, int);
	void fill(int, int);
	RandomIntGaussians(int, int, int);
	RandomIntGaussians(int, int, int, int);
	float intToFloat(uint16_t);
	uint16_t floatToInt(float);
	void getRandomIntGaussians(int, int, int, uint16_t*);
	~RandomIntGaussians();

	const static int numGaussiansPerSeed = 32;
};

class Beta
{
	public:
	float alpha, beta;
};

class PwUnif{
	public:
	float unifPct;
	float p1, p2;
	float lb, ub;
	float changePt;
	float logChangePt;
	float log1minusChangePt;
	float logp1, logp2;
};

class TwoPwUnif{
	public:
	float p1, p2, p3;
	float lb, ub;
	float changePt1, changePt2;
};

class KPwUnif{

	float *base;
	public:
	float *ps;
	float *logps;
	float lb, ub;
	float *changePts;
	float *logChangePts;
	float *log1minusChangePts;
	int k;

	void print();
//	KPwUnif(float *, const int, const int, const
//	float, const float);
	KPwUnif(float *, const int, const int, const
	float, const float, const float);
//	~KPwUnif();
};

class PwLinear{

	float *base;
	public:
	float *as, *bs;
	int k;
	float lb, ub;
	float *changePts;

	PwLinear(float *, const int, const int, const float, const
	float, const float);

};

class Unif{
	public:
	float lb, ub;
};

class BoolVector
{
	public:
		uint8_t *a;
		int length; // measured in bits, must be a multiple of 8.
		int allocSize; // also measured in bits.

		BoolVector():a(NULL), length(0), allocSize(0) {}
		~BoolVector();
};

class IntVector
{
	public:
		int *a;
		int length; 
		int allocSize;

		IntVector():a(NULL), length(0), allocSize(0){}
		~IntVector();
};

class MinwisePerms
{
	public:
	int *randoms_base;
	int *a, *b;
	static const int largePrime = 2147483647;
	
	int numHashes, nHashes;

	int initRandom;

	int hash(int key, int nHash) const;

	MinwisePerms(int, int);
	~MinwisePerms();
};

class MinwiseIntSketches
{
	public:

	int appendHashes; 
	int appendAllocSize;

	int *hashes;

	int nVectors, nPoints;
	IntVector* extraSketches;
	int* minSketches;
	int minSketchSize;

	int *mins;

//	MinwisePerms *mp;

	void buildMinSketches(const RectMatrix*, const MinwisePerms *);
	void appendToSketch(const RectMatrix*, const int, const
	MinwisePerms*, const int );
	void appendToSketch(const RectMatrix*, const int, const
	MinwisePerms*);
//	void sketch(const int*, int, int, int, const MinwisePerms*, int*, int*);
	
	MinwiseIntSketches(int, int, int, int);
	~MinwiseIntSketches();
};

struct bssOptionsDef{
	float cosThreshold;
	float angleThreshold;
	float cosDelta;
	float eps1;
	float eps2;
	int numSamples;
	int onDemandSketches;
	int gaussians;
	int useIntGaussians;
	int uniformPrior;
	int betaPrior;
	int pwLinearPrior;
	int num_os;
	int availableMemoryInMB;
	int availableRandomGaussiansMemoryInMB;

	int binaryCosine;
	int jaccard;
	int useMinwiseBitSketches;

	int backoffPrior;
	float backoffPct;

	int additionalPruning;
	int doublePruning; // compute partial sums also?
	int useWithLSH;
	float lsh_recall;

	int pruneOnly;

	int lsh_useExact;

	int hash_initRandomSeed;

	int MIN_MATCHES_STEP;
	int MIN_MIN_HASHES;
	int HASHES_STEP;
	int APPEND_HASHES;
//	int sketchesMemoryInMB; 
};

typedef struct bssOptionsDef BssOptions;


class CosineSketches
{
	public:

	int appendHashes; 
	int appendAllocSize;
	//const static int appendHashes = 64;
	//const static int appendAllocSize = appendHashes*3;

	int gaussians;

	int *minHashes, *mins;
	MinwisePerms *integerToBitPerms;

	wgttype *sums;
	uint8_t *tempFeatureBits;
	float* tempRandomGaussians;
	uint16_t* tempRandomIntGaussians;
	void *tempMurmurHash;

	int nPoints;
	BoolVector *extraSketches;
	uint8_t *minSketches; 
	int minSketchSize; // must be a multiple of 8.

	RandomIntGaussians *rig;

	int hashIntegerToBit(int, int);
	void appendToSketch(const RectMatrix*, const int,
	MinwisePerms*, const int);
	void appendToSketch_mp(const RectMatrix*, const int, void*,
	const int);

	void setupMinwiseBitSketches(int, int);

	void buildMinSketches(const RectMatrix*, RandomIntGaussians *);

	void buildMinSketches(const RectMatrix*, MinwisePerms *);

	void buildBatch(const RectMatrix*, RandomIntGaussians*,
	uint8_t*, int, int);
	void buildMinSketchesInBatches(const RectMatrix*, RandomIntGaussians*);

	void (CosineSketches::*old_appendToSketchPtr) (const RectMatrix*,
	const int);

	void (CosineSketches::*appendToSketchPtr) (const RectMatrix*,
	const int, void*, const int);

	void appendToSketch_rig(const RectMatrix*, const int);
	void appendToSketch_rb(const RectMatrix*, const int);
	void appendToSketch_rg(const RectMatrix*, const int);
	void appendToSketch_prg(const RectMatrix*, const int);

	void appendToSketch_rig2(const RectMatrix*, const int, void*,
	const int);

	void appendToSketch(const RectMatrix*, const int,
	RandomIntGaussians*, const int);

	void appendToSketch(const RectMatrix*, const int,
	RandomIntGaussians*);

	CosineSketches(int, int, int, int);
	void printSketch(int, int, int);
	~CosineSketches();
};

class BayesLSH{
	public:

	static const int BITS_MIN_MATCHES_STEP = 32;
	static const int INTS_MIN_MATCHES_STEP = 32;
	static const int BITS_MAX_MAX_HASHES = 1024*12;
	static const int INTS_MAX_MAX_HASHES = 1024*2;

	
	BssOptions opts;
	int *minMatches;
	long *numPruned_stages;
	long *numEarlyConcentrated_stages;

	float **posteriorModeCache;
	float *posteriorModeCache_base;

	int8_t **concentratedEnough;
	int8_t *concentratedEnough_base;

	float minMatchesThreshold;

	int stopAppending;
	int consumedSketchMemoryInMB;

	int maxSketchSize, minSketchSize, maxHashesAtOnce;
	int MIN_MIN_HASHES, HASHES_STEP, MIN_MATCHES_STEP,
	APPEND_HASHES;
	long numNotConcentrated;

	int minMinSketchSize;

	int doneGeneratingMinSketches;

	CosineSketches *cs;
	RandomIntGaussians *rig;
	RectMatrix* fvs;
	MinwiseIntSketches *mis;
	MinwisePerms *mp;

	Unif uprior;

	void *prior;
	float (*posteriorMode) (void *, int, int, float**, const int,
	const int*);
	float (*posteriorCdf) (void *, float, float, float);

	int (*worstCaseHashesToPrune)(void *, BssOptions, float
	(*)(void*, float, float, float));
	float (*confidenceOfConcentration)(float, const int, const
	int, float (*)(void*, float, float, float), void*, float);
	int8_t isConcentratedEnough(const int, const int);

	float processPair_cos(const int, const int);
	float processPair_jac(const int, const int);
	float (BayesLSH::*processPair)(const int, const int);

	void setMinMinSketchSize(int);
	int getDefaultMinSketchSize();
	int getMaxMaxSketchSize();

	void setupFunctionPointers();
	void setupParameters();
	void setupCaches();
	void setupSketches();

	void setupUnifPrior();
	void setupPwUnifPrior(PwUnif*);
	void setupBetaPrior(Beta*);
	void setupKPwUnifPrior(KPwUnif*);
	void setupPwLinearPrior(PwLinear *);
	
	void setup();

	void printPruneStatistics();

	void jaccardLSH_setupSketchesFirst();

	BayesLSH(RectMatrix *fvs, BssOptions opts): fvs(fvs),
	opts(opts), minMinSketchSize(0),
	doneGeneratingMinSketches(0),
	MIN_MIN_HASHES(opts.MIN_MIN_HASHES) {}
	~BayesLSH();

};


class LSH{
	int nPoints, nVectors, nDimensions;

	idxtype *xadj,  *adjncy, *invIndexAdjncy,
	*invIndex_starts, *invIndex_ends, *inDegrees;
	wgttype *adjwgt, *invIndexWgts;

	CosineSketches *cs;
	MinwiseIntSketches *mis;
	MinwisePerms *mp;

	RectMatrix *fvs;
	int nHashes;
	float cosThreshold;
	long numDotProducts;

	int cardinalityOfHash;
	int sizeOfBlock;
	int maxHashes;
	int numValuesPerBlock;
	int numBlocks;
	float recall;
	int* invArrays_base;
	int** invArrays;
	int* invArraysXadj_base;
	int** invArraysXadjs;

	BssOptions opts;

	BayesLSH *blsh;

	void addPairsForVertex(const int, Matrix*, uint8_t*);
	Matrix* simPairs();
	void buildLSHIndexes();
//	Matrix* processIndexes();
	Matrix* processIndexes( );
	void freeLSHIndexes();

	void (LSH::*processPair)(const int, const int, Matrix*,
	RectMatrix*);

	void processPair_blshWrapper(const int, const int, Matrix*,
	RectMatrix*);
	void processPair_jac_exact(const int, const int, Matrix*,
	RectMatrix*);
	void processPair_cos_exact(const int, const int, Matrix*,
	RectMatrix*);
	void processPair_jac(const int, const int, Matrix*,
	RectMatrix*);
	void processPair_cos(const int, const int, Matrix*,
	RectMatrix*);
	
	public:

	LSH(RectMatrix* fvs, int nHashes, float ct) : fvs(fvs),
	nHashes(nHashes), cosThreshold(ct) {}

	LSH(RectMatrix* fvs, float threshold, int maxHashes, int sb, float rec) :
	fvs(fvs), cosThreshold(threshold), maxHashes(maxHashes), recall(rec){}

	LSH(RectMatrix* fvs, BssOptions opts, int maxHashes,
	float rec): fvs(fvs), opts(opts), cosThreshold(opts.cosThreshold),
	maxHashes(maxHashes), recall(rec), nPoints(fvs->nVectors),
	nVectors(fvs->nVectors), nDimensions(fvs->nDimensions),
	adjncy(fvs->adjncy), xadj(fvs->xadj) {}

	Matrix* runRecall();
//	void runRecall();
//	Matrix* run();

//	Matrix* lshPlusBayesLSH(RectMatrix*, BssOptions);
//	Matrix* lshPlusBayesLSH();
	Matrix* bayesLSH();

	static int getNumBlocksForRecall(float, float, int);
	static const float marginPct = 0.9;
};

class AllPairsRealCosine{

	int nPoints;
	int nDimensions;

	idxtype *xadj, *xadj_ends, *adjncy, *invIndexAdjncy,
	*invIndex_starts, *invIndex_ends;
	wgttype *adjwgt, *invIndexWgts;

	long numDotProducts;
	long numTotalCandidates;

	BayesLSH *blsh;

	Matrix* simPairs();

	void addPairsForVertex(const int, Matrix*, uint8_t*);
	void pruneAddPairsForVertex(const int, Matrix*, uint8_t*,int*);
	void doublePruneAddPairsForVertex(const int, Matrix*,
	uint8_t*, int*, float*);

	void (AllPairsRealCosine::*processPair)(const int, const int,
	Matrix*);

	void processPair_blshWrapper(const int, const int, Matrix*);

//	Matrix *fvs, *reordered ;
	RectMatrix *fvs, *reordered ;
	Beta betaPrior;
	PwUnif puPrior;
	Unif unifPrior;
	PwLinear *plprior;
	KPwUnif *kpprior;
	
	timer priorTimer, skTimer, simTimer;

	BssOptions opts;

	wgttype* maxInWeights, *maxOutWeights;
	int* inDegrees;


	public:

	AllPairsRealCosine(RectMatrix *fvs, BssOptions opts) : fvs(fvs),
	opts(opts){}

	Matrix* bayesLSH();
};

class AllPairsBinary{
	int nPoints, nVectors;
	int nDimensions;

	idxtype *xadj, *xadj_ends, *adjncy, *invIndexAdjncy,
	*invIndex_starts, *invIndex_ends, *inDegrees, *outDegrees,
	*ascOutDegreeOrder;
	float threshold;

	long numDotProducts;
	long numTotalCandidates;

	BayesLSH *blsh;

	Matrix* simPairs();
	void pruneAddPairsForVertex(const int, Matrix*, uint8_t*,
	int*);
	void addPairsForVertex(const int, Matrix*, uint8_t*);
	void doublePruneAddPairsForVertex(int, Matrix*, uint8_t*,
	int*, float*);

	void (AllPairsBinary::*processPair)(const int, const int,
	Matrix*);

	void processPair_blshWrapper(const int, const int, Matrix*);
	void processPair_noBlsh(const int, const int, Matrix*);

	RectMatrix *fvs, *reordered ;
	Beta betaPrior;
	PwUnif puPrior;
	Unif unifPrior;
	KPwUnif *kprior;
	PwLinear *plprior;

	timer priorTimer, skTimer, simTimer;
	BssOptions opts;

	public:

	Matrix* bayesLSH();
	AllPairsBinary(RectMatrix *fvs, BssOptions opts) : fvs(fvs),
	opts(opts){}

};

typedef unsigned int UInt32;
//typedef int UInt32;

/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
struct graphdef {
  idxtype *gdata, *rdata;	/* Memory pools for graph and refinement data.
                                   This is where memory is allocated and used
                                   the rest of the fields in this structure */

  int nvtxs, nedges;		/* The # of vertices and edges in the graph */
  idxtype *xadj;		/* Pointers to the locally stored vertices */
  idxtype *vwgt;		/* Vertex weights */
  idxtype *vsize;		/* Vertex sizes for min-volume formulation */
  idxtype *adjncy;		/* Array that stores the adjacency lists of nvtxs */
  idxtype *adjwgt;		/* Array that stores the weights of the adjacency lists */

  idxtype *adjwgtsum;		/* The sum of the adjacency weight of each vertex */

  idxtype *label;

  idxtype *cmap;
  
  /* Venu: my addition
  * Refine map: maps vertices in coarse graph to vertices in
   * refined graph. Two maps needed as 2 vertices are mapped to
   * one vertex in the coarse graph */
  idxtype *rmap1;
  idxtype *rmap2;
  idxtype *numDescendants;
  /* indicates if this graph is the original
  (i.e. the most refined) graph*/
  int isOrgGraph; 
  int isDirected; 
  wgttype *pagerank;

  /* Partition parameters */
  int mincut, minvol;
  idxtype *where, *pwgts;
  int nbnd;
  idxtype *bndptr, *bndind;

  /* Additional info needed by the MOC routines */
  int ncon;			/* The # of constrains */ 
  float *nvwgt;			/* Normalized vertex weights */
  float *npwgts;		/* The normalized partition weights */

  struct graphdef *coarser, *finer;
};

typedef struct graphdef GraphType;

struct floatintDef{
	float f;
	int i;
};
typedef struct floatintDef FloatInt;

typedef struct {
	int i, j;
}IntInt;




