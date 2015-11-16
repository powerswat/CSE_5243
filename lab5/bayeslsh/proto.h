/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.1.1.1 2011-12-14 18:42:37 venu Exp $
 *
 */

/* debug.c */
void fitToPiecewiseUniform_knownKnee(const float *, const
int, const float, const float, const float, float * const);
void fitToPiecewiseUniform(float*, int, float, float, float*,
float*);
void convertCosineForHashing(float *, int, float *);
void fitToBeta(float*, int, float*, float*, float lowerlimit=0,
float upperlimit=1);
float* getOrderStatistics(float*, int, int);
int float_compare(const void*, const void*);

void generateRandoms(int , int* );

float *allPairsSampleCandidates_wrapper(GraphType *, int, float);

void getInAndOutDegrees(int, idxtype*, idxtype*, idxtype**,
idxtype**);

/* mclutils.c */
int checkSorted(int, idxtype*, idxtype*);
Matrix* removeSelfLoops(Matrix*, int);
void freeMatrix(Matrix*);
Matrix* allocMatrix(int, int, int, int, int);
void dumpMatrix(Matrix*);
void sortAdjLists(int, idxtype*, idxtype*);
void sortAdjLists(int, idxtype*, idxtype*, wgttype*);
Matrix* permuteRowsAndColumns(const Matrix*, const idxtype*,
const idxtype*);

/* mclbase.c */
int doTheyIntersect(idxtype*, int, idxtype*, int, int sort=0);
wgttype dotProduct(RectMatrix*, int, int, int sort=0);
wgttype dotProduct(idxtype*, wgttype*, int, idxtype*, wgttype*,
				int, int sort =0);
float jaccard(RectMatrix*, int, int, int sort=0);
float cosine_binary(RectMatrix*, int, int, int sort=0);
int sizeOfSetIntersect(idxtype*, int, idxtype*, int, int sort=0);
Matrix* add(Matrix*, Matrix*);

/* myqsort.c */
void iidxsort(int, idxtype *);
void iintsort(int, int *);

/* timing.c */
double seconds(void);

/* util.c */
void errexit(char *,...);
#ifndef DMALLOC
int *imalloc(int, const char *);
idxtype *idxmalloc(int, const char *);
float *fmalloc(int, const char *);
int *ismalloc(int, int, const char *);
idxtype *idxsmalloc(int, idxtype, const char *);
//void *GKmalloc(int, char *);
void *GKmalloc(long, const char *);
#endif
void GKfree(void **,...); 
int *iset(int n, int val, int *x);
idxtype *idxset(int n, idxtype val, idxtype *x);
float *sset(int n, float val, float *x);
int iamax(int, int *);
int idxamax(int, idxtype *);
int idxamax_strd(int, idxtype *, int);
int samax(int, float *);
int samax2(int, float *);
int idxamin(int, idxtype *);
int samin(int, float *);
int idxsum(int, idxtype *);
int idxsum_strd(int, idxtype *, int);
void idxadd(int, idxtype *, idxtype *);
int charsum(int, const char *);
int isum(int, int *);
float ssum(int, float *);
float ssum_strd(int n, float *x, int);
void sscale(int n, float, float *x);
float snorm2(int, float *);
float sdot(int n, float *, float *);
void saxpy(int, float, float *, int, float *, int);
void ParallelQSortFloatsInts(wgttype*, idxtype*, int, int);
void ParallelQSortIntsUsingScores(idxtype*, idxtype*, idxtype*,
			int, int);
void ParallelQSort(idxtype*,wgttype*,int,int);
void ParallelQSortInts(idxtype*,idxtype*,int,int);
void QSortIntsUsingInts(idxtype*, idxtype*, int, int);
void ParallelQSortLongs(long*,wgttype*,int,int);
int bsearch_insertPos(idxtype*, int, int, int);
void RandomPermute(int, idxtype *, int);
void permuteDegreeOrder(int, idxtype*, idxtype*);
wgttype RandomSelect(wgttype*, int, int, int);
idxtype RandomSelectInts(idxtype*, int, int, int);
double drand48();
void srand48(long);
int ispow2(int);
void InitRandom(int);
int log2(int);

/* io.c */
void WriteBinaryMatrixInBinaryFormat(const char*, int, int*, int*);
void WriteMatrixInBinaryFormat(const char*, int, int*, int*,
float*);
void writeFloats(float*, int, char*);
void ReadMatrixFromBinary(Matrix*, char*, wgttype threshold=0);
void ReadBinaryMatrixFromBinary(Matrix*, char*);
void ReadMatrix(Matrix *, char*, wgttype threshold=0);
void readMemberships(char*, int, idxtype*, idxtype**, idxtype**,
				idxtype*);
idxtype* getNodesToComponentMap(Matrix*, int*, wgttype);
idxtype* compSizeDistribution(GraphType*, int*);
int isGraphConnected(GraphType*);
void WriteRMap(const char*, idxtype*, int);
void printHistogram(idxtype*, int, FILE*);
int readClustering(const char *, int *, int);
void ReadGraph(GraphType *, const char *, int *, int, int);
//void ReadTxtGraph(GraphType *, char *, int *, int);
void WritePartition(const char *, idxtype *, int, int);
void my_WritePartition(const char *, idxtype *, int, float);
void my_WritePartitionAddOne(const char *, idxtype *, int);
void my_WriteMappedPartition(const char *, idxtype *, idxtype *, int);
void WriteMeshPartition(const char *, int, int, idxtype *, int, idxtype *);
void WritePermutation(const char *, idxtype *, int);
int CheckGraph(GraphType *);
idxtype *ReadMesh(const char *, int *, int *, int *);
void WriteGraph(const char *, int, idxtype *, idxtype *);
void WriteTxtGraph(const char *, int, idxtype *, idxtype *);
void WriteMatrix(const char *, int, idxtype*, idxtype*, wgttype*);
void WriteGraphWithWts(const char *, int, idxtype*, idxtype*,
idxtype*);
void WriteMappedTxtGraphWithWts(const char *, int,
		idxtype*,idxtype*,idxtype*,idxtype*, int);


RectMatrix* allPairsReorderData(int, RectMatrix*, int**, float**,
float**, int**, int**, int**, int**);
RectMatrix* reorderBinaryData(const RectMatrix*, idxtype*,
idxtype*, idxtype*);
RectMatrix* reorderDimensionsAndVectors(const RectMatrix*,
	idxtype*, wgttype*, wgttype*, idxtype*);
RectMatrix* convertMatrixToRectMatrix(Matrix*);
void freeRectMatrix(RectMatrix*);
RectMatrix* allocRectMatrix(int, int, int, int binary=0);

Matrix* finalizeResultMatrix(Matrix*);

RectMatrix* transpose(RectMatrix*);
Matrix* getTranspose(Matrix*);
void addToMatrix(Matrix*, const int, const int, const wgttype);
wgttype dotProduct_partial(const int nvtxs, const idxtype *xadj, const idxtype *xadj_ends, 
		const idxtype *adjncy, const wgttype *adjwgt, int a, int
		b);
int compareFloatInt( const void *a, const void *b);
int addUniformSamples(float **sims, int nsims, float lb,
	float ub, float pct);
