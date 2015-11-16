
#include <all.h>

/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

//  sprintf(out2, "Error! %s", out1);

  printf("Error! %s", out1);
  fflush(stdout);

  abort();
}



#ifndef DMALLOC
/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
int *imalloc(int n, const char *msg)
{
  if (n == 0)
    return NULL;

  return (int *)GKmalloc(sizeof(int)*n, msg);
}


/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
idxtype *idxmalloc(int n, const char *msg)
{
  if (n == 0)
    return NULL;

  return (idxtype *)GKmalloc(sizeof(idxtype)*n, msg);
}


/*************************************************************************
* The following function allocates an array of float 
**************************************************************************/
float *fmalloc(int n, const char *msg)
{
  if (n == 0)
    return NULL;

  return (float *)GKmalloc(sizeof(float)*n, msg);
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *ismalloc(int n, int ival, const char *msg)
{
  if (n == 0)
    return NULL;

  return iset(n, ival, (int *)GKmalloc(sizeof(int)*n, msg));
}



/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
idxtype *idxsmalloc(int n, idxtype ival, const char *msg)
{
  if (n == 0)
    return NULL;

  return idxset(n, ival, (idxtype *)GKmalloc(sizeof(idxtype)*n, msg));
}

/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *GKmalloc(long nbytes, const char *msg)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL) 
  {
    printf("In GKmalloc for long\n");
    errexit("***Memory allocation failed for %s. Requested size: %ld bytes", msg, nbytes);
  }

  return ptr;
}

#endif

/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void GKfree(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
    free(*ptr1);
  *ptr1 = NULL;

  va_start(plist, ptr1);

  /* while ((int)(ptr = va_arg(plist, void **)) != -1) { */
  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
      free(*ptr);
    *ptr = NULL;
  }

  va_end(plist);
}            



int float_compare(const void *a, const void *b)
{
	float c = *((float*)a) - *((float*)b);
	return ( c > 0 ) ? 1 : ( (c < 0) ? -1 : 0); 
}

int compareFloatInt( const void *a, const void *b)
{
	FloatInt *aa = (FloatInt*)a;
	FloatInt *bb = (FloatInt*)b;
	if ( aa->f < bb->f )
		return -1;
	else if ( aa->f == bb->f )
		return 0;
	else
		return 1;
}

wgttype dotProduct_partial(const int nvtxs, const idxtype *xadj, const idxtype *xadj_ends, 
		const idxtype *adjncy, const wgttype *adjwgt, int a, int b)
{
	wgttype sum = 0;
	if ( xadj_ends[a] <= xadj[a] || xadj_ends[b] <= xadj[b] )
	{
		// one of the two is a null vector.
		return sum;
	}
	int a_counter = xadj[a];
	int b_counter = xadj[b];
	for ( ; a_counter < xadj_ends[a] && b_counter < xadj_ends[b]; )
	{
		if ( adjncy[a_counter] < adjncy[b_counter] )
			a_counter++;
		else if ( adjncy[b_counter] < adjncy[a_counter] )
			b_counter++;
		else 
		{
			sum += adjwgt[a_counter] * adjwgt[b_counter];
			a_counter++;
			b_counter++;
		}
	}
	return sum;

}

void generateRandoms(int numToGenerate, int* randoms )
{
	int i;
	int max = largePrime;
	int minRandomNumber = largePrime/100;
	for ( i=0; i<numToGenerate; i++ )
	{
		randoms[i] = (rand()) % max;
	//	if ( primes[i] == 0 )
		if ( randoms[i] < minRandomNumber )
			i--;
//		primes[i] = primes[i] % max;
	}
}

int checkSorted(int nvtxs, idxtype* xadj, idxtype* adjncy)
{
	for ( int i=0; i<nvtxs; i++ )
	{
		for ( int j=xadj[i]; j<xadj[i+1]-1; j++ )
		{
			if ( adjncy[j] > adjncy[j+1] )
				return 0;	
		}
	}
	return 1;
}


void sortAdjLists(int nvtxs, idxtype* xadj, idxtype *adjncy)
{
	for ( int i=0; i<nvtxs; i++ )
	{
		iidxsort(xadj[i+1]-xadj[i], adjncy+xadj[i]);
	}
}

void sortAdjLists(int nvtxs, idxtype* xadj, idxtype* adjncy, wgttype* adjwgt)
{
	int i,j;
	for(i=0;i<nvtxs;i++)
	{
		if ( xadj[i] > xadj[i+1] )
		{
			printf("Yikes! something wrong with xadjs.");
			printf("xadj of %d < xadj of %d\n", (i+1), i);
			abort();
		}
		ParallelQSort(adjncy,adjwgt,xadj[i],xadj[i+1]-1);
	}

}

void n_sortAdjLists(int nvtxs, idxtype* xadj, idxtype* adjncy, wgttype* adjwgt)
{
	int i,j;
	for(i=0;i<nvtxs;i++)
	{
		ParallelQSort(adjncy,adjwgt,xadj[i],xadj[i+1]-1);
	}
}

/* NOTE: I assume there are no self loops already. */
Matrix* addSelfLoops(Matrix* M)
{
	Matrix* ret=allocMatrix(M->nvtxs, M->nnz+M->nvtxs, 0, 0, 0);
	int i=0;

	ret->xadj[0]=0;
	for ( i=0; i<M->nvtxs; i++ )
	{
		int j,jj;
		/* NOTE: I assume the adjacency lists are sorted in increasing
		 * order! */
		int addedLoop=0;
		for ( j=M->xadj[i], jj=ret->xadj[i]; j<M->xadj[i+1]; j++)
		{
			if ( M->adjncy[j] > i && !addedLoop )
			{
				ret->adjncy[jj] = i;
				ret->adjwgt[jj++] = 1;
				addedLoop=1;
			}
			ret->adjncy[jj] = M->adjncy[j];
			ret->adjwgt[jj++] = M->adjwgt[j];
		}
		if ( !addedLoop )
		{
			ret->adjncy[jj] = i;
			ret->adjwgt[jj++] = 1;
			addedLoop=1;
		}
		ret->xadj[i+1] = jj;
	}

	return ret;
}

void removeSelfLoops(GraphType* g)
{
	int i,j,k;

	j=0;
	for ( i=0; i<g->nvtxs; i++)
	{
		k = g->xadj[i];
		g->xadj[i] = j;
		for ( ; k<g->xadj[i+1]; k++ )
		{
			if ( g->adjncy[k] == i )
				continue;
			g->adjncy[j] = g->adjncy[k];
			if ( g->adjwgt != NULL )
				g->adjwgt[j] = g->adjwgt[k];
			j++;
		}
	}
	g->xadj[g->nvtxs] = j;
	g->nedges = j;
}

Matrix* removeSelfLoops(Matrix* M, int normalizeColumns)
{
	int nvtxs=M->nvtxs,nnz=M->nnz,i, ret_counter;
	Matrix* ret=allocMatrix(nvtxs, nnz, 0, 0, 0);
	wgttype sum;

	ret->xadj[0]=0;
	ret_counter=0;
	for ( i=0; i<nvtxs; i++ )
	{
		int j;
		sum=0;
		for( j=M->xadj[i]; j<M->xadj[i+1]; j++ )
		{
			if ( M->adjncy[j] == i )
				continue;
			ret->adjncy[ret_counter]=M->adjncy[j];
			ret->adjwgt[ret_counter++]=M->adjwgt[j];
			sum+=M->adjwgt[j];
		}
		ret->xadj[i+1]=ret_counter;

		if ( normalizeColumns > 0 && sum > 0 )
		{
			for( j=ret->xadj[i]; j<ret->xadj[i+1]; j++ )
				ret->adjwgt[j]/=sum;
		}
	}

	return ret;
}

void addToMatrix(Matrix* m, const int vId, const int nbr, 
	const wgttype value)
{
	if ( m->xadj[vId+1]+1 > m->currentSize )
	{
		int new_currentSize = m->currentSize +
		m->sizeIncrement;
		m->adjncy = (idxtype*)	realloc( m->adjncy, 
					sizeof(idxtype)*new_currentSize);
		m->adjwgt = (wgttype*)	realloc( m->adjwgt,
					sizeof(wgttype)*new_currentSize);
		if ( m->adjncy == NULL || m->adjwgt == NULL )
		{
			printf("Could not allocate");
			printf(" %d ints or floats!\n",
					new_currentSize);
			abort();							
		}
		m->currentSize = new_currentSize;
	}
	m->adjncy[m->xadj[vId+1]] = nbr;
	m->adjwgt[m->xadj[vId+1]] = value;
	m->xadj[vId+1]++;
}

void freeMatrix(Matrix* a)
{
	if ( a == NULL )
		return;

	if ( a->adjncy != NULL )
	{
		free(a->adjncy);
		a->adjncy = NULL;
		free(a->xadj);
		a->xadj = NULL;
		free(a->adjwgt);
		a->adjwgt = NULL;
	}
	if ( a->adjwgtsum != NULL )
	{
		free(a->adjwgtsum);
		a->adjwgtsum = NULL;
	}
	if ( a->maxwgt != NULL )
	{
		free(a->maxwgt);
		a->maxwgt = NULL;
	}
	if ( a->attractors != NULL )
	{
		free(a->attractors);
		a->attractors = NULL;
	}
	if ( a->rmap != NULL )
	{
		free(a->rmap);
		a->rmap = NULL;
	}

	free(a);
	a=NULL;
}

RectMatrix* convertMatrixToRectMatrix(Matrix *m)
{
	RectMatrix *ret;
	ret = (RectMatrix*) malloc(sizeof(RectMatrix));

	ret->xadj = m->xadj;
	ret->adjncy = m->adjncy;
	ret->adjwgt = m->adjwgt;
	ret->nVectors = m->nvtxs;
	ret->currentSize = ret->xadj[ret->nVectors];
	
	int ndims = 0;
	int smallest = -1;
	for ( int i=0; i<ret->xadj[ret->nVectors]; i++ )
	{
		if ( smallest == -1 || ret->adjncy[i] < smallest )
			smallest = ret->adjncy[i];
		if ( ret->adjncy[i] > ndims )
			ndims = ret->adjncy[i];
	}

	ret->nDimensions = ndims + 1;
//	printf("Smallest dimension: %d\n", smallest);
//	printf("nDimensions: %d\n", ret->nDimensions);

	return ret;
}

Matrix* getTranspose(Matrix* M)
{
	if ( M->nvtxs == 0 || M->xadj[M->nvtxs] == 0 )
	{
		printf("Matrix getTranspose invoked with empty matrix\n");
		return M;
	}

	int nvtxs = M->nvtxs;
	Matrix* ret=allocMatrix(M->nvtxs, M->xadj[nvtxs], 0, 0, 0);
	idxtype *inDegrees = idxmalloc(M->nvtxs,
	"getTranspose2:inDegrees");
	for ( int i=0; i<nvtxs; i++ )
		inDegrees[i] = 0;
	for ( int i=0; i<M->xadj[nvtxs]; i++ )
	{
/*		if ( M->adjncy[i] < 0 || M->adjncy[i] > M->nvtxs-1 )
		{
			printf("Yikes! Mistake in input, edge to %d\n",
			M->adjncy[i]);
			abort();
		}
*/		inDegrees[M->adjncy[i]]++;
	}

	ret->xadj[0]=0;
	for ( int i = 0; i < nvtxs; i++)
		ret->xadj[i+1] = ret->xadj[i] + inDegrees[i];
	
	idxtype *counters = inDegrees;
	for ( int i = 0; i < nvtxs; i++)
		counters[i] = ret->xadj[i];

	for ( int i = 0; i < M->nvtxs; i++ )
	{
		for ( int j = M->xadj[i]; j < M->xadj[i+1]; j++ )
		{
			int k = M->adjncy[j];
			ret->adjncy[counters[k]] = i;
			ret->adjwgt[counters[k]] = M->adjwgt[j];
			counters[k]++;
		}
	}
	return ret;
}

RectMatrix* transpose(RectMatrix* m)
{
	RectMatrix* ret;
	ret = (RectMatrix*) malloc(sizeof(RectMatrix));

	ret->nVectors = m->nDimensions;
	ret->nDimensions = m->nVectors;
	ret->xadj = (int*)malloc(sizeof(int)*(ret->nVectors+1));
	ret->currentSize = m->xadj[m->nVectors];
	ret->adjncy =
	(int*)malloc(sizeof(int)*((long)ret->currentSize));
	if ( m->adjwgt != NULL )
		ret->adjwgt =
		(float*)malloc(sizeof(float)*((long)ret->currentSize));
	
	int* inDegrees = (int*)calloc(m->nDimensions, sizeof(int));
	for ( int i=0; i<m->xadj[m->nVectors]; i++ )
		inDegrees[m->adjncy[i]]++;
	
	ret->xadj[0]=0;
	for ( int i=1; i<=ret->nVectors; i++ )
		ret->xadj[i] = ret->xadj[i-1] + inDegrees[i-1];
	
//	printf("m->nVectors:%d\n", m->nVectors);
//	printf("m->nDimensions:%d\n", m->nDimensions);
	for ( int k=0; k<m->nVectors; k++ )
	{
		for ( int i=m->xadj[k]; i<m->xadj[k+1]; i++ )
		{
			int j =
			ret->xadj[m->adjncy[i]+1]-inDegrees[m->adjncy[i]];
			inDegrees[m->adjncy[i]]--;
			ret->adjncy[j] = k;
			if ( m->adjwgt != NULL )
				ret->adjwgt[j] = m->adjwgt[i];
		}
	}
	
	free(inDegrees);

	return ret;
}

RectMatrix* allocRectMatrix(int nVectors, int nDimensions, int
nnz, int binary)
{
	RectMatrix* ret = (RectMatrix*) malloc(sizeof(RectMatrix));

	ret->nVectors = nVectors;
	ret->nDimensions = nDimensions;

	ret->xadj = (idxtype*) malloc(sizeof(idxtype) * (nVectors+1));
	ret->adjncy = (idxtype*) malloc(sizeof(idxtype) * nnz);
	if ( !binary )
		ret->adjwgt = (wgttype*) malloc(sizeof(wgttype) * nnz);
	else
		ret->adjwgt = NULL;

	ret->currentSize = nnz;

	return ret;
}

void freeRectMatrix(RectMatrix *a)
{
	if ( a == NULL )
		return;

	if ( a->adjncy != NULL )
	{
		free(a->adjncy);
		a->adjncy = NULL;
		free(a->xadj);
		a->xadj = NULL;
		if ( a->adjwgt != NULL )
		{
			free(a->adjwgt);
			a->adjwgt = NULL;
		}
	}

	free(a);
	a=NULL;
}

Matrix* allocMatrix(int nvtxs, int nedges, int allocSum, int
						allocMax, int allocAttractor)
{
	Matrix* a;
	a=(Matrix*)malloc(sizeof(Matrix));
	a->nvtxs=nvtxs;
	a->nnz=nedges;
	a->xadj = (idxtype*) malloc( sizeof(idxtype)*(nvtxs+1) ); 
//	a->xadj = (idxtype*) GKmalloc( sizeof(idxtype)*(nvtxs+1) , 
//					"allocMatrix:xadj" );
	a->adjncy=(idxtype*) GKmalloc(sizeof(idxtype)*nedges , 
						"allocMatrix:adjncy" );
	a->adjwgt=(wgttype*)GKmalloc( sizeof(wgttype)*nedges , 
						"allocMatrix:adjwgt" );
	if ( allocSum )
		a->adjwgtsum=(wgttype*)malloc(sizeof(wgttype)*nvtxs);
	else
		a->adjwgtsum=NULL;
		
 	if ( allocSum )
		a->maxwgt=(wgttype*)malloc(sizeof(wgttype)*nvtxs);
	else
		a->maxwgt=NULL;

	if ( allocAttractor )
	{
		a->attractors=(idxtype*)malloc(sizeof(idxtype)*nvtxs);
//		printf("Allocating attractors\n");
	}
	else
		a->attractors=NULL;

	a->rmap=NULL;
	a->currentSize = nedges;
	return a;
}

void dumpMatrix(Matrix* a)
{
	FILE* fp=stderr;
	int i,j;
	for( i=0; i<a->nvtxs;i++)
	{
		fprintf(fp,"%d",i);
		for(j=a->xadj[i];j<a->xadj[i+1];j++)
			fprintf(fp," %d:%1.4f", a->adjncy[j],a->adjwgt[j]);
		fprintf(fp,"\n");
	}
}

void getPermutedGraph(idxtype* perm, idxtype* revPerm, int nvtxs, 
		int nedges, idxtype* xadj, idxtype* adjncy, idxtype* adjwgt, 
		idxtype** p_xadj, idxtype** p_adjncy, idxtype** p_adjwgt)
{
	*p_xadj = idxmalloc(nvtxs+1, "getPermutedGraph:p_xadj");
	*p_adjncy = idxmalloc(nedges, "getPermutedGraph:p_adjncy");
	if ( adjwgt != NULL )
		*p_adjwgt = idxmalloc(nedges, "getPermutedGraph:p_adjwgt");
	else
		*p_adjwgt = NULL;

	int i;
	(*p_xadj)[0]=0;
	for ( i=0; i<nvtxs; i++ )
	{
		int orgI = revPerm[i];
		int j;
		(*p_xadj)[i+1] = (*p_xadj)[i] + ( xadj[orgI+1] -
								xadj[orgI] );
		for ( j=(*p_xadj)[i]; j<(*p_xadj)[i+1]; j++ )
		{
			int orgJ = xadj[orgI] + j - (*p_xadj)[i];
			(*p_adjncy)[j] = perm[adjncy[orgJ]];
			if ( adjwgt != NULL )
				(*p_adjwgt)[j] = adjwgt[orgJ];
		}
		
		if ( adjwgt != NULL )
		{
			ParallelQSortInts( (*p_adjncy), (*p_adjwgt), (*p_xadj)[i],
			(*p_xadj)[i+1]-1 );						
		}
		else
		{
			iidxsort( (*p_xadj)[i+1]-(*p_xadj)[i],
						(*p_adjncy)+(*p_xadj)[i] );	
		}
	}

	return;
}

// this actually undoes the permutation
Matrix* permuteRowsAndColumns(Matrix* M, idxtype* perm)
{
	Matrix *ret = allocMatrix(M->nvtxs, M->nnz, 0, 0, 0);
	int i;
	idxtype* revPerm = idxmalloc(M->nvtxs,
						"permuteRowsAndColumns:revPerm");

	for ( i=0; i<M->nvtxs; i++ )
	{
		revPerm[perm[i]] = i;
	}

	ret->xadj[0]=0;
	for ( i=0; i<M->nvtxs; i++ )
	{
		int orgI = revPerm[i];
		int j;
		ret->xadj[i+1] = ret->xadj[i] + ( M->xadj[orgI+1] -
								M->xadj[orgI] );
		for ( j=ret->xadj[i]; j<ret->xadj[i+1]; j++ )
		{
			int orgJ = M->xadj[orgI] + j - ret->xadj[i];
			ret->adjncy[j] = perm[M->adjncy[orgJ]];
			ret->adjwgt[j] = M->adjwgt[orgJ];
		}
								
		ParallelQSort( ret->adjncy, ret->adjwgt, ret->xadj[i],
		ret->xadj[i+1]-1 );						
	}

	ret->nnz = M->nnz;

	free ( revPerm );

	return ret;
}

// this actually undoes the permutation
Matrix* permuteRowsAndColumns(const Matrix *M, const idxtype*
rowPerm, const idxtype* colPerm)
{
	Matrix *ret = allocMatrix(M->nvtxs, M->nnz, 0, 0, 0);
	int i;
	idxtype* revRowPerm = idxmalloc(M->nvtxs,
						"permuteRowsAndColumns:revRowPerm");
//	idxtype* revColPerm = idxmalloc(M->nvtxs, 
//						"permuteRowsAndColumns:revColPerm");

	for ( i=0; i<M->nvtxs; i++ )
	{
		revRowPerm[rowPerm[i]] = i;
//		revColPerm[colPerm[i]] = i;
	}

	ret->xadj[0]=0;
	for ( i=0; i<M->nvtxs; i++ )
	{
		int orgI = revRowPerm[i];
		int j;
		ret->xadj[i+1] = ret->xadj[i] + ( M->xadj[orgI+1] -
								M->xadj[orgI] );
		for ( j=ret->xadj[i]; j<ret->xadj[i+1]; j++ )
		{
			int orgJ = M->xadj[orgI] + j - ret->xadj[i];
			ret->adjncy[j] = colPerm[M->adjncy[orgJ]];
			ret->adjwgt[j] = M->adjwgt[orgJ];
		}
								
		ParallelQSort( ret->adjncy, ret->adjwgt, ret->xadj[i],
		ret->xadj[i+1]-1 );						
	}

	ret->nnz = M->nnz;

//	GKfree ( (void**)&revRowPerm, (void**)&revColPerm, LTERM );
	GKfree ( (void**)&revRowPerm, LTERM );

	return ret;
	
}

/*************************************************************************
* This function returns the seconds
**************************************************************************/
double seconds(void)
{
  return((double) clock()/CLOCKS_PER_SEC);
}


/*************************************************************************
* This function returns true if the a is a power of 2
**************************************************************************/
int ispow2(int a)
{
  for (; a%2 != 1; a = a>>1);
  return (a > 1 ? 0 : 1);
}


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
void InitRandom(int seed)
{
  if (seed == -1) {
#ifndef __VC__
    srand48(7654321L);  
#endif
    srand(4321);  
  }
  else {
#ifndef __VC__
    srand48(seed);  
#endif
    srand(seed);  
  }
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
int log2(int a)
{
  int i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}

void swapwgttype(wgttype* a, wgttype* b)
{
	wgttype t;
	t=*a;
	*a=*b;
	*b=t;
}

void swapidxtype(idxtype* a, idxtype* b)
{
	idxtype t;
	t=*a;
	*a=*b;
	*b=t;
}

int doTheyIntersect(idxtype* a, int a_length, idxtype* b, int
b_length, int sort)
{
	if ( sort > 0 )
	{
		iidxsort(a_length, a);
		iidxsort(b_length, b);
	}

	for ( int i=0, j=0; i<a_length && j < b_length; )
	{
		if ( a[i] == b[j] )
		{
			return 1;
		}
		else if ( a[i] > b[j] )
		{
			j++;
		}
		else 
			i++;
	}

	return 0;
}

wgttype dotProduct(RectMatrix *m, int i, int j, int sort)
{
	return dotProduct(m->adjncy+m->xadj[i], m->adjwgt+m->xadj[i],
	m->xadj[i+1]-m->xadj[i], m->adjncy+m->xadj[j],
	m->adjwgt+m->xadj[j], m->xadj[j+1]-m->xadj[j], sort);
}

float cosine_binary(RectMatrix *m, int i, int j, int sort)
{
	int l1 = m->xadj[i+1]-m->xadj[i];
	int l2 = m->xadj[j+1]-m->xadj[j];
	int intersect = sizeOfSetIntersect(m->adjncy+m->xadj[i],l1,
	m->adjncy+m->xadj[j], l2, sort);
	return ((float) (intersect*1.0)/sqrt(1.0*l1*l2));
}

wgttype dotProduct(idxtype* a, wgttype *a_wgts, int a_length,
idxtype* b, wgttype *b_wgts, int b_length, int sort)
{
	if ( sort > 0 )
	{
		ParallelQSort(a, a_wgts, 0, a_length-1);
		ParallelQSort(b, b_wgts, 0, b_length-1);
	}

	wgttype result= 0;
	for ( int i=0, j=0; i<a_length && j < b_length; )
	{
		if ( a[i] == b[j] )
		{
			result += a_wgts[i]*b_wgts[j];
			i++; j++;
		}
		else if ( a[i] > b[j] )
		{
			j++;
		}
		else 
			i++;
	}
	return result;

}

float jaccard(RectMatrix *m, int i, int j, int sort)
{
	int l1 = m->xadj[i+1] - m->xadj[i];
	int l2 = m->xadj[j+1] - m->xadj[j];
	float num=sizeOfSetIntersect(m->adjncy+m->xadj[i], l1,
	m->adjncy+m->xadj[j], l2, sort);
	return num/(l1+l2-num);
}

int sizeOfSetIntersect(idxtype* a, int a_length, idxtype* b, int
b_length, int sort)
{
	if ( sort > 0 )
	{
		iidxsort(a_length, a);
		iidxsort(b_length, b);
	}

	int numMatches = 0;
	for ( int i=0, j=0; i<a_length && j < b_length; )
	{
		if ( a[i] == b[j] )
		{
			numMatches++;
			i++; j++;
		}
		else if ( a[i] > b[j] )
		{
			j++;
		}
		else 
			i++;
	}
	return numMatches;
}


/* This function can be used only when the adjacency lists for
 * each vertex are sorted! Beware!
 **/
Matrix* add(Matrix* M0, Matrix* M1)
{
	if ( checkSorted(M0->nvtxs, M0->xadj, M0->adjncy) == 0 || 
			checkSorted(M1->nvtxs, M1->xadj, M1->adjncy) == 0)
	{
		printf("Yikes! Input to matrix add is not sorted!\n");
		fflush(stdout);
		return NULL;
	}

	Matrix* ret;
	int i,j,k; 
	long ret_size=(M0->nnz > M1->nnz)? M0->nnz: M1->nnz ;
	int init_ret_size = ret_size;
	wgttype t;

	ret=allocMatrix(M0->nvtxs, init_ret_size ,0,0,0);
	ret->xadj[0]=0;
	for(i=0;i<M0->nvtxs;i++)
	{
		int l=ret->xadj[i];
		for(j=M0->xadj[i],k=M1->xadj[i];
			 j<M0->xadj[i+1] && k<M1->xadj[i+1];)
		{
			wgttype a=M0->adjwgt[j],b=M1->adjwgt[k];

			if ( l+2 >= ret_size )
			{
				ret_size+=ret_size;
				ret->adjncy=(idxtype*)realloc(ret->adjncy,
								(ret_size)*sizeof(idxtype));
				ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
								(ret_size)*sizeof(wgttype));
				if ( ret->adjncy == NULL || ret->adjwgt == NULL )
				{
					printf("Could not allocate");
					printf(" %ld bytes!\n",ret_size*sizeof(idxtype));
					abort();
				}
			}

			if (M0->adjncy[j]==M1->adjncy[k])
			{
				ret->adjwgt[l]=a+b;
				ret->adjncy[l]=M0->adjncy[j];
				j++;k++;l++;
			}
			else 
			{
				if ( M0->adjncy[j] < M1->adjncy[k] )
				{
					ret->adjwgt[l]=a;
					ret->adjncy[l]=M0->adjncy[j];
					j++;l++;
				}
				else
				{
					ret->adjwgt[l]=b;
					ret->adjncy[l]=M1->adjncy[k];
					k++;l++;
				}
			}
		}
		if ( l+(M0->xadj[i+1]-j)+(M1->xadj[i+1]-k) >= ret_size )	
		{
			ret_size+=ret_size;
			ret->adjncy=(idxtype*)realloc(ret->adjncy,
							(ret_size)*sizeof(idxtype));
			ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
							(ret_size)*sizeof(wgttype));
			if ( ret->adjncy == NULL || ret->adjwgt == NULL )
			{
				printf("Could not allocate");
				printf(" %ld bytes!\n",ret_size*sizeof(idxtype));
				abort();
			}
		}

		for(;j<M0->xadj[i+1];j++)
		{
			ret->adjncy[l]=M0->adjncy[j];
			ret->adjwgt[l++]=M0->adjwgt[j];
		}
		for(;k<M1->xadj[i+1];k++)
		{
			ret->adjncy[l]=M1->adjncy[k];
			ret->adjwgt[l++]=M1->adjwgt[k];
		}

		ret->xadj[i+1]=l;
	}

	ret->nnz=ret->xadj[ret->nvtxs];

	return ret;
}

void convertCosineForHashing(float *sims, int n, float *out)
{
	for ( int i=0; i<n; i++ )
	//	out[i] = convertCosineForHashing(sims[i]);
		out[i] = 1.0-acos(sims[i])/PI;
}

void fitToPiecewiseUniform_knownKnee(const float * const sims,
const int numSamples, const float lb, const float ub, const float
knee, float * const pct)
{
	int numLessThanKnee = 0;
	for ( int i=0; i<numSamples; i++ )
	{
		if ( sims[i] < knee )
			numLessThanKnee++;
	}
	*pct = numLessThanKnee*1.0/numSamples;
}

void fitToPiecewiseUniform(float *orderStats, int nOrderStats,
float lb, float ub, float *knee, float *pct)
{
	// need to choose a point t in [lb,ub], 
	// such that having a piece-wise uniform distribution, with
	// one uniform distribution on [lb,t] and another one on
	// [t,ub] is the best division. 

	// i=0 is the minimum.
	float bestLL = 0;
	float bestT = -1;
	float bestCdf = 0;
	for ( int i=1; i<nOrderStats; i++ )
	{
		float t = orderStats[i];
		if ( t-lb < 1.0e-04 || ub-t < 1.0e-04 )
			continue;
		float cdf = i*1.0/(nOrderStats-1);
		float ll = cdf * log(cdf/(t-lb)) + (1.0-cdf) * 
					log((1.0-cdf)/(ub-t));
		if ( bestT == -1 || ll > bestLL )
		{
			bestT = t;
			bestLL = ll;
			bestCdf = cdf;
		}
	}
	
	*knee = bestT;
	*pct = bestCdf;
	
}

/* fit using method-of-moments */
void fitToBeta(float* samples, int numSamples, float *alpha,
float *beta, float lowerlimit, float upperlimit)
{
	float sum = 0, sumOfSquares = 0;

	for ( int i=0; i < numSamples; i++ )
	{
		sum += samples[i];
		sumOfSquares += samples[i]*samples[i];
	}

	float mean = sum/numSamples;
	float variance = sumOfSquares/numSamples - mean*mean;
	if ( lowerlimit != 0 || upperlimit != 1 )
	{
		mean = (mean-lowerlimit)/(upperlimit-lowerlimit);
		variance = variance / 
				(upperlimit-lowerlimit)*(upperlimit-lowerlimit);
	}

//	printf("fitting to beta: mean: %.6f, variance:%.6f\n", mean, variance);
	float k = ((mean*(1.0f-mean))/variance - 1.0f); 
	*alpha = (float) mean * k;
	*beta = (float) (1.0-mean)*k;

//	printf("beta: k: %.4f, alpha: %.4f, beta: %.4f\n", k, *alpha,
//		*beta);

	return;
}

float* getOrderStatistics(float *sims, int numSims, int
numOrderStatistics)
{
	float* os =
		(float*)malloc(sizeof(float)*(numOrderStatistics+1));

	qsort( sims, numSims, sizeof(float), float_compare);
	for ( int i=0; i<numOrderStatistics+1; i++)
	{
		int k = (int)round(i*numSims*1.0/numOrderStatistics);
		if ( k >= numSims )
			k--;
		os[i] = sims[k];
	}

	return os;
}

/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
int *iset(int n, int val, int *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
idxtype *idxset(int n, idxtype val, idxtype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}

wgttype* wgtset(int n, wgttype val, wgttype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
float *sset(int n, float val, float *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int iamax(int n, int *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax(int n, idxtype *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax_strd(int n, idxtype *x, int incx)
{
  int i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int samax(int n, float *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the almost maximum element in a vector
**************************************************************************/
int samax2(int n, float *x)
{
  int i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int idxamin(int n, idxtype *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int samin(int n, float *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


