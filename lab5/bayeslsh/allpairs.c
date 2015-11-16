#include <all.h>

long AllPairsSimSearch::binaryAddPairsForVertex(int vId, Matrix *result)
{
	long numDotProducts = 0;	
	
	wgttype remScore = outDegrees[vId];
	
	wgttype minSize;
	if ( binaryCosine)
		minSize = (wgttype) ((double)threshold*threshold*outDegrees[vId]);
	else
		minSize = (wgttype) threshold * outDegrees[vId];

	result->xadj[vId+1] = result->xadj[vId];

	int numCandidates = 0;
	for ( int i=xadj[vId+1]-1; i >= xadj[vId]; i-- )
	{
		int k = adjncy[i];

		// remove those vectors from the inverted index that do
		// not have sufficient size.
		for ( int j = invIndex_starts[k]; j < invIndex_ends[k]; j++ )
		{
			int uId = invIndexAdjncy[j];
			if ( outDegrees[uId] >= minSize )
				break;
			invIndex_starts[k]++;
		}

		for ( int j = invIndex_starts[k]; j < invIndex_ends[k];	j++ )
		{
			int uId = invIndexAdjncy[j];

			if ( remScore >= minSize || temp[uId] > 0 )
			{
				if ( temp[uId] <= 0 )
				{
					candidates[numCandidates++] = uId;
				}
				temp[uId]++; 
			}
		}
		remScore--; 
	}

	numTotalCandidates2 += numCandidates;

	for ( int i=0; i < numCandidates; i++ )
	{
		int c = candidates[i];

		int maxIntersect = (int) round(temp[c]) + xadj_ends[c]-xadj[c];
		double upperBound;
		double binaryCosineDen = 0;

		if ( binaryCosine )
		{
			binaryCosineDen = sqrt(outDegrees[vId]*1.0*outDegrees[c]);
			upperBound = maxIntersect*1.0/binaryCosineDen;
		}
		else
		{
			upperBound = maxIntersect*1.0/(outDegrees[vId] +
			outDegrees[c] - maxIntersect);
		}

		if ( upperBound < threshold )
		{
			temp[c] = 0;
			continue;
		}
	
		double sum = temp[c] + sizeOfSetIntersect(adjncy +
		xadj[vId], xadj_ends[vId]-xadj[vId], adjncy+xadj[c],
		xadj_ends[c]-xadj[c]);
		float sim;
		if ( binaryCosine )
		{
			sim = sum/binaryCosineDen;
		}
		else
		{
			sim = sum/(outDegrees[vId] + outDegrees[c] - sum);
		}
			
		numDotProducts++;
		if ( sim >= threshold )
		{
			addToMatrix(result, vId, c, sim);
		}

		temp[c] = 0;
	}

	return numDotProducts;
}

long AllPairsSimSearch::addPairsForVertex(int vId, Matrix *result)
{
	long numDotProducts = 0;	
	
	wgttype remScore = 0;
	for ( int i=xadj[vId]; i < xadj[vId+1]; i++ )
		remScore += adjwgt[i] * maxInWeights[adjncy[i]];
	wgttype minSize = (wgttype) (threshold / maxOutWeights[vId]);

	result->xadj[vId+1] = result->xadj[vId];
	if ( remScore < threshold || minSize > nDimensions)
	{
		return numDotProducts;
	}

	int numCandidates = 0;
	for ( int i=xadj[vId+1]-1; i >= xadj[vId]; i-- )
	{
		int k = adjncy[i];

		// remove those vectors from the inverted index that do
		// not have sufficient size.
		for ( int j = invIndex_starts[k]; j < invIndex_ends[k]; j++ )
		{
			int uId = invIndexAdjncy[j];
			if ( xadj[uId+1]-xadj[uId] >= minSize )
				break;
			invIndex_starts[k]++;
		}
		
		for ( int j = invIndex_starts[k]; j < invIndex_ends[k];	j++ )
		{
			int uId = invIndexAdjncy[j];

			if ( remScore >= threshold || temp[uId] > 0 )
			{
				wgttype add = invIndexWgts[j] * adjwgt[i];
				if ( temp[uId] <= 0 && add > 0)
				{
					candidates[numCandidates++] = uId;
				}
				temp[uId] += add;
			}
		}
		remScore -= adjwgt[i] * maxInWeights[k];
	}

	numTotalCandidates2 += numCandidates;

	for ( int i=0; i < numCandidates; i++ )
	{
		int c = candidates[i];

		wgttype upperbound = amin(xadj_ends[c]-xadj[c],
		xadj_ends[vId]-xadj[vId]) * maxOutWeights[vId] *
		maxOutWeights[c] + temp[c];
		if ( upperbound < threshold )
		{
			temp[c] = 0;
			continue;
		}
	
		wgttype sum = temp[c] + 
			dotProduct_partial(nVectors, xadj, xadj_ends, adjncy, adjwgt,
			vId, c);
		numDotProducts++;
		if ( sum >= threshold )
		{
			addToMatrix(result, vId, c, sum);
		}
		temp[c] = 0;
	}

	return numDotProducts;
}

void AllPairsSimSearch::initSimPairs()
{
	xadj = reordered->xadj;
	adjncy = reordered->adjncy;
	adjwgt = reordered->adjwgt;

	threshold = opts.threshold;

	long size = sizeof(idxtype) * ((long) xadj[nVectors]);
	invIndexAdjncy = (idxtype*) malloc(size);

	if ( adjwgt != NULL )
	{
		size = sizeof(wgttype) * ((long) xadj[nVectors]);
		invIndexWgts = (wgttype*) malloc(size);
	}

	size = sizeof(idxtype) * nDimensions;
	invIndex_starts = (idxtype*) malloc(size);
	invIndex_ends = (idxtype*) malloc(size);

	invIndex_starts[0] = 0;
	for ( int i = 1; i < nDimensions; i++ )
		invIndex_starts[i] = invIndex_starts[i-1] +	inDegrees[i-1];
	for ( int i = 0; i < nDimensions; i++ )
		invIndex_ends[i] = invIndex_starts[i];

	xadj_ends = idxmalloc(nVectors, "xadj ends");
	for ( int i=0; i < nVectors; i++ )
		xadj_ends[i] = xadj[i+1];

	temp = (wgttype*)calloc(nVectors,sizeof(wgttype));
	candidates = idxmalloc(nVectors, "candidates");

}

void AllPairsSimSearch::finalizeSimPairs()
{
	free(invIndexAdjncy);
	if ( adjwgt != NULL )
		free(invIndexWgts);
	free(invIndex_starts);
	free(invIndex_ends);
	free(xadj_ends);

	free(temp);
	free(candidates);
}

Matrix* AllPairsSimSearch::binarySimPairs()
{
	int init_ret_size = nVectors * 50;
	Matrix *result = allocMatrix(nVectors, init_ret_size, 0, 0, 0);
	result->sizeIncrement = nVectors*20;
	result->xadj[0] = 0;

	timer simtimer;
	cleartimer(simtimer);
	starttimer(simtimer);
	long numDotProducts = 0;
	numTotalCandidates2 = 0;
	long totalIndexed = 0;
//	printf("0..");
//	fflush(stdout);
	printf("Progress Indicator (no. of records):");
	for( int i = 0; i < nVectors; i++ )
	{
		float b = 0;
		numDotProducts += binaryAddPairsForVertex(i, result );
		int j;
		for ( j = xadj[i]; j < xadj_ends[i]; j++ )
		{
			b++;
			if ( b/outDegrees[i] >= threshold )
				break;
		}
		// truncate the rest of the vector, since we're going to
		// index it.
		xadj_ends[i] = j;

		for ( ; j < xadj[i+1]; j++ )
		{
			invIndexAdjncy[invIndex_ends[adjncy[j]]] = i;
			invIndex_ends[adjncy[j]]++;
			totalIndexed++;
		}

		if ( i % 20000 == 0 )
		{
			printf("%d..",i);
			fflush(stdout);
		}
	}
	printf("\n");
	stoptimer(simtimer);
//	printf("Sim. time: %.3f\n", gettimer(simtimer));
	printf("Num. total candidates: %ld\n", numTotalCandidates2);
	printf("Num. dot products: %ld\n", numDotProducts);
//	printf("Total indexed: %ld\n", totalIndexed);

	printf("Num. similar pairs: %d\n",
	result->xadj[result->nvtxs]);
	result->nnz = result->xadj[result->nvtxs];
	fflush(stdout);

	return result;
}

Matrix* AllPairsSimSearch::simPairs_reordered()
{
	int init_ret_size = nVectors * 50;
	Matrix *result = allocMatrix(nVectors, init_ret_size, 0, 0, 0);
	result->sizeIncrement = nVectors*20;
	result->xadj[0] = 0;

	timer simtimer;
	cleartimer(simtimer);
	starttimer(simtimer);
	long numDotProducts = 0;
	numTotalCandidates2 = 0;
//	printf("0..");
//	fflush(stdout);
	printf("Progress Indicator (no. of records):");
	for( int i = 0; i < nVectors; i++ )
	{
		wgttype b = 0;
		numDotProducts += addPairsForVertex(i, result );
		int j;
		wgttype new_maxOutWeight = 0;
		for ( j = xadj[i]; j < xadj_ends[i]; j++ )
		{
			wgttype min = maxInWeights[adjncy[j]];
			if ( min > maxOutWeights[i] )
				min = maxOutWeights[i];
			b = b + min * adjwgt[j];
			if ( b >= threshold )
			{
				break;
			}
			if ( adjwgt[j] > new_maxOutWeight )
				new_maxOutWeight = adjwgt[j];
		}
		// truncate the rest of the vector, since we're going to
		// index it.
		xadj_ends[i] = j;
		// update maxOutWeights[i]
		maxOutWeights[i] = new_maxOutWeight;
		for ( ; j < xadj[i+1]; j++ )
		{
			invIndexAdjncy[invIndex_ends[adjncy[j]]] = i;
			invIndexWgts[invIndex_ends[adjncy[j]]] = adjwgt[j];
			invIndex_ends[adjncy[j]]++;
		}

		if ( i % 20000 == 0 )
		{
			printf("%d..",i);
			fflush(stdout);
		}
	}
	printf("\n");
	stoptimer(simtimer);
//	printf("Sim. time: %.3f\n", gettimer(simtimer));
	printf("Num. dot products: %ld\n", numDotProducts);
	printf("Num. total candidates: %ld\n", numTotalCandidates2);
	printf("Num. similar pairs: %d\n",
	result->xadj[result->nvtxs]);
	result->nnz = result->xadj[result->nvtxs];
	fflush(stdout);

	return result;	
}

RectMatrix* allPairsReorderData(int binary, RectMatrix
*featureVectors, int **inDegrees, float **maxOutWeights, float
**maxInWeights, int **descMaxOutWeightOrder, int
**ascOutDegreeOrder, int **outDegrees, int **permOrder)
{
	int nVectors = featureVectors->nVectors;
	int nDimensions = featureVectors->nDimensions;

	(*inDegrees) = (int*) calloc(nDimensions, sizeof(int) );

	if ( !binary )
	{
		(*maxOutWeights) = (float*) calloc(nVectors, sizeof(float));
		(*maxInWeights) = (float*) calloc(nDimensions, sizeof(float));
		(*descMaxOutWeightOrder) = (int*) malloc(sizeof(int) * nVectors);
	}
	else
	{
		*ascOutDegreeOrder = (int*) malloc(sizeof(int)*nVectors);
		*outDegrees = (int*) malloc(sizeof(int)*nVectors);
	}

	timer reorderTimer;
	cleartimer(reorderTimer);
	starttimer(reorderTimer);

	RectMatrix * reordered;
	if ( !binary )
	{
		reordered = reorderDimensionsAndVectors(featureVectors,
		*inDegrees, *maxInWeights, *maxOutWeights,
		*descMaxOutWeightOrder);
		*permOrder = *descMaxOutWeightOrder;
	}
	else
	{
		reordered = reorderBinaryData(featureVectors, *inDegrees,
					*outDegrees, *ascOutDegreeOrder);
		*permOrder = *ascOutDegreeOrder;
	}

	stoptimer(reorderTimer);
	printf("Time for re-ordering: %.3f\n",
	gettimer(reorderTimer));

	return reordered;
}

Matrix* AllPairsSimSearch::run()
{
	threshold = opts.threshold;
	nVectors = featureVectors->nVectors;
	nDimensions = featureVectors->nDimensions;

	if ( binary == 2 )
	{
//		printf("Binary, cosine similarity\n");
		binaryCosine = 1;
		binary = 1;
	}
	else if ( binary == 1 )
	{
//		printf("Binary, jaccard similarity\n");
		binaryCosine = 0;
	}

	timer totalTimer;
	cleartimer(totalTimer);
	starttimer(totalTimer);

	int *permOrder;
	reordered = allPairsReorderData(binary, featureVectors, 
 	 &inDegrees,&maxOutWeights, &maxInWeights, &descMaxOutWeightOrder,
	 &ascOutDegreeOrder, &outDegrees, &permOrder);

//	reordered = featureVectors;

	freeRectMatrix(featureVectors);

	Matrix *simMatrix;
	initSimPairs();
	if ( !binary )
		simMatrix = simPairs_reordered();
	else
		simMatrix = binarySimPairs();
	finalizeSimPairs();

	timer endTime;
	cleartimer(endTime);
	starttimer(endTime);

	freeRectMatrix(reordered);

	if ( !binary )
	{
		GKfree( (void**)&maxInWeights,
		(void**)&maxOutWeights, (void**)&inDegrees, LTERM);
	}
	else
	{
		GKfree( (void**)&inDegrees, (void**)&outDegrees, LTERM);
	}

	if ( simMatrix->nvtxs == 0 ||
	simMatrix->xadj[simMatrix->nvtxs] == 0 )
	{
		stoptimer(endTime);
		printf("end time: %.3f\n", gettimer(endTime));
		fflush(stdout);
		return simMatrix;
	}

	sortAdjLists(simMatrix->nvtxs, simMatrix->xadj,
			simMatrix->adjncy, simMatrix->adjwgt);
	Matrix* result = simMatrix;
	Matrix* revMatrix = getTranspose(simMatrix);
	result = add(simMatrix, revMatrix);
	freeMatrix(revMatrix);
	freeMatrix(simMatrix);

	Matrix*	result_perm = permuteRowsAndColumns(result, permOrder,permOrder);
	GKfree((void**)&permOrder, LTERM);
	
	freeMatrix(result);

	stoptimer(endTime);
//	printf("end time: %.3f\n", gettimer(endTime));

	stoptimer(totalTimer);
	printf("-----------------------------\n");
	printf("Time for total computation: %f\n",
			gettimer(totalTimer));
	printf("-----------------------------\n");

	fflush(stdout);
	return result_perm;

}

int compareIntInt(const void *a, const void *b)
{
	IntInt *aa = (IntInt*)a;
	IntInt *bb = (IntInt*)b;
	if ( aa->j < bb->j )
		return -1;
	else if ( aa->j == bb->j )
		return 0;
	else
		return 1;
}

RectMatrix* reorderBinaryData(const RectMatrix* featureVectors,
idxtype* inDegrees, idxtype* outDegrees, idxtype*
ascOutDegreeOrder)
{
	int nVectors, nDimensions;
	nVectors = featureVectors->nVectors;
	nDimensions = featureVectors->nDimensions;
	idxtype* xadj = featureVectors->xadj;
	idxtype* adjncy = featureVectors->adjncy;

	RectMatrix *reordered;
	
	for ( int i=0; i<xadj[nVectors]; i++ )
		inDegrees[adjncy[i]]++;
	
	IntInt* negInsAndOrder =
	(IntInt*)malloc(sizeof(IntInt)*nDimensions);

	for ( int i=0; i < nDimensions ; i++ )
	{
		negInsAndOrder[i].i = i;
		negInsAndOrder[i].j = 0 - inDegrees[i];
	}

	qsort((void *)negInsAndOrder, nDimensions, sizeof(IntInt), 
		compareIntInt);

	for ( int i=0; i < nDimensions; i++ )
	{
		inDegrees[i] = (0 - negInsAndOrder[i].j);
	}

	idxtype* descInDegreeOrder_rev = (idxtype*) 
								malloc(sizeof(idxtype)*nDimensions);
	for ( int i=0; i < nDimensions; i++ )
	{
		descInDegreeOrder_rev[negInsAndOrder[i].i] = i;
	}	
	GKfree( (void**)&negInsAndOrder, LTERM);

	for ( int i=0; i < xadj[nVectors]; i++ )
		adjncy[i] = descInDegreeOrder_rev[adjncy[i]];

	GKfree( (void**)&descInDegreeOrder_rev , LTERM);

	// Ok, by now, we have re-ordered the dimension ids and
	// changed adjncy as well.

	IntInt* outsAndOrder =
		(IntInt*)malloc(sizeof(IntInt)*nVectors);
	for ( int i=0; i<nVectors; i++ )
	{
		outsAndOrder[i].i = i;
		outsAndOrder[i].j = xadj[i+1]-xadj[i];
	}

	qsort((void*)outsAndOrder, nVectors, sizeof(IntInt),
	compareIntInt);

	for ( int i=0; i<nVectors; i++ )
	{
		outDegrees[i] = outsAndOrder[i].j;
		ascOutDegreeOrder[i] = outsAndOrder[i].i;
	}

	GKfree( (void**) &outsAndOrder, LTERM);	
		
	reordered = allocRectMatrix(nVectors, nDimensions,
	xadj[nVectors], 1);

/*	reordered = (RectMatrix*) malloc(sizeof(RectMatrix));
	reordered->xadj = (int*)malloc(sizeof(int)*(nVectors+1));
	reordered->adjncy = (int*)malloc(sizeof(int)*xadj[nVectors]);
*/
	reordered->xadj[0] = 0;
	for ( int i=0; i<nVectors; i++ )
	{
		int oldi = ascOutDegreeOrder[i];
		int j = reordered->xadj[i];
		for ( int oldj = xadj[oldi]; oldj < xadj[oldi+1]; oldj++, j++)
		{
			reordered->adjncy[j] = adjncy[oldj];
		}
		reordered->xadj[i+1] = j;
	}

	sortAdjLists(nVectors, reordered->xadj, reordered->adjncy);

	return reordered;
}

RectMatrix* reorderDimensionsAndVectors(const RectMatrix*
featureVectors, idxtype *inDegrees, wgttype *maxInWeights,
wgttype *maxOutWeights, idxtype* descMaxOutWeightOrder)
{
	int nVectors, nDimensions;
	nVectors = featureVectors->nVectors;
	nDimensions = featureVectors->nDimensions;
	idxtype* xadj = featureVectors->xadj;
	idxtype* adjncy = featureVectors->adjncy;
	wgttype* weights = featureVectors->adjwgt;

	RectMatrix *reordered;

	FloatInt* negInsAndOrder = (FloatInt*) 
							malloc(sizeof(FloatInt)*nDimensions);
	
	for ( int i=0; i<xadj[nVectors]; i++ )
		inDegrees[adjncy[i]]++;

	for ( int i=0; i<nDimensions; i++ )
	{
		negInsAndOrder[i].i = i;
		negInsAndOrder[i].f = 0 - inDegrees[i];
	}

	qsort((void *)negInsAndOrder, nDimensions, sizeof(FloatInt), 
			compareFloatInt);

	for ( int i=0; i < nDimensions; i++ )
	{
		inDegrees[i] = (idxtype) (0 - negInsAndOrder[i].f);
	}

	idxtype* descInDegreeOrder_rev = (idxtype*) 
								malloc(sizeof(idxtype)*nDimensions);
	for ( int i=0; i < nDimensions; i++ )
	{
		descInDegreeOrder_rev[negInsAndOrder[i].i] = i;
	}	
	GKfree( (void**)&negInsAndOrder, LTERM);

	for ( int i=0; i < xadj[nVectors]; i++ )
		adjncy[i] = descInDegreeOrder_rev[adjncy[i]];

	GKfree( (void**)&descInDegreeOrder_rev , LTERM);

	// Ok, by now, we have re-ordered the dimension ids and
	// changed adjncy as well.

	// Now, let's fill up maxOutWeights and maxInWeights.
	for ( int i=0; i<nVectors; i++)
	{
		maxOutWeights[i] = 0;
		for ( int j=xadj[i]; j<xadj[i+1]; j++)
		{
			if ( weights[j] > maxOutWeights[i] )
				maxOutWeights[i] = weights[j];
			if ( weights[j] > maxInWeights[adjncy[j]] )
				maxInWeights[adjncy[j]] = weights[j];
		}
	}

	FloatInt* descMaxOutAndOrder = (FloatInt*) 
			malloc(sizeof(FloatInt) * nVectors);
	
	for ( int i=0; i<nVectors; i++ )
	{
		descMaxOutAndOrder[i].f = 0 - maxOutWeights[i];
		descMaxOutAndOrder[i].i = i;
	}

	qsort( descMaxOutAndOrder, nVectors, sizeof(FloatInt), 
				compareFloatInt);
	
	for ( int i = 0; i<nVectors; i++ )
	{
		maxOutWeights[i] = 0 - descMaxOutAndOrder[i].f;
		descMaxOutWeightOrder[i] = descMaxOutAndOrder[i].i;
	}

	GKfree( (void**) &descMaxOutAndOrder, LTERM);

	reordered = allocRectMatrix(nVectors, nDimensions, xadj[nVectors]); 

	reordered->xadj[0] = 0;
	for ( int i=0; i<nVectors; i++ )
	{
		int oldi = descMaxOutWeightOrder[i];
		//int oldi = descMaxOutAndOrder[i].i;
		int j = reordered->xadj[i];
		for ( int oldj = xadj[oldi]; oldj < xadj[oldi+1]; oldj++, j++)
		{
			reordered->adjwgt[j] = weights[oldj];
			reordered->adjncy[j] = adjncy[oldj];
		}
		reordered->xadj[i+1] = j;
	}

	sortAdjLists(nVectors, reordered->xadj, reordered->adjncy,
		reordered->adjwgt);

	return reordered;
}





