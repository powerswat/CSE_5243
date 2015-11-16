#include <all.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>

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

float* allPairsBinarySampleCandidates(RectMatrix* reordered, int* inDegrees,
int* outDegrees, float threshold, int* numSamples, int
binaryCosine)
{
	int nVectors = reordered->nVectors;
	int nDimensions = reordered->nDimensions;
	idxtype *xadj = reordered->xadj;
	idxtype *adjncy = reordered->adjncy;

	float* sims = (float*) malloc(sizeof(float)*(*numSamples));
	int numRejectedSamples = 0;
	int nrs1=0, nrs2=0, nrs3=0;
	idxtype* indexedFeatures;
	int actualSamples = 0;
	long maxRejectSamples = 10000000;
	for ( int i = 0; i<(*numSamples); )
	{
		if ( numRejectedSamples > maxRejectSamples )
			break;

		int v1 = (rand())%(nVectors);
		int v2 = (rand())%(nVectors);

		if ( v1 == v2 )
			continue; // reject sample

		if ( v1 > v2 )
		{
			int tmp;
			SWAP(v1, v2, tmp);
		}

		int d1 = xadj[v1+1]-xadj[v1];
		int d2 = xadj[v2+1]-xadj[v2];

		int doIntersect = doTheyIntersect(adjncy+xadj[v1], d1,
			adjncy+xadj[v2], d2);
		if ( doIntersect <= 0 )
		{
			numRejectedSamples++;
			nrs1++;
			continue; // reject sample.
		}

		float b=0;
		int j;
		for ( j=xadj[v1]; j<xadj[v1+1]; j++ )
		{
			b++;
			if ( b/d1 >= threshold )
				break;
		}

		if ( j == xadj[v1+1] )
		{
			numRejectedSamples++;
			nrs2++;
			continue;
		}

		int numIndexedFeatures = xadj[v1+1] - j;
		indexedFeatures = (idxtype*) 
					malloc(sizeof(idxtype)*numIndexedFeatures);
		for ( int k = 0; k<numIndexedFeatures; k++ )
		{
			indexedFeatures[k] = adjncy[j+k];
		}
	
		doIntersect = doTheyIntersect(indexedFeatures,
			numIndexedFeatures, adjncy+xadj[v2], d2);
		free(indexedFeatures);
		if ( doIntersect <= 0 )
		{
			nrs3++;
			numRejectedSamples++;
			continue;
		}

		// d2 always >= d1
		if ( (binaryCosine && d1/sqrt(1.0*d1*d2) < threshold) ||
			(!binaryCosine && d1*1.0/d2 < threshold) )
		{
			numRejectedSamples++;
			continue;
		}

		sims[i] = sizeOfSetIntersect(adjncy+xadj[v1],
		xadj[v1+1]-xadj[v1], adjncy+xadj[v2],
		xadj[v2+1]-xadj[v2]);
		sims[i] = sims[i]/(d1+d2-sims[i]);
		i++;
		actualSamples++;
	}

	if ( actualSamples < (*numSamples) )
	{
		sims = (float*)realloc(sims,actualSamples*sizeof(float)); 
		*numSamples = actualSamples;
	}

	#ifndef NDEBUG
	printf("Sampled %d similarities\n", actualSamples);
	printf("Total no. of rejected samples: %d\n", numRejectedSamples);
	printf("nrs1: %d, nrs2: %d, nrs3: %d\n", nrs1, nrs2, nrs3);
	fflush(stdout);
	#endif

	return sims;
}

void AllPairsBinary::doublePruneAddPairsForVertex(int vId, Matrix*
result, uint8_t* done, int* candidates, float* temp)
{
	wgttype remScore = outDegrees[vId];
	wgttype minSize;
	if ( opts.binaryCosine)
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

	numTotalCandidates += numCandidates;
	for ( int i=0 ; i < numCandidates; i++)
	{
		int c = candidates[i];

		int maxIntersect = (int) round(temp[c]) + xadj_ends[c]-xadj[c];
		double upperBound;
		double binaryCosineDen = 0;

		if ( opts.binaryCosine )
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

		numDotProducts++;
		(this->*processPair)(vId, c, result);
		temp[c] = 0;
	}

}

void AllPairsBinary::pruneAddPairsForVertex(int vId, Matrix*
result, uint8_t* done, int* candidates)
{
	wgttype remScore = outDegrees[vId];
	wgttype minSize = (wgttype) threshold * outDegrees[vId];

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

			if ( remScore >= minSize && done[uId] <= 0 )
			{
				candidates[numCandidates++] = uId;
				done[uId] = 1; 
			}
		}
		remScore--; 
	}
	
	numTotalCandidates += numCandidates;

	for ( int i=0 ; i < numCandidates; i++)
	{
		numDotProducts++;
		int c = candidates[i];
		(this->*processPair)(vId, c, result);
		done[c] = 0;
	}

}

void AllPairsBinary::processPair_blshWrapper(const int vId, const
int uId, Matrix *result)
{
	float s = (blsh->*(blsh->processPair))(vId, uId);
	if ( s > 0 )
		addToMatrix(result, vId, uId, s);
}

void AllPairsBinary::addPairsForVertex(const idxtype vId, Matrix*
result, uint8_t* done)
{
	result->xadj[vId+1] = result->xadj[vId];

	for ( int i=xadj[vId+1]-1; i >= xadj[vId]; i-- )
	{
		int k = adjncy[i];

		for ( int j = invIndex_starts[k]; j < invIndex_ends[k];	j++ )
		{
			int uId = invIndexAdjncy[j];

			if ( done[uId] )
				continue;
			done[uId] = 1;

			numTotalCandidates++;

			(this->*processPair)(vId, uId, result);
		}
	}
	memset( (void*)done, 0, nPoints);

}

Matrix* AllPairsBinary::simPairs()
{
	xadj = reordered->xadj;
	adjncy = reordered->adjncy;
	
	invIndexAdjncy = idxmalloc(xadj[nPoints],"Inverted index adjncy");
	invIndex_starts = idxmalloc(nDimensions, "Inverted index xadj starts");
	invIndex_ends = idxmalloc(nDimensions, "Inverted index xadj ends");

	invIndex_starts[0] = 0;
	for ( int i = 1; i < nDimensions; i++ )
		invIndex_starts[i] = invIndex_starts[i-1] +	inDegrees[i-1];
	for ( int i = 0; i < nDimensions; i++ )
		invIndex_ends[i] = invIndex_starts[i];

	xadj_ends = idxmalloc(nPoints, "xadj ends");
	for ( int i=0; i < nPoints; i++ )
		xadj_ends[i] = xadj[i+1];
	
	int init_ret_size = nPoints * 50;
	Matrix *result = allocMatrix(nPoints, init_ret_size, 0, 0, 0);
	result->sizeIncrement = nPoints*20;
	result->xadj[0] = 0;

	uint8_t* done = (uint8_t*)calloc(nPoints, sizeof(uint8_t));
	int* candidates;
	if ( opts.additionalPruning || opts.doublePruning )
	{
		candidates = (int*)malloc(sizeof(int)*nPoints);
	}
	float *temp;
	if ( opts.doublePruning )
	{
		temp = (float*)calloc(nPoints, sizeof(float));
	}

	cleartimer(simTimer);
	starttimer(simTimer);
	numDotProducts = 0;
	numTotalCandidates = 0;

	printf("Progress Indicator (no. of records):");
	for( int i = 0; i < nVectors; i++ )
	{
		float b = 0;
		if ( opts.doublePruning )
			doublePruneAddPairsForVertex(i, result, done,
			candidates, temp);
		else if ( opts.additionalPruning )
			pruneAddPairsForVertex(i, result, done, candidates);
		else
			addPairsForVertex(i, result, done);
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
		}

		if ( i % 20000 == 0 )
		{
			printf("%d..",i);
			fflush(stdout);
		}
	}
	printf("\n");
	stoptimer(simTimer);
//	printf("Sim. time: %.3f\n", gettimer(simTimer));
	printf("Num. candidates supplied to BayesLSH: %ld\n", numDotProducts);
	printf("Num. similar pairs: %d\n",result->xadj[result->nvtxs]);
	result->nnz = result->xadj[result->nvtxs];

	if ( opts.doublePruning )
		free(temp);
	if ( opts.doublePruning || opts.additionalPruning )
		free(candidates);
	free(done);
	free(invIndexAdjncy);
	free(invIndex_starts);
	free(invIndex_ends);
	free(xadj_ends);

	return result;

}

void setupPwUnifPrior(float *sims, BssOptions opts, PwUnif *puPrior)
{
	timer priorTimer;
	cleartimer(priorTimer);
	starttimer(priorTimer);

	puPrior->lb = 0;
	if ( opts.useMinwiseBitSketches || opts.binaryCosine )
		puPrior->lb = 0.5;
	puPrior->ub = 1.0;
	puPrior->unifPct = opts.backoffPct;

	float knee, pct;

	int num_os = opts.num_os;
	float *os = getOrderStatistics(sims, opts.numSamples, num_os);
	fitToPiecewiseUniform(os, num_os+1, puPrior->lb,
	puPrior->ub, &knee, &pct);
	GKfree( (void**)&os, LTERM);

	puPrior->p1 = pct/(knee-puPrior->lb);
	puPrior->p2 = (1.0-pct)/(1-knee);
	puPrior->changePt = knee;
	puPrior->logChangePt = log(puPrior->changePt);
	puPrior->log1minusChangePt = log(1.0-puPrior->changePt);
	puPrior->logp1 = log(puPrior->p1);
	puPrior->logp2 = log(puPrior->p2);

	stoptimer(priorTimer);
	#ifndef NDEBUG
	printf("Time to sample prior: %.3f\n", gettimer(priorTimer));
	printf("Piece-wise uniform change pt: %f, Percentile: %.3f\n", 
			knee, pct);
	fflush(stdout);
	#endif

}

float* allPairsBinarySampleCandidates_wrapper(RectMatrix *fvs, int* inDegrees, int*
outDegrees, BssOptions *opts)
{
	timer priorTimer;
	cleartimer(priorTimer);
	starttimer(priorTimer);

	#ifndef NDEBUG
	printf("Going to sample similarities\n");
	fflush(stdout);
	#endif

	float *sims = allPairsBinarySampleCandidates(fvs, inDegrees,
	outDegrees, opts->cosThreshold, &(opts->numSamples),
	opts->binaryCosine);

	if ( opts->backoffPrior && opts->backoffPct > 0 )
	{
		float ub = 1.0;
		float lb = 0;

		#ifndef NDEBUG
		printf("Going to backoff prior by %f uniform fraction\n",
		opts->backoffPct);
		#endif
		opts->numSamples = addUniformSamples(&sims,
		opts->numSamples, lb, ub, opts->backoffPct );
	}

	if ( opts->useMinwiseBitSketches )
	{
		for ( int i=0; i<opts->numSamples; i++ )
			sims[i] = convertJaccardForHashing(sims[i]);
	}
	else if ( opts->binaryCosine )
	{
		for ( int i=0; i<opts->numSamples; i++ )
			sims[i] = convertCosineForHashing(sims[i]);
	}

	stoptimer(priorTimer);
	#ifndef NDEBUG
	printf("Time to sample prior: %.3f\n", gettimer(priorTimer));
	fflush(stdout);
	#endif

	return sims;
}

void setupBetaPrior(RectMatrix *fvs, int* inDegrees, int*
outDegrees, BssOptions opts, Beta *bPrior)
{
	float *sims = allPairsBinarySampleCandidates_wrapper(fvs, inDegrees,
	outDegrees, &opts);

	fitToBeta(sims, opts.numSamples, &(bPrior->alpha),
	&(bPrior->beta), 0, 1);
	#ifndef NDEBUG
	printf("beta.alpha: %f, beta.beta: %f\n", bPrior->alpha,
	bPrior->beta);
	#endif
}

Matrix* AllPairsBinary::bayesLSH()
{
	nPoints = fvs->nVectors;
	nVectors = fvs->nVectors;
	nDimensions = fvs->nDimensions; 
	RectMatrix *featureVectors = fvs;
	threshold = opts.cosThreshold;

	inDegrees = (idxtype*) calloc(nDimensions, sizeof(idxtype));
	ascOutDegreeOrder = (idxtype*) malloc(sizeof(idxtype)*nVectors);
	outDegrees = (idxtype*) malloc(sizeof(idxtype)*nVectors);

	timer totalTimer;
	cleartimer(totalTimer);
	starttimer(totalTimer);

	idxtype* permOrder;

	reordered = reorderBinaryData(featureVectors, inDegrees,
	outDegrees, ascOutDegreeOrder);
	permOrder = ascOutDegreeOrder;

	freeRectMatrix(featureVectors);

	printf("Finished reordering data\n");
	fflush(stdout);

	gsl_set_error_handler_off();

	blsh = new BayesLSH(reordered, opts);
	if ( !opts.uniformPrior )
	{
		if ( opts.pwLinearPrior )
		{
			float *sims = allPairsBinarySampleCandidates_wrapper(reordered,
			inDegrees, outDegrees, &opts);

			int k=5;
			float lb = opts.binaryCosine ? 0.5 : 0;
			plprior = new PwLinear(sims, opts.numSamples, k, lb,
			1,
			convertCosineForHashing(opts.cosThreshold-opts.cosDelta));
			blsh->setupPwLinearPrior(plprior);
		}
		else if ( !opts.betaPrior )
		{
			float *sims = allPairsBinarySampleCandidates_wrapper(reordered,
			inDegrees, outDegrees, &opts);

		/*	setupPwUnifPrior(sims, opts, &puPrior);
			blsh->setupPwUnifPrior(&puPrior);
		*/
			int k=5;
			float lb= opts.binaryCosine ? 0.5 : 0 ;
			kprior = new KPwUnif(sims, opts.numSamples, k, lb,
			1,convertCosineForHashing(opts.cosThreshold-opts.cosDelta));
			kprior->print();
			blsh->setupKPwUnifPrior(kprior);
		
		}
		else
		{
			if ( opts.binaryCosine )
			{
				printf("Sorry, beta prior cannot be used for binary cosine\n");
				exit(0);
			}
			Beta *bPrior = new Beta();
			setupBetaPrior(reordered, inDegrees, outDegrees,
			opts, bPrior);
			blsh->setupBetaPrior(bPrior);
		}
	}
	else
	{
		blsh->setupUnifPrior();
	}
	blsh->setup();

	processPair = &AllPairsBinary::processPair_blshWrapper;
	Matrix* simMatrix = simPairs();

	#ifndef NDEBUG
	blsh->printPruneStatistics();
	#endif

	delete blsh;

	Matrix *result = finalizeResultMatrix(simMatrix);
	if ( result->xadj[result->nvtxs] > 0 )
	{
		Matrix* result_perm = permuteRowsAndColumns(result,
		permOrder, permOrder);
		freeMatrix(result);

		GKfree( (void**)&permOrder, LTERM);

		result = result_perm;
	}

	stoptimer(totalTimer);
	printf("-----------------------------\n");
	printf("Time for total computation: %f\n",
			gettimer(totalTimer));
	printf("-----------------------------\n");

	return result;
}



