#include <all.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>

static inline float convertCosineForHashing(float sim)
{
	return 1.0-acos(sim)/PI;
}

static inline float convertHashingForCosine(float x)
{
	return cos(PI*(1.0-x));
}

/* This function samples pairs that will be generated as
 * candidates. 
 * Rejection sampling: generate random samples, and keep only
 * those samples which go unpruned.
 */
float* allPairsSampleCandidates(RectMatrix *m, wgttype *maxInWeights, wgttype
*maxOutWeights, idxtype *inDegrees, wgttype cosThreshold, int*
numSamples)
{
	timer priorTimer;
	cleartimer(priorTimer);
	starttimer(priorTimer);

	int nVectors = m->nVectors;
	int nDimensions = m->nDimensions;
	idxtype *xadj = m->xadj;
	idxtype *adjncy = m->adjncy;
	wgttype *adjwgt = m->adjwgt;

	float* sims = (float*) malloc(sizeof(float)*(*numSamples));
	int numRejectedSamples = 0;
	int nrs1=0, nrs2=0, nrs3=0;
	idxtype* indexedFeatures;
	long maxRejectSamples = 10000000;
	int actualSamples=0;
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

		if ( cosThreshold > 0 ) 
		{
			/* v1 will be processed before v2. For this pair to be
			 * not pruned, v1 must be added to the indexes of at
			 * least one of the dimensions that both v1 and v2 share.
			 * Therefore, make up a new array containing only those
			 * features to whose indexes v1 is added to.
			 */
			wgttype b = 0;
			int j;
			for ( j = xadj[v1]; j < xadj[v1+1]; j++ )
			{
				wgttype min = maxInWeights[adjncy[j]];
				if ( min > maxOutWeights[v1] )
					min = maxOutWeights[v1];
				b = b + min * adjwgt[j];
				if ( b >= cosThreshold )
				{
					break;
				}
			}
			/* now, j to xadj_ends[i] are features to whose indexes v1 is
			 * going to be added. Put this in an array */
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
		}

		sims[i] = dotProduct(adjncy+xadj[v1], adjwgt+xadj[v1],
		xadj[v1+1]-xadj[v1], adjncy+xadj[v2], adjwgt+xadj[v2],
		xadj[v2+1]-xadj[v2]);

		sims[i] = convertCosineForHashing(sims[i]);
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
	stoptimer(priorTimer);
	printf("Time to sample similarities: %.3f\n", gettimer(priorTimer));
	fflush(stdout);
	#endif

	return sims;
}

int addUniformSamples(float **sims, int nsims, float lb,
float ub, float pct)
{
	int extra = (int) round(nsims * pct);
	if ( nsims == 0 )
	{
		extra = 100;
		#ifndef NDEBUG
		printf("Going to add %d uniform samples\n", extra);
		#endif
	}
	*sims = (float*) realloc(*sims, (nsims+extra)*sizeof(float));

	gsl_rng *r;
	r = gsl_rng_alloc(gsl_rng_taus2);

	for ( int i=nsims; i<nsims + extra; i++ )
	{
		(*sims)[i] = gsl_rng_uniform(r);
		(*sims)[i] = lb + ((*sims)[i])/(ub-lb);
	}
	gsl_rng_free(r);

	return nsims+extra;
}

float* allPairsSampleCandidates(Matrix *m, wgttype *maxInWeights, wgttype
*maxOutWeights, idxtype *inDegrees, wgttype cosThreshold, int*
numSamples)
{
	RectMatrix *rm = convertMatrixToRectMatrix(m);	
	return allPairsSampleCandidates(rm, maxInWeights, maxOutWeights,
	inDegrees, cosThreshold, numSamples);
}

void AllPairsRealCosine::doublePruneAddPairsForVertex(const idxtype vId, Matrix*
result, uint8_t* done, int *candidates, float *temp)
{
	result->xadj[vId+1] = result->xadj[vId];

	wgttype remScore = 0;
	for ( int i=xadj[vId]; i < xadj[vId+1]; i++ )
		remScore += adjwgt[i] * maxInWeights[adjncy[i]];
	wgttype minSize = (wgttype) (opts.cosThreshold / maxOutWeights[vId]);

	result->xadj[vId+1] = result->xadj[vId];
	if ( remScore < opts.cosThreshold || minSize > nDimensions)
	{
		return;
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

			if ( remScore >= opts.cosThreshold || temp[uId] > 0 )
			{
				wgttype add = invIndexWgts[j] * adjwgt[i];
				if ( temp[uId] <=0 && add > 0 )
					candidates[numCandidates++] = uId;
				temp[uId] += add;
			}
		}
		remScore -= adjwgt[i] * maxInWeights[k];
	}

	numTotalCandidates += numCandidates;

	for ( int i=0 ; i < numCandidates; i++)
	{
		idxtype uId = candidates[i];	
		int c = uId;

		wgttype upperbound = amin(xadj_ends[c]-xadj[c],
		xadj_ends[vId]-xadj[vId]) * maxOutWeights[vId] *
		maxOutWeights[c] + temp[c];
		if ( upperbound < opts.cosThreshold )
		{
			temp[c] = 0;
			continue;
		}
		numDotProducts++;
		(this->*processPair)(vId, uId, result);

		temp[c] = 0;
	}
}

void AllPairsRealCosine::pruneAddPairsForVertex(const idxtype vId, Matrix*
result, uint8_t* done, int *candidates)
{
	result->xadj[vId+1] = result->xadj[vId];

	wgttype remScore = 0;
	for ( int i=xadj[vId]; i < xadj[vId+1]; i++ )
		remScore += adjwgt[i] * maxInWeights[adjncy[i]];
	wgttype minSize = (wgttype) (opts.cosThreshold / maxOutWeights[vId]);

	result->xadj[vId+1] = result->xadj[vId];
	if ( remScore < opts.cosThreshold || minSize > nDimensions)
	{
		return;
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

			if ( remScore >= opts.cosThreshold && done[uId] <= 0 )
			{
				candidates[numCandidates++] = uId;
				done[uId] = 1;
			}
		}
		remScore -= adjwgt[i] * maxInWeights[k];
	}

	numTotalCandidates += numCandidates;

	for ( int i=0 ; i < numCandidates; i++)
	{
		idxtype uId = candidates[i];	
		numDotProducts++;
		(this->*processPair)(vId, uId, result);
		done[uId] = 0;
	}
		
}	

void AllPairsRealCosine::processPair_blshWrapper(const int vId, const
int uId, Matrix* result)
{
	float s = (blsh->*(blsh->processPair))(vId, uId);
	if ( s > 0 )
		addToMatrix(result, vId, uId, s);
}

void AllPairsRealCosine::addPairsForVertex(const idxtype vId, Matrix*
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
			numDotProducts++;

			(this->*processPair)(vId, uId, result);
		}
	}
	memset( (void*)done, 0, nPoints);
}

Matrix* AllPairsRealCosine::simPairs()
{
	xadj = reordered->xadj;
	adjncy = reordered->adjncy;
	adjwgt = reordered->adjwgt;
	
	invIndexAdjncy = idxmalloc(xadj[nPoints],"Inverted index adjncy");
	invIndexWgts =	(wgttype*) malloc(sizeof(wgttype)*xadj[nPoints]);
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
	for( int i = 0; i < nPoints; i++ )
	{
		wgttype b = 0;

		if ( opts.doublePruning )
			doublePruneAddPairsForVertex(i, result, done,
			candidates, temp);
		else if ( opts.additionalPruning )
			pruneAddPairsForVertex(i, result, done, candidates);
		else
			addPairsForVertex(i, result, done);

		int j;
		wgttype new_maxOutWeight = 0;
		for ( j = xadj[i]; j < xadj_ends[i]; j++ )
		{
			wgttype min = maxInWeights[adjncy[j]];
			if ( min > maxOutWeights[i] )
				min = maxOutWeights[i];
			b = b + min * adjwgt[j];
			if ( b >= opts.cosThreshold )
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
	free(invIndexWgts);
	free(invIndex_starts);
	free(invIndex_ends);
	free(xadj_ends);

	return result;
}

Matrix* AllPairsRealCosine::bayesLSH()
{
	nPoints = fvs->nVectors;
	nDimensions = fvs->nDimensions; 
	RectMatrix *featureVectors = fvs;

	inDegrees = (idxtype*) calloc(nDimensions, sizeof(idxtype));
	maxOutWeights = (wgttype*) calloc(nPoints, sizeof(wgttype));
	maxInWeights = (wgttype*) calloc(nDimensions, sizeof(wgttype));
	idxtype* descMaxOutWeightOrder = (idxtype*) 
				malloc( sizeof(idxtype) * nPoints);

	timer totalTimer;
	cleartimer(totalTimer);
	starttimer(totalTimer);

	reordered =
		reorderDimensionsAndVectors(featureVectors, inDegrees,
		maxInWeights, maxOutWeights, descMaxOutWeightOrder);

	freeRectMatrix(featureVectors);

	printf("Finished reordering data\n");
	fflush(stdout);

	gsl_set_error_handler_off();

	blsh = new BayesLSH(reordered, opts);
	if ( !opts.uniformPrior  )
	{
		float *sims = allPairsSampleCandidates(reordered,
		maxInWeights, maxOutWeights, inDegrees,
		opts.cosThreshold, &(opts.numSamples) );

		if ( opts.backoffPrior && opts.backoffPct > 0 )
		{
			#ifndef NDEBUG
			printf("Going to backoff prior by %f uniform fraction\n",
			opts.backoffPct);
			#endif
			opts.numSamples = addUniformSamples(&sims,
			opts.numSamples, 0.5, 1, opts.backoffPct );
		}

		if ( opts.pwLinearPrior )
		{
			int k=5;
			plprior = new PwLinear(sims, opts.numSamples, k, 0.5,
			1,
			convertCosineForHashing(opts.cosThreshold-opts.cosDelta));
			blsh->setupPwLinearPrior(plprior);
		}
		else
		{
			int k=5;
			kpprior = new KPwUnif(sims, opts.numSamples, k, 0.5,
			1,
			convertCosineForHashing(opts.cosThreshold-opts.cosDelta));
			kpprior->print();
			blsh->setupKPwUnifPrior(kpprior);
		}
	}
	else
		blsh->setupUnifPrior();
	
	blsh->setup();
	processPair = &AllPairsRealCosine::processPair_blshWrapper;

	Matrix* simMatrix = simPairs();

	#ifndef NDEBUG
	blsh->printPruneStatistics();
	#endif

	delete blsh;

	Matrix *result = finalizeResultMatrix(simMatrix);
	if ( result->xadj[result->nvtxs] > 0 )
	{
		Matrix* result_perm = permuteRowsAndColumns(result,
		descMaxOutWeightOrder, descMaxOutWeightOrder);
		freeMatrix(result);

		GKfree((void**)&descMaxOutWeightOrder, LTERM);
		
		result = result_perm;
	}

	stoptimer(totalTimer);
	printf("-----------------------------\n");
	printf("Time for total computation: %f\n",
			gettimer(totalTimer));
	printf("-----------------------------\n");
	fflush(stdout);

	return result;

}

void setupPwUnifPrior(RectMatrix *fvs, float
*maxOutWeights, float *maxInWeights, int *inDegrees, BssOptions
opts, PwUnif *puPrior)
{
	timer priorTimer;
	cleartimer(priorTimer);
	starttimer(priorTimer);

	float *sims = allPairsSampleCandidates(fvs, maxInWeights,
	maxOutWeights, inDegrees, opts.cosThreshold,
	&(opts.numSamples) );

	puPrior->lb = 0.5;
	puPrior->ub = 1.0;
	puPrior->unifPct = opts.backoffPct;

	if ( opts.backoffPrior && opts.backoffPct > 0 )
	{
		#ifndef NDEBUG
		printf("Going to backoff prior by %f uniform fraction\n",
		puPrior->unifPct);
		#endif
		opts.numSamples = addUniformSamples(&sims,
		opts.numSamples, puPrior->lb, puPrior->ub, puPrior->unifPct );
	}

	int num_os = opts.num_os;
	float *os = getOrderStatistics(sims, opts.numSamples, num_os);
	float knee, pct;
	fitToPiecewiseUniform(os, num_os+1, 0.5, 1, &knee, &pct);

	GKfree((void**)&sims, (void**)&os, LTERM);

	puPrior->p1 = pct/(knee-0.5);
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


