#include <all.h>
#include <assert.h>

static inline uint32_t popcount ( uint32_t v )
{
	v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
	uint32_t c = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count

	return c;
}

static inline float convertJaccardForHashing(float sim)
{
	return (1.0+sim)/2;
}

static inline float convertHashingForJaccard(float sim)
{
	return (2*sim-1);
}

static inline float convertHashingForCosine(float x)
{
	return cos(PI*(1.0-x));
}

static inline float convertCosineForHashing(float sim)
{
	return 1.0-acos(sim)/PI;
}

void LSH::processPair_jac(const int vId, const int uId,
Matrix *result, RectMatrix *fvs)
{
	if ( fvs->xadj[vId+1] == fvs->xadj[vId] || fvs->xadj[uId+1]
	== fvs->xadj[uId] )
		return;

	int* uw = (mis)->minSketches + ((long) uId)*((mis)->minSketchSize);
	int* vw = (mis)->minSketches + ((long) vId)*((mis)->minSketchSize);
	
	int matches = 0;
	int nTrials;
	for (nTrials=0; nTrials < (mis)->minSketchSize; nTrials ++ )
	{
		if (*(vw++) == *(uw++)) 
			matches++;  
	}

	float sim = 1.0*matches / nTrials;

	if ( sim >= (LSH::marginPct)*(cosThreshold) )
		addToMatrix(result, vId, uId, sim);
}

void LSH::processPair_blshWrapper(const int vId, const int
uId, Matrix *result, RectMatrix *fvs)
{
	float s = (blsh->*(blsh->processPair))(vId, uId);
	if ( s > 0 )
		addToMatrix(result, vId, uId, s);
}

void LSH::processPair_jac_exact(const int vId, const int uId,
Matrix *result, RectMatrix *fvs)
{
	float sim = jaccard(fvs, vId, uId);
	if ( sim > cosThreshold )
		addToMatrix(result, vId, uId, sim);
}

void LSH::processPair_cos_exact(const int vId, const int uId,
Matrix *result, RectMatrix *fvs)
{
	float sim;
	if ( opts.binaryCosine )
		sim = cosine_binary(fvs, vId, uId);
	else
		sim = dotProduct(fvs, vId, uId);
	if ( sim > cosThreshold )
		addToMatrix(result, vId, uId, sim);
}

void LSH::processPair_cos(const int vId, const int uId,
Matrix *result, RectMatrix *fvs)
{
	if ( fvs->xadj[vId+1] == fvs->xadj[vId] || fvs->xadj[uId+1]
	== fvs->xadj[uId] )
		return;

	uint32_t *uw = (uint32_t*) ( (uint8_t*)(cs)->minSketches + 
					((long) uId)*((cs)->minSketchSize/8) );
	uint32_t *vw = (uint32_t*) ( (uint8_t*)(cs)->minSketches + 
					((long) vId)*((cs)->minSketchSize/8) );

	int matches = 0;
	int nTrials;
	for (nTrials=0; nTrials < (cs)->minSketchSize;
	nTrials += 32 )
	{
		matches += 32 - (int) popcount(*(vw++) ^ *(uw++));
	}

	float sim = 1.0*matches / nTrials;
	sim = convertHashingForCosine(sim);
	
	if ( sim >= (LSH::marginPct)*(cosThreshold) )
		addToMatrix(result, vId, uId, sim);
}

void LSH::freeLSHIndexes()
{
	free(invArrays_base);
	free(invArrays);
	free(invArraysXadj_base);
	free(invArraysXadjs);
}

float* sampleLSHCandidates(RectMatrix *m, int
*invArrayXadj, int *invArrayAdjncy, int hashCardinality, float
(*similarity)(RectMatrix*, int, int, int), int numSamples)
{
	long numCandidates = 0;
	for ( int i=0; i<hashCardinality; i++ )
	{
		long numEquals = invArrayXadj[i+1] - invArrayXadj[i];
		numCandidates += numEquals*((numEquals-1)/2);
	}

	double samplePct = numSamples*1.0/numCandidates;
	//printf("samplePct: %f\n", samplePct);
	//fflush(stdout);
	if ( samplePct > 1 )
	{
		samplePct = 1.0;
	}

	float* sims = (float*) malloc(sizeof(float)*numSamples);
	int sims_counter = 0;
	int done = 0;

	// when samplePct is very small, the naive approach can be
	// very slow. Hence, this one, refer "Efficient generation of
	// large random networks" by Bataglej and Brandes
	for ( int i=0, xadjIndex=0, vIdIndex=0, uIdIndex=1; 
				i<numSamples; i++ )
	{
		double r = (rand())*1.0/RAND_MAX;
		long numToSkip = 1 + (long)
		floor(log(1-r)/log(1-samplePct));

		int curr_len = invArrayXadj[xadjIndex+1] -
		invArrayXadj[xadjIndex];

		for ( long j=0; j<numToSkip; )
		{
			if ( j + (curr_len - uIdIndex - 1) > numToSkip )
			{
				uIdIndex += numToSkip - j;
				break;
			}

			j += curr_len - uIdIndex;
			if ( vIdIndex+2 >= curr_len )
			{
				if ( xadjIndex + 1 >= hashCardinality )
				{
					done = 1;
					break;
				}

				do
				{
					xadjIndex++;
					curr_len =
					invArrayXadj[xadjIndex+1]-invArrayXadj[xadjIndex];
				} while ( curr_len < 2 && xadjIndex <
				hashCardinality );

				if ( xadjIndex >= hashCardinality )
				{
					done = 1; break;
				}

				vIdIndex=0;
			}
			else
			{
				vIdIndex++;
			}

			uIdIndex = vIdIndex+1;
		}

		if ( done )
			break;

		assert(xadjIndex < hashCardinality);
		assert(uIdIndex < curr_len);
		assert(vIdIndex < curr_len);
		assert(uIdIndex > vIdIndex);

		int vId =
		invArrayAdjncy[invArrayXadj[xadjIndex]+vIdIndex];
		int uId = 
		invArrayAdjncy[invArrayXadj[xadjIndex]+uIdIndex];
		if ( m->xadj[vId+1]-m->xadj[vId] < 1 ||
		m->xadj[uId+1]-m->xadj[uId] < 1 )
			continue;
		sims[sims_counter++] = similarity(m, vId, uId, 0);
	}

	#ifndef NDEBUG
	printf("Got %d samples\n", sims_counter);
	#endif

	if ( sims_counter < numSamples )
	{
		sims = (float *)realloc(sims, sims_counter*sizeof(float));
	}

	numSamples = sims_counter;

	return sims;
}

void setupBetaPrior_lsh(RectMatrix *fvs, int
*invArrayXadj, int *invArrayAdjncy, int hashCardinality, float
(*similarity)(RectMatrix*, int, int, int), BssOptions opts, Beta *bPrior)
{
	timer sampleTimer;

	cleartimer(sampleTimer);
	starttimer(sampleTimer);

	float *sims = sampleLSHCandidates(fvs, invArrayXadj,
	invArrayAdjncy, hashCardinality, similarity, opts.numSamples);

	if ( opts.backoffPrior && opts.backoffPct > 0 )
	{
		#ifndef NDEBUG
		printf("Going to backoff prior by %f uniform fraction\n",
		opts.backoffPct);
		#endif
		opts.numSamples = addUniformSamples(&sims,
		opts.numSamples, 0, 1, opts.backoffPct );
	}

	fitToBeta(sims, opts.numSamples, &(bPrior->alpha),
	&(bPrior->beta), 0, 1);

	stoptimer(sampleTimer);

	#ifndef NDEBUG
	printf("beta.alpha: %f, beta.beta: %f\n", bPrior->alpha,
	bPrior->beta);
	printf("Time for sampling: %f\n", gettimer(sampleTimer));
	fflush(stdout);
	#endif
}


Matrix* LSH::processIndexes()
{
	int init_ret_size = nPoints * 50;
	Matrix *result = allocMatrix(nPoints, init_ret_size, 0, 0, 0);
	result->sizeIncrement = nPoints*20;
	result->xadj[0] = 0;

	for ( int block=0; block<numBlocks; block++ )
	{
		sortAdjLists(numValuesPerBlock, invArraysXadjs[block],
		invArrays[block]);
	}

	uint8_t* done = (uint8_t*)calloc(nPoints, sizeof(uint8_t));
	int* currCandidates = (int*)malloc(nPoints*sizeof(int));

	uint8_t* csk;
	int* jsk;
	if ( !opts.jaccard || opts.useMinwiseBitSketches )
		csk = (uint8_t*) cs->minSketches;
	else
		jsk = mis->minSketches;
	long numCandidates = 0;
	printf("Progress Indicator (no. of records):");
	for ( long vId=0; vId<nVectors; vId++ )
	{
		result->xadj[vId+1] = result->xadj[vId];
		int numCurrCandidates = 0;
		for ( int block=0; block<numBlocks; block++ )
		{
			int *xadj = invArraysXadjs[block];
			int *adjncy = invArrays[block];
			// generate candidate neighbors from this block.
			int hv;
			if ( !opts.jaccard || opts.useMinwiseBitSketches )
			{
				hv = *(csk + vId*(cs->minSketchSize/8) +
				block*sizeOfBlock/8);
			}
			else
			{
				hv = *(jsk + vId*mis->minSketchSize +
				block*sizeOfBlock);
			}
			for ( int j=xadj[hv];j<xadj[hv+1]; j++ )
			{
				int uId = adjncy[j];
				if ( uId >= vId )
					break;
				if ( done[uId] > 0 )
					continue;

				currCandidates[numCurrCandidates++] = uId;
				done[uId] = 1;

				numCandidates++;

				(this->*processPair)(vId, uId, result, fvs);
			}
		}

		for ( int i=0; i<numCurrCandidates; i++ )
			done[currCandidates[i]] = 0;

		if ( vId % 20000 == 0 )
		{
			printf("%ld..", vId); fflush(stdout);
		}
	}
	printf("\n");

	free(done);
	free(currCandidates);

	printf("Num. candidates supplied to BayesLSH: %ld\n",
	numCandidates);
	return result;
}

void allocLSHIndexes(const int cardinalityOfHash, const int
sizeOfBlock, const int numBlocks, const int nVectors, int** invArrays_base, int***
invArrays, int** invArraysXadj_base, int*** invArraysXadjs)
{
	int	numValuesPerBlock = (int) pow(cardinalityOfHash, sizeOfBlock);
	int numInvArrays = numBlocks;
	long invArrays_base_size = numInvArrays*((long)nVectors);
	(*invArrays_base) =
	(int*)malloc(sizeof(int)*invArrays_base_size);
	(*invArrays) = (int**)malloc(sizeof(int*)*numInvArrays);
	for ( int i=0; i<numInvArrays; i++ )
		(*invArrays)[i] = (*invArrays_base) + i*nVectors;
	
	long iaXadj_base_size = sizeof(int) * numInvArrays * (long)
	(numValuesPerBlock+1);
	(*invArraysXadj_base) = (int*)malloc(iaXadj_base_size);
	(*invArraysXadjs) = (int**)malloc(sizeof(int*)*numInvArrays);
	for ( int i=0; i<numInvArrays; i++ )
	{
		(*invArraysXadjs)[i] = (*invArraysXadj_base) +
		i*(numValuesPerBlock+1);
	}
}

void LSH::buildLSHIndexes()
{
	assert ( ((!opts.jaccard || opts.useMinwiseBitSketches) &&
	sizeOfBlock == 8) || (opts.jaccard
	&& sizeOfBlock==1) ); 
/*	if ( ! ( ((!opts.jaccard || opts.useMinwiseBitSketches) &&
	sizeOfBlock == 8) || (opts.jaccard
	&& sizeOfBlock==1) ) )
	{
		// for now, only this will work
		printf("Can't deal with sizeOfBlock %d\n", sizeOfBlock);
		exit(0);
	}
*/
	int cardinalityOfHash;
	if ( opts.jaccard && !opts.useMinwiseBitSketches )
		cardinalityOfHash = nDimensions;
	else
		cardinalityOfHash = 2;
	numValuesPerBlock = (int) pow(cardinalityOfHash, sizeOfBlock);

	allocLSHIndexes(cardinalityOfHash, sizeOfBlock, numBlocks,
	nVectors, &invArrays_base, &invArrays, &invArraysXadj_base,
	&invArraysXadjs);

	int* countsOfBlockValues = (int*)calloc(numValuesPerBlock, sizeof(int));

	uint8_t* csk;
	int* jsk;
	if ( !opts.jaccard || opts.useMinwiseBitSketches )
		csk = (uint8_t*) cs->minSketches;
	else
		jsk = mis->minSketches;

	// initialize invArraysXadjs
	for ( int block=0; block<numBlocks; block++ )
	{
		// first go through all values in this block, and update
		// countsOfBlockValues
		for ( long i=0; i<nVectors; i++ )
		{
			int hv;
			if ( !opts.jaccard || opts.useMinwiseBitSketches )
			{
				hv = *(csk + i*(cs->minSketchSize/8) +
				block*(sizeOfBlock/8));
			}
			else
			{
				hv = *(jsk + i*mis->minSketchSize +
				block*sizeOfBlock);
			}
			countsOfBlockValues[hv]++;
		}
		
		invArraysXadjs[block][0] = 0;
		for ( int value=1; value<=numValuesPerBlock; value++ )
		{
			invArraysXadjs[block][value] =
			invArraysXadjs[block][value-1] +
			countsOfBlockValues[value-1];
		}

		for ( long i=0; i<nVectors; i++ )
		{
			int hv;
			if ( !opts.jaccard || opts.useMinwiseBitSketches )
			{
				hv = *(csk + i*(cs->minSketchSize/8) +
				block*(sizeOfBlock/8));
			}
			else
			{
				hv = *(jsk + i*mis->minSketchSize +
				block*sizeOfBlock);
			}
			
			int j = invArraysXadjs[block][hv+1] - countsOfBlockValues[hv];
			invArrays[block][j] = i;
			countsOfBlockValues[hv]--;
		}

	}
		
	free(countsOfBlockValues);
	#ifndef NDEBUG
	printf("Done building inverted indexes\n");
	fflush(stdout);
	#endif
}

int LSH::getNumBlocksForRecall(float recall, float
threshold, int sizeOfBlock)
{
	return (int) ceil(log(1-recall)/log(1-pow(threshold, sizeOfBlock)));

}

Matrix* LSH::runRecall()
{

	timer totalTimer, skTimer;
	cleartimer(totalTimer);
	starttimer(totalTimer);

	float threshold = cosThreshold;
	if ( !opts.jaccard )
	{
		threshold = convertCosineForHashing(cosThreshold);
		sizeOfBlock = 8;
	}
	else
		sizeOfBlock = 1;

	if ( opts.useMinwiseBitSketches )
	{
		opts.useMinwiseBitSketches = 0;
	}

	numBlocks = getNumBlocksForRecall(recall, threshold, sizeOfBlock);
	#ifndef NDEBUG
	printf("numBlocks %d, for recall %f, threshold %f, sizeofBlock %d\n",
	numBlocks, recall, threshold,sizeOfBlock);
	fflush(stdout);
	#endif

	nHashes = sizeOfBlock * numBlocks;
	if ( !opts.lsh_useExact )
	{
		if ( nHashes > maxHashes )
		{
			nHashes = sizeOfBlock * (int)
			round(maxHashes*1.0/sizeOfBlock);
		}
		if ( nHashes < maxHashes )
		{
			nHashes = maxHashes;
		}
	}

	cleartimer(skTimer);
	starttimer(skTimer);
	if ( !opts.jaccard )
	{
		RandomIntGaussians *rig = new
		RandomIntGaussians(nDimensions, nHashes,
		opts.availableRandomGaussiansMemoryInMB,
		opts.hash_initRandomSeed);
		cs = new CosineSketches(nPoints, nHashes, 8, 8);
		cs->buildMinSketchesInBatches(fvs, rig);
		delete rig;
	}
	else
	{
		mp = new MinwisePerms(nHashes,opts.hash_initRandomSeed );
		mis = new MinwiseIntSketches(nVectors, nHashes, 1, 1);
		mis->buildMinSketches(fvs, mp);
		delete mp;
	}
	stoptimer(skTimer);
	#ifndef NDEBUG
	printf("Time to build %d*%d minSketches: %.3f\n",
	nPoints, nHashes, gettimer(skTimer));
	fflush(stdout);
	#endif

	buildLSHIndexes();
	Matrix *simMatrix;

	if ( opts.jaccard )
	{
		if ( opts.lsh_useExact )
			processPair = &LSH::processPair_jac_exact;
		else
			processPair = &LSH::processPair_jac;
	}
	else
	{
		if ( opts.lsh_useExact )
			processPair = &LSH::processPair_cos_exact;
		else
			processPair = &LSH::processPair_cos;
	}

	simMatrix = processIndexes();
		
	printf("Num. similar pairs: %d\n", simMatrix->xadj[simMatrix->nvtxs]);
	freeLSHIndexes();
	Matrix *result = finalizeResultMatrix(simMatrix);

	if ( !opts.jaccard )
	{
		delete cs;
	}
	else
	{
		delete mis;
	}

	stoptimer(totalTimer);
	printf("-----------------------------\n");
	printf("Time for total computation: %f\n",
			gettimer(totalTimer));
	printf("-----------------------------\n");

	return result;
}

Matrix* finalizeResultMatrix(Matrix* simMatrix)
{
	Matrix* result;

	if ( simMatrix->nvtxs == 0 ||
	simMatrix->xadj[simMatrix->nvtxs] == 0 )
	{
		result = simMatrix;
		return result;
	}

	sortAdjLists(simMatrix->nvtxs, simMatrix->xadj,
			simMatrix->adjncy, simMatrix->adjwgt);
	result = simMatrix;
	Matrix* revMatrix = getTranspose(simMatrix);
	result = add(simMatrix, revMatrix);
	freeMatrix(revMatrix);
	freeMatrix(simMatrix);

	return result;
}

Matrix* LSH::bayesLSH()
{
	timer totalTimer;
	cleartimer(totalTimer);
	starttimer(totalTimer);

	float threshold = cosThreshold;
	if ( !opts.jaccard && !opts.useMinwiseBitSketches )
	{
		threshold = convertCosineForHashing(cosThreshold);
		sizeOfBlock = 8;
	}
	else if ( opts.useMinwiseBitSketches )
	{
		threshold = convertJaccardForHashing(cosThreshold);
		sizeOfBlock = 8;
	}
	else
		sizeOfBlock = 1;

	numBlocks = getNumBlocksForRecall(recall, threshold, sizeOfBlock);
	#ifndef NDEBUG
	printf("numBlocks %d, for recall %f, threshold %f, sizeofBlock %d\n",
	numBlocks, recall, threshold,sizeOfBlock);
	fflush(stdout);
	#endif

	gsl_set_error_handler_off();

	blsh = new BayesLSH(fvs, opts);
	blsh->setMinMinSketchSize(numBlocks*sizeOfBlock);

	int doneBuildingLSHIndexes = 0;
	if ( opts.jaccard && !opts.useMinwiseBitSketches &&
	!opts.uniformPrior && opts.betaPrior )
	{
		blsh->jaccardLSH_setupSketchesFirst();
		mis = blsh->mis;
		buildLSHIndexes();
		doneBuildingLSHIndexes = 1;

		Beta *bPrior = new Beta();
		setupBetaPrior_lsh(fvs, invArraysXadjs[0], invArrays[0],
		fvs->nDimensions, &jaccard, opts, bPrior);
		blsh->setupBetaPrior(bPrior);
		blsh->setup();
	}
	else
	{
		blsh->setupUnifPrior();
		blsh->setup();
	}

	if ( opts.jaccard && !opts.useMinwiseBitSketches )
		mis = blsh->mis;
	else
		cs = blsh->cs;

	processPair = &LSH::processPair_blshWrapper;

	if ( !doneBuildingLSHIndexes )
		buildLSHIndexes();

	Matrix* simMatrix;
	simMatrix = processIndexes();

	printf("Num. similar pairs: %d\n", simMatrix->xadj[simMatrix->nvtxs]);
	freeLSHIndexes();

	#ifndef NDEBUG
	printf("Num. not concentrated: %ld\n",
	blsh->numNotConcentrated);

	blsh->printPruneStatistics();
	#endif

	delete blsh;

	Matrix *result = finalizeResultMatrix(simMatrix);

	stoptimer(totalTimer);
	printf("-----------------------------\n");
	printf("Time for total computation: %f\n",
			gettimer(totalTimer));
	printf("-----------------------------\n");

	return result;
	
}


