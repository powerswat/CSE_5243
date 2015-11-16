/*
 *
 * $Id: bayeslsh_main.c,v 1.1.1.1 2011-12-14 18:42:37 venu Exp $
 *
 */

#include <all.h>

void print_help(const char * program_name)
{
	fprintf(stderr, "Usage: %s <InputFile> [-t threshold] [-s <similarity] ", program_name);
	fprintf(stderr, "[-c <candidateGenerator>] [-l <BayesLSH_Lite?>] "); 
	fprintf(stderr, "[-e <epsilon>] [-g <gamma>] [-d <delta>] ");
	//fprintf(stderr, "[-n <noGaussians>] [-u <uniformPrior>] [-o outputFile]\n");
	fprintf(stderr, "[-r <randomSeed?>]");
	fprintf(stderr, " [-f <inputFileFormat>] [-o outputFile]\n");
	fprintf(stderr, "-t, -s, -c and -o are mandatory arguments\n");
	fprintf(stderr, "similarity: 0 - real, cosine; 1 - binary,jaccard, 2 - binary,cosine\n");
	fprintf(stderr, "candidateGenerator: 0 - AllPairs, 1 - LSH\n");
	fprintf(stderr, "inputFileFormat: 0 - Metis, 1 - binary\n");
	fprintf(stderr, "BayesLSH_Lite: default is 0\n");
	fprintf(stderr, "randomSeed: default is 0\n");
	fprintf(stderr, "epsilon, gamma: default is 0.03\n");
	fprintf(stderr, "delta: default is 0.05\n");
}

static inline float convertCosineForHashing(float sim)
{
	return 1.0-acos(sim)/PI;
}

static inline float convertJaccardForHashing(float sim)
{
	return (1.0+sim)/2;
}

int main(int argc, char *argv[])
{
	GraphType graph;
	char filename[256];
	char outputFile[256];
	int outputFileGiven = 0, similarityGiven = 0,
	thresholdGiven=0, generatorGiven=0;
	
	int inputFormatBinary = 0;

	float fudgeFactor = 0.00001;

	BssOptions opts;
	opts.numSamples = 10000;
	opts.num_os = 200;
	opts.cosThreshold = 0.01;
	opts.onDemandSketches = 1;
	opts.gaussians = 1;
	opts.useIntGaussians = 1;
	opts.uniformPrior = 1;
	opts.eps1 = 0.03;
	opts.eps2 = 0.03;
	opts.cosDelta = 0.05;
	opts.availableMemoryInMB = 15*(1024/2); // 
	opts.availableRandomGaussiansMemoryInMB = 9*(1024/2);
	opts.MIN_MIN_HASHES = 128;
	opts.HASHES_STEP = 256;
	opts.APPEND_HASHES = 512;
	opts.additionalPruning = 1;
	opts.backoffPrior = 1;
	opts.backoffPct = 0.5;
	opts.jaccard = 0;
	opts.binaryCosine = 0;
	opts.betaPrior = 0;
	opts.pwLinearPrior = 0;
	opts.useWithLSH = 0;
	opts.lsh_recall = 0.98;
	opts.doublePruning = 1;
	opts.pruneOnly = 0;
	opts.hash_initRandomSeed = 0;

	opts.useMinwiseBitSketches = 0;

	int graphInput = 0;

	if ( argc < 2 )
	{
		print_help(argv[0]);
		exit(0);
	}

	for (argv++; *argv != NULL; argv++)
	{
		if ((*argv)[0] == '-')
		{
			switch ((*argv)[1])
			{
			case 'f':
			case 'F':
			  inputFormatBinary = atoi(*(++argv));
			  break;
			case 'o':
			case 'O':
			  strcpy(outputFile,*(++argv));
			  outputFileGiven=1;
			  break;
			case 'r':
			case 'R':
				opts.hash_initRandomSeed = atoi(*(++argv));
				break;
			case 'l':
			case 'L':
				opts.pruneOnly = atoi(*(++argv));
				break;
			case 's':
			case 'S':
				similarityGiven = 1;
				opts.jaccard = atoi(*(++argv));
				if ( opts.jaccard == 2 )
				{
					opts.jaccard = 0;
					opts.binaryCosine = 1;
				}
				break;
			case 'g':
			case 'G':
				opts.eps2 = atof(*(++argv));
				break;
			case 'e':
			case 'E':
				opts.eps1 = atof(*(++argv));
				break;
			case 'd':
			case 'D':
				opts.cosDelta = atof(*(++argv));
				break;
			case 't':
			case 'T':
				opts.cosThreshold = atof(*(++argv));
				opts.cosThreshold -= fudgeFactor;
				thresholdGiven = 1;
				break;
			case 'c':
			case 'C':
				opts.useWithLSH = atoi(*(++argv));
				generatorGiven = 1;
				break;
			default:
			  printf("Invalid switch %s\n", *argv);
			  print_help(argv[0]);
			  exit(0);
			}
		}
	    else
		{
	      strcpy(filename, *argv);
	    }
	}

	if ( !outputFileGiven )
	{
		fprintf(stderr, "Ouput file not given!\n");
		exit(0);
	}
	if ( !similarityGiven )
	{
		fprintf(stderr, "Similarity type not given!\n");
		exit(0);
	}
	if ( !thresholdGiven )
	{
		fprintf(stderr, "Threshold not given!\n");
		exit(0);
	}
	if ( !generatorGiven )
	{
		fprintf(stderr, "Candidate generation algorithm not given!\n");
		exit(0);
	}


	if ( opts.jaccard && !opts.binaryCosine )
	{
		opts.betaPrior = 1;
		opts.uniformPrior = 0;
		opts.backoffPrior = 1;
		opts.backoffPct = 0.1;
		opts.HASHES_STEP = 32;
		opts.APPEND_HASHES = 64;
		opts.MIN_MIN_HASHES = 64;
	}

	if ( opts.uniformPrior == 3 )
	{
		opts.uniformPrior = 0;
		opts.pwLinearPrior = 1;
	}

	if ( !opts.jaccard || opts.useMinwiseBitSketches )
		opts.MIN_MATCHES_STEP = BayesLSH::BITS_MIN_MATCHES_STEP;
	else
	{
		opts.MIN_MATCHES_STEP = BayesLSH::INTS_MIN_MATCHES_STEP;
		if ( opts.pruneOnly )
			opts.MIN_MATCHES_STEP = 16;
	}

	if ( !opts.jaccard )
		opts.angleThreshold = convertCosineForHashing( opts.cosThreshold);
	else
		opts.angleThreshold =
		convertJaccardForHashing(opts.cosThreshold);

	if ( opts.binaryCosine )
	{
		opts.doublePruning = 1;
		opts.additionalPruning = 0;
	}
	else if ( !opts.jaccard )
		opts.doublePruning = 0;

	printf("Input arguments\n");
	printf("----------------\n");
	printf("Input file: %s\n", filename);
	printf("Algorithm: ");
	opts.useWithLSH ? (printf("LSH + ")) : (printf("AllPairs + ")) ;
	opts.pruneOnly ? (printf("BayesLSH-Lite")) : (printf("BayesLSH"));
	printf("\nSimilarity: %s\n", (opts.jaccard ? "Jaccard" : "Cosine"));
	printf("Similarity threshold: %f\n", opts.cosThreshold);
	printf("epsilon : %f\n", opts.eps1);
	printf("gamma : %f\n", opts.eps2);
	printf("delta : %f\n", opts.cosDelta);
	printf("Randomly initialize seed for hash functions? %s\n",
			(opts.hash_initRandomSeed ? "Yes" : "No"));
	printf("Input File Format: %s\n", (inputFormatBinary ?
			"Binary" : "Metis"));
	printf("Output file: %s\n", outputFile);
	printf("------------------\n\n");
	fflush(stdout);

	Matrix* ret;
	RectMatrix *inputFvs;
	InitRandom(0);

	if ( inputFormatBinary )
	{
		Matrix fvs;
		if ( opts.jaccard || opts.binaryCosine )
			ReadBinaryMatrixFromBinary(&fvs, filename);
		else
			ReadMatrixFromBinary(&fvs, filename);

		inputFvs = convertMatrixToRectMatrix(&fvs);
	}
	else
	{
		if ( opts.jaccard || opts.binaryCosine )
		{
			int wgtflag, addSelfLoop=0, txtFormat=0;
			ReadGraph(&graph, filename, &wgtflag, addSelfLoop, txtFormat);
	//		InitRandom(time(NULL));

			sortAdjLists(graph.nvtxs, graph.xadj, graph.adjncy);

			Matrix fvs;
			fvs.nvtxs = graph.nvtxs;
			fvs.xadj = graph.xadj;
			fvs.adjncy = graph.adjncy;
			fvs.adjwgt = NULL;
			fvs.nnz = graph.nedges;
			inputFvs = convertMatrixToRectMatrix(&fvs);
		}
		else
		{
			Matrix featureVectors;
			int inputFv_threshold = 0;
			ReadMatrix(&featureVectors, filename, inputFv_threshold);
			
			inputFvs = convertMatrixToRectMatrix(&featureVectors);

		}
	}

	if ( !opts.useWithLSH )
	{
		if ( opts.jaccard || opts.binaryCosine )
		{
			AllPairsBinary bjs(inputFvs, opts);
		//	ret = bjs.run();
			ret = bjs.bayesLSH();
		}
		else
		{
			AllPairsRealCosine bss(inputFvs, opts);
		//	ret = bss.run();
			ret = bss.bayesLSH();
		}
	}
	else
	{
		LSH plsh(inputFvs, opts, 0, opts.lsh_recall);	
	//	ret = plsh.lshPlusBayesLSH();
		ret = plsh.bayesLSH();
	}


	WriteMatrix(outputFile, ret->nvtxs, ret->xadj,
		ret->adjncy, ret->adjwgt);

	return 0;

}

