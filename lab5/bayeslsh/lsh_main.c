
#include <all.h>

main(int argc, char *argv[])
{
	char filename[256], outputFile[256];
	int outputFileGiven = 0, similarityGiven = 0, thresholdGiven
	= 0;
	float recall = 0.97;
	int sizeOfBlock = 8;
	BssOptions opts;
	opts.jaccard = 0;
	int numHashes = 1024;
	int inputFormatBinary = 0;
	opts.lsh_useExact = 1;
	opts.binaryCosine = 0;

	float fudgeFactor = 0.00001;

	opts.availableMemoryInMB = 15*(1024/2); // 
	opts.availableRandomGaussiansMemoryInMB = 5*(1024/2);

	if ( argc < 2) 
	{
		fprintf(stderr, "Usage: %s <InputFile> [-t <threshold>]",
					argv[0]);
		fprintf(stderr, " [-s <similarity>] [-a <approximate?>]");
		fprintf(stderr, " [-e <expectedRecall>] [-h <numHashes>] [-r <initRandomSeed>]");
		fprintf(stderr, " [-f <inputFileFormat>] [-o <outputFile>]\n");
		fprintf(stderr, "-t, -s and -o are mandatory arguments\n");
		fprintf(stderr, "similarity: 0 - real, cosine; 1 - binary,jaccard, 2 - binary,cosine\n");
		fprintf(stderr, "inputFileFormat: 0 - Metis, 1 - binary\n");
		fprintf(stderr, "approximate: default is 0. If you set approximate 1, use -h to indicate number of hashes\n");
		fprintf(stderr, "randomSeed: default is 0\n");
		fprintf(stderr, "expectedRecall: default is 0.97\n");
		exit(0);
	}

	for (argv++; *argv != NULL; argv++)
	{
	    if ((*argv)[0] == '-')
		{
	      	switch ((*argv)[1])
			{
			case 'a':
			case 'A':
			  opts.lsh_useExact = atoi(*(++argv));
			  if ( opts.lsh_useExact )
			  	opts.lsh_useExact = 0;
			  else
			  	opts.lsh_useExact = 1;
			  break;
			case 'r':
			case 'R':
				opts.hash_initRandomSeed = atoi(*(++argv));
				break;
			case 'f':
			case 'F':
			  inputFormatBinary = atoi(*(++argv));
			  break;
			case 'e':
			case 'E':
			  recall = atof(*(++argv));
			  break;
			case 't':
			case 'T':
			  opts.cosThreshold = atof(*(++argv));
			  opts.cosThreshold -= fudgeFactor;
			  thresholdGiven = 1;
			  break;
			case 'h':
			case 'H':
			  numHashes = atoi(*(++argv));
			  break;
			case 's':
			case 'S':
			  opts.jaccard = atoi(*(++argv));
			  if ( opts.jaccard == 2 )
			  {
			  	opts.binaryCosine = 1;
				opts.jaccard = 0;
			  }
			  similarityGiven = 1;
			  break;
			case 'o':
			case 'O':
			  strcpy(outputFile,*(++argv));
			  outputFileGiven=1;
			  break;
			default:
			  printf("Invalid switch %s\n", *argv);
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
		printf("Please specify output file using -o option\n");
		return 1;
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

	if ( opts.cosThreshold <= 0 )
	{
		printf("Please specify a threshold > 0 using -t option\n");
		return 1;
	}

	int nnumHashes;
	if ( opts.jaccard )
		nnumHashes = 8 * (int) round(numHashes/8);
	else
		nnumHashes = 32 * (int) round(numHashes/32);
	if ( nnumHashes != numHashes )
	{
		printf("numHashes has been rounded out to %d\n",
		nnumHashes);
		numHashes = nnumHashes;
	}

	printf("Input arguments\n");
	printf("-------------\n");
	printf("Input file: %s\n", filename);
	printf("Similarity: %s\n", (opts.jaccard ? "Jaccard" : "Cosine"));
	printf("Similarity threshold: %f\n", opts.cosThreshold);
	printf("Approximate similarity using hashes? %s\n",
				(opts.lsh_useExact ? "No" : "Yes"));
	printf("If approximate, number of hashes for approximation: %d\n", numHashes);
	printf("Expected Recall: %f\n", recall);
	printf("-------------\n\n");

	
	Matrix featureVectors;
	RectMatrix *fvs;
	if ( inputFormatBinary )
	{
		if ( opts.jaccard || opts.binaryCosine)
			ReadBinaryMatrixFromBinary(&featureVectors,	filename);
		else
			ReadMatrixFromBinary(&featureVectors, filename);
	}
	else
	{
		if ( !opts.jaccard )
		{
			int inputFv_threshold = 0;
			ReadMatrix(&featureVectors, filename, inputFv_threshold);
		}
		else
		{
			GraphType g;
			int addSelfLoop=0;
			int wgtflag=0;
			ReadGraph(&g, filename, &wgtflag, addSelfLoop, 0);
			if ( g.adjwgt != NULL )
				free(g.adjwgt);

			featureVectors.nvtxs = g.nvtxs;
			featureVectors.xadj = g.xadj;
			featureVectors.adjncy = g.adjncy;
			featureVectors.adjwgt = NULL;
		}
	}
	
	fvs = convertMatrixToRectMatrix(&featureVectors);
	LSH plsh(fvs, opts, numHashes, recall);
	Matrix *ret = 	plsh.runRecall();

	printf("Nnz in similarity matrix: %d\n", ret->xadj[ret->nvtxs]);
	
	WriteMatrix(outputFile, ret->nvtxs, ret->xadj, ret->adjncy,
		ret->adjwgt);

}
