
#include <all.h>

main(int argc, char *argv[])
{
	char filename[256], outputFile[256];
	float threshold = 0;
	int binary = 0, inputFormatBinary = 0 ;
	float fudgeFactor = 0.00001;

	int outputFileGiven = 0, similarityGiven = 0,
	thresholdGiven=0, generatorGiven=0;

	if ( argc < 2) 
	{
		fprintf(stderr, "Usage: %s <InputFile> [-t <threshold>]",
		argv[0]);
		fprintf(stderr, " [-s <similarity>] [-f <inputFileFormat>] [-o <outputFile>]\n");
		fprintf(stderr, "-t, -s and -o are mandatory arguments\n");
		fprintf(stderr, "similarity: 0 - real, cosine; 1 - binary,jaccard, 2 - binary,cosine\n");
		fprintf(stderr, "inputFileFormat: 0 - Metis, 1 - binary\n");
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
			case 's':
			case 'S':
			  binary = atoi(*(++argv));
			  similarityGiven = 1;
			  break;
			case 't':
			case 'T':
			  threshold = atof(*(++argv));
			  threshold -= fudgeFactor;
			  thresholdGiven = 1;
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

	if ( threshold <= 0 )
	{
		printf("Please specify a threshold > 0 using -t option\n");
		return 1;
	}

	printf("Input arguments:\n");
	printf("-------------------\n");
	printf("Input file: %s\n", filename);
	printf("Similarity: %s\n", ( binary==1 ? "Jaccard" : "Cosine"));
	printf("Similarity Threshold: %f\n", threshold);
	printf("Output file: %s\n", outputFile);
	printf("-------------------\n");
	fflush(stdout);

	Matrix featureVectors;
	if ( inputFormatBinary )
	{
		if ( binary )
			ReadBinaryMatrixFromBinary(&featureVectors,	filename);
		else
			ReadMatrixFromBinary(&featureVectors, filename);
	}
	else
	{
		if ( binary )
		{
			GraphType g;
			int wgtflag, txtFormat=0, addSelfLoop=0;
			ReadGraph(&g, filename, &wgtflag, addSelfLoop,
			txtFormat);
			if ( g.adjwgt != NULL )
			{
				printf("Warning: graph contains weights, but ");
				printf("they will be ignored\n");
			}
			featureVectors.nvtxs = g.nvtxs;
			featureVectors.nnz = g.nedges;
			featureVectors.xadj = g.xadj;
			featureVectors.adjncy = g.adjncy;
			featureVectors.adjwgt = NULL;
		}
		else
		{
			int inputFv_threshold = 0;
			ReadMatrix(&featureVectors, filename, inputFv_threshold);
		}
	}
	
	RectMatrix* fvs = convertMatrixToRectMatrix(&featureVectors);
	AllPairsOpts bo;
	bo.threshold = threshold;

	AllPairsSimSearch bss(fvs, bo, binary);
	Matrix *ret = bss.run();
	printf("Nnz in similarity matrix: %d\n", ret->xadj[ret->nvtxs]);
	
	WriteMatrix(outputFile, ret->nvtxs, ret->xadj, ret->adjncy,
		ret->adjwgt);
}
