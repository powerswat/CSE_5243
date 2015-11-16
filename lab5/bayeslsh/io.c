/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * io.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 * $Id: io.c,v 1.1.1.1 2011-12-14 18:42:37 venu Exp $
 *
 */

#include <all.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


void ReadBinaryMatrixFromBinary(Matrix *m, char *filename)
{
// Format is as follows:
// for each record:
// 	<recordId> <recordLength> <feature1> .. <featureN>
// each is an int
// Crucially, there's no information on how many records there
// are!

  FILE *fpin;
  if ((fpin = fopen(filename, "rb")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  struct stat st;
  stat(filename, &st);
  long size = st.st_size;
 // m->nnz = (size - sizeof(int)*2*m->nvtxs)/(sizeof(int));
  m->nnz = (size)/(sizeof(int));
//  printf("Upper-bound on nnz in %s as %d\n", filename, m->nnz);

  // Going to make the pessimistic assumption that there are as
  // many records as there are nnz
  m->xadj = idxmalloc(m->nnz+1, "ReadMatrix:xadj");
  m->adjncy = idxmalloc(m->nnz, "ReadMatrix: adjncy");
  m->nvtxs = 0;
  m->xadj[0] = 0;
  m->adjwgt = NULL;

  int rid;
  while ( fread(&rid, sizeof(int), 1, fpin) == 1 )
  {
  	int len, f;
  	fread(&len, sizeof(int), 1, fpin);
	fread( m->adjncy+m->xadj[m->nvtxs], sizeof(int), len, fpin);
	m->xadj[m->nvtxs+1] = m->xadj[m->nvtxs] + len;
	m->nvtxs++;

  }

  m->xadj = (int*) realloc(m->xadj, (m->nvtxs+1)*sizeof(int));
  if ( m->xadj[m->nvtxs] != m->nnz )
  {
	m->nnz = m->xadj[m->nvtxs];
	m->adjncy = (int*)realloc(m->adjncy, sizeof(int)*m->nnz);
  }

//  printf("Going to decrement all feature ids by 1\n");
  for ( int i=0; i<m->nnz; i++ )
  	m->adjncy[i]--;

  printf("Successfully read input\n");
  printf("Num. records: %d, total nnz: %d\n", m->nvtxs, m->nnz);
  fflush(stdout);

  if ( !checkSorted(m->nvtxs, m->xadj, m->adjncy) )
	  sortAdjLists(m->nvtxs, m->xadj, m->adjncy);

  fclose(fpin);
}

void ReadMatrixFromBinary(Matrix *m, char *filename, wgttype
threshold)
{
// The format is as follows:
// <numRecords> (int - 4 bytes>
// for each record:
// 		<recordLength> (int - 4 bytes) 
// 		<feature1> .. <featureN> (N ints)
// 		<featureWt1> .. <featureWtN> (N floats)

  FILE *fpin;
  if ((fpin = fopen(filename, "rb")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fread(&m->nvtxs, sizeof(int), 1, fpin);
		
//  printf("Number of records in %s is %d\n", filename, m->nvtxs);

  struct stat st;
  stat(filename, &st);
  long size = st.st_size;
  m->nnz = (size - sizeof(int)*m->nvtxs-1)
  			/(sizeof(int)+sizeof(float));
 // printf("Determined nnz in %s as %d\n", filename, m->nnz);

  m->xadj = idxmalloc(m->nvtxs+1, "ReadMatrix:xadj");
  m->adjncy = idxmalloc(m->nnz, "ReadMatrix:adjncy");
  long s = sizeof(float) * ((long) m->nnz);
  m->adjwgt = (float*) malloc(s);
  if ( m->adjwgt == NULL)
  {
  	printf("Yikes! Could not alloc m->adjwgt %ld floats\n", s);
	exit(0);
  }
  	
  m->xadj[0] = 0;
  for ( int i=0; i<m->nvtxs; i++ )
  {
  	int len;
	fread(&len, sizeof(int), 1, fpin);
	fread( m->adjncy+m->xadj[i], sizeof(int), len, fpin);
	fread( m->adjwgt+m->xadj[i], sizeof(float), len, fpin);
  	m->xadj[i+1] = m->xadj[i] + len; 
  }

  if ( m->xadj[m->nvtxs] != m->nnz )
  {
//  	printf("Yikes! promised %d nnz, but I got %d nnz\n", m->nnz,
//	m->xadj[m->nvtxs]);
	m->nnz = m->xadj[m->nvtxs];
	m->adjncy = (int*)realloc(m->adjncy, sizeof(int)*m->nnz);
	m->adjwgt = (float*)realloc(m->adjwgt, sizeof(float)*m->nnz);
  }

//  printf("Going to decrement all feature ids by 1\n");
  for ( int i=0; i<m->nnz; i++ )
  	m->adjncy[i]--;

//  printf("m->adjwgt[10]:%f\n", m->adjwgt[10]);

  printf("Successfully read input\n");
  printf("Num. records: %d, total nnz: %d\n", m->nvtxs, m->nnz);
  fflush(stdout);

  if ( !checkSorted(m->nvtxs, m->xadj, m->adjncy) )
  	sortAdjLists(m->nvtxs, m->xadj, m->adjncy, m->adjwgt);

  fclose(fpin);

}

void ReadGraph(GraphType *graph, const char *filename, int *wgtflag,
	int addSelfLoop, int txtFormat)
{
  int i, j, k, l, fmt, readew, readvw, ncon, edge, ewgt,
  actualnvtxs=0, noOfSelfLoops=0, allocedEdges; 
  
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  line = (char *)malloc(sizeof(char)*(MAXLINE+1));

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  if (feof(fpin)) {
    graph->nvtxs = 0;
    free(line);
    return;
  }

  fmt = ncon = 0;
  sscanf(line, "%d %d %d %d", &(graph->nvtxs), &(graph->nedges), &fmt, &ncon);
  
  readew = (fmt%10 > 0);
  readvw = ((fmt/10)%10 > 0);
  if (fmt >= 100) {
    printf("Cannot read this type of file format!");
    exit(0);
  }

  *wgtflag = 0;
  if (readew)
    *wgtflag += 1;
  if (readvw)
    *wgtflag += 2;

  if (ncon > 0 && !readvw) {
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("You specified ncon=%d, but the fmt parameter does not specify vertex weights\n", ncon);
    printf("Make sure that the fmt parameter is set to either 10 or 11.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  graph->nedges *=2;
  /* Venu: my addition */
  if ( addSelfLoop > 0 )
  	graph->nedges += graph->nvtxs;

  ncon = graph->ncon = (ncon == 0 ? 1 : ncon);

  /*printf("%d %d %d %d %d [%d %d]\n", fmt, fmt%10, (fmt/10)%10, ncon, graph->ncon, readew, readvw);*/

  if (graph->nvtxs > MAXIDX) 
    errexit("\nThe matrix is too big: %d [%d %d]\n", graph->nvtxs, MAXIDX, sizeof(idxtype));

  xadj = graph->xadj = idxsmalloc(graph->nvtxs+1, 0, "ReadGraph: xadj");
  adjncy = graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: adjncy");

  vwgt = graph->vwgt = (readvw ? idxmalloc(ncon*graph->nvtxs, "ReadGraph: vwgt") : NULL);
  adjwgt = graph->adjwgt = (readew ? idxmalloc(graph->nedges, "ReadGraph: adjwgt") : NULL);

  allocedEdges = graph->nedges;
  noOfSelfLoops = 0;

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0, actualnvtxs=0; i<graph->nvtxs;i++) 
  {
	int vertexid;
    
	do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));

	if ( feof(fpin) )
		break;
	
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    if (readvw) {
      for (l=0; l<ncon; l++) {
        vwgt[i*ncon+l] = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }
    }

	if ( txtFormat > 0 )
	{
		vertexid=(int)strtol(oldstr,&newstr,10)-1;
		oldstr=newstr;

		/* add missing nodes */
		while( i < vertexid )
		{
			if ( addSelfLoop > 0 )
			{
			/* self loop for missing nodes */
				xadj[i+1]=xadj[i]+1;
				adjncy[xadj[i]]=i;
				k++;
				if ( readew )
					adjwgt[xadj[i]]=1;
			}
			else
				xadj[i+1]=xadj[i];

			i++;
		}
		actualnvtxs++;
	}
	if ( addSelfLoop > 0 )
	{
		adjncy[k]=i; 
		/* For now, assigning an edge weight of 1 */
		// this weight is modified later on, being set to the
		// max. weight out of the node.
		if ( readew )
			adjwgt[k] = 1; 
		k++;
	}  

	idxtype maxwgt = 0;
	// keep track of max. edge weight to this node
    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      if (readew) {
        ewgt = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }

      if (edge < 0)
        break;

	  /* Venu: my addition. Since we've already added self-loop
	   * above, disregard self loop here. */
	  if ( edge == i && addSelfLoop > 0 )
	  {
		graph->nedges--;
		noOfSelfLoops++;
		continue;
	  }
	  // if addSelfLoop <= 0, do not include self loops.
	  else if ( edge == i )
	  {
	  	graph->nedges--;
		noOfSelfLoops++;
		continue;
	  }

	  if ( k >= allocedEdges )
	  {
	  	printf("No. of edges more than %d\n", allocedEdges );
		printf("Reading list of vertex %d\n",i+1);
		abort();
	  }

      adjncy[k] = edge;
      if (readew) 
	  {
        adjwgt[k] = ewgt;
		if ( ewgt > maxwgt )
			maxwgt = ewgt;
	  }
      k++;
    } 
    xadj[i+1] = k;
	if ( addSelfLoop > 0 && readew )
	{
		// weight on the self-loop is set to the max. edge weight
		// out of the node.
		adjwgt[xadj[i]] = maxwgt;
	}
  }


  /* add missing nodes at the end of the file */
  while( i < graph->nvtxs )
  {
  	/* self loop for missing nodes */
  	if ( addSelfLoop > 0 )
  	{
  		xadj[i+1]=xadj[i]+1;
  		adjncy[xadj[i]]=i;
  		k++;
  		if ( readew )
  			adjwgt[xadj[i]]=1;
  	}
  	else
  		xadj[i+1]=xadj[i];
  	i++;
  }

  if ( addSelfLoop == 0  && noOfSelfLoops > 0 )
  	printf("removed %d self loops\n", noOfSelfLoops);

  if ( graph->nedges != k )
  {
	  printf("expected no. of edges:%d, xadj[nvtxs]:%d\n", 
  			graph->nedges, k);
  }
	
  if (k > allocedEdges) {
	  if ( addSelfLoop > 0 )
	  {
	  	k -= graph->nvtxs;
		graph->nedges -= graph->nvtxs;
	  }
	if ( txtFormat <= 0 )
	{
		k /= 2;
		graph->nedges /= 2;
	}
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("In the first line of the file, you specified that the graph contained\n%d edges. However, I found %d edges in the file.\n", graph->nedges, k);
    if (2*k == graph->nedges) {
      printf("\n *> I detected that you specified twice the number of edges that you have in\n");
      printf("    the file. Remember that the number of edges specified in the first line\n");
      printf("    counts each edge between vertices v and u only once.\n\n");
    }
    printf("Please specify the correct number of edges in the first line of the file.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  graph->nedges = k;

  free(line);
}

void ReadMatrix(Matrix *m, char *filename, wgttype threshold)
{
  idxtype *xadj, *adjncy;
  wgttype *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  line = (char *)malloc(sizeof(char)*(MAXLINE+1));

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  int fmt=0, readew = 0;
  sscanf(line, "%d %d %d", &(m->nvtxs), &(m->nnz), &fmt);

  if ( fmt > 0 )
  	readew = 1;

  m->nnz *= 2;
  m->nnz++;
  m->xadj = idxmalloc(m->nvtxs+1, "ReadMatrix:xadj");
  m->adjncy = idxmalloc(m->nnz, "ReadMatrix:adjncy");
  if ( readew )
  {
	  m->adjwgt = (wgttype*)GKmalloc(sizeof(wgttype)*m->nnz, "ReadMatrix:adjwgt");
  }
  else
  	m->adjwgt = NULL;

  long k=0, numPruned = 0;
  int edge,i;
  wgttype wt;
  for (m->xadj[0]=0, k=0, i=0 ; i<m->nvtxs; i++) 
  {
    
	do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));

	if ( feof(fpin) )
		break;
	
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

   
    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;
	  
	  if ( readew ) 
	  {
	      wt = (wgttype)strtod(oldstr, &newstr);
    	  oldstr = newstr;
	  }

      if (edge < 0)
        break;

	  if ( readew && wt < threshold )
	  {
	  	numPruned++;
		continue;
	  }

	  if ( k >= m->nnz )
	  {
	  	printf("No. of edges more than %d\n", m->nnz);
		printf("Reading list of vertex %d\n",i+1);
		abort();
	  }

      m->adjncy[k] = edge;
	  if ( readew )
	      m->adjwgt[k] = wt;
      k++;
    } 
    m->xadj[i+1] = k;
  }

  if ( threshold > 0 || numPruned > 0 )
  {
  	printf("Pruned %ld edges with wt less than %f\n", numPruned,
	threshold);
  }
  if ( k != m->nnz )
  {
//  	printf("Promised edges: %d, Actual edges: %ld\n", m->nnz, k);
	m->adjncy = (idxtype*) realloc(m->adjncy, sizeof(idxtype)*k);
	if ( readew )
		m->adjwgt = (wgttype*) realloc(m->adjwgt, sizeof(wgttype)*k);
	m->nnz = k;
  }

  printf("Successfully read input\n");
  printf("Num. records: %d, Nnz: %d\n", m->nvtxs, m->nnz);
  fflush(stdout);

  if ( !checkSorted(m->nvtxs, m->xadj, m->adjncy) )
	  sortAdjLists(m->nvtxs, m->xadj, m->adjncy, m->adjwgt);

  fclose(fpin);
  
}

void WriteMatrixInBinaryFormat(const char *filename, int
nvtxs, int *xadj, int *adjncy, float *adjwgt)
{
  int i, j;
  FILE *fpout;

  if ((fpout = fopen(filename, "wb")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  // first add 1 to all features
  for ( i=0; i<xadj[nvtxs]; i++ )
  	adjncy[i]++;

  for ( i=0; i<nvtxs; i++ )
  {
	int len = xadj[i+1] - xadj[i];
	fwrite( &len, sizeof(int), 1 , fpout);

	if ( len > 0 )
	{
		fwrite(adjncy+xadj[i], sizeof(int), len, fpout);
		fwrite(adjwgt+xadj[i], sizeof(float), len, fpout);
	}
  }

  fclose(fpout);
}

void WriteBinaryMatrixInBinaryFormat(const char *filename, int
nvtxs, int *xadj, int *adjncy)
{
  int i, j;
  FILE *fpout;

  if ((fpout = fopen(filename, "wb")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  // first add 1 to all features
  for ( i=0; i<xadj[nvtxs]; i++ )
  	adjncy[i]++;

  for ( i=0; i<nvtxs; i++ )
  {
  	int vId = i+1;
	fwrite( &vId, sizeof(int), 1, fpout);

	int len = xadj[i+1] - xadj[i];
	fwrite( &len, sizeof(int), 1 , fpout);

	if ( len > 0 )
		fwrite(adjncy+xadj[i], sizeof(int), len, fpout);
  }

  fclose(fpout);

}

void WriteMatrix(const char *filename, int nvtxs, idxtype *xadj,
idxtype *adjncy, wgttype *adjwgt)
{
  int i, j;
  FILE *fpout;

  printf("Will write to file %s\n",filename);

  if ( adjwgt == NULL )
  {
  	printf("WriteMatrix: adjwgt is null\n");
	abort();
  }

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d 1", nvtxs, xadj[nvtxs]/2);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
	if ( xadj[i] == xadj[i+1] )
	{
/*		printf("No neighbours for node %d\n", i);
		abort(); */
		continue;
	}
	fprintf(fpout, "%d %.6f", adjncy[xadj[i]]+1,adjwgt[xadj[i]]);
    for (j=xadj[i]+1; j<xadj[i+1]; j++)
      fprintf(fpout, " %d %.6f", adjncy[j]+1, adjwgt[j]);
  }
  fprintf(fpout, "\n");

  fclose(fpout);
}


