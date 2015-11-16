/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * myqsort.c
 * 
 * This file contains a fast idxtype increasing qsort algorithm.
 * Addopted from TeX
 * 
 * Started 10/18/96
 * George
 * 
 * $Id: myqsort.c,v 1.1.1.1 2011-12-14 18:42:37 venu Exp $
 */

#include <all.h>			/* only for type declarations */

#define		THRESH		1	/* threshold for insertion */
#define		MTHRESH		6	/* threshold for median */


static void siqst(idxtype *, idxtype *);
static void iiqst(int *, int *);

/*************************************************************************
* Entry point of idxtype increasing sort
**************************************************************************/
void iidxsort(int n, idxtype *base)
{
  register idxtype *i;
  register idxtype *j;
  register idxtype *lo;
  register idxtype *hi;
  register idxtype *min;
  register idxtype c;
  idxtype *max;

  if (n <= 1)
    return;

  max = base + n;

  if (n >= THRESH) {
    siqst(base, max);
    hi = base + THRESH;
  }
  else 
    hi = max;

  for (j = lo = base; lo++ < hi;) {
    if (*j > *lo)
      j = lo;
  }
  if (j != base) { /* swap j into place */
    c = *base;
    *base = *j;
    *j = c;
  }

  for (min = base; (hi = min += 1) < max;) {
    while (*(--hi) > *min);
    if ((hi += 1) != min) {
      for (lo = min + 1; --lo >= min;) {
	c = *lo;
	for (i = j = lo; (j -= 1) >= hi; i = j)
	   *i = *j;
	*i = c;
      }
    }
  }
}

static void siqst(idxtype *base, idxtype *max)
{
  register idxtype *i;
  register idxtype *j;
  register idxtype *jj;
  register idxtype *mid;
  register int ii;
  register idxtype c;
  idxtype *tmp;
  int lo;
  int hi;

  lo = max - base;		/* number of elements as idxtype */
  do {
    mid = base + ((unsigned) lo>>1);
    if (lo >= MTHRESH) {
      j = (*base > *mid ? base : mid);
      tmp = max - 1;
      if (*j > *tmp) {
        j = (j == base ? mid : base); /* switch to first loser */
        if (*j < *tmp)
          j = tmp;
      }

      if (j != mid) {  /* SWAP */ 
        c = *mid;
        *mid = *j;
        *j = c;
      }
    }

    /* Semi-standard quicksort partitioning/swapping */
    for (i = base, j = max - 1;;) {
      while (i < mid && *i <= *mid)
        i++;
      while (j > mid) {
        if (*mid <= *j) {
          j--;
          continue;
        }
        tmp = i + 1;	/* value of i after swap */
        if (i == mid) 	/* j <-> mid, new mid is j */
          mid = jj = j;
        else 		/* i <-> j */
          jj = j--;
        goto swap;
      }

      if (i == mid) 
	break;
      else {		/* i <-> mid, new mid is i */
        jj = mid;
        tmp = mid = i;	/* value of i after swap */
        j--;
      }
swap:
      c = *i;
      *i = *jj;
      *jj = c;
      i = tmp;
    }

    i = (j = mid) + 1;
    if ((lo = j - base) <= (hi = max - i)) {
      if (lo >= THRESH)
        siqst(base, j);
      base = i;
      lo = hi;
    }
    else {
      if (hi >= THRESH)
        siqst(i, max);
      max = j;
    }
  } while (lo >= THRESH);
}





/*************************************************************************
* Entry point of int increasing sort
**************************************************************************/
void iintsort(int n, int *base)
{
  register int *i;
  register int *j;
  register int *lo;
  register int *hi;
  register int *min;
  register int c;
  int *max;

  if (n <= 1)
    return;

  max = base + n;

  if (n >= THRESH) {
    iiqst(base, max);
    hi = base + THRESH;
  }
  else 
    hi = max;

  for (j = lo = base; lo++ < hi;) {
    if (*j > *lo)
      j = lo;
  }
  if (j != base) { /* swap j into place */
    c = *base;
    *base = *j;
    *j = c;
  }

  for (min = base; (hi = min += 1) < max;) {
    while (*(--hi) > *min);
    if ((hi += 1) != min) {
      for (lo = min + 1; --lo >= min;) {
	c = *lo;
	for (i = j = lo; (j -= 1) >= hi; i = j)
	   *i = *j;
	*i = c;
      }
    }
  }
}


static void iiqst(int *base, int *max)
{
  register int *i;
  register int *j;
  register int *jj;
  register int *mid;
  register int ii;
  register int c;
  int *tmp;
  int lo;
  int hi;

  lo = max - base;		/* number of elements as ints */
  do {
    mid = base + ((unsigned) lo>>1);
    if (lo >= MTHRESH) {
      j = (*base > *mid ? base : mid);
      tmp = max - 1;
      if (*j > *tmp) {
        j = (j == base ? mid : base); /* switch to first loser */
        if (*j < *tmp)
          j = tmp;
      }

      if (j != mid) {  /* SWAP */ 
        c = *mid;
        *mid = *j;
        *j = c;
      }
    }

    /* Semi-standard quicksort partitioning/swapping */
    for (i = base, j = max - 1;;) {
      while (i < mid && *i <= *mid)
        i++;
      while (j > mid) {
        if (*mid <= *j) {
          j--;
          continue;
        }
        tmp = i + 1;	/* value of i after swap */
        if (i == mid) 	/* j <-> mid, new mid is j */
          mid = jj = j;
        else 		/* i <-> j */
          jj = j--;
        goto swap;
      }

      if (i == mid) 
	break;
      else {		/* i <-> mid, new mid is i */
        jj = mid;
        tmp = mid = i;	/* value of i after swap */
        j--;
      }
swap:
      c = *i;
      *i = *jj;
      *jj = c;
      i = tmp;
    }

    i = (j = mid) + 1;
    if ((lo = j - base) <= (hi = max - i)) {
      if (lo >= THRESH)
        iiqst(base, j);
      base = i;
      lo = hi;
    }
    else {
      if (hi >= THRESH)
        iiqst(i, max);
      max = j;
    }
  } while (lo >= THRESH);
}

int RandomPartitionInts(idxtype* a, int start, int end)
{
	int n=end-start+1,i,j;
	int tmp,x;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
		SWAP(a[i],a[end],tmp);

	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmp);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmp);

	return i+1;
}

int RandomPartition(wgttype* a, int start, int end)
{
	int n=end-start+1,i,j;
	wgttype tmp,x;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
		SWAP(a[i],a[end],tmp);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmp);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmp);

	return i+1;
}

int ParallelRandomPartitionInts(idxtype* a, idxtype* b, int start,
									int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	idxtype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

void ParallelQSortInts(idxtype *a, idxtype *b, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartitionInts(a,b,start,end);
		ParallelQSortInts(a,b,start,q-1);
		ParallelQSortInts(a,b,q+1,end);
	}
}

int ParallelRandomPartitionIntsUsingScores(idxtype* a, idxtype*
				b, idxtype *scores,	int start, int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	idxtype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

	//x=a[end];
	x = scores[a[end]];
	i=start-1;
	for(j=start;j<end;j++)
	{
		//if (a[j] <= x)
		if ( scores[a[j]] <= x )
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

// we will use scores as the comparer to sort a and b together. 
void ParallelQSortIntsUsingScores(idxtype *a, idxtype *b, idxtype
*scores, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartitionIntsUsingScores(a,b,scores,start,end);
		ParallelQSortIntsUsingScores(a,b,scores,start,q-1);
		ParallelQSortIntsUsingScores(a,b,scores,q+1,end);
	}
}

int RandomPartitionIntsUsingInts(idxtype* a, const idxtype *scores,
				int start, int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	idxtype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
	}

	//x=a[end];
	x = scores[a[end]];
	i=start-1;
	for(j=start;j<end;j++)
	{
		//if (a[j] <= x)
		if ( scores[a[j]] <= x )
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);

	return i+1;
}

// we will use scores as the comparer to sort a. 
void QSortIntsUsingInts(idxtype *a, idxtype *scores, int start,
int end)
{
	int q;
	if ( start < end )
	{
		q=RandomPartitionIntsUsingInts(a,scores,start,end);
		QSortIntsUsingInts(a,scores,start,q-1);
		QSortIntsUsingInts(a,scores,q+1,end);
	}
}

/* This routine partitions both the arrays a and b, according to
 * the order imposed by a. */
int ParallelRandomPartitionFloatsInts(wgttype* a, idxtype* b, int start,
int end)
{
	int n=end-start+1,i,j;
	int tmpidx;
	wgttype tmpwgt, x;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
	{
		SWAP(a[i],a[end],tmpwgt);
		SWAP(b[i],b[end],tmpidx);
	}

	x=a[end];
	i=start-1;
	for(j=start;j<end;j++)
	{
		if (a[j] <= x)
		{
			i++;
			if ( i == j )
				continue;
			SWAP(a[i],a[j],tmpwgt);
			SWAP(b[i],b[j],tmpidx);
		}
	}
	SWAP(a[i+1],a[end],tmpwgt);
	SWAP(b[i+1],b[end],tmpidx);

	return i+1;
}

/* sorts both a and b, according to the order imposed by a.*/
void ParallelQSortFloatsInts(wgttype *a, idxtype *b, int start,
int end)
{
	static int count=0;
	int q;
	count++;
	if ( count % 10000 == 0 )
	{
		printf("sortFloatsInts called %d times\n", count);
		fflush(stdout);
	}
/*	if ( start < end )
	{
		q=ParallelRandomPartitionFloatsInts(a,b,start,end);
		ParallelQSortFloatsInts(a,b,start,q-1);
		ParallelQSortFloatsInts(a,b,q+1,end);
	} */
	while ( start < end )
	{
		q=ParallelRandomPartitionFloatsInts(a,b,start,end);
		if ( q-start > end-q )
		{
			ParallelQSortFloatsInts(a,b,start,q-1);
			start = q+1;
		}
		else
		{
			ParallelQSortFloatsInts(a,b,q+1, end);
			end = q-1;
		}
	}
}

/* This routine partitions both the arrays a and b, according to
 * the order imposed by a. */
int ParallelRandomPartitionLongs(long* a, wgttype* b, int start,
int end)
{
	int n=end-start+1,i,j;
	long tmpidx,x;
	wgttype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

/* Sorts the arrays a and b, according to the order imposed by a.
   [start,end] is the inclusive range of the two arrays i.e.
   the length of the two arrays is end-start+1, not end-start. */
void ParallelQSortLongs(long *a, wgttype *b, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartitionLongs(a,b,start,end);
		ParallelQSortLongs(a,b,start,q-1);
		ParallelQSortLongs(a,b,q+1,end);
	}
}

/* This routine partitions both the arrays a and b, according to
 * the order imposed by a. */
int ParallelRandomPartition(idxtype* a, wgttype* b, int start,
int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	wgttype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

/* Sorts the arrays a and b, according to the order imposed by a.
   [start,end] is the inclusive range of the two arrays i.e.
   the length of the two arrays is end-start+1, not end-start. */
void ParallelQSort(idxtype *a, wgttype *b, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartition(a,b,start,end);
		ParallelQSort(a,b,start,q-1);
		ParallelQSort(a,b,q+1,end);
	}
}

idxtype RandomSelectInts(idxtype *a, int start, int end, int i)
{
	int q,k,j;
	if (start==end)
		return a[start];
	
	q = RandomPartitionInts(a,start,end);

//	printf("After rand partition, q:%d, end:%d\n",q,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
	abort();
*/
	k=q-start+1;
	if (k==i)
		return a[q];
	if (i < k )
		return RandomSelectInts(a,start,q-1,i);
	else
		return RandomSelectInts(a,q+1,end,i-k);
}

wgttype RandomSelect(wgttype *a, int start, int end, int i)
{
	int q,k,j;
	if (start==end)
		return a[start];
	
	q=RandomPartition(a,start,end);

//	printf("After rand partition, q:%d, end:%d\n",q,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
	abort();
*/
	k=q-start+1;
	if (k==i)
		return a[q];
	if (i < k )
		return RandomSelect(a,start,q-1,i);
	else
		return RandomSelect(a,q+1,end,i-k);
}


