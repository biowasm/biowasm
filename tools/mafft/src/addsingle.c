#include "mltaln.h"

#define SMALLMEMORY 1

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

static int treein;
static int topin;
static int treeout;
static int distout;
static int noalign;
static int multidist;
static int maxdist = 2; // scale -> 2bai
static int allowlongadds;
static int addtotop;
static int keeplength;
static int diffout; // Incomplete. Not used
static int ndeleted;
static int mapout;
static int smoothing;
static double hitout;

static int tuplesize;

#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define D6LENFACA 0.01
#define D6LENFACB 2500
#define D6LENFACC 2500
#define D6LENFACD 0.1
#define D10LENFACA 0.01
#define D10LENFACB 1000000
#define D10LENFACC 1000000
#define D10LENFACD 0.0

typedef struct _thread_arg
{
	int njob; 
	int nadd; 
	int *nlen; 
	int *follows; 
	char **name; 
	char **seq; 
	LocalHom **localhomtable; 
	double **iscore; 
	double **nscore; 
	int *istherenewgap; 
	int **newgaplist;
	RNApair ***singlerna; 
	double *eff_kozo_mapped; 
	int alloclen;
	Treedep *dep;
	int  ***topol;
	double  **len;
	Addtree *addtree;
	GapPos **deletelist;
	GapPos **difflist;
#ifdef enablemultithread
	int *iaddshare;
	int thread_no;
	pthread_mutex_t *mutex_counter;
#endif
} thread_arg_t;


#ifdef enablemultithread
typedef struct _gaplist2alnxthread_arg
{
//	int thread_no;
	int ncycle;
	int *jobpospt;
	int tmpseqlen;
	int lenfull;
	char **seq;
	int *newgaplist;
	int *posmap;
	pthread_mutex_t *mutex;
} gaplist2alnxthread_arg_t;

typedef struct _distancematrixthread_arg
{
	int thread_no;
	int njob;
	int norg;
	int *jobpospt;
	int **pointt;
	int *nogaplen;
	double **imtx;
	double **nmtx;
	double *selfscore;
	pthread_mutex_t *mutex;
} distancematrixthread_arg_t;

typedef struct _jobtable2d
{
    int i;  
    int j;  
} Jobtable2d;

typedef struct _dndprethread_arg
{
	int njob;
	int thread_no;
	double *selfscore;
	double **mtx;
	char **seq;
	Jobtable2d *jobpospt;
	pthread_mutex_t *mutex;
} dndprethread_arg_t;

#endif

typedef struct _blocktorealign
{
	int start;
	int end;
	int nnewres;
} Blocktorealign;

static void cnctintvec( int *res, int *o1, int *o2 )
{
	while( *o1 != -1 ) *res++ = *o1++;
	while( *o2 != -1 ) *res++ = *o2++;
	*res = -1;
}

static void countnewres( int len, Blocktorealign *realign, int *posmap, int *gaplist )
{
	int i, regstart, regend, len1;
	regstart = 0;
	len1 = len+1;
	for( i=0; i<len1; i++ )
	{

		if( realign[i].nnewres || gaplist[i] )
		{
			regend = posmap[i]-1;
			realign[i].start = regstart;
			realign[i].end = regend;
		}
		if( gaplist[i] )
		{
			realign[i].nnewres++;
//			fprintf( stderr, "hit? reg = %d-%d\n", regstart, regend );
		}
		regstart = posmap[i]+1;
	}
}
static void fillgap( char *s, int len )
{
	int orilen = strlen( s );
	s += orilen;
	len -= orilen;
	while( len-- )
		*s++ = '-';
	*s = 0;
}

static int lencomp( const void *a, const void *b ) // osoikamo
{
	char **ast = (char **)a;
	char **bst = (char **)b;
	int lena = strlen( *ast );
	int lenb = strlen( *bst );
//	fprintf( stderr, "a=%s, b=%s\n", *ast, *bst );
//	fprintf( stderr, "lena=%d, lenb=%d\n", lena, lenb );
	if( lena > lenb ) return -1;
	else if( lena < lenb ) return 1;
	else return( 0 );
}

static int dorealignment_tree( Blocktorealign *block, char **fullseq, int *fullseqlenpt, int norg, int ***topol, int *follows )
{
	int i, j, k, posinold, newlen, *nmem;
	int n0, n1, localloclen, nhit, hit1, hit2;
	int *pickhistory;
	int nprof1, nprof2, pos, zure;
	char **prof1, **prof2;
	int *iinf0, *iinf1;
	int *group, *nearest, *g2n, ngroup;
	char ***mem;
	static char **tmpaln0 = NULL;
	static char **tmpaln1 = NULL;
	static char **tmpseq;
	int ***topolpick;
	int *tmpint;
	int *intptr, *intptrx;
	char *tmpseq0, *cptr, **cptrptr;


	localloclen = 4 * ( block->end - block->start + 1 );	 // ookisugi?
	tmpaln0 = AllocateCharMtx( njob, localloclen );
	tmpaln1 = AllocateCharMtx( njob, localloclen );
	tmpseq = AllocateCharMtx( 1, *fullseqlenpt * 4 );
	iinf0 = AllocateIntVec( njob );
	iinf1 = AllocateIntVec( njob );
	nearest = AllocateIntVec( njob ); // oosugi

	posinold = block->start;

	n0 = 0;
	n1 = 0;
	for( i=0; i<njob; i++ )
	{
		strncpy( tmpseq[0], fullseq[i] + block->start, block->end - block->start + 1 );
		tmpseq[0][block->end - block->start + 1] = 0;
		commongappick( 1, tmpseq );
		if( tmpseq[0][0] != 0 )
		{
			if( i < norg )
			{
				fprintf( stderr, "BUG!!!!\n" );
				exit( 1 );
			}
			strcpy( tmpaln0[n0], tmpseq[0] );
			iinf0[n0] = i;
			nearest[n0] = follows[i-norg];
			n0++;
		}
		else
		{
			strcpy( tmpaln1[n0], "" );
			iinf1[n1] = i;
			n1++;
		}
	}
	mem = AllocateCharCub( n0, n0+1, 0 ); // oosugi
	nmem = AllocateIntVec( n0 ); // oosugi
	g2n = AllocateIntVec( n0 ); // oosugi
	group = AllocateIntVec( n0 ); // oosugi
	for( i=0; i<n0; i++ ) mem[i][0] = NULL;
	for( i=0; i<n0; i++ ) nmem[i] = 0;
	ngroup = 0;
	for( i=0; i<n0; i++ )
	{
		for( j=0; j<i; j++ ) if( nearest[j] == nearest[i] ) break;
		if( j == i ) group[i] = ngroup++;
		else group[i] = group[j];

		for( j=0; mem[group[i]][j]; j++ )
			;
		mem[group[i]][j] = tmpaln0[i];
		mem[group[i]][j+1] = NULL;
		nmem[group[i]]++;
		g2n[group[i]] = nearest[i];
//		fprintf( stderr, "%d -> %d -> group%d\n", i, nearest[i], group[i] );
//		fprintf( stderr, "mem[%d][%d] = %s\n", group[i], j, mem[group[i]][j] );
	}

	for( i=0; i<ngroup; i++ )
	{
//		fprintf( stderr, "before sort:\n" );
//		for( j=0; j<nmem[i]; j++ ) fprintf( stderr, "%s\n", mem[i][j] );
//		fprintf( stderr, "\n" );
		qsort( mem[i], nmem[i], sizeof( char * ), lencomp );
//		fprintf( stderr, "after sort:\n" );
//		for( j=0; j<nmem[i]; j++ ) fprintf( stderr, "%s\n", mem[i][j] );
//		fprintf( stderr, "\n" );
	}

#if 0
	for( i=1; i<n0; i++ )
	{
		profilealignment( 1, n1, i, tmpaln0+i, tmpaln1, tmpaln0, localloclen, alg );
	}
	newlen = strlen( tmpaln0[0] );
	for( i=0; i<n1; i++ ) eq2dash( tmpaln1[i] );
#else
//	newlen = 0;
	for( i=0; i<ngroup; i++ )
	{
//		for( k=0; mem[i][k]; k++ ) fprintf( stderr, "mem[%d][%d] = %s\n", i, j, mem[i][k] );

		for( j=1; j<nmem[i]; j++ )
		{
			profilealignment2( 1, j, mem[i]+j, mem[i], localloclen, alg );
		}
//		for( j=0; j<nmem[i]; j++ ) fprintf( stderr, "j=%d, %s\n", j, mem[i][j] );

#if 0 // iru
		if( ( j = strlen( mem[i][0] ) ) > newlen ) newlen = j;
		for( j=0; j<=i; j++ )
		{
			for( k=0; mem[j][k]; k++ )
				fillgap( mem[j][k], newlen );
		}
#endif

	}
#if 0
	fprintf( stderr, "After ingroupalignment (original order):\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "%s\n", tmpaln0[i] );
#endif
#endif

	topolpick = AllocateIntCub( ngroup, 2, ngroup );
	pickhistory = AllocateIntVec( ngroup );
	tmpint = AllocateIntVec( 2 );
	prof1 = AllocateCharMtx( n0, 0 );
	prof2 = AllocateCharMtx( n0, 0 );
	for( i=0; i<ngroup; i++ )
	{
		topolpick[i][0][0] = -1;
		topolpick[i][1][0] = -1;
		pickhistory[i] = -1;
	}

	nhit = 0;
	for( i=0; i<norg-1; i++ )
	{
		for( intptr=topol[i][0]; *intptr>-1; intptr++ )
		{
			for( intptrx=g2n,k=0; k<ngroup; intptrx++,k++ ) 
			{
				if( *intptr == *intptrx )
				{
					hit1 = k;
					goto exitloop1;
				}
			}
		}
		continue; 
		exitloop1:
//		fprintf( stderr, "hit1! group%d -> %d\n", k, topol[i][0][j] );

		for( intptr=topol[i][1]; *intptr>-1; intptr++ )
		{
			for( intptrx=g2n,k=0; k<ngroup; intptrx++,k++ ) 
			{
				if( *intptr == *intptrx )
				{
					hit2 = k;
					goto exitloop2;
				}
			}
		}
		continue; 
		exitloop2:

		if( pickhistory[hit1] == -1 )
		{
			topolpick[nhit][0][0] = hit1;
			topolpick[nhit][0][1] = -1;
		}
		else
		{
			intcpy( topolpick[nhit][0], topolpick[pickhistory[hit1]][0] );
			intcat( topolpick[nhit][0], topolpick[pickhistory[hit1]][1] );
		}
		if( pickhistory[hit2] == -1 )
		{
			topolpick[nhit][1][0] = hit2;
			topolpick[nhit][1][1] = -1;
		}
		else
		{
			intcpy( topolpick[nhit][1], topolpick[pickhistory[hit2]][0] );
			intcat( topolpick[nhit][1], topolpick[pickhistory[hit2]][1] );
		}

		pickhistory[hit1] = nhit;
		pickhistory[hit2] = nhit;
		nhit++;
//		g2n[hit1] = -1;
//		g2n[hit2] = -1;

//		fprintf( stderr, "hit2! group%d -> %d\n", k, topol[i][1][j] );

#if 0
		fprintf( stderr, "\nHIT!!! \n" );
		fprintf( stderr, "\nSTEP %d\n", i );
		for( j=0; topol[i][0][j]>-1; j++ ) fprintf( stderr, "%3d ", topol[i][0][j] );
		fprintf( stderr, "\n" );
		for( j=0; topol[i][1][j]>-1; j++ ) fprintf( stderr, "%3d ", topol[i][1][j] );
		fprintf( stderr, "\n" );
#endif
	}

	for( i=0; i<ngroup-1; i++ )
	{
#if 0
		fprintf( stderr, "\npickSTEP %d\n", i );
		for( j=0; topolpick[i][0][j]>-1; j++ ) fprintf( stderr, "%3d ", topolpick[i][0][j] );
		fprintf( stderr, "\n" );
		for( j=0; topolpick[i][1][j]>-1; j++ ) fprintf( stderr, "%3d ", topolpick[i][1][j] );
		fprintf( stderr, "\n" );
#endif

		pos = 0;
//		for( j=0; topolpick[i][0][j]>-1; j++ ) for( k=0; (cptr=mem[topolpick[i][0][j]][k]); k++ ) prof1[pos++] = cptr;
		for( intptr=topolpick[i][0]; *intptr>-1; intptr++ ) 
			for( cptrptr=mem[*intptr]; (cptr=*cptrptr); cptrptr++ ) 
				prof1[pos++] = cptr;
		nprof1 = pos;
		pos = 0;
//		for( j=0; topolpick[i][1][j]>-1; j++ ) for( k=0; (cptr=mem[topolpick[i][1][j]][k]); k++ ) prof2[pos++] = cptr;
		for( intptr=topolpick[i][1]; *intptr>-1; intptr++ ) 
			for( cptrptr=mem[*intptr]; (cptr=*cptrptr); cptrptr++ ) 
				prof2[pos++] = cptr;
		nprof2 = pos;


		profilealignment2( nprof1, nprof2, prof1, prof2, localloclen, alg );
#if 0
		for( j=0; j<nprof1; j++ ) fprintf( stderr, "prof1[%d] = %s\n", j, prof1[j] );
		for( j=0; j<nprof2; j++ ) fprintf( stderr, "prof2[%d] = %s\n", j, prof2[j] );
		fprintf( stderr, "done.\n" );
#endif
	}
	newlen = strlen( tmpaln0[0] );
	for( j=0; j<n1; j++ ) fillgap( tmpaln1[j], newlen );

#if 0
	fprintf( stderr, "After rerealignment (original order):\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "%s\n", tmpaln0[i] );
#endif

//	newlen = strlen( tmpaln0[0] );
	zure = ( block->end - block->start + 1 - newlen );
//	fprintf( stderr, "zure = %d, localloclen=%d, newlen=%d\n", zure, localloclen, newlen );


	if( *fullseqlenpt < strlen( fullseq[0] ) - (block->end-block->start+1)  + newlen + 1 )
	{
		*fullseqlenpt = strlen( fullseq[0] ) * 2;
		fprintf( stderr, "reallocating..." );
		for( i=0; i<njob; i++ )
		{
			fullseq[i] = realloc( fullseq[i], *fullseqlenpt * sizeof( char ) );
			if( !fullseq[i] )
			{
				fprintf( stderr, "Cannot reallocate seq[][]\n" );
				exit( 1 );
			}
		}
		fprintf( stderr, "done.\n" );
	}


	tmpseq0 = tmpseq[0];
	posinold = block->end+1;
	for( i=0; i<n0; i++ ) 
	{
		strncpy( tmpseq0, tmpaln0[i], newlen );
		strcpy( tmpseq0+newlen, fullseq[iinf0[i]] + posinold );
		strcpy( fullseq[iinf0[i]]+block->start, tmpseq0 );
	}
	for( i=0; i<n1; i++ ) 
	{
//		eq2dash( tmpaln1[i] );
		strncpy( tmpseq0, tmpaln1[i], newlen );
		strcpy( tmpseq0+newlen, fullseq[iinf1[i]] + posinold );
		strcpy( fullseq[iinf1[i]]+block->start, tmpseq0 );
	}
	FreeCharMtx( tmpaln0 );
	FreeCharMtx( tmpaln1 );
	FreeCharMtx( tmpseq );
	for( i=0; i<n0; i++ ) 
	{
//		for( j=0; j<njob; j++ ) free( mem[i][j] );
		free( mem[i] );
	}
	free( mem );
	free( nmem );
	free( iinf0 );
	free( iinf1 );
	free( group );
	free( g2n );
	free( prof1 );
	free( prof2 );
	free( nearest );
	FreeIntCub( topolpick );
	free( pickhistory );
	free( tmpint );

	return( zure );
}


#if 0
static int dorealignment( Blocktorealign *block, char **fullseq, int alloclen, int fullseqlen, int norg )
{
	int i, posinnew, posinold, newlen;
	int n0, n1;
	int zure;
	static int *iinf0, *iinf1;
	static char **tmpaln0 = NULL;
	static char **tmpaln1 = NULL;
	static char **tmpseq;
	char *opt, *npt;

	if( tmpaln0 == NULL )
	{
		tmpaln0 = AllocateCharMtx( njob, alloclen );
		tmpaln1 = AllocateCharMtx( njob, alloclen );
		tmpseq = AllocateCharMtx( 1, fullseqlen );
		iinf0 = AllocateIntVec( njob );
		iinf1 = AllocateIntVec( njob );
	}
	posinold = block->start;


	n0 = 0;
	n1 = 0;
	for( i=0; i<njob; i++ )
	{
		strncpy( tmpseq[0], fullseq[i] + block->start, block->end - block->start + 1 );
		tmpseq[0][block->end - block->start + 1] = 0;
		commongappick( 1, tmpseq );
//		if( strlen( tmpseq[0] ) > 0 )
		if( tmpseq[0][0] != 0 )
		{
			if( i < norg )
			{
				fprintf( stderr, "BUG!!!!\n" );
				exit( 1 );
			}
			strcpy( tmpaln0[n0], tmpseq[0] );
			iinf0[n0] = i;
			n0++;
		}
		else
		{
			strcpy( tmpaln1[n0], "" );
			iinf1[n1] = i;
			n1++;
		}
	}


	for( i=1; i<n0; i++ )
	{
		profilealignment( 1, n1, i, tmpaln0+i, tmpaln1, tmpaln0, alloclen, alg ); // n1 ha allgap
	}

#if 1
	fprintf( stderr, "After realignment:\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "%s\n", tmpaln0[i] );
//	for( i=0; i<n1; i++ ) fprintf( stderr, "%s\n", tmpaln1[i] );
#endif

	newlen = strlen( tmpaln0[0] );
	for( i=0; i<n0; i++ ) strncpy( fullseq[iinf0[i]]+block->start, tmpaln0[i], newlen );
	for( i=0; i<n1; i++ ) 
	{
		eq2dash( tmpaln1[i] );
		strncpy( fullseq[iinf1[i]] + block->start, tmpaln1[i], newlen );
	}

	posinold = block->end+1;
	posinnew = block->start + newlen;


	zure = ( block->end - block->start + 1 - strlen( tmpaln0[0] ) );

	for( i=0; i<njob; i++ ) 
	{
#if 0
		strcpy( fullseq[i]+posinnew, fullseq[i]+posinold ); // ??
#else
		opt = fullseq[i] + posinold;
		npt = fullseq[i] + posinnew;
		while( ( *npt++ = *opt++ ) );
		*npt = 0;
#endif
	}

	return( zure );
}
#endif

static void adjustposmap( int len, int *posmap, int *gaplist )
{
	int *newposmap;
	int *mpt1, *mpt2;
	int lenbk, zure;
	newposmap = calloc( len+2, sizeof( int ) );
	lenbk = len;
	zure = 0;
	mpt1 = newposmap;
	mpt2 = posmap;

#if 0
	int i;
	fprintf( stderr, "posmapa = " );
	for( i=0; i<len+2; i++ )
	{
		fprintf( stderr, "%3d ", posmap[i] );
	}
	fprintf( stderr, "\n" );
#endif

	while( len-- ) 
	{
		zure += *gaplist++;
		*mpt1++ = *mpt2++ + zure;
	}
	zure += *gaplist++;
	*mpt1 = *mpt2 + zure;

	mpt1 = newposmap;
	mpt2 = posmap;
	len = lenbk;
	while( len-- ) *mpt2++ = *mpt1++;
	*mpt2 = *mpt1;
	free( newposmap );
#if 0
	fprintf( stderr, "posmapm = " );
	for( i=0; i<lenbk+2; i++ )
	{
		fprintf( stderr, "%3d ", posmap[i] );
	}
	fprintf( stderr, "\n" );
#endif
}

static int insertgapsbyotherfragments_compact( int len, char *a, char *s, int *l, int *p )
{
	int gaplen;
	int i, pos, pi;
	int prevp = -1;
	int realignment = 0;
//	fprintf( stderr, "### insertgapsbyotherfragments\n" );
	for( i=0; i<len; i++ )
	{
		gaplen = l[i];
		pi = p[i];
		pos = prevp + 1;
//		fprintf( stderr, "gaplen = %d\n", gaplen );
		while( gaplen-- )
		{
			pos++;
			*a++ = *s++;
		}
//		fprintf( stderr, "pos = %d, pi = %d\n", pos, pi );
		while( pos++ < pi )
		{
			*a++ = '=';
			realignment = 1;
		}
		*a++ = *s++;
		prevp = pi;
	}
	gaplen = l[i];
	pi = p[i];
	pos = prevp + 1;
	while( gaplen-- )
	{
		pos++;
		*a++ = *s++;
	}
	while( pos++ < pi )
	{
		*a++ = '=';
		realignment = 1;
	}
	*a = 0;
	return( realignment );
}

void makegaplistcompact( int len, int *p, int *c, int *l )
{
	int i;
	int pg;
	int prep = -1;
	for( i=0; i<len+2; i++ )
	{
		if( ( pg = p[i]-prep-1) > 0 && l[i] > 0 )
		{
			if( pg < l[i] ) 
			{
				c[i] = l[i] - pg;
			}
			else
			{
				c[i] = 0;
			}
		}
		else
		{
			c[i] = l[i];
		}
		prep = p[i];
	}
}


void gaplist2alnx( int len, char *a, char *s, int *l, int *p, int lenlimit )
{
	int gaplen;
	int pos, pi, posl;
	int prevp = -1;
	int reslen = 0;
	char *sp;
//	char *abk = a;

#if 0
	int i;
	char *abk = a;
	fprintf( stderr, "s = %s\n", s );
	fprintf( stderr, "posmap  = " );
	for( i=0; i<len+2; i++ )
	{
		fprintf( stderr, "%3d ", p[i] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "gaplist = " );
	for( i=0; i<len+2; i++ )
	{
		fprintf( stderr, "%3d ", l[i] );
	}
	fprintf( stderr, "\n" );
#endif
	while( len-- )
	{
		gaplen = *l++;
		pi = *p++;

		if( (reslen+=gaplen) > lenlimit )
		{
			fprintf( stderr, "Length over. Please recompile!\n" );
			exit( 1 );
		}
		while( gaplen-- ) *a++ = '-';

		pos = prevp + 1;
		sp = s + pos;
		if( ( posl = pi - pos ) )
		{
			if( ( reslen += posl ) > lenlimit )
			{
				fprintf( stderr, "Length over. Please recompile\n" );
				exit( 1 );
			}
			while( posl-- ) *a++ = *sp++;
		}

		if( reslen++ > lenlimit )
		{
			fprintf( stderr, "Length over. Please recompile\n" );
			exit( 1 );
		}
		*a++ = *sp;
		prevp = pi;
	}

	gaplen = *l;
	pi = *p;
	if( (reslen+=gaplen) > lenlimit )
	{
		fprintf( stderr, "Length over. Please recompile\n" );
		exit( 1 );
	}
	while( gaplen-- ) *a++ = '-';

	pos = prevp + 1;
	sp = s + pos;
	if( ( posl = pi - pos ) )
	{
		if( ( reslen += posl ) > lenlimit )
		{
			fprintf( stderr, "Length over. Please recompile\n" );
			exit( 1 );
		}
		while( posl-- ) *a++ = *sp++;
	}
	*a = 0;
//	fprintf( stderr, "reslen = %d, strlen(a) = %d\n", reslen, strlen( abk ) );
//	fprintf( stderr, "a = %s\n", abk );
}

static void makenewgaplist( int *l, char *a )
{
	while( 1 )
	{
		while( *a == '=' )
		{
			a++;
			(*l)++;
//			fprintf( stderr, "a[] (i) = %s, *l=%d\n", a, *(l) );
		}
		*++l = 0;
		if( *a == 0 ) break;
		a++;
	}
	*l = -1;
}


void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	codonpos = 0;
	codonscore = 0;
	outnumber = 0;
	scoreout = 0;
	treein = 0;
	topin = 0;
	rnaprediction = 'm';
	rnakozo = 0;
	nevermemsave = 0;
	inputfile = NULL;
	addfile = NULL;
	addprofile = 1;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0; // chuui
	force_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'X';
	sueff_global = 0.1;
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNAppenalty_ex = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	TMorJTT = JTT;
	consweight_multi = 1.0;
	consweight_rna = 0.0;
	nadd = 0;
	multidist = 0;
	tuplesize = -1;
	legacygapcost = 0;
	allowlongadds = 0;
	addtotop = 0;
	keeplength = 0;
	diffout = 0;
	mapout = 0;
	smoothing = 0;
	distout = 0;
	hitout = 0.0;
	nwildcard = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
                    goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'e':
					RNApthr = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'o':
					RNAppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					--argc;
					goto nextoption;
				case 'r':
					consweight_rna = atof( *++argv );
					rnakozo = 1;
					--argc;
					goto nextoption;
				case 'c':
					consweight_multi = atof( *++argv );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 'R':
					codonpos = 1;
					break;
				case 'S':
					codonscore = 1;
					break;
#if 0
				case 'R':
					rnaprediction = 'r';
					break;
				case 's':
					RNAscoremtx = 'r';
					break;
#endif
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
				case 'y':
					distout = 1;
					break;
				case '^':
					hitout = atof( *++argv );
					--argc;
					goto nextoption;
				case 't':
					treeout = 1;
					break;
				case 'T':
					noalign = 1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
#if 1
				case 'O':
					outgap = 0;
					break;
#else
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
//				case 'S':
//					scoreout = 1;
//					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
				case 's':
					treemethod = 's';
					break;
#endif
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
				case 'E':
					treemethod = 'E';
					break;
				case 'q':
					treemethod = 'q';
					break;
				case 'n' :
					outnumber = 1;
					break;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'Q':
					alg = 'Q';
					break;
#endif
				case 'H':
					alg = 'H';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'B':  // hitsuyou! memopt -M -B no tame
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					force_fft = 1;
					use_fft = 1;
					break;
				case 'U':
					treein = 1;
					break;
				case 's':
					addtotop = 2;
					break;
				case 'x':
					addtotop = 1;
					break;
				case 'V':
					allowlongadds = 1;
					break;
				case 'p':
					smoothing = 1;
					break;
#if 0
				case 'V':
					topin = 1;
					break;
#endif
				case 'u':
					tbrweight = 0;
					weight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					multidist = 1;
					break;
				case 'W':
					tuplesize = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'd':
					disp = 1;
					break;
#endif
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
#if 0
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
#endif
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'Z':
					checkC = 1;
					break;
#endif
				case 'L':
					legacygapcost = 1;
					break;
				case 'Y':
					keeplength = 1;
					break;
				case 'z':
					mapout = 2;
					break;
				case 'Z':
					mapout = 1;
					break;
				case '%':
					diffout = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
	if( alg == 'C' && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : C, o\n" );
		exit( 1 );
	}
}


static double treebase( int nseq, int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, double *effarr, int *alloclen, LocalHom **localhomtable, RNApair ***singlerna, double *effarr_kozo )
{

	int i, l, m;
	int len1nocommongap, len2nocommongap;
	int len1, len2;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	double *effarr1_kozo = NULL;
	double *effarr2_kozo = NULL;
	LocalHom ***localhomshrink = NULL;
	int *fftlog;
	int m1, m2;
	int *gaplen;
	int *gapmap;
	int *alreadyaligned;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	RNApair ***grouprna1, ***grouprna2;

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( nseq, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( nseq, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	fftlog = AllocateIntVec( nseq );
	effarr1 = AllocateDoubleVec( nseq );
	effarr2 = AllocateDoubleVec( nseq );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	alreadyaligned = AllocateIntVec( nseq );
	if( constraint )
	{
		localhomshrink = (LocalHom ***)calloc( nseq, sizeof( LocalHom ** ) );
#if SMALLMEMORY
		if( multidist )
		{
			for( i=0; i<nseq; i++) localhomshrink[i] = (LocalHom **)calloc( 1, sizeof( LocalHom *) );
		}
		else
#endif
		{
			for( i=0; i<nseq; i++) localhomshrink[i] = (LocalHom **)calloc( nseq, sizeof( LocalHom *) );
		}
	}
	effarr1_kozo = AllocateDoubleVec( nseq ); //tsuneni allocate sareru.
	effarr2_kozo = AllocateDoubleVec( nseq ); //tsuneni allocate sareru.
	for( i=0; i<nseq; i++ ) effarr1_kozo[i] = 0.0;
	for( i=0; i<nseq; i++ ) effarr2_kozo[i] = 0.0;

	gaplen = AllocateIntVec( *alloclen+10 ); // maikai shokika
	gapmap = AllocateIntVec( *alloclen+10 ); // maikai shokika
	for( i=0; i<nseq-1; i++ ) alreadyaligned[i] = 1;
	alreadyaligned[nseq-1] = 0;

	for( l=0; l<nseq; l++ ) fftlog[l] = 1;


	if( constraint ) 
	{
#if SMALLMEMORY
		if( multidist )
			dontcalcimportance_firstone( nseq, effarr, aseq, localhomtable );
		else
			calcimportance( nseq, effarr, aseq, localhomtable );
#else
		calcimportance( nseq, effarr, aseq, localhomtable );
#endif
	}

	tscore = 0.0;
	for( l=0; l<nseq-1; l++ )
	{
		if( mergeoralign[l] == 'n' )
		{
//			fprintf( stderr, "SKIP!\n" );
#if 0
			free( topol[l][0] );
			free( topol[l][1] );
			free( topol[l] );
#endif
			continue;
		}

		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
        len1 = strlen( aseq[m1] );
        len2 = strlen( aseq[m2] );
        if( *alloclen < len1 + len2 )
        {
#if 0
			fprintf( stderr, "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, nseq, *alloclen + 10  );
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				fprintf( stderr, "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				fprintf( stderr, "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			fprintf( stderr, "done. *alloclen = %d\n", *alloclen );
#else
			fprintf( stderr, "Length over!\n" );
			exit( 1 );
#endif
		}

		if( effarr_kozo )
		{
			clus1 = fastconjuction_noname_kozo( topol[l][0], aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
			clus2 = fastconjuction_noname_kozo( topol[l][1], aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
		}
		else
		{
			clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0, NULL );
			clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0, NULL );
		}

		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' )
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";


		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap );
			commongappick( clus2, mseq2 );
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap );
			commongappick( clus1, mseq1 );
			len1nocommongap = strlen( mseq1[0] );
		}
		

//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );
//
#if 1
//		fprintf( stderr, "\rSTEP % 5d /%d ", l+1, nseq-1 );
//		fflush( stderr );
#else
		fprintf( stdout, "STEP %d /%d\n", l+1, nseq-1 );
		fprintf( stderr, "STEP %d /%d\n", l+1, nseq-1 );
		fprintf( stderr, "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
#endif



//		for( i=0; i<clus1; i++ ) fprintf( stderr, "## STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );

		if( constraint )
		{
#if SMALLMEMORY
			if( multidist )
			{
				fastshrinklocalhom_one( topol[l][0], topol[l][1], nseq-1, localhomtable, localhomshrink );
			}
			else
#endif
			{
				fastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
			}

//			msfastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			fprintf( stdout, "localhomshrink =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
//			weightimportance4( clus1, clus2, effarr1, effarr2, localhomshrink );
//			fprintf( stderr, "after weight =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
		}
		if( rnakozo && rnaprediction == 'm' )
		{
			makegrouprna( grouprna1, singlerna, topol[l][0] );
			makegrouprna( grouprna2, singlerna, topol[l][1] );
		}


/*
		fprintf( stderr, "before align all\n" );
		display( aseq, nseq );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/


		if( !nevermemsave && ( constraint != 2  && alg != 'M'  && ( len1 > 50000 || len2 > 50000 ) ) )
		{
			fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL; // 2013/Jul17
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}


		if( alg == 'M' ) // hoka no thread ga M ni shitakamo shirenainode
		{
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}


//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000 ); // v6.708
//		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
//		fprintf( stderr, "f=%d, clus1=%d, fftlog[m1]=%d, clus2=%d, fftlog[m2]=%d\n", ffttry, clus1, fftlog[m1], clus2, fftlog[m2] );

		if( constraint == 2 )
		{
			if( codonpos )
			{
				reporterr( "\n\nThe --codonpos option is supported only for --6merpair --addfragments at this point. Add the --6merpair flag, for now.\n\n\n" );
				exit( 1 );
			}
			if( alg == 'M' )
			{
				fprintf( stderr, "\n\nMemory saving mode is not supported.\n\n" );
				exit( 1 );
			}
			fprintf( stderr, "c" );
			if( alg == 'A' )
			{
				imp_match_init_strict( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, NULL, 1, topol[l][0], topol[l][1], NULL, NULL, NULL, -1, 0 );
				if( rnakozo ) imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = A__align( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, constraint, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1, NULL, NULL, NULL, 0.0, 0.0 ); // cpmxchild0 tsukaeru??
			}
			else if( alg == 'Q' )
			{
				fprintf( stderr, "Q has been disabled.\n" );
				exit( 1 );
			}
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			fprintf( stderr, "f" );
			if( alg == 'M' )
			{
				fprintf( stderr, "m" );
				pscore = Falign_udpari_long( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
				pscore = Falign( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
		}
		else
		{
			fprintf( stderr, "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					fprintf( stderr, "m" );
					pscore = MSalignmm( n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, NULL, NULL, NULL, 0.0, 0.0 );
					break;
				case( 'A' ):
					pscore = A__align( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1, NULL, NULL, NULL, 0.0, 0.0 ); // cpmxchild0 tsukaeru??
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}


		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		fprintf( stderr, "aseq[last] = %s\n", aseq[nseq-1] );

#if SCOREOUT
		fprintf( stderr, "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
#if 0 // New gaps = '='
		fprintf( stderr, "Original msa\n" );
		for( i=0; i<clus1; i++ )
			fprintf( stderr, "%s\n", mseq1[i] );
		fprintf( stderr, "Query\n" );
		for( i=0; i<clus2; i++ )
			fprintf( stderr, "%s\n", mseq2[i] );
#endif

//		writePre( nseq, name, nlen, aseq, 0 );

		if( disp ) display( aseq, nseq );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
//			if( deleteadditionalinsertions ) ndeleted += deletenewinsertions( clus2, clus1, mseq2, mseq1, deleterecord );
			adjustgapmap( strlen( mseq2[0] )-len2nocommongap+len2, gapmap, mseq2[0] );
			restorecommongaps( nseq, nseq-(clus1+clus2), aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			findnewgaps( clus2, 0, mseq2, gaplen );
			insertnewgaps( nseq, alreadyaligned, aseq, topol[l][1], topol[l][0], gaplen, gapmap, *alloclen, alg, '-' );
//			for( i=0; i<nseq; i++ ) eq2dash( aseq[i] );
			for( i=0; (m=topol[l][0][i])>-1; i++ ) alreadyaligned[m] = 1;
		}
		if( mergeoralign[l] == '2' )
		{
//			if( deleteadditionalinsertions ) ndeleted += deletenewinsertions( clus1, clus2, mseq1, mseq2, deleterecord );
			adjustgapmap( strlen( mseq1[0] )-len1nocommongap+len1, gapmap, mseq1[0] );
			restorecommongaps( nseq, nseq-(clus1+clus2), aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			findnewgaps( clus1, 0, mseq1, gaplen );
			insertnewgaps( nseq, alreadyaligned, aseq, topol[l][0], topol[l][1], gaplen, gapmap, *alloclen, alg, '-' );
//			for( i=0; i<nseq; i++ ) eq2dash( aseq[i] );
			for( i=0; (m=topol[l][1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

#if 0
		free( topol[l][0] );
		free( topol[l][1] );
		free( topol[l] );
#endif
	}

//for( i=0; i<nseq; nseq++ )
//reporterr( "In treebase() before deletenewinsertions, %s\n", aseq[i] );
//	if( keeplength ) ndeleted += deletenewinsertions_withoutusingequal( nseq-1, 1, 0, aseq, aseq+nseq-1, NULL, deletemapiadd, deletelagiadd, deletelistiadd );

#if SCOREOUT
	fprintf( stderr, "totalscore = %10.2f\n\n", tscore );
#endif
	free( gaplen );
	free( gapmap );
	if( rnakozo && rnaprediction == 'm' )
	{
		free( grouprna1 );
		free( grouprna2 );
	}
	free( fftlog ); // iranai
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( alreadyaligned );
	if( constraint )
	{
		for( i=0; i<nseq; i++ ) free( localhomshrink[i] ); // ??
		free( localhomshrink );
	}
	free( effarr1_kozo );
	free( effarr2_kozo );


	return( pscore );
}




static void mtxcpy( int norg, int njobc, double ***iscorec, double **iscore )
{
	int i, nlim, n;
	double *fpt, *fptc;
	
	*iscorec = AllocateFloatHalfMtx( njobc );
	nlim = norg-1;
	for( i=0; i<nlim; i++ )
	{
		fptc = (*iscorec)[i]+1;
		fpt  = iscore[i]+1;
		n = norg-i-1;
		while( n-- )
			*fptc++ = *fpt++;
//		for( j=i+1; j<norg; j++ )	
//			(*iscorec)[i][j-i] = iscore[i][j-i];
	}
}


static void	*addsinglethread( void *arg )
	{
		thread_arg_t *targ = (thread_arg_t *)arg;
		int *nlenc = NULL;
		char **namec = NULL;
		Treedep *depc = NULL;
		char **mseq1 = NULL, **mseq2 = NULL;
		double **iscorec;
//		double **iscorecbk; // to speedup
		double *effc = NULL;
		int ***topolc = NULL;
		double **lenc = NULL;
		LocalHom **localhomtablec = NULL;
		int *memlist0 = NULL;
		int *memlist1 = NULL;
		int *addmem = NULL;
		int njobc, norg;
		char **bseq = NULL;
		int i, j, k, m, iadd, rep, neighbor;
		char *mergeoralign = NULL;
		int *nogaplenjusttodecideaddhereornot = NULL;
		char *tmpseq = NULL;

#ifdef enablemultithread
		int thread_no = targ->thread_no;
		int *iaddshare = targ->iaddshare; 
#endif
		int njob = targ->njob; 
		int *follows = targ->follows;
		int nadd = targ->nadd;
		int *nlen = targ->nlen;
		char **name = targ->name;
		char **seq = targ->seq;
		LocalHom **localhomtable = targ->localhomtable;
		double **iscore = targ->iscore;
		double **nscore = targ->nscore;
		int *istherenewgap = targ->istherenewgap;
		int **newgaplist = targ->newgaplist;
		RNApair ***singlerna = targ->singlerna;
		double *eff_kozo_mapped = targ->eff_kozo_mapped;
		int alloclen = targ->alloclen;
		Treedep *dep = targ->dep;
		int ***topol = targ->topol;
		double **len = targ->len;
		Addtree *addtree = targ->addtree;
		GapPos **deletelist = targ->deletelist;
		GapPos **difflist = targ->difflist;
		double pscore;
		int *alnleninnode = NULL;
		char *targetseq;



//		fprintf( stderr, "\nPreparing thread %d\n", thread_no );
		norg = njob - nadd;
		njobc = norg+1;

#if 0
		alnleninnode = AllocateIntVec( norg );
		addmem = AllocateIntVec( nadd+1 );
		depc = (Treedep *)calloc( njobc, sizeof( Treedep ) );
		mseq1 = AllocateCharMtx( njob, 0 );
		mseq2 = AllocateCharMtx( njob, 0 );
		bseq = AllocateCharMtx( njobc, alloclen );
		namec = AllocateCharMtx( njob, 0 );
		nlenc = AllocateIntVec( njob );
		mergeoralign = AllocateCharVec( njob );
		nogaplenjusttodecideaddhereornot = AllocateIntVec( njobc );
		tmpseq = calloc( alloclen, sizeof( char ) );
#else
		alnleninnode = AllocateIntVec( norg );
		addmem = AllocateIntVec( nadd+1 );
		depc = (Treedep *)calloc( njobc, sizeof( Treedep ) );
		mseq1 = AllocateCharMtx( njobc, 0 );
		mseq2 = AllocateCharMtx( njobc, 0 );
		bseq = AllocateCharMtx( njobc, alloclen );
		namec = AllocateCharMtx( njobc, 0 );
		nlenc = AllocateIntVec( njobc );
		mergeoralign = AllocateCharVec( njobc );
		nogaplenjusttodecideaddhereornot = AllocateIntVec( njobc );
		tmpseq = calloc( alloclen, sizeof( char ) );
#endif

		if( allowlongadds ) // hontou ha iranai.
		{
			for( i=0; i<njobc; i++ ) nogaplenjusttodecideaddhereornot[i] = 0;
		}
		else
		{
			for( i=0; i<norg; i++ ) 
			{
				gappick0( tmpseq, seq[i] );
				nogaplenjusttodecideaddhereornot[i] = strlen( tmpseq );
			}
		}

		for( i=0; i<norg; i++ ) strcpy( bseq[i], seq[i] );
		if( norg == 1 )
		{
			alnleninnode[0] = strlen( bseq[0] );
		}
		else
		{
			for( i=norg-2; i>=0; i-- ) 
//			for( i=norg-2; i; i-- )  // BUG!!!!
			{
//				reporterr( "\nstep %d\n", i );
				k = 0;
				for( j=0; (m=topol[i][0][j])!=-1; j++ ) 
				{
					mseq1[k++] = bseq[m];
//					reporterr( "%d ", m );
				}
				for( j=0; (m=topol[i][1][j])!=-1; j++ ) 
				{
					mseq1[k++] = bseq[m];
//					reporterr( "%d ", m );
				}
//				reporterr( "\n" );
				commongappick( k, mseq1 );
				alnleninnode[i] = strlen( mseq1[0] );
//				fprintf( stderr, "alnleninnode[%d] = %d\n", i, alnleninnode[i] );
			}
		}
//		for( i=0; i<norg-1; i++ )
//			fprintf( stderr, "alnleninnode[%d] = %d\n", i, alnleninnode[i] );


		if( constraint )
		{
			localhomtablec = (LocalHom **)calloc( njobc, sizeof( LocalHom *) ); // motto chiisaku dekiru.
#if SMALLMEMORY
			if( multidist )
			{
				for( i=0; i<njobc; i++) localhomtablec[i] = (LocalHom *)calloc( 1, sizeof( LocalHom ) ); // motto chiisaku dekiru.
			}
			else
#endif
			{
				for( i=0; i<njobc; i++) localhomtablec[i] = (LocalHom *)calloc( njobc, sizeof( LocalHom ) ); // motto chiisaku dekiru.
				for( i=0; i<norg; i++ ) for( j=0; j<norg; j++ ) localhomtablec[i][j] = localhomtable[i][j]; // iru!
			}
		}


		topolc = AllocateIntCub( njobc, 2, 0 );
		lenc = AllocateFloatMtx( njobc, 2 );
		effc = AllocateDoubleVec( njobc );
//		for( i=0; i<norg; i++ ) nlenc[i] = strlen( seq[i] );
		for( i=0; i<norg; i++ ) nlenc[i] = nlen[i];
		for( i=0; i<norg; i++ ) namec[i] = name[i];
		memlist0 = AllocateIntVec( norg+1 );
		memlist1 = AllocateIntVec( 2 );
		for( i=0; i<norg; i++ ) memlist0[i] = i;
		memlist0[norg] = -1;

//		fprintf( stderr, "\ndone. %d\n", thread_no );

//		mtxcpy( norg, norg, &iscorecbk, iscore ); // to speedup?


		iadd = -1;
		while( 1 )
		{
#ifdef enablemultithread
			if( nthread )
			{
				pthread_mutex_lock( targ->mutex_counter );
				iadd = *iaddshare;
				if( iadd == nadd )
				{
					pthread_mutex_unlock( targ->mutex_counter );
					break;
				}
				if( iadd < 500 )
					fprintf( stderr, "\rSTEP %d / %d (thread %d)                    \r", iadd, nadd, thread_no );
				else if( iadd % 100 == 0 )
					fprintf( stderr, "\nSTEP %d / %d (thread %d)                    \n", iadd, nadd, thread_no );
				++(*iaddshare);
				targetseq = seq[norg+iadd];
				pthread_mutex_unlock( targ->mutex_counter );
			}
			else
#endif
			{
				iadd++;
				if( iadd == nadd ) break;
				targetseq = seq[norg+iadd];
				if( iadd < 500 )
					fprintf( stderr, "\rSTEP %d / %d                    \r", iadd, nadd );
				else if( iadd % 100 == 0 )
					fprintf( stderr, "\nSTEP %d / %d                    \n", iadd, nadd );
			}

			for( i=0; i<norg; i++ ) strcpy( bseq[i], seq[i] );
//			gappick0( bseq[norg], seq[norg+iadd] );
			gappick0( bseq[norg], targetseq );

			if( allowlongadds ) // missed in v7.220
				nogaplenjusttodecideaddhereornot[norg] = 0;
			else
				nogaplenjusttodecideaddhereornot[norg] = strlen( bseq[norg] );

			mtxcpy( norg, njobc, &iscorec, iscore );
	
			if( multidist || tuplesize > 0 )
			{
				for( i=0; i<norg; i++ ) iscorec[i][norg-i] = nscore[i][iadd];
			}
			else
			{
				for( i=0; i<norg; i++ ) iscorec[i][norg-i] = iscore[i][norg+iadd-i];
			}


#if 0
			for( i=0; i<njobc-1; i++ )
			{
				fprintf( stderr, "i=%d\n", i );
				for( j=i+1; j<njobc; j++ ) 
				{
					fprintf( stderr, "%d-%d, %f\n", i, j, iscorec[i][j-i] );
				}
			}
#endif
			nlenc[norg] = nlen[norg+iadd];
			namec[norg] = name[norg+iadd];
			if( constraint) 
			{
				for( i=0; i<norg; i++ )
				{
#if SMALLMEMORY
					if( multidist )
					{
						localhomtablec[i][0] = localhomtable[i][iadd];
//						localhomtablec[norg][i] = localhomtable[norg+iadd][i];
					}
					else
#endif
					{
						localhomtablec[i][norg] = localhomtable[i][norg+iadd];
						localhomtablec[norg][i] = localhomtable[norg+iadd][i];
					}
				}
//				localhomtablec[norg][norg] = localhomtable[norg+iadd][norg+iadd]; // iranai!!
			}
	
//			fprintf( stderr, "Constructing a UPGMA tree %d ... ", iadd );
//			fflush( stderr );
	

//			if( iadd == 0 )
//			{
//			}
//			fixed_musclesupg_double_realloc_nobk_halfmtx( njobc, iscorec, topolc, lenc, depc, 0, 1 );

			if( addtotop == 2 ) // root
				neighbor = addonetip2fixedpos( njobc, topolc, lenc, iscorec, topol, len, dep, treeout, addtree, iadd, name, alnleninnode, nogaplenjusttodecideaddhereornot, noalign, 'r' );
			else if( addtotop == 1 ) // top
				neighbor = addonetip2fixedpos( njobc, topolc, lenc, iscorec, topol, len, dep, treeout, addtree, iadd, name, alnleninnode, nogaplenjusttodecideaddhereornot, noalign, 't' );
			else
				neighbor = addonetip( njobc, topolc, lenc, iscorec, topol, len, dep, treeout, addtree, iadd, name, alnleninnode, nogaplenjusttodecideaddhereornot, noalign );


			FreeFloatHalfMtx( iscorec, njobc );

	
			if( tbrweight )
			{
				weight = 3; 
				counteff_simple_double_nostatic( njobc, topolc, lenc, effc );
			}
			else
			{
				for( i=0; i<njobc; i++ ) effc[i] = 1.0;
			}
		
//			FreeFloatMtx( lenc );

			if( noalign ) // nen no tame weight wo keisan.
			{
//				FreeFloatHalfMtx( iscorec, njobc ); // saki ni continue suru baai ha fukkatsu.
				continue;
			}

//			reporterr( "iadd = %d\n", iadd );
		
#if 0
			for( i=0; i<njobc-1; i++ )
			{
				fprintf( stderr, "\n step %d\n", i );
				fprintf( stderr, "topol[%d] = \n", i );
				for( j=0; topolc[i][0][j]!=-1; j++ ) fprintf( stderr, "%d ", topolc[i][0][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "len=%f\n", lenc [i][0] );
				for( j=0; topolc[i][1][j]!=-1; j++ ) fprintf( stderr, "%d ", topolc[i][1][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "len=%f\n", lenc [i][1] );
			}

			fprintf( stderr, "\nneighbor = %d, iadd = %d\n", neighbor, iadd );
#endif
			follows[iadd] = neighbor;

			for( i=0; i<njobc-1; i++ ) mergeoralign[i] = 'n';
			for( j=njobc-1; j<njobc; j++ )
			{
				addmem[0] = j;
				addmem[1] = -1;
				for( i=0; i<njobc-1; i++ )
				{
					if( samemembern( topolc[i][0], addmem, 1 ) ) // arieru
					{
//						fprintf( stderr, "HIT!\n" );
						if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
						else mergeoralign[i] = '1';
					}
					else if( samemembern( topolc[i][1], addmem, 1 ) )
					{
//						fprintf( stderr, "HIT!\n" );
						if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
						else mergeoralign[i] = '2';
					}
				}
			}
	
//			for( i=0; i<1; i++ ) addmem[i] = njobc-1+i;
			addmem[0] = njobc-1;
			addmem[1] = -1;
			for( i=0; i<njobc-1; i++ )
			{
				if( includemember( topolc[i][0], addmem ) && includemember( topolc[i][1], addmem ) )
				{
					mergeoralign[i] = 'w';
				}
				else if( includemember( topolc[i][0], addmem ) )
				{
					mergeoralign[i] = '1';
//					fprintf( stderr, "HIT 1! iadd=%d", iadd );
				}
				else if( includemember( topolc[i][1], addmem ) )
				{
					mergeoralign[i] = '2';
//					fprintf( stderr, "HIT 2! iadd=%d", iadd );
				}
			}
#if 0
			for( i=0; i<njob-1; i++ )
			{
				fprintf( stderr, "mem0 = " );
				for( j=0; topol[i][0][j]>-1; j++ )	fprintf( stderr, "%d ", topol[i][0][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "mem1 = " );
				for( j=0; topol[i][1][j]>-1; j++ )	fprintf( stderr, "%d ", topol[i][1][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "i=%d, mergeoralign[] = %c\n", i, mergeoralign[i] );
			}
#endif


#if 0
			for( i=0; i<norg; i++ ) fprintf( stderr, "seq[%d, iadd=%d] = \n%s\n", i, iadd, seq[i] );
			fprintf( stderr, "gapmapS (iadd=%d) = \n", iadd );
			for( i=0; i<lennocommongap; i++ ) fprintf( stderr, "%d\n", gapmapS[i] );
#endif


//			fprintf( stderr, "Progressive alignment ... \r" );
		
#if 0
			pthread_mutex_lock( targ->mutex_counter );
			fprintf( stdout, "\nmergeoralign (iadd=%d) = ", iadd );
			for( i=0; i<njobc-1; i++ ) fprintf( stdout, "%c", mergeoralign[i] );
			fprintf( stdout, "\n" );
			pthread_mutex_unlock( targ->mutex_counter );
#endif

			singlerna = NULL;
			pscore = treebase( njobc, nlenc, bseq, 1, mergeoralign, mseq1, mseq2, topolc, effc, &alloclen, localhomtablec, singlerna, eff_kozo_mapped );
#if 0
			pthread_mutex_lock( targ->mutex_counter );
//			fprintf( stdout, "res (iadd=%d) = %s, pscore=%f\n", iadd, bseq[norg], pscore );
//			fprintf( stdout, "effc (iadd=%d) = ", iadd );
//			for( i=0; i<njobc; i++ ) fprintf( stdout, "%f ", effc[i] );
//			fprintf( stdout, "\n" );
			pthread_mutex_unlock( targ->mutex_counter );
#endif
	
	
#if 0
			fprintf( trap_g, "done.\n" );
			fclose( trap_g );
#endif
//			fprintf( stdout, "\n>seq[%d, iadd=%d] = \n%s\n", norg+iadd, iadd, seq[norg+iadd] );
//			fprintf( stdout, "\n>bseq[%d, iadd=%d] = \n%s\n", norg, iadd, bseq[norg] );
	
//			strcpy( seq[norg+iadd], bseq[norg] );


			if( diffout )
			{
				reporterr( "Not yet written\n" );
				exit( 1 );
//				deletenewinsertions_difflist( norg, 1, bseq, bseq+norg, difflist+iadd );
//				strcpy( targetseq, bseq[norg] );
//				i = norg; // no new gap!!
			}
			else if( keeplength )
			{
//				reporterr( "deletelist = %p\n", deletelist );
//				reporterr( "deletelist+iadd = %p\n", deletelist+iadd );
				ndeleted += deletenewinsertions_whole_eq( norg, 1, bseq, bseq+norg, deletelist+iadd );
//				for( i=0; i<norg+1; i++ ) reporterr( ">\n%s\n", bseq[i] );
				strcpy( targetseq, bseq[norg] );
				i = norg; // no new gap!!
			}
			else
			{
				strcpy( targetseq, bseq[norg] );
				rep = -1;
				for( i=0; i<norg; i++ )
				{
//					fprintf( stderr, "Checking %d/%d\n", i, norg );
					if( strchr( bseq[i], '=' ) ) break;
				}
			}

			if( i == norg ) 
				istherenewgap[iadd] = 0;
			else
			{
				rep = i;
				istherenewgap[iadd] = 1;


				makenewgaplist( newgaplist[iadd], bseq[rep] );
//				for( i=0; newgaplist[iadd][i]!=-1; i++ ) fprintf( stderr, "%d: %d\n", i, newgaplist[iadd][i] );
			}
			eq2dash( targetseq );

		}


#if 1
		if( constraint && localhomtablec )
		{
			for( i=0; i<njobc; i++ ) free( localhomtablec[i] );
			free( localhomtablec );
			localhomtablec = NULL;
		}
		if( mergeoralign ) free( mergeoralign );  mergeoralign = NULL;
		if( nogaplenjusttodecideaddhereornot ) free( nogaplenjusttodecideaddhereornot ); nogaplenjusttodecideaddhereornot = NULL;
		if( alnleninnode ) free( alnleninnode ); alnleninnode = NULL;
		if( tmpseq ) free( tmpseq ); tmpseq = NULL;
		if( bseq ) FreeCharMtx( bseq ); bseq = NULL;
		if( namec ) free( namec ); namec = NULL;
		if( nlenc ) free( nlenc  ); nlenc = NULL;
		if( depc ) free( depc ); depc = NULL;
		if( topolc ) FreeIntCub( topolc ); topolc = NULL;
		if( lenc ) FreeFloatMtx( lenc ); lenc = NULL;
		if( effc ) FreeDoubleVec( effc ); effc = NULL;
		if( memlist0 ) free( memlist0 ); memlist0 = NULL;
		if( memlist1 ) free( memlist1 ); memlist1 = NULL;
		if( addmem ) free( addmem ); addmem = NULL;
		if( mseq1 ) free( mseq1 ); mseq1 = NULL;
		if( mseq2 ) free( mseq2 ); mseq2 = NULL;
		Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
		Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
		A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
		if( commonIP ) FreeIntMtx( commonIP );
		commonIP = NULL;
		commonAlloc1 = commonAlloc2 = 0;
#endif
//		FreeFloatHalfMtx( iscorecbk, norg );

		return( NULL );
	}

static int nunknown = 0;

void seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < tuplesize )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void makecompositiontable_p( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		table[point]++;
#if 0 // kakunin shinai
		if( (unsigned int)table[point]++ >= INT_MAX )
		{
			reporterr( "Overflow. table[point]=%d>INT_MAX(%d).\n", table[point], INT_MAX );
			exit( 1 );
		}
#endif
	}
}


void makepointtable_nuc_dectet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *262144;
	point += *n++ * 65536;
	point += *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ *262144;
		point *= 4;
		point += *n++;
		*pointt++ = point;

	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc_octet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 16384;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

#ifdef enablemultithread

void *dndprethread( void *arg )
{
	dndprethread_arg_t *targ = (dndprethread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	double *selfscore = targ->selfscore;
	double **mtx = targ->mtx;
	char **seq = targ->seq;
	Jobtable2d *jobpospt = targ->jobpospt;

	int i, j;
	double ssi, ssj, bunbo;
	double mtxv;

	if( njob == 1 ) return( NULL );
	
	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
//		fprintf( stderr, "\n i=%d, j=%d before check\n", i, j );
		if( j == njob )
		{
//			fprintf( stderr, "\n j = %d, i = %d, njob = %d\n", j, i, njob );
			fprintf( stderr, "%4d/%4d (thread %4d), dndpre\r", i+1, njob, thread_no );
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
//				fprintf( stderr, "\n i=%d, njob-1=%d\n", i, njob-1 );
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
		}
//		fprintf( stderr, "\n i=%d, j=%d after check\n", i, j );
		jobpospt->j = j;
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex );

		ssi = selfscore[i];
		ssj = selfscore[j];

		bunbo = MIN( ssi, ssj );
		if( bunbo == 0.0 )
			mtxv = maxdist;
		else
			mtxv = maxdist * ( 1.0 - (double)naivepairscore11( seq[i], seq[j], penalty * 10  ) / bunbo );
#if 1
		if( mtxv < 0.0 )
		{
			fprintf( stderr, "WARNING: distance %d-%d is strange, %f.\n", i, j, mtxv );
			mtxv = 0.0;
//			exit( 1 ); // 2016/Aug/3
		}
		if( mtxv > 9.9 )
		{
			fprintf( stderr, "WARNING: distance %d-%d is strange, %f.\n", i, j, mtxv );
			mtxv = 9.9;
//			exit( 1 ); // 2016/Aug/3
		}
#else // CHUUI!!!  2012/05/16
		if( mtxv > 2.0 )
		{
			mtxv = 2.0;
		}
		if( mtxv < 0.0 )
		{
			fprintf( stderr, "Distance %d-%d is strange, %f.\n", i, j, mtxv );
			exit( 1 );
		}
#endif
		mtx[i][j-i] = mtxv;
	}
}

static void *gaplist2alnxthread( void *arg )
{
	gaplist2alnxthread_arg_t *targ = (gaplist2alnxthread_arg_t *)arg;
//	int thread_no = targ->thread_no;
	int ncycle = targ->ncycle;
	char **seq = targ->seq;
	int *newgaplist = targ->newgaplist;
	int *posmap = targ->posmap;
	int *jobpospt = targ->jobpospt;
	int tmpseqlen = targ->tmpseqlen;
	int lenfull = targ->lenfull;
	char *tmpseq1;
	int i;

	tmpseq1 = AllocateCharVec( tmpseqlen );

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = *jobpospt;
		if( i == ncycle )
		{
			pthread_mutex_unlock( targ->mutex );
			free( tmpseq1 );
			return( NULL );
		}
		*jobpospt = i+1;
		pthread_mutex_unlock( targ->mutex );

 		gaplist2alnx( lenfull, tmpseq1, seq[i], newgaplist, posmap, tmpseqlen  );
//		fprintf( stderr, ">%s (iadd=%d)\n%s\n", name[i], iadd, tmpseq1 );
		strcpy( seq[i], tmpseq1 );
	}
}

static void *distancematrixthread( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int norg = targ->norg;
	int *jobpospt = targ->jobpospt;
	int **pointt = targ->pointt;
	double **imtx = targ->imtx;
	double **nmtx = targ->nmtx;
	double *selfscore = targ->selfscore;
	int *nogaplen = targ->nogaplen;

	double lenfac, bunbo, longer, shorter, mtxv;
	int *table1;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = *jobpospt;
		if( i == norg )
		{
			pthread_mutex_unlock( targ->mutex );
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		*jobpospt = i+1;
		pthread_mutex_unlock( targ->mutex );

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			fprintf( stderr, "\r% 5d / %d (thread %4d)", i+1, norg, thread_no );
		}
		makecompositiontable_p( table1, pointt[i] );
	
		for( j=i+1; j<njob; j++ ) 
		{
			mtxv = (double)commonsextet_p( table1, pointt[j] );
			if( nogaplen[i] > nogaplen[j] )
			{
				longer=(double)nogaplen[i];
				shorter=(double)nogaplen[j];
			}
			else
			{
				longer=(double)nogaplen[j];
				shorter=(double)nogaplen[i];
			}
			lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
			bunbo = MIN( selfscore[i], selfscore[j] );

			if( j < norg )
			{
				if( bunbo == 0.0 )
					imtx[i][j-i] = maxdist;
				else
					imtx[i][j-i] = maxdist * ( 1.0 - mtxv / bunbo ) * lenfac;

			}
			else
			{
				if( bunbo == 0.0 )
					nmtx[i][j-norg] = maxdist;
				else
					nmtx[i][j-norg] = maxdist * ( 1.0 - mtxv / bunbo ) * lenfac;
			}
		} 
		free( table1 );

//		for( j=i+1; j<norg; j++ ) 
//			imtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] );
//		for( j=norg; j<njob; j++ ) 
//			nmtx[i][j-norg] = (double)commonsextet_p( table1, pointt[j] );
//		free( table1 );
	}
}
#endif


void ktupledistancematrix( int nseq, int norg, int nlenmax, char **seq, char **name, double **imtx, double **nmtx )
{
	char *tmpseq;
	int *grpseq;
	int **pointt;
	int i, j;
	int *nogaplen;
	int *table1;
	double lenfac, bunbo, longer, shorter, mtxv;
	double *selfscore;
	selfscore = AllocateFloatVec( nseq );

	fprintf( stderr, "\n\nMaking a distance matrix ..\n" );
	fflush( stderr );

	tmpseq = AllocateCharVec( nlenmax+1 );
	grpseq = AllocateIntVec( nlenmax+1 );
	pointt = AllocateIntMtx( nseq, nlenmax+1 );
	nogaplen = AllocateIntVec( nseq ); 

	if( dorp == 'd' ) tsize = (int)pow( 4, tuplesize );
	else              tsize = (int)pow( 6, 6 );

	if( dorp == 'd' && tuplesize == 6 )
	{
		lenfaca = D6LENFACA;
		lenfacb = D6LENFACB;
		lenfacc = D6LENFACC;
		lenfacd = D6LENFACD;
	}
	else if( dorp == 'd' && tuplesize == 10 )
	{
		lenfaca = D10LENFACA;
		lenfacb = D10LENFACB;
		lenfacc = D10LENFACC;
		lenfacd = D10LENFACD;
	}
	else    
	{
		lenfaca = PLENFACA;
		lenfacb = PLENFACB;
		lenfacc = PLENFACC;
		lenfacd = PLENFACD;
	}

	maxl = 0;
	for( i=0; i<nseq; i++ ) 
	{
		gappick0( tmpseq, seq[i] );
		nogaplen[i] = strlen( tmpseq );
		if( nogaplen[i] < 6 )
		{
//			fprintf( stderr, "Seq %d, too short, %d characters\n", i+1, nogaplen[i] );
//			fprintf( stderr, "Please use mafft-ginsi, mafft-linsi or mafft-ginsi\n\n\n" );
//			exit( 1 );
		}
		if( nogaplen[i] > maxl ) maxl = nogaplen[i];
		if( dorp == 'd' ) /* nuc */
		{
			seq_grp_nuc( grpseq, tmpseq );
//			makepointtable_nuc( pointt[i], grpseq );
//			makepointtable_nuc_octet( pointt[i], grpseq );
			if( tuplesize == 10 )
				makepointtable_nuc_dectet( pointt[i], grpseq );
			else if( tuplesize == 6 )
				makepointtable_nuc( pointt[i], grpseq );
			else
			{
				fprintf( stderr, "tuplesize=%d: not supported\n", tuplesize );
				exit( 1 );
			}
		}
		else                 /* amino */
		{
			seq_grp( grpseq, tmpseq );
			makepointtable( pointt[i], grpseq );
		}

	}
	if( nunknown ) fprintf( stderr, "\nThere are %d ambiguous characters\n", nunknown );

	for( i=0; i<nseq; i++ ) // serial de jubun
	{
		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		makecompositiontable_p( table1, pointt[i] );

		selfscore[i] = (double)commonsextet_p( table1, pointt[i] );
		free( table1 );
	}

#ifdef enablemultithread
	if( nthread > 0 )
	{
		distancematrixthread_arg_t *targ; 
		int jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex;

		jobpos = 0; 
		targ = calloc( nthread, sizeof( distancematrixthread_arg_t ) ); 
		handle = calloc( nthread, sizeof( pthread_t ) ); 
		pthread_mutex_init( &mutex, NULL );

		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].njob = nseq;
			targ[i].norg = norg;
			targ[i].jobpospt = &jobpos;
			targ[i].pointt = pointt;
			targ[i].imtx = imtx;
			targ[i].nmtx = nmtx;
			targ[i].selfscore = selfscore;
			targ[i].nogaplen = nogaplen;
			targ[i].mutex = &mutex;

			pthread_create( handle+i, NULL, distancematrixthread, (void *)(targ+i) );
		}
	
		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex );
		free( handle );
		free( targ );

	}
	else
#endif
	{
		for( i=0; i<norg; i++ )
		{
			table1 = (int *)calloc( tsize, sizeof( int ) );
			if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
			if( i % 100 == 0 )
			{
				fprintf( stderr, "\r% 5d / %d", i+1, norg );
				fflush( stderr );
			}
			makecompositiontable_p( table1, pointt[i] );
	
			for( j=i+1; j<nseq; j++ ) 
			{
				mtxv = (double)commonsextet_p( table1, pointt[j] );
				if( nogaplen[i] > nogaplen[j] )
				{
					longer=(double)nogaplen[i];
					shorter=(double)nogaplen[j];
				}
				else
				{
					longer=(double)nogaplen[j];
					shorter=(double)nogaplen[i];
				}
				lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
				bunbo = MIN( selfscore[i], selfscore[j] );

				if( j < norg )
				{
					if( bunbo == 0.0 )
						imtx[i][j-i] = maxdist;
					else
						imtx[i][j-i] = maxdist * ( 1.0 - mtxv / bunbo ) * lenfac;
				}
				else
				{
					if( bunbo == 0.0 )
						nmtx[i][j-norg] = maxdist;
					else
						nmtx[i][j-norg] = maxdist * ( 1.0 - mtxv / bunbo ) * lenfac;
				}
			} 
			free( table1 );
		}
	}

	fprintf( stderr, "\ndone.\n\n" );
	fflush( stderr );

	free( grpseq );
	free( tmpseq );
	FreeIntMtx( pointt );
    free( nogaplen );
    free( selfscore );


	if( hitout<0.0 )
	{
		fprintf( stdout, "Threshold=%f\n\n", -hitout );
		for( i=0; i<norg; i++ )
		{
			for( j=norg; j<nseq; j++ ) 
			{
				if( nmtx[i][j-norg] < -hitout )
					break;
			}
			if( j<nseq ) 
			{
				fprintf( stdout, "%s may be similar to:\n", name[i]+1 );
				for( j=norg; j<nseq; j++ ) 
				{
					if( nmtx[i][j-norg] < -hitout )
						fprintf( stdout, "    %s, %f\n", name[j]+1, nmtx[i][j-norg] );
				}
				fprintf( stdout, "\n" );
			}
		}
		exit( 1 );
	}
	if( hitout>0.0 )
	{
		fprintf( stdout, "Threshold=%f\n\n", hitout );
		for( i=norg; i<nseq; i++ )
		{
			for( j=0; j<norg; j++ ) 
			{
				if( nmtx[j][i-norg] < hitout )
					break;
			}
			if( j<norg ) 
			{
				fprintf( stdout, "%s may be similar to:\n", name[i]+1 );
				for( j=0; j<norg; j++ ) 
				{
					if( nmtx[j][i-norg] < hitout )
						fprintf( stdout, "    %s, %f\n", name[j]+1, nmtx[j][i-norg] );
				}
				fprintf( stdout, "\n" );
			}
		}
		exit( 1 );
	}

#if 0 // writehat2 wo kakinaosu
	if( distout )
	{
		hat2p = fopen( "hat2", "w" );
		WriteFloatHat2_pointer_halfmtx( hat2p, nseq, name, mtx );
		fclose( hat2p );
	}
#endif
}

void dndpre( int nseq, char **seq, double **mtx ) // not used yet
{
	int i, j, ilim;
	double *selfscore;
	double mtxv;
	double ssi, ssj, bunbo;

	selfscore = AllocateFloatVec( nseq );

	for( i=0; i<nseq; i++ )
	{
		selfscore[i] = (double)naivepairscore11( seq[i], seq[i], 0 );
	}
#ifdef enablemultithread
	if( nthread > 0 )
	{
		dndprethread_arg_t *targ;
		Jobtable2d jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex;

		jobpos.i = 0;
		jobpos.j = 0;

		targ = calloc( nthread, sizeof( dndprethread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex, NULL );

		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].njob = nseq;
			targ[i].selfscore = selfscore;
			targ[i].mtx = mtx;
			targ[i].seq = seq;
			targ[i].jobpospt = &jobpos;
			targ[i].mutex = &mutex;

			pthread_create( handle+i, NULL, dndprethread, (void *)(targ+i) );
		}

		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex );

	}
	else
#endif
	{
		ilim = nseq-1;
		for( i=0; i<ilim; i++ )
		{
			ssi = selfscore[i];
			fprintf( stderr, "%4d/%4d\r", i+1, nseq );

			for( j=i+1; j<nseq; j++ )
			{
				ssj = selfscore[j];
				bunbo = MIN( ssi, ssj );
				if( bunbo == 0.0 )
					mtxv = maxdist;
				else
					mtxv = maxdist * ( 1.0 - (double)naivepairscore11( seq[i], seq[j], penalty * 10 ) / bunbo );

#if 1
				if( mtxv < 0.0 )
				{
					fprintf( stderr, "WARNING: distance %d-%d is strange, %f.\n", i, j, mtxv );
					mtxv = 0.0;
//					exit( 1 ); // 2016/Aug/3
				}
				if( mtxv > 9.9 )
				{
					fprintf( stderr, "WARNING: distance %d-%d is strange, %f.\n", i, j, mtxv );
					mtxv = 9.9;
//					exit( 1 ); // 2016/Aug/3
				}
#else // CHUUI!!!  2012/05/16
				if( mtxv > 2.0 )
				{
					mtxv = 2.0;
				}
				if( mtxv < 0.0 )
				{
					fprintf( stderr, "Distance %d-%d is strange, %f.\n", i, j, mtxv );
					exit( 1 );
				}
#endif
				mtx[i][j-i] = mtxv;
			}
		}
	}
	
#if TEST
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) 
		fprintf( stdout, "i=%d, j=%d, mtx[][] = %f\n", i, j, mtx[i][j] );
#endif
	free( selfscore );

}

static int searchlet( char *p1, char *p2 )
{
	char *p;
	for( p=p1; p<=p2; p++ )
	{
		if( *p != '-' )
		{
			return( p-p1 );
		}
	}
	return( -1 );
}

static void smoothing2( int njob, int nadd, int lenfull, char **seq, Blocktorealign *realign )
{
	int i, j, norg = njob-nadd;

	reporterr( "Smoothing2\n" );
	for( i=2; i<lenfull+1; i++ )
	{
		int postoshiftfrom;
		int shiftonadd;
		int postoshiftto;
		if( realign[i-2].nnewres && realign[i-1].nnewres == 0 && realign[i].nnewres ) 
		{
			postoshiftto = realign[i-2].start;
			postoshiftfrom = realign[i-2].end + 1;
			if( postoshiftfrom != realign[i].start -2 )
			{
				reporterr( "Unexpected pattern??? i=%d, realign[i-1].end=%d, realign[i].start=%d\n", i, realign[i-1].end, realign[i].start );
				exit( 1 );
			}

			realign[i].nnewres += realign[i-2].nnewres; 
			realign[i].start = realign[i-2].start+2; 
			realign[i-2].nnewres = 0;
			realign[i-2].start = realign[i-2].end = 0;


//			reporterr( "SHIFT %d -> %d\n", postoshiftfrom, postoshiftto );
			for( j=0; j<norg; j++ )
			{
				if( seq[j][postoshiftto] != '-' )
				{
					reporterr( "Unexpected pattern 2. i=%d, postoshiftto=%d, postoshiftfrom=%d, seq[%d][%d]=%c\n", i, postoshiftto, postoshiftfrom, j, postoshiftto, seq[j][postoshiftto] );
//					reporterr( seq[j] );
					exit( 1 );
				}
				seq[j][postoshiftto] = seq[j][postoshiftfrom];
				seq[j][postoshiftfrom] = '-';

				if( seq[j][postoshiftto+1] != '-' )
				{
					reporterr( "Unexpected pattern 3???\n" );
					exit( 1 );
				}
				seq[j][postoshiftto+1] = seq[j][postoshiftfrom+1];
				seq[j][postoshiftfrom+1] = '-';
			}
			for( j=norg; j<njob; j++ )
			{
				if( seq[j][postoshiftto] == '-' )
				{
					shiftonadd = searchlet( seq[j]+postoshiftto, seq[j]+postoshiftfrom );
					if( shiftonadd != -1 )
					{
						seq[j][postoshiftto] = seq[j][postoshiftto+shiftonadd];
						seq[j][postoshiftto+shiftonadd] = '-';
					}
				}
				if( seq[j][postoshiftto+1] == '-' )
				{
					shiftonadd = searchlet( seq[j]+postoshiftto+1, seq[j]+postoshiftfrom+1 );
					if( shiftonadd != -1 )
					{
						seq[j][postoshiftto+1] = seq[j][postoshiftto+1+shiftonadd];
						seq[j][postoshiftto+1+shiftonadd] = '-';
					}
				}
			}
		}
	}
//	for( i=0; i<lenfull+1; i++ ) fprintf( stderr, "i=%d, nnewres=%d, start=%d, end=%d\n", i, realign[i].nnewres, realign[i].start, realign[i].end );
}
static void smoothing1( int njob, int nadd, int lenfull, char **seq, Blocktorealign *realign )
{
	int i, j, norg = njob-nadd;

	reporterr( "Smoothing1\n" );
	for( i=1; i<lenfull+1; i++ )
	{
		int postoshiftfrom;
		int shiftonadd;
		int postoshiftto;
		if( realign[i-1].nnewres && realign[i].nnewres ) 
		{
			postoshiftto = realign[i-1].start;
			postoshiftfrom = realign[i-1].end + 1;
			if( postoshiftfrom != realign[i].start -1 )
			{
				reporterr( "Unexpected pattern??? i=%d, realign[i-1].end=%d, realign[i].start=%d\n", i, realign[i-1].end, realign[i].start );
				exit( 1 );
			}

			realign[i].nnewres += realign[i-1].nnewres; 
			realign[i].start = realign[i-1].start+1; 
			realign[i-1].nnewres = 0;
			realign[i-1].start = realign[i-1].end = 0;


//			reporterr( "SHIFT %d -> %d\n", postoshiftfrom, postoshiftto );
			for( j=0; j<norg; j++ )
			{
				if( seq[j][postoshiftto] != '-' )
				{
					reporterr( "Unexpected pattern 2???\n" );
					exit( 1 );
				}
				seq[j][postoshiftto] = seq[j][postoshiftfrom];
				seq[j][postoshiftfrom] = '-';
			}
			for( j=norg; j<njob; j++ )
			{
				if( seq[j][postoshiftto] == '-' )
				{
					shiftonadd = searchlet( seq[j]+postoshiftto, seq[j]+postoshiftfrom );
					if( shiftonadd != -1 )
					{
						seq[j][postoshiftto] = seq[j][postoshiftto+shiftonadd];
						seq[j][postoshiftto+shiftonadd] = '-';
					}
				}
			}
		}
	}
//	for( i=0; i<lenfull+1; i++ ) fprintf( stderr, "i=%d, nnewres=%d, start=%d, end=%d\n", i, realign[i].nnewres, realign[i].start, realign[i].end );
}
	
static double ratioambiguouscolumn_simple( char **s, int n, int l )
{
	int i, j;
	int nunambiguous;
	int nambiguous;
	if( dorp == 'd' )
		nunambiguous = 4;
	else
		nunambiguous = 6;

	nambiguous = 0;
	for( i=0; i<l; i++ )
	{
		for( j=0; j<n; j++ )
			if( amino_grp[(int)s[j][i]] < nunambiguous )
				break;
		if( j == n )
			nambiguous += 1;
	}

	reporterr( "nambiguous = %d / %d, ratio=%f\n", nambiguous, l, (double)nambiguous / l );
	return( (double)nambiguous / l );
}

int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static char **name, **seq;
	static char  **tmpseq;
	static char  *tmpseq1;
//	static char  *check1, *check2;
	static double **iscore, **iscore_kozo;
	static double *eff_kozo, *eff_kozo_mapped = NULL;
	int i, j, f, ien;
	int iadd;
	static int ***topol_kozo;
	Treedep *dep;
	static double **len_kozo;
	FILE *prep;
	FILE *infp;
	FILE *hat2p;
	int alignmentlength;
	char c;
	int alloclen, fullseqlen, tmplen;
	LocalHom **localhomtable = NULL;
	static char *kozoarivec;
	int nkozo;
	int njobc, norg, lenfull;
	int **newgaplist_o;
	int *newgaplist_compact;
	int **follower;
	int *follows;
	int *istherenewgap;
	int zure;
	int *posmap;
	int *ordertable;
	FILE *orderfp;
	int tmpseqlen;
	Blocktorealign *realign;
	RNApair ***singlerna;
	int ***topol;
	double **len;
	double **iscoreo, **nscore;
	FILE *fp;
	GapPos **deletelist = NULL;
	GapPos **difflist = NULL;
	char **addbk = NULL;
	char *originalgaps = NULL;
	Addtree *addtree;


	arguments( argc, argv );
#ifndef enablemultithread
	nthread = 0;
#endif


	if( fastathreshold < 0.0001 ) constraint = 0;

	if( inputfile )
	{
		infp = fopen( inputfile, "r" );
		if( !infp ) 
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else    
		infp = stdin;

	getnumlen( infp );
	rewind( infp );


	nkozo = 0;

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	norg = njob-nadd;
	njobc = norg+1;
	fprintf( stderr, "norg = %d\n", norg );
	fprintf( stderr, "njobc = %d\n", njobc );
//	if( norg > 1000 || nadd > 1000 ) use_fft = 0;
//	if( norg > 1000 ) use_fft = 0;


	fullseqlen = alloclen = nlenmax*4+1; //chuui!
	seq = AllocateCharMtx( njob, alloclen );

	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );

	ndeleted = 0;


	if( multidist || tuplesize > 0 )
	{
		iscore = AllocateFloatHalfMtx( norg );
		nscore = AllocateFloatMtx( norg, nadd );
	}
	else
	{
		iscore = AllocateFloatHalfMtx( njob );
		nscore = NULL;
	}

	kozoarivec = AllocateCharVec( njob );


	ordertable = AllocateIntVec( norg+1 );


	if( constraint )
	{
#if SMALLMEMORY
		if( multidist )
		{
			localhomtable = (LocalHom **)calloc( norg, sizeof( LocalHom *) );
			for( i=0; i<norg; i++)
			{
				localhomtable[i] = (LocalHom *)calloc( nadd, sizeof( LocalHom ) );
				for( j=0; j<nadd; j++)
				{
					localhomtable[i][j].start1 = -1;
					localhomtable[i][j].end1 = -1;
					localhomtable[i][j].start2 = -1;
					localhomtable[i][j].end2 = -1;
					localhomtable[i][j].overlapaa = -1.0;
					localhomtable[i][j].opt = -1.0;
					localhomtable[i][j].importance = -1.0;
					localhomtable[i][j].next = NULL;
					localhomtable[i][j].korh = 'h';
				}
			}
//			localhomtable = (LocalHom **)calloc( norg+nadd, sizeof( LocalHom *) );
//			for( i=norg; i<norg+nadd; i++) // hontou ha iranai
//			{
//				localhomtable[i] = (LocalHom *)calloc( norg, sizeof( LocalHom ) );
//				for( j=0; j<norg; j++)
//				{
//					localhomtable[i][j].start1 = -1;
//					localhomtable[i][j].end1 = -1;
//					localhomtable[i][j].start2 = -1;
//					localhomtable[i][j].end2 = -1;
//					localhomtable[i][j].overlapaa = -1.0;
//					localhomtable[i][j].opt = -1.0;
//					localhomtable[i][j].importance = -1.0;
//					localhomtable[i][j].next = NULL;
//					localhomtable[i][j].korh = 'h';
//				}
//			}
		}
		else
#endif
		{
			localhomtable = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
			for( i=0; i<njob; i++)
			{
				localhomtable[i] = (LocalHom *)calloc( njob, sizeof( LocalHom ) );
				for( j=0; j<njob; j++)
				{
					localhomtable[i][j].start1 = -1;
					localhomtable[i][j].end1 = -1;
					localhomtable[i][j].start2 = -1;
					localhomtable[i][j].end2 = -1;
					localhomtable[i][j].overlapaa = -1.0;
					localhomtable[i][j].opt = -1.0;
					localhomtable[i][j].importance = -1.0;
					localhomtable[i][j].next = NULL;
					localhomtable[i][j].korh = 'h';
				}
			}
		}

		fprintf( stderr, "Loading 'hat3' ... " );
		prep = fopen( "hat3", "r" );
		if( prep == NULL ) ErrorExit( "Make hat3." );
#if SMALLMEMORY
		if( multidist )
		{
//			readlocalhomtable_two( prep, norg, nadd, localhomtable, localhomtable+norg, kozoarivec );
			readlocalhomtable_one( prep, norg, nadd, localhomtable, kozoarivec );
		}
		else
#endif
		{
			readlocalhomtable( prep, njob, localhomtable, kozoarivec );
		}

		fclose( prep );
		fprintf( stderr, "\ndone.\n" );


		nkozo = 0;
		for( i=0; i<njob; i++ ) 
		{
//			fprintf( stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i] );
			if( kozoarivec[i] ) nkozo++;
		}
		if( nkozo )
		{
			topol_kozo = AllocateIntCub( nkozo, 2, 0 );
			len_kozo = AllocateFloatMtx( nkozo, 2 );
			iscore_kozo = AllocateFloatHalfMtx( nkozo );
			eff_kozo = AllocateDoubleVec( nkozo );
			eff_kozo_mapped = AllocateDoubleVec( njob );
		}


#if SMALLMEMORY
//		outlocalhom_part( localhomtable, norg, nadd );
#else
//		outlocalhom( localhomtable, njob );
#endif

#if 0
		fprintf( stderr, "Extending localhom ... " );
		extendlocalhom2( njob, localhomtable );
		fprintf( stderr, "done.\n" );
#endif
	}

#if 0
	readData( infp, name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
	fclose( infp );
#endif

	constants( njob, seq );

#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

//	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

	alignmentlength = strlen( seq[0] );
	for( i=0; i<norg; i++ )
	{
		if( hitout == 0.0 && alignmentlength != strlen( seq[i] ) )
		{
			fprintf( stderr, "#################################################################################\n" );
			fprintf( stderr, "# ERROR!                                                                        #\n" );
			fprintf( stderr, "# The original%4d sequences must be aligned                                    #\n", njob-nadd );
			fprintf( stderr, "#################################################################################\n" );
			exit( 1 );
		}
	}
	if( addprofile )
	{
		fprintf( stderr, "Not supported!\n" );
		exit( 1 );
	}

	if( tuplesize > 0 ) // if mtx is internally computed
	{
		if( multidist == 1 )
		{
			ktupledistancematrix( njob, norg, nlenmax, seq, name, iscore, nscore ); // iscore ha muda.

//			hat2p = fopen( "hat2-1", "w" );
//			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, iscore );
//			fclose( hat2p );

			dndpre( norg, seq, iscore );
//			fprintf( stderr, "Loading 'hat2i' (aligned sequences) ... " );
//			prep = fopen( "hat2i", "r" );
//			if( prep == NULL ) ErrorExit( "Make hat2i." );
//			readhat2_doublehalf_pointer( prep, njob-nadd, name, iscore );
//			fclose( prep );
//			fprintf( stderr, "done.\n" );

//			hat2p = fopen( "hat2-2", "w" );
//			WriteFloatHat2_pointer_halfmtx( hat2p, norg, name, iscore );
//			fclose( hat2p );
		}
		else
		{
			ktupledistancematrix( njob, norg, nlenmax, seq, name, iscore, nscore );
		}
	}
	else
	{
		if( multidist == 1 )
		{
			fprintf( stderr, "Loading 'hat2n' (aligned sequences - new sequences) ... " );
			prep = fopen( "hat2n", "r" );
			if( prep == NULL ) ErrorExit( "Make hat2n." );
			readhat2_doublehalf_part_pointer( prep, njob, nadd, name, nscore );
			fclose( prep );
			fprintf( stderr, "done.\n" );
		
			fprintf( stderr, "Loading 'hat2i' (aligned sequences) ... " );
			prep = fopen( "hat2i", "r" );
			if( prep == NULL ) ErrorExit( "Make hat2i." );
			readhat2_doublehalf_pointer( prep, njob-nadd, name, iscore );
			fclose( prep );
			fprintf( stderr, "done.\n" );
		}
		else
		{
			fprintf( stderr, "Loading 'hat2' ... " );
			prep = fopen( "hat2", "r" );
			if( prep == NULL ) ErrorExit( "Make hat2." );
			readhat2_doublehalf_pointer( prep, njob, name, iscore );
			fclose( prep );
			fprintf( stderr, "done.\n" );
		}
	}

#if 1
	if( distout )
	{
		fprintf( stderr, "Writing distances between new sequences and existing msa.\n" );
		hat2p = fopen( "hat2", "w" );
		if( multidist || tuplesize > 0 )
		{
			for( iadd=0; iadd<nadd; iadd++ ) 
			{
				fprintf( hat2p, "Distance between new sequence %d and %d sequences in existing msa\n", iadd+1, norg );
				for( i=0; i<norg; i++ ) 
				{
					fprintf( hat2p, "%5.3f ", nscore[i][iadd] );
					if( (i+1) % 12 == 0 ) fprintf( hat2p, "\n" );
				}
				fprintf( hat2p, "\n\n" );
			}
		}
		else
		{
			for( iadd=0; iadd<nadd; iadd++ ) 
			{
				fprintf( hat2p, "Distance between new sequence %d and %d sequences in existing msa\n", iadd+1, norg );
				for( i=0; i<norg; i++ ) 
				{
					fprintf( hat2p, "%5.3f ", iscore[i][norg+iadd-i] );
					if( (i+1) % 12 == 0 ) fprintf( hat2p, "\n" );
				}
				fprintf( hat2p, "\n\n" );
			}
		}
		fclose( hat2p );
//		exit( 1 );
//		hat2p = fopen( "hat2", "w" );
//		WriteFloatHat2_pointer_halfmtx( hat2p, norg, name, iscore );
//		fclose( hat2p );
//		exit( 1 );
	}
#endif


	singlerna = NULL;


	if( diffout )
	{
		lenfull = strlen( seq[0] );
		originalgaps = (char *)calloc( lenfull+1, sizeof( char ) );
		recordoriginalgaps( originalgaps, norg, seq );

		difflist = (GapPos **)calloc( nadd+1, sizeof( GapPos * ) );
		for( i=0; i<nadd; i++ )
		{
			difflist[i] = calloc( 1, sizeof( GapPos ) );
			difflist[i][0].pos = -1;
			difflist[i][0].len = 0;
		}
		difflist[nadd] = NULL;
		deletelist = NULL;
	}

	else if( keeplength )
	{
		lenfull = strlen( seq[0] );
		originalgaps = (char *)calloc( lenfull+1, sizeof( char ) );
		recordoriginalgaps( originalgaps, norg, seq );


		deletelist = (GapPos **)calloc( nadd+1, sizeof( GapPos * ) );
		for( i=0; i<nadd; i++ )
		{
			deletelist[i] = calloc( 1, sizeof( GapPos ) );
			deletelist[i][0].pos = -1;
			deletelist[i][0].len = 0;
		}
		deletelist[nadd] = NULL;
		difflist = NULL;
	}
	else
	{
		originalgaps = NULL;
		deletelist = NULL;
		difflist = NULL;
	}

	commongappick( norg, seq );
	lenfull = strlen( seq[0] );

	if( 0 && ratioambiguouscolumn_simple( seq, norg, lenfull ) > 0.03 )
	{
		fprintf( stderr, "################################################################################\n" );
		fprintf( stderr, "# \n" );
		fprintf( stderr, "# The reference MSA has >3%% ambiguous columns.\n" );
		fprintf( stderr, "# Please prepare a better reference.\n" );
		fprintf( stderr, "# \n" );
		fprintf( stderr, "################################################################################\n" );
		exit( 1 );
	}

	if( keeplength && mapout )
	{
		addbk = (char **)calloc( nadd+1, sizeof( char * ) );
		for( i=0; i<nadd; i++ )
		{
			ien = strlen( seq[norg+i] );
			addbk[i] = (char *)calloc( ien + 1, sizeof( char ) );
			gappick0( addbk[i], seq[norg+i] );
		}
		addbk[nadd] = NULL;
	}
	else
	{
		addbk = NULL;
	}

	

//	newgaplist_o = AllocateIntMtx( nadd, alloclen ); //ookisugi
	newgaplist_o = AllocateIntMtx( nadd, lenfull*2 );
	newgaplist_compact = AllocateIntVec( lenfull*2 );
	istherenewgap = AllocateIntVec( nadd );
	follower = AllocateIntMtx( norg, 1 );
	for( i=0; i<norg; i++ ) follower[i][0] = -1;
	follows = AllocateIntVec( nadd );

	dep = (Treedep *)calloc( norg, sizeof( Treedep ) );
	topol = AllocateIntCub( norg, 2, 0 );
	len = AllocateFloatMtx( norg, 2 );
//	iscoreo = AllocateFloatHalfMtx( norg );
	mtxcpy( norg, norg, &iscoreo, iscore );

	if( treeout )
	{
		addtree = (Addtree *)calloc( nadd, sizeof( Addtree ) );
		if( !addtree )
		{
			fprintf( stderr, "Cannot allocate addtree\n" );
			exit( 1 );
		}
	}


//	nlim = norg-1;
//	for( i=0; i<nlim; i++ )
//	{
//		fptc = iscoreo[i]+1;
//		fpt  = iscore[i]+1;
//		j = norg-i-1;
//		while( j-- )
//			*fptc++ = *fpt++;
////	for( j=i+1; j<norg; j++ )	
////		iscoreo[i][j-i] = iscore[i][j-i];
//	}

//	fprintf( stderr, "building a tree.." );
	if( treein )
	{
		reporterr(       "Loading a tree ... " );
		loadtop( norg, iscoreo, topol, len, name, NULL, dep ); // nogaplen?
		reporterr( "\ndone.\n\n" );
	}
	else if( treeout )
		fixed_musclesupg_double_realloc_nobk_halfmtx_treeout( norg, iscoreo, topol, len, name, nlen, dep, 1 );
	else
		fixed_musclesupg_double_realloc_nobk_halfmtx( norg, iscoreo, topol, len, dep, 0, 1 );
//	fprintf( stderr, "done.\n" );

	if( norg > 1 ) 
		cnctintvec( ordertable, topol[norg-2][0], topol[norg-2][1] );
	else
	{
		ordertable[0] = 0; ordertable[1] = -1;
	}
	FreeFloatHalfMtx( iscoreo, norg );

#ifdef enablemultithread
	if( nthread )
	{
		pthread_t *handle;
		pthread_mutex_t mutex_counter;
		thread_arg_t *targ;
		int *iaddsharept;

		targ = calloc( nthread, sizeof( thread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex_counter, NULL );
		iaddsharept = calloc( 1, sizeof(int) );
		*iaddsharept = 0;

		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].follows = follows;
			targ[i].njob = njob; 
			targ[i].nadd = nadd; 
			targ[i].nlen = nlen; 
			targ[i].name = name;
			targ[i].seq = seq;
			targ[i].localhomtable = localhomtable;
			targ[i].iscore = iscore;
			targ[i].nscore = nscore;
			targ[i].istherenewgap = istherenewgap;
			targ[i].newgaplist = newgaplist_o;
			targ[i].singlerna = singlerna;
			targ[i].eff_kozo_mapped = eff_kozo_mapped;
			targ[i].alloclen = alloclen;
			targ[i].iaddshare = iaddsharept;
			targ[i].dep = dep;
			targ[i].topol = topol;
			targ[i].len = len;
			targ[i].addtree = addtree;
			targ[i].deletelist = deletelist;
			targ[i].difflist = difflist;
			targ[i].mutex_counter = &mutex_counter;
			pthread_create( handle+i, NULL, addsinglethread, (void *)(targ+i) );
		}
		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex_counter );
		free( handle );
		free( targ );
		free( iaddsharept );
	}
	else
#endif
	{
		thread_arg_t *targ;
		targ = calloc( 1, sizeof( thread_arg_t ) );
		targ[0].follows = follows;
		targ[0].njob = njob; 
		targ[0].nadd = nadd; 
		targ[0].nlen = nlen; 
		targ[0].name = name;
		targ[0].seq = seq;
		targ[0].localhomtable = localhomtable;
		targ[0].iscore = iscore;
		targ[0].nscore = nscore;
		targ[0].istherenewgap = istherenewgap;
		targ[0].newgaplist = newgaplist_o;
		targ[0].singlerna = singlerna;
		targ[0].eff_kozo_mapped = eff_kozo_mapped;
		targ[0].alloclen = alloclen;
		targ[0].dep = dep;
		targ[0].topol = topol;
		targ[0].len = len;
		targ[0].addtree = addtree;
		targ[0].deletelist = deletelist;
		targ[0].difflist = difflist;
		addsinglethread( targ );
		free( targ );
	}
	free( dep );
	FreeFloatMtx( len );
	if( multidist || tuplesize > 0 )
	{
		FreeFloatHalfMtx( iscore, norg );
		FreeFloatMtx( nscore );
	}
	else
	{
		FreeFloatHalfMtx( iscore, njob );
	}


//	for( i=0; i<nadd; i++ ) fprintf( stdout, ">%s (%d) \n%s\n", name[norg+i], norg+i, seq[norg+i] );

	if( treeout )
	{
		fp = fopen( "infile.tree", "a" );
		if( fp == 0 )
		{
			fprintf( stderr, "File error!\n" );
			exit( 1 );
		}
		for( i=0; i<nadd; i++ )
		{
			fprintf( fp, "\n" );
			fprintf( fp, "%8d: %s\n", norg+i+1, name[norg+i]+1 );
			fprintf( fp, "          nearest sequence: %d\n", addtree[i].nearest + 1 );
			fprintf( fp, "          approximate distance: %f\n", addtree[i].dist1 );
			fprintf( fp, "          sister group: %s\n", addtree[i].neighbors );
			fprintf( fp, "          approximate distance: %f\n", addtree[i].dist2 );
			free( addtree[i].neighbors );
		}
		fclose( fp );
		free( addtree );
	}

	for( iadd=0; iadd<nadd; iadd++ )
	{
		f = follows[iadd];
		for( i=0; follower[f][i]!=-1; i++ )
			;
		if( !(follower[f] = realloc( follower[f], (i+2)*sizeof(int) ) ) )
		{
			fprintf( stderr, "Cannot reallocate follower[]" );
			exit( 1 );
		}
		follower[f][i] = iadd;
		follower[f][i+1] = -1;
#if 0
		fprintf( stderr, "\nfollowers of %d = ", f );
		for( i=0; follower[f][i]!=-1; i++ )
			fprintf( stderr, "%d ", follower[f][i]  );
		fprintf( stderr, "\n" );
#endif
	}

	orderfp = fopen( "order", "w" );
	if( !orderfp )
	{
		fprintf( stderr, "Cannot open 'order'\n" );
		exit( 1 );
	}
	for( i=0; ordertable[i]!=-1; i++ )
	{
		fprintf( orderfp, "%d\n", ordertable[i] );
//		for( j=0; follower[i][j]!=-1; j++ )
//			fprintf( orderfp, "%d\n", follower[i][j]+norg );
		for( j=0; follower[ordertable[i]][j]!=-1; j++ )
			fprintf( orderfp, "%d\n", follower[ordertable[i]][j]+norg );
//			fprintf( orderfp, "%d -> %d\n", follower[i][j]+norg, i );
	}
	fclose( orderfp );

	posmap = AllocateIntVec( lenfull+2 );
	realign = calloc( lenfull+2, sizeof( Blocktorealign ) );
	for( i=0; i<lenfull+1; i++ ) posmap[i] = i;
	for( i=0; i<lenfull+1; i++ )
	{
		realign[i].nnewres = 0;
		realign[i].start = 0;
		realign[i].end = 0;
	}

	fprintf( stderr, "\n\nCombining ..\n" );
	fflush( stderr );
	tmpseqlen = alloclen * 100;
	tmpseq = AllocateCharMtx( 1, tmpseqlen );


//	check1 = AllocateCharVec( tmpseqlen );
//	check2 = AllocateCharVec( tmpseqlen );
//	gappick0( check2, seq[0] );
	for( iadd=0; iadd<nadd; iadd++ )
	{
//		fprintf( stderr, "%d / %d\r", iadd, nadd );
		fflush( stderr );

//		fprintf( stderr, "\niadd == %d\n", iadd );
		makegaplistcompact( lenfull, posmap, newgaplist_compact, newgaplist_o[iadd] );
		if( iadd == 0 || istherenewgap[iadd] )
		{
			tmpseq1 = tmpseq[0];
// 			gaplist2alnx( lenfull, tmpseq1, seq[0], newgaplist_o[iadd], posmap, tmpseqlen );
 			gaplist2alnx( lenfull, tmpseq1, seq[0], newgaplist_compact, posmap, tmpseqlen );
//			fprintf( stderr, "len = %d ? %d\n", strlen( tmpseq1 ), alloclen );
			if( ( tmplen = strlen( tmpseq1 ) ) >= fullseqlen )
			{
				fullseqlen = tmplen * 2+1;
//				fprintf( stderr, "Length over!\n" );
//				fprintf( stderr, "strlen(tmpseq1)=%d\n", (int)strlen( tmpseq1 ) );
				fprintf( stderr, "reallocating..." );
//				fprintf( stderr, "alloclen=%d\n", alloclen );
//				fprintf( stderr, "Please recompile!\n" );
//				exit( 1 );
				for( i=0; i<njob; i++ )
				{
					seq[i] = realloc( seq[i], fullseqlen * sizeof( char ) );
					if( !seq[i] )
					{
						fprintf( stderr, "Cannot reallocate seq[][]\n" );
						exit( 1 );
					}
				}
				fprintf( stderr, "done.\n" );
			}
			strcpy( seq[0], tmpseq1 );

			ien = norg+iadd;
#ifdef enablemultithread
			if( nthread > 0 && ien > 500 )
			{
				gaplist2alnxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				fprintf( stderr, "%d / %d (threads %d-%d)\r", iadd, nadd, 0, nthread );

				targ = calloc( nthread, sizeof( gaplist2alnxthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
				jobpos = 1;
				for( i=0; i<nthread; i++ )
				{
//					targ[i].thread_no = i;
					targ[i].ncycle = ien;
					targ[i].jobpospt = &jobpos;
					targ[i].tmpseqlen = tmpseqlen;
					targ[i].lenfull = lenfull;
					targ[i].seq = seq;
//					targ[i].newgaplist = newgaplist_o[iadd];
					targ[i].newgaplist = newgaplist_compact;
					targ[i].posmap = posmap;
					targ[i].mutex = &mutex;

					pthread_create( handle+i, NULL, gaplist2alnxthread, (void *)(targ+i) );
				}
				for( i=0; i<nthread; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
				fprintf( stderr, "%d / %d\r", iadd, nadd );
				for( i=1; i<ien; i++ )
				{
					tmpseq1 = tmpseq[0];
					if( i == 1 ) fprintf( stderr, " %d / %d\r", iadd, nadd );
// 					gaplist2alnx( lenfull, tmpseq1, seq[i], newgaplist_o[iadd], posmap, tmpseqlen  );
 					gaplist2alnx( lenfull, tmpseq1, seq[i], newgaplist_compact, posmap, tmpseqlen  );
//					fprintf( stderr, ">%s (iadd=%d)\n%s\n", name[i], iadd, tmpseq1 );
					strcpy( seq[i], tmpseq1 );
				}
			}
		}
		tmpseq1 = tmpseq[0];
//		insertgapsbyotherfragments_simple( lenfull, tmpseq1, seq[norg+iadd], newgaplist_o[iadd], posmap );
		insertgapsbyotherfragments_compact( lenfull, tmpseq1, seq[norg+iadd], newgaplist_o[iadd], posmap );
//		fprintf( stderr, "%d = %s\n", iadd, tmpseq1 );
		eq2dash( tmpseq1 );
		strcpy( seq[norg+iadd], tmpseq1 );

//		adjustposmap( lenfull, posmap, newgaplist_o[iadd] );
		adjustposmap( lenfull, posmap, newgaplist_compact );
		countnewres( lenfull, realign, posmap, newgaplist_o[iadd] ); // muda?
//		countnewres( lenfull, realign, posmap, newgaplist_compact ); // muda?

	}
	fprintf( stderr, "\r   done.                      \n\n" );

#if 0
	for( i=0; i<njob; i++ )
	{
		fprintf( stderr, ">%s\n", name[i] );
		fprintf( stderr, "%s\n", seq[i] );
	}
#endif

#if 0
	fprintf( stderr, "realign[].nnewres = " );
	for( i=0; i<lenfull+1; i++ )
	{
		fprintf( stderr, "%d ", realign[i].nnewres );
	}
	fprintf( stderr, "\n" );
	for( i=0; i<lenfull+1; i++ ) fprintf( stderr, "i=%d, nnewres=%d, start=%d, end=%d\n", i, realign[i].nnewres, realign[i].start, realign[i].end );
#endif

	
	if( smoothing )
	{
//		for( i=0; i<lenfull+1; i++ ) fprintf( stderr, "i=%d, nnewres=%d, start=%d, end=%d\n", i, realign[i].nnewres, realign[i].start, realign[i].end );
		smoothing1( njob, nadd, lenfull, seq, realign );
//		for( i=0; i<lenfull+1; i++ ) fprintf( stderr, "i=%d, nnewres=%d, start=%d, end=%d\n", i, realign[i].nnewres, realign[i].start, realign[i].end );
		smoothing2( njob, nadd, lenfull, seq, realign );
	}

	for( i=0; i<lenfull+1; i++ )
	{
		if( realign[i].nnewres > 1 ) 
		{
//			fprintf( stderr, "i=%d: %d-%d\n", i, realign[i].start, realign[i].end );
			fprintf( stderr, "\rRealigning %d/%d           \r", i, lenfull );
//			zure = dorealignment_compact( realign+i, seq, &fullseqlen, norg );
//			zure = dorealignment_order( realign+i, seq, &fullseqlen, norg, ordertable, follows );
			zure = dorealignment_tree( realign+i, seq, &fullseqlen, norg, topol, follows );
#if 0
			gappick0( check1, seq[0] );
			fprintf( stderr, "check1 = %s\n", check1 );
			if( strcmp( check1, check2 ) )
			{
				fprintf( stderr, "CHANGED!!!!!\n" );
				exit( 1 );
			}
#endif
			for( j=i+1; j<lenfull+1; j++ )
			{
				if( realign[j].nnewres )
				{
					realign[j].start -= zure;
					realign[j].end -= zure;
				}
			}
		}
	}
	FreeIntCub( topol );
	fprintf( stderr, "\r   done.                      \n\n" );
	fflush( stderr );


	if( keeplength )
	{
		FILE *dlf;
		restoreoriginalgaps( njob, seq, originalgaps );

		dlf = fopen( "_deletelist", "w" );
		for( i=0; i<nadd; i++ )
		{
			if( deletelist[i] )
				for( j=0; deletelist[i][j].pos!=-1; j++ )
					fprintf( dlf, "%d %d %d\n", njob-nadd+i, deletelist[i][j].pos, deletelist[i][j].len ); // 0origin
		}
		fclose( dlf );

		if( mapout )
		{
			dlf = fopen( "_deletemap", "w" );
			if( mapout == 1 )
				reconstructdeletemap( nadd, addbk, deletelist, seq+njob-nadd, dlf, name+njob-nadd );
			else
				reconstructdeletemap_compact( nadd, addbk, deletelist, seq+njob-nadd, dlf, name+njob-nadd );
			fclose( dlf );
		}
	}



	if( commonIP ) FreeIntMtx( commonIP );
	commonIP = NULL;
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	FreeIntMtx( newgaplist_o );
	FreeIntVec( newgaplist_compact );
	FreeIntVec( posmap );
	free( realign );
	free( istherenewgap );
	FreeIntMtx( follower );
	free( follows );
	free( ordertable );
	free( kozoarivec );
	free( nlen );
	FreeCharMtx( tmpseq );
	freeconstants();
	if( addbk ) FreeCharMtx( addbk ); addbk = NULL;
	if( deletelist ) 
	{
		for( i=0; deletelist[i] != NULL; i++ ) free( deletelist[i] );
		free( deletelist );
		deletelist = NULL;
	}
	if( difflist ) 
	{
		for( i=0; difflist[i] != NULL; i++ ) free( difflist[i] );
		free( difflist );
		difflist = NULL;
	}
	if( originalgaps ) free( originalgaps ); originalgaps = NULL;

	writeData_pointer( prep_g, njob, name, nlen, seq );
#if 0
	writeData( stdout, njob, name, nlen, bseq );
	writePre( njob, name, nlen, bseq, !contin );
	writeData_pointer( prep_g, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif

#if SMALLMEMORY
	if( multidist )
	{
//		if( constraint ) FreeLocalHomTable_two( localhomtable, norg, nadd );
		if( constraint ) FreeLocalHomTable_one( localhomtable, norg, nadd );
	}
	else
#endif
	{
		if( constraint ) FreeLocalHomTable( localhomtable, njob );
	}

	FreeCharMtx( seq );
	FreeCharMtx( name );
	closeFiles();
	commonsextet_p( NULL, NULL );

	SHOWVERSION;
	if( ndeleted > 0 )
	{
		reporterr( "\nTo keep the alignment length, %d letters were DELETED.\n", ndeleted );
		if( mapout )
			reporterr( "The deleted letters are shown in the (filename).map file.\n" );
		else
			reporterr( "To know the positions of deleted letters, rerun the same command with the --mapout option.\n" );
	}
	return( 0 );
}
