#include "mltaln.h"

static void strncpy0( char *s1, char *s2, int n )
{
	while( n-- ) *s1++ = *s2++;
	*s1 = 0;
}

#if 0
static void strncpy0x( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0b0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0b1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0b2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0n0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0n1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0n2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0a0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0a1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0a2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0o0( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0o1( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
static void strncpy0o2( char *s1, char *s2, int n ) { while( n-- ) *s1++ = *s2++; *s1 = 0; }
#endif

static void eqpick( char *aseq, char *seq )
{
	for( ; *seq != 0; seq++ )
	{
		if( *seq != '=' )
			*aseq++ = *seq;
	}
	*aseq = 0;

}

void profilealignment2( int n0, int n2, char **aln0, char **aln2, int alloclen, char alg ) // n1 ha allgap
{
	int i, newlen;
	double *effarr0, *effarr2;
	int *allgap0, *allgap2;
	double dumdb;
	int alcount0, alcount2;

	if( aln0[0][1] == 0 && aln2[0][1] == 0 ) return; // --allowshift no tokiha...
//	reporterr( "profilealignment!\n" );

	commongappick( n0, aln0 );
	commongappick( n2, aln2 );

	effarr0 = AllocateDoubleVec( n0 );
	effarr2 = AllocateDoubleVec( n2 );
	allgap0 = AllocateIntVec( n0 );
	allgap2 = AllocateIntVec( n2 );

#if 1 // new weight 2015/Jun
	alcount0 = 0;
	for( i=0; i<n0; i++ ) 
	{
		if( isallgap( aln0[i] ) ) allgap0[i] = 1;
		else
		{
			alcount0++;
			allgap0[i] = 0;
		}
	}

	alcount2 = 0;
	for( i=0; i<n2; i++ ) 
	{
		if( isallgap( aln2[i] ) ) allgap2[i] = 1;
		else
		{
			alcount2++;
			allgap2[i] = 0;
		}
	}

	for( i=0; i<n0; i++ ) if( !allgap0[i] ) effarr0[i] = 1.0 / (double)( alcount0 );
	for( i=0; i<n2; i++ ) if( !allgap2[i] ) effarr2[i] = 1.0 / (double)( alcount2 );
#else
	eff = 1.0 / (double)n0; for( i=0; i<n0; i++ ) effarr0[i] = eff;
	eff = 1.0 / (double)n2; for( i=0; i<n2; i++ ) effarr2[i] = eff;
#endif

	newgapstr = "-";
	if( alg == 'M' )
		MSalignmm( n_dis_consweight_multi, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, 1, 1, NULL, NULL, NULL, 0.0, 0.0 ); //outgap=1, 2014/Dec/1
	else
		A__align( n_dis_consweight_multi, penalty, penalty_ex, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, 1, 1, -1, -1, NULL, NULL, NULL, 0.0, 0.0 ); //outgap=1, 2014/Dec/1

	newlen = strlen( aln0[0] );

#if 0 // tabun hitsuyou
	for( j=0; j<newlen; j++ )
	{
//		fprintf( stderr, "j=%d\n", j );
		for( i=0; i<n0; i++ )
		{
			if( aln0[i][j] != '-' ) break;
		}
		if( i == n0 ) 
		{
			for( i=0; i<n1; i++ ) 
			{
				if( aln1[i][j] != '-' ) break;
			}
		}
		else i = -1;

		if( i == n1 ) 
		{
			for( i=0; i<n1; i++ ) aln1[i][j] = '=';
		}
	}
	fprintf( stderr, "in profilealignment,\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "\n>aln0[%d] = \n%s\n", i, aln0[i] );
	for( i=0; i<n1; i++ ) fprintf( stderr, "\n>aln1[%d] = \n%s\n", i, aln1[i] );
	for( i=0; i<n2; i++ ) fprintf( stderr, "\n>aln2[%d] = \n%s\n", i, aln2[i] );
#endif

	free( effarr0 );
	free( effarr2 );
	free( allgap0 );
	free( allgap2 );
}

static void profilealignment( int n0, int n1, int n2, char **aln0, char **aln1, char **aln2, int alloclen, char alg ) // n1 ha allgap
{
	int i, j, newlen;
	double *effarr0 = NULL, *effarr2 = NULL;
	int *allgap0 = NULL, *allgap2 = NULL;
	double dumdb;
	int alcount0, alcount2;
	char *cptr;
	
//	effarr0 = AllocateDoubleVec( n0 );
//	effarr2 = AllocateDoubleVec( n2 );
//	allgap0 = AllocateIntVec( n0 );
//	allgap2 = AllocateIntVec( n2 );
//
	if( aln0[0][1] == 0 && aln2[0][1] == 0 ) return; // --allowshift no tokiha...

//	reporterr( "In profilealignment(), strlen( aln0[0] ) %d\n", strlen( aln0[0] ) );
//	reporterr( "In profilealignment(), strlen( aln2[0] ) %d\n", strlen( aln2[0] ) );

	commongappick( n0, aln0 );
	commongappick( n2, aln2 );

//	reporterr( "after commongappick, strlen( aln0[0] ) %d\n", strlen( aln0[0] ) );
//	reporterr( "after commongappick, strlen( aln2[0] ) %d\n", strlen( aln2[0] ) );

//	reporterr( "\n\n\n" );

	if( aln2[0][0] == 0 )
	{
		newlen = j = strlen( aln0[0] );
		cptr = aln2[0];
		while( j-- ) *cptr++ = '-';
		*cptr = 0;

		cptr = aln2[0];
		for( i=1; i<n2; i++ ) strcpy( aln2[i], cptr );
		return;
	}

#if 1
	effarr0 = (double *)malloc( n0 * sizeof( double ) );
	effarr2 = (double *)malloc( n2 * sizeof( double ) );
	allgap0 = (int *)malloc( n0 * sizeof( int ) );
	allgap2 = (int *)malloc( n2 * sizeof( int ) );
#else
	effarr0 = (double *)calloc( n0, sizeof( double ) );
	effarr2 = (double *)calloc( n2, sizeof( double ) );
	allgap0 = (int *)calloc( n0, sizeof( int ) );
	allgap2 = (int *)calloc( n2, sizeof( int ) );
#endif

#if 1 // new weight 2015/Jun
	alcount0 = 0;
	for( i=0; i<n0; i++ ) 
	{
		if( isallgap( aln0[i] ) ) allgap0[i] = 1;
		else
		{
			alcount0++;
			allgap0[i] = 0;
		}
	}

	alcount2 = 0;
	for( i=0; i<n2; i++ ) 
	{
		if( isallgap( aln2[i] ) ) allgap2[i] = 1;
		else
		{
			alcount2++;
			allgap2[i] = 0;
		}
	}

	for( i=0; i<n0; i++ ) if( !allgap0[i] ) effarr0[i] = 1.0 / (double)( alcount0 ); else effarr0[i] = 0.0; // malloc / alloca no baai
	for( i=0; i<n2; i++ ) if( !allgap2[i] ) effarr2[i] = 1.0 / (double)( alcount2 ); else effarr2[i] = 0.0; // malloc / alloca no baai
#else
	for( i=0; i<n0; i++ ) if( !allgap0[i] ) effarr0[i] = 1.0 / (double)( alcount0 ); // calloc no baai
	for( i=0; i<n2; i++ ) if( !allgap2[i] ) effarr2[i] = 1.0 / (double)( alcount2 ); // calloc no baai
	eff = 1.0 / (double)n0; for( i=0; i<n0; i++ ) effarr0[i] = eff;
	eff = 1.0 / (double)n2; for( i=0; i<n2; i++ ) effarr2[i] = eff;
#endif

	newgapstr = "-";
	if( alg == 'M' )
		MSalignmm( n_dis_consweight_multi, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, 1, 1, NULL, NULL, NULL, 0.0, 0.0 ); //outgap=1, 2014/Dec/1
	else
		A__align( n_dis_consweight_multi, penalty, penalty_ex, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, 1, 1, -1, -1, NULL, NULL, NULL, 0.0, 0.0 ); //outgap=1, 2014/Dec/1

	newlen = strlen( aln0[0] );

	for( i=0; i<newlen; i++ ) aln1[0][i] = '-';
	aln1[0][i] = 0;
	for( i=1; i<n1; i++ ) strcpy( aln1[i], aln1[0] );

	for( j=0; j<newlen; j++ )
	{
//		fprintf( stderr, "j=%d\n", j );
		for( i=0; i<n0; i++ )
		{
			if( aln0[i][j] != '-' ) break;
		}
		if( i == n0 ) 
		{
			for( i=0; i<n1; i++ ) 
			{
				if( aln1[i][j] != '-' ) break;
			}
		}
		else i = -1;

		if( i == n1 ) 
		{
			for( i=0; i<n1; i++ ) aln1[i][j] = '=';
		}
	}
#if 0
	fprintf( stderr, "in profilealignment, before commongappick\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "\n>aln0[%d] = %s\n", i, aln0[i] );
	for( i=0; i<n1; i++ ) fprintf( stderr, "\n>aln1[%d] = %s\n", i, aln1[i] );
	for( i=0; i<n2; i++ ) fprintf( stderr, "\n>aln2[%d] = %s\n", i, aln2[i] );
#endif

#if 0
	fprintf( stderr, "in profilealignment, after commongappick\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "\n>aln0[%d] = %s\n", i, aln0[i] );
	for( i=0; i<n1; i++ ) fprintf( stderr, "\n>aln1[%d] = %s\n", i, aln1[i] );
	for( i=0; i<n2; i++ ) fprintf( stderr, "\n>aln2[%d] = %s\n", i, aln2[i] );
#endif


	free( effarr0 );
	free( effarr2 );
	free( allgap0 );
	free( allgap2 );
}

void eq2dashmatomete( char **s, int n )
{
	int i, j;
	char sj;

	for( j=0; (sj=s[0][j]); j++ )
	{
		if( sj == '=' )
		{
			for( i=0; i<n; i++ )
			{
				s[i][j] = '-';
			}
		}
	}
}

void eq2dashmatometehayaku( char **s, int n )
{
	int i, j, c;
	int *tobechanged;
	int len = strlen( s[0] );

	tobechanged = calloc( len+1, sizeof( int ) ); // len+1, 2017/Nov/15
	c = 0;
	for( j=0; j<len; j++ )
	{
		if( s[0][j] == '=' ) tobechanged[c++] = j;
	}
	tobechanged[c] = -1;

	for( i=0; i<n; i++ )
	{
		for( c=0; (j=tobechanged[c])!=-1; c++ )
			s[i][j] = '-';
	}
	free( tobechanged );
}

void eq2dash( char *s )
{
	while( *s )
	{
		if( *s == '=' ) 
		{
			*s = '-';
		}
		s++;
	}
}

static void plus2gapchar( char *s, char gapchar )
{
	while( *s )
	{
		if( *s == '+' ) 
		{
			*s = gapchar;
		}
		s++;
	}
}

void findnewgaps( int n, int rep, char **seq, int *gaplen )
{
	int i, pos, len, len1;

	len = strlen( seq[0] );	
//	for( i=0; i<len; i++ ) gaplen[i] = 0; // calloc de shokika sareteirukara hontou ha iranai
	len1 = len + 1;
	for( i=0; i<len1; i++ ) gaplen[i] = 0; // realloc de shokika sareteirukara iru!
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		if( seq[rep][i] == '=' ) 
		{
			if( disp ) fprintf( stderr, "Newgap! pos = %d\n", pos );
			gaplen[pos]++;
		}
		else
			pos++;
	}

#if 0
	if( disp )
	{
		fprintf( stderr, "\ngaplen[] in findnewgaps() = \n" );
		for(i=0; i<pos; i++ ) fprintf( stderr, "%d ", gaplen[i] );
		fprintf( stderr, "\n" );
		reporterr( "pos=%d\n", pos );
	}
#endif
}

void findcommongaps( int n, char **seq, int *gapmap )
{
	int i, j, pos, len, len1;
	len = strlen( seq[0] );	
	len1 = len+1;

//	fprintf( stderr, "seq[0] = %s\n", seq[0] );
	for( i=0; i<len1; i++ ) gapmap[i] = 0;
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		for( j=0; j<n; j++ )
			if( seq[j][i] != '-' ) break;

		if( j == n ) gapmap[pos]++;
		else
			pos++;
	}
#if 0
	for( i=0; i<pos; i++ )
	{
		fprintf( stderr, "vec[%d] = %d\n", i, gapmap[i] );
	}
#endif
}

void adjustgapmap( int newlen, int *gapmap, char *seq )
{
	int j;
	int pos;
	int newlen1 = newlen+1;
	int *tmpmap;


	tmpmap = AllocateIntVec( newlen+2 );
	j = 0;
	pos = 0;
	while( *seq )
	{
//		fprintf( stderr, "j=%d *seq = %c\n", j, *seq );
		if( *seq++ == '=' )
			tmpmap[j++] = 0;
		else
		{
			tmpmap[j++] = gapmap[pos++];
		}
	}
	tmpmap[j++] = gapmap[pos];

	for(j=0; j<newlen1; j++)
		gapmap[j] = tmpmap[j];

	free( tmpmap );

#if 0
	reporterr( "gapmap in adjustgapmap() =\n" );
	for(j=0; j<newlen1; j++) reporterr( "%d ", gapmap[j] );
	reporterr( "length = %d\n", newlen );
#endif
}


static int countnogaplen( int *gaplen, int *term )
{
	int v = 0;
	while( gaplen < term )
	{
		if( *gaplen++ == 0 ) v++;
		else break;
	}
	return( v );
}

static int countgapmap( int *gapmap, int *term )
{
	int v = 0;
	while( gapmap < term )
	{
//		reporterr( "*gapmap = %d\n", *gapmap );
		if( *gapmap++ == 0 ) v++;
		else break;
	}
	return( v );
}

void insertnewgaps( int njob, int *alreadyaligned, char **seq, int *ex1, int *ex2, int *gaplen, int *gapmap, int alloclen, char alg, char gapchar )
{
	int *mar;
	char *gaps;
	char *cptr;
	int i, j, k, len, rep, len0, lp, blocklen;
	char **mseq2, **mseq0, **mseq1;
	char **aseq, *newchar;
	int ngroup2, ngroup0, ngroup1;
	int *list0, *list1, *list2;
	int posin12, gapshift, newpos;
	int mlen1, mlen0, mlen2;


	mar = calloc( njob, sizeof( int ) );
	list0 = calloc( njob, sizeof( int ) );
	list1 = calloc( njob, sizeof( int ) );
	list2 = calloc( njob, sizeof( int ) );

	for( i=0; i<njob; i++ ) mar[i] = 0;
	for( i=0; i<njob; i++ ) 
	{
		if( alreadyaligned[i]==0 ) mar[i] = 3;
	}
	for( i=0; (k=ex1[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex1[i] );
	}
	for( i=0; (k=ex2[i])>-1; i++ ) 
	{
		mar[k] = 2;
//		fprintf( stderr, "excluding %d\n", ex2[i] );
	}

	ngroup2 = ngroup1 = ngroup0 = 0;
	for( i=0; i<njob; i++ )
	{
		if( mar[i] == 2 ) 
		{
			list2[ngroup2] = i;
			ngroup2++;
		}
		if( mar[i] == 1 ) 
		{
			list1[ngroup1] = i;
			ngroup1++;
		}
		if( mar[i] == 0 ) 
		{
			list0[ngroup0] = i;
//			fprintf( stderr, "inserting new gaps to %d\n", i );
			ngroup0++;
		}
	}
	list0[ngroup0] = list1[ngroup1] = list2[ngroup2] = -1;
	if( ngroup0 == 0 )
	{
//		fprintf( stderr, "Nothing to do\n" );
		free( mar );
		free( list0 );
		free( list1 );
		free( list2 );
		return;
	}

	for( i=0; i<njob; i++ ) if( mar[i] == 0 ) break;
	rep = i;
	len = strlen( seq[rep] );
	len0 = len+1;

//
//	if( i == njob )
//	{
////		fprintf( stderr, "Nothing to do\n" );
//		free( mar );
//		return;
//	}

	mseq2 = AllocateCharMtx( ngroup2, alloclen );
	mseq1 = AllocateCharMtx( ngroup1, alloclen );
	mseq0 = AllocateCharMtx( ngroup0, alloclen );
	aseq = AllocateCharMtx( njob, alloclen );
	gaps = calloc( alloclen, sizeof( char ) );

	for( i=0; i<njob; i++ ) aseq[i][0] = 0;
	newpos = 0;
	posin12 = 0;
#if 0
	fprintf( stderr, "\ngaplen[] = \n" );
	for(i=0; i<len0; i++ ) fprintf( stderr, "%d", gaplen[i] );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\ngapmap[] = \n" );
	for(i=0; i<len0; i++ ) fprintf( stderr, "%d", gapmap[i] );
	fprintf( stderr, "\n" );
#endif

	for( j=0; j<len0; j++ )
	{
//		fprintf( stderr, "\nj=%d, gaplen[%d]=%d\n", j, j, gaplen[j] );
		if( gaplen[j] )
		{
//			fprintf( stderr, "j=%d GAP!\n", j );
			for( i=0; i<ngroup0; i++ ) mseq0[i][0] = 0;
			for( i=0; i<ngroup1; i++ ) mseq1[i][0] = 0;
			for( i=0; i<ngroup2; i++ ) mseq2[i][0] = 0;
			mlen0 = mlen1 = mlen2 = 0;

			gapshift = gaplen[j];
			cptr = gaps;
			while( gapshift-- ) *cptr++ = gapchar;
			*cptr = 0;
			gapshift = gaplen[j];

			for( i=0; i<ngroup0; i++ ) strncpy0( mseq0[i]+mlen0, gaps, gapshift );
			for( i=0; i<ngroup1; i++ ) strncpy0( mseq1[i]+mlen1, seq[list1[i]]+posin12, gapshift );
			for( i=0; i<ngroup2; i++ ) strncpy0( mseq2[i]+mlen2, seq[list2[i]]+posin12, gapshift );
			posin12 += gapshift;
			mlen0 += gapshift;
			mlen1 += gapshift;
			mlen2 += gapshift;

			gapshift = gapmap[posin12];
//			fprintf( stderr, "gapmap[%d] kouho = %d\n", posin12, gapmap[posin12] );


			for( i=0; i<ngroup0; i++ ) strncpy0( mseq0[i]+mlen0, seq[list0[i]]+j, gapshift );
			for( i=0; i<ngroup1; i++ ) strncpy0( mseq1[i]+mlen1, seq[list1[i]]+posin12, gapshift );
			for( i=0; i<ngroup2; i++ ) strncpy0( mseq2[i]+mlen2, seq[list2[i]]+posin12, gapshift );
			mlen0 += gapshift;
			mlen1 += gapshift;
			mlen2 += gapshift;
#if 0
			for( i=0; i<ngroup0; i++ ) fprintf( stderr, "### mseq0[%d] = %s\n", i, mseq0[i] );
			for( i=0; i<ngroup1; i++ ) fprintf( stderr, "### mseq1[%d] = %s\n", i, mseq1[i] );
			for( i=0; i<ngroup2; i++ ) fprintf( stderr, "### mseq2[%d] = %s\n", i, mseq2[i] );
#endif

			if( gapshift ) 
			{
//				reporterr( "profilealignment (j=%d)!!!\n", j );

				profilealignment( ngroup0, ngroup1, ngroup2, mseq0, mseq1, mseq2, alloclen, alg );
			}

			j += gapshift;
			posin12 += gapshift;

			newpos = strlen( aseq[rep] ); // kufuu?
			for( i=0; i<ngroup0; i++ ) strcpy( aseq[list0[i]]+newpos, mseq0[i] );
			for( i=0; i<ngroup1; i++ ) strcpy( aseq[list1[i]]+newpos, mseq1[i] );
			for( i=0; i<ngroup2; i++ ) strcpy( aseq[list2[i]]+newpos, mseq2[i] );

//			fprintf( stderr, "gapshift = %d\n", gapshift );
		}
		blocklen = 1 + countnogaplen( gaplen+j+1, gaplen+len0 );
//		fprintf( stderr, "\nj=%d, blocklen=%d, len0=%d\n", j, blocklen, len0 );
//		blocklen = 1;

		newpos = strlen( aseq[rep] );

#if 0
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos] = seq[list0[i]][j];
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos] = seq[list1[i]][posin12];
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos] = seq[list2[i]][posin12];
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos+1] = 0;
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos+1] = 0;
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos+1] = 0;
#else

		for( i=0; i<ngroup0; i++ )
		{
			lp = list0[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+j, blocklen );
		}
		for( i=0; i<ngroup1; i++ )
		{
			lp = list1[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+posin12, blocklen );
		}
		for( i=0; i<ngroup2; i++ )
		{
			lp = list2[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+posin12, blocklen );
		}
//		fprintf( stderr, "### aseq[l0] = %s\n", aseq[list0[0]] );
//		fprintf( stderr, "### aseq[l1] = %s\n", aseq[list1[0]] );
//		fprintf( stderr, "### aseq[l2] = %s\n", aseq[list2[0]] );
#endif

//		fprintf( stderr, "j=%d -> %d\n", j, j+blocklen-1 );
		j += (blocklen-1);


		posin12 += (blocklen-1);


		posin12++;
	}
#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<ngroup0; i++ ) fprintf( stderr, " seq[l0i] = %s\n", seq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) fprintf( stderr, " seq[l1i] = %s\n", seq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) fprintf( stderr, " seq[l2i] = %s\n", seq[list2[i]] );
	fprintf( stderr, "=====>\n" );
	for( i=0; i<ngroup0; i++ ) fprintf( stderr, "aseq[l0i] = %s\n", aseq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) fprintf( stderr, "aseq[l1i] = %s\n", aseq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) fprintf( stderr, "aseq[l2i] = %s\n", aseq[list2[i]] );
#endif

//	for( i=0; i<njob; i++ ) if( mar[i] != 3 ) strcpy( seq[i], aseq[i] );
	for( i=0; i<ngroup0; i++ ) strcpy( seq[list0[i]], aseq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) strcpy( seq[list1[i]], aseq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) strcpy( seq[list2[i]], aseq[list2[i]] );


	free( mar );
	free( gaps );
	free( list0 );
	free( list1 );
	free( list2 );
	FreeCharMtx( mseq2 );
	FreeCharMtx( mseq1 ); // ? added 2012/02/12
	FreeCharMtx( mseq0 );
	FreeCharMtx( aseq ); // ? added 2012/02/12
}

void insertnewgaps_bothorders( int njob, int *alreadyaligned, char **seq, int *ex1, int *ex2, int *gaplen, int *gapmap, int gapmaplen, int alloclen, char alg, char gapchar )
{
//	int disp = 0;
	int *mar;
	char *gaps;
	char *cptr;
	int i, j, k, len, rep, len0, lp, blocklen, blockmap;
	char **mseq2, **mseq0, **mseq1;
	char **aseq, *newchar;
	int ngroup2, ngroup0, ngroup1;
	int *list0, *list1, *list2;
	int posin12, gapshifta, gapshiftn, gapshiftb, newpos;
	int mlen1, mlen0, mlen2;
	int jinc;

#if 0 // for debug
	int bug = 0;
	char *nogapseq1, *nogapseq2;
	nogapseq1 = calloc( alloclen, sizeof( char ) );
	nogapseq2 = calloc( alloclen, sizeof( char ) );
#endif // for debug

	mar = calloc( njob, sizeof( int ) );
	list0 = calloc( njob, sizeof( int ) );
	list1 = calloc( njob, sizeof( int ) );
	list2 = calloc( njob, sizeof( int ) );

	for( i=0; i<njob; i++ ) mar[i] = 0;
	for( i=0; i<njob; i++ ) 
	{
		if( alreadyaligned[i]==0 ) mar[i] = 3;
	}
	for( i=0; (k=ex1[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex1[i] );
	}
	for( i=0; (k=ex2[i])>-1; i++ ) 
	{
		mar[k] = 2;
//		fprintf( stderr, "excluding %d\n", ex2[i] );
	}

	ngroup2 = ngroup1 = ngroup0 = 0;
	for( i=0; i<njob; i++ )
	{
		if( mar[i] == 2 ) 
		{
			list2[ngroup2] = i;
			ngroup2++;
		}
		if( mar[i] == 1 ) 
		{
			list1[ngroup1] = i;
			ngroup1++;
		}
		if( mar[i] == 0 ) 
		{
			list0[ngroup0] = i;
//			fprintf( stderr, "inserting new gaps to %d\n", i );
			ngroup0++;
		}
	}
	list0[ngroup0] = list1[ngroup1] = list2[ngroup2] = -1;
	if( ngroup0 == 0 )
	{
//		fprintf( stderr, "Nothing to do\n" );
		free( mar );
		free( list0 );
		free( list1 );
		free( list2 );
		return;
	}

	for( i=0; i<njob; i++ ) if( mar[i] == 0 ) break;
	rep = i;
	len = strlen( seq[rep] );
	len0 = len+1;

//	reporterr( "alloclen = %d\n", alloclen );
//	reporterr( "len0 = %d\n", strlen( seq[list0[0]] ) );
//	reporterr( "len1 = %d\n", strlen( seq[list1[0]] ) );
//	reporterr( "gapmaplen = %d\n", gapmaplen );
//	reporterr( "ng0, ng1, ng2 = %d, %d, %d\n", ngroup0, ngroup1, ngroup2 );

//
//	if( i == njob )
//	{
////		fprintf( stderr, "Nothing to do\n" );
//		free( mar );
//		return;
//	}

	mseq2 = AllocateCharMtx( ngroup2, alloclen );
	mseq1 = AllocateCharMtx( ngroup1, alloclen );
	mseq0 = AllocateCharMtx( ngroup0, alloclen );
	aseq = AllocateCharMtx( njob, alloclen );
	gaps = calloc( alloclen, sizeof( char ) );

	for( i=0; i<njob; i++ ) aseq[i][0] = 0;
	newpos = 0;
	posin12 = 0;

#if 0
	if( disp )
	{
		int p;
		fprintf( stderr, "\n" );
		fprintf( stderr, "len0 = %d\n", len0 );
		fprintf( stderr, "\ngaplen[] = \n" );
		reporterr( "seq0[0]  = %s\n", seq[list0[0]] );
		reporterr( "seq1[1]  = %s\n", seq[list1[0]] );
		reporterr( "seq2[2]  = %s\n", seq[list2[0]] );
		reporterr( "seq[rep] = %s\n", seq[rep] );
		for(i=0,p=0; i<gapmaplen; i++,p++ )
		{
			fprintf( stderr, "gaplen %d:%d %-*.*s\n", p, gaplen[p], gaplen[p]+1, gaplen[p]+1, seq[list1[0]]+i );
			i += gaplen[p];
		}
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "\ngapmap[] = \n" );
	for(i=0; i<gapmaplen; i++ ) fprintf( stderr, "gapmap %d:%d %-*.*s\n", i, gapmap[i], gapmap[i]+1, gapmap[i]+1, seq[list1[0]]+i );

	reporterr( "seq1 = \n" );
	reporterr( "%s\n", seq[list1[0]] );
#endif

	

	for( j=0; j<len0; j++ )
	{
//		fprintf( stderr, "\nj=%d, gaplen[%d]=%d\n", j, j, gaplen[j] );
//		if( gaplen[j] || gapmap[posin12] )
		if( gaplen[j] || gapmap[posin12] )
		{
//			fprintf( stderr, "j=%d GAP!\n", j );
			for( i=0; i<ngroup0; i++ ) mseq0[i][0] = 0;
			for( i=0; i<ngroup1; i++ ) mseq1[i][0] = 0;
			for( i=0; i<ngroup2; i++ ) mseq2[i][0] = 0;
			mlen0 = mlen1 = mlen2 = 0;


			gapshiftb = gapmap[posin12];
//			fprintf( stderr, "\ngapmap[%d] kouho = %d, posint12 = %d\n", posin12, gapmap[posin12], posin12 );

			if( gapshiftb ) // koko ga osoi!
			{
				for( i=0; i<ngroup0; i++ ) strncpy0( mseq0[i]+mlen0, seq[list0[i]]+j, gapshiftb ); // tokuni koko!
				for( i=0; i<ngroup1; i++ ) strncpy0( mseq1[i]+mlen1, seq[list1[i]]+posin12, gapshiftb );
				for( i=0; i<ngroup2; i++ ) strncpy0( mseq2[i]+mlen2, seq[list2[i]]+posin12, gapshiftb );
				mlen0 += gapshiftb;
				mlen1 += gapshiftb;
				mlen2 += gapshiftb;
#if 0
				reporterr( "\n\n" );
				for( i=0; i<1; i++ ) fprintf( stderr, "##b mseq0[%d] = %s\n", i, mseq0[i] );
				for( i=0; i<1; i++ ) fprintf( stderr, "##b mseq1[%d] = %s\n", i, mseq1[i] );
				for( i=0; i<1; i++ ) fprintf( stderr, "##b mseq2[%d] = %s\n", i, mseq2[i] );
#endif
				posin12 += gapshiftb;
			}
#if 0 // nen no tame
			for( i=0, jinc=0; i<gapshiftb; i++ ) jinc += 1+gaplen[j+i];
			if( jinc != gapshiftb )
			{
				reporterr( "\n#################!\n" );
				reporterr( "# Unexpected gap pattern!\n" );
				reporterr( "# There are overlapped %d gaps in gaplen[] and gapmap[]. j=%d-%d, posin12=%d-%d\n", jinc, j-gapshiftb-jinc, j, posin12-gapshiftb, posin12 );
				reporterr( "\n#################!\n" );
				exit( 1 );
			}
			j += gapshiftb;
#else
			j += gapshiftb;
#endif

#if 0
			if( disp && gapshiftb )
			{
				reporterr( "after gapshiftb, j=%d, posin12=%d\n", j, posin12 );
				reporterr( "mseq0[0] = %s\n", mseq0[0] );
				reporterr( "mseq1[0] = %s\n", mseq1[0] );
				reporterr( "mseq2[0] = %s\n", mseq2[0] );
			}
#endif
	
//			fprintf( stderr, "gaplen[%d]=%d, posin12 = %d\n", j, gaplen[j], posin12 );

			while( 1 )
			{
				gapshiftn = gaplen[j];
				if( gapshiftn )
				{
					cptr = gaps;
					while( gapshiftn-- ) *cptr++ = gapchar;
					*cptr = 0;
					gapshiftn = gaplen[j];
		
					for( i=0; i<ngroup0; i++ ) strncpy0( mseq0[i]+mlen0, gaps, gapshiftn );
					for( i=0; i<ngroup1; i++ ) strncpy0( mseq1[i]+mlen1, seq[list1[i]]+posin12, gapshiftn );
					for( i=0; i<ngroup2; i++ ) strncpy0( mseq2[i]+mlen2, seq[list2[i]]+posin12, gapshiftn );
					posin12 += gapshiftn;
					mlen0 += gapshiftn;
					mlen1 += gapshiftn;
					mlen2 += gapshiftn;
#if 0
					for( i=0; i<1; i++ ) fprintf( stderr, "##n mseq0[%d] = %s\n", i, mseq0[i] );
					for( i=0; i<1; i++ ) fprintf( stderr, "##n mseq1[%d] = %s\n", i, mseq1[i] );
					for( i=0; i<1; i++ ) fprintf( stderr, "##n mseq2[%d] = %s\n", i, mseq2[i] );
#endif
				}
	
#if 0
				if( disp && gapshiftn )
				{
					reporterr( "after gapshiftn (j=%d, gaplen[j]=%d, posin12=%d, gapshiftn=%d)\n", j, gaplen[j], posin12-gapshiftn, gapshiftn );
					reporterr( "mseq0[0] = %s\n", mseq0[0] );
					reporterr( "mseq1[0] = %s\n", mseq1[0] );
					reporterr( "mseq2[0] = %s\n", mseq2[0] );
				}
#endif
	
	
				gapshifta = gapmap[posin12];
//				fprintf( stderr, "gapmap[%d] kouho = %d, posin12 = %d\n", posin12, gapmap[posin12], posin12 );
	
				if( gapshifta )
				{
					for( i=0; i<ngroup0; i++ ) strncpy0( mseq0[i]+mlen0, seq[list0[i]]+j, gapshifta );
					for( i=0; i<ngroup1; i++ ) strncpy0( mseq1[i]+mlen1, seq[list1[i]]+posin12, gapshifta );
					for( i=0; i<ngroup2; i++ ) strncpy0( mseq2[i]+mlen2, seq[list2[i]]+posin12, gapshifta );
					mlen0 += gapshifta;
					mlen1 += gapshifta;
					mlen2 += gapshifta;
					posin12 += gapshifta;

#if 0
					for( i=0; i<1; i++ ) fprintf( stderr, "##a mseq0[%d] = %s\n", i, mseq0[i] );
					for( i=0; i<1; i++ ) fprintf( stderr, "##a mseq1[%d] = %s\n", i, mseq1[i] );
					for( i=0; i<1; i++ ) fprintf( stderr, "##a mseq2[%d] = %s\n", i, mseq2[i] );
#endif
				}
#if 0
				j += gapshifta; // BUG!!
#else
//				for( i=1, jinc=1; i<gapshifta; i++ ) jinc += 1+gaplen[j+i];
				for( i=1, jinc=0; i<=gapshifta; i++ ) jinc += 1+gaplen[j+i];
//				j += jinc;
				j += gapshifta;
	
				if( jinc == gapshifta ) break;

				reporterr( "(a) There are overlapped %d gaps in gaplist[] and gapmap[]. j=%d-%d, posin12=%d-%d, jinc=%d, gapshifta=%d\n", jinc, j-gapshifta-jinc, j, posin12-gapshifta, posin12, jinc, gapshifta );
#endif
			}

#if 0
			if( disp && gapshifta )
			{
				reporterr( "after gapshifta, j=%d, posin12=%d\n", j, posin12 );
				reporterr( "mseq0[0] = %s\n", mseq0[0] );
				reporterr( "mseq1[0] = %s\n", mseq1[0] );
				reporterr( "mseq2[0] = %s\n", mseq2[0] );
			}
#endif


			if( gapshiftb + gapshifta ) 
			{
#if 0
				for( i=0; i<1; i++ ) fprintf( stderr, "### mseq0[%d] = %s\n", i, mseq0[i] );
				for( i=0; i<1; i++ ) fprintf( stderr, "### mseq1[%d] = %s\n", i, mseq1[i] );
				for( i=0; i<1; i++ ) fprintf( stderr, "### mseq2[%d] = %s\n", i, mseq2[i] );
#endif
//				if( disp ) reporterr( "profilealignment (j=%d)!!!\n", j );

				profilealignment( ngroup0, ngroup1, ngroup2, mseq0, mseq1, mseq2, alloclen, alg );
			}


			newpos = strlen( aseq[rep] ); // kufuu?
			for( i=0; i<ngroup0; i++ ) strcpy( aseq[list0[i]]+newpos, mseq0[i] );
			for( i=0; i<ngroup1; i++ ) strcpy( aseq[list1[i]]+newpos, mseq1[i] );
			for( i=0; i<ngroup2; i++ ) strcpy( aseq[list2[i]]+newpos, mseq2[i] );

#if 0
			if( disp )
			{
				reporterr( "after profilealignment\n" );
				reporterr( "mseq0[0] = %s\n", mseq0[0] );
				reporterr( "mseq1[0] = %s\n", mseq1[0] );
				reporterr( "mseq2[0] = %s\n", mseq2[0] );

				gappick0equalminustmptmptmp( nogapseq1, aseq[list1[0]] );
				gappick0equalminustmptmptmp( nogapseq2, seq[list1[0]] );
	
				reporterr( "aseq[list1[0].nogap = %s\n", nogapseq1 );
				reporterr( " seq[list1[0].nogap = %s\n", nogapseq2 );
			}
#endif

//			fprintf( stderr, "gapshift = %d\n", gapshift );
		}
		newpos = strlen( aseq[rep] );
		blocklen = 1 + countnogaplen( gaplen+j+1, gaplen+len0 );
//		fprintf( stderr, "\nj=%d, blocklen=%d, len0=%d\n", j, blocklen, len0 );

		blockmap = 1 + countgapmap( gapmap+posin12+1, gapmap+gapmaplen );
//		fprintf( stderr, "posin12=%d, blockmap=%d, len0=%d\n", posin12, blockmap, len0 );

//		if( disp ) reporterr( "newpos = %d, blocklen = %d, blockmap = %d, j=%d, posin12=%d\n", newpos, blocklen, blockmap, j, posin12 );

		if( blockmap < blocklen ) blocklen = blockmap;




#if 0
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos] = seq[list0[i]][j];
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos] = seq[list1[i]][posin12];
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos] = seq[list2[i]][posin12];
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos+1] = 0;
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos+1] = 0;
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos+1] = 0;
#else

//		if( j >= len0 ) break; // iru?


		for( i=0; i<ngroup0; i++ )
		{
			lp = list0[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+j, blocklen );
//			*(newchar+blocklen) = 0; iranai
		}
		for( i=0; i<ngroup1; i++ )
		{
			lp = list1[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+posin12, blocklen );
//			*(newchar+blocklen) = 0; iranai
		}
		for( i=0; i<ngroup2; i++ )
		{
			lp = list2[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+posin12, blocklen );
//			*(newchar+blocklen) = 0; iranai
		}


//			reporterr( "adding %c to aseq[list1[0]]\n", seq[list1[0]][posin12] );


//		for( i=0; i<ngroup0; i++ ) fprintf( stderr, "### aseq0[%d] = %s\n", i, aseq[list0[i]] );
//		for( i=0; i<ngroup1; i++ ) fprintf( stderr, "### aseq1[%d] = %s\n", i, aseq[list1[i]] );
//		for( i=0; i<ngroup2; i++ ) fprintf( stderr, "### aseq2[%d] = %s\n", i, aseq[list2[i]] );
#endif

//		fprintf( stderr, "j=%d -> %d\n", j, j+blocklen-1 );

		j += (blocklen-1);
//		j += gaplen[j];


		posin12 += (blocklen-1); // sono aida ni gapmap wo miotosu?


		posin12++;

#if 0
		if( disp )
		{
			gappick0equalminustmptmptmp( nogapseq1, aseq[list1[0]] );
			gappick0equalminustmptmptmp( nogapseq2, seq[list1[0]] );

			reporterr( "aseq[list1[0].nogap = %s\n", nogapseq1 );
			reporterr( " seq[list1[0].nogap = %s\n", nogapseq2 );
			reporterr( "" );
//			reporterr( "seq[list1[0]] = %s\n", seq[list1[0]] );
//			reporterr( "seq[list2[0]] = %s\n", seq[list2[0]] );
		}
#endif
	}
#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<ngroup0; i++ ) fprintf( stderr, " seq[l0i] = \n%s\n", seq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) fprintf( stderr, " seq[l1i] = \n%s\n", seq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) fprintf( stderr, " seq[l2i] = \n%s\n", seq[list2[i]] );
	reporterr( "0         1         2         3         4         5         6         7       \n" );
	reporterr( "012345678901234567890123456789012345678901234567890123456789012345678901234567\n" );
	fprintf( stderr, "=====>\n" );
	for( i=0; i<ngroup0; i++ ) fprintf( stderr, "aseq[l0i] = \n%s\n", aseq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) fprintf( stderr, "aseq[l1i] = \n%s\n", aseq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) fprintf( stderr, "aseq[l2i] = \n%s\n", aseq[list2[i]] );
#endif

#if 0
	reporterr( "list0[0]=%d\n", list0[0] );
	reporterr( "list0[%d-1]=%d\n", ngroup0, list0[ngroup0-1] );
	reporterr( "seq0[list[0]]=%s\n", seq[list0[0]] );
	reporterr( "aseq0[list[0]]=%s\n", aseq[list0[0]] );
#endif



#if 0
	for( i=0; i<ngroup0; i++ )
	{
		if( strlen( aseq[list0[0]] ) != strlen( aseq[list0[i]] ) )
		{
			reporterr( "Length error! len[0] = %d, but len[%d] = %d\n", strlen( aseq[list0[0]] ), list0[i], strlen( aseq[list0[i]] ) );
			bug = 1;
			break;
		}
	}
	for( i=0; i<ngroup0; i++ )
	{
		gappick0equalminustmptmptmp( nogapseq1, aseq[list0[i]] );
		gappick0equalminustmptmptmp( nogapseq2, seq[list0[i]] );
		if( strcmp( nogapseq1, nogapseq2 ) ) bug = 1;
	}
	for( i=0; i<ngroup1; i++ )
	{
		gappick0equalminustmptmptmp( nogapseq1, aseq[list1[i]] );
		gappick0equalminustmptmptmp( nogapseq2, seq[list1[i]] );
		if( strcmp( nogapseq1, nogapseq2 ) ) bug = 1;
	}
	for( i=0; i<ngroup2; i++ )
	{
		gappick0equalminustmptmptmp( nogapseq1, aseq[list2[i]] );
		gappick0equalminustmptmptmp( nogapseq2, seq[list2[i]] );
		if( strcmp( nogapseq1, nogapseq2 ) ) bug = 1;
	}

	free( nogapseq1 );
	free( nogapseq2 );
	
	if( bug )
	{
		reporterr( "ERROR!!!!!!!\n" );
		reporterr( ">aseq1[%d], len = %d\n%s\n", list1[0], strlen( aseq[list1[0]] ), aseq[list0[0]] );
		reporterr( ">seq1[%d], len = %d\n%s\n", list1[0], strlen( seq[list1[0]] ), seq[list0[i]] );
		exit( 1 );


		for( i=0; i<ngroup0; i++ ) reporterr( ">aseq0[%d], len = %d\n%s\n", list0[i], strlen( aseq[list0[i]] ), aseq[list0[i]] );
		for( i=0; i<ngroup1; i++ ) reporterr( ">aseq1[%d], len = %d\n%s\n", list1[i], strlen( aseq[list1[i]] ), aseq[list1[i]] );
		for( i=0; i<ngroup2; i++ ) reporterr( ">aseq2[%d], len = %d\n%s\n", list2[i], strlen( aseq[list2[i]] ), aseq[list2[i]] );

		for( i=0; i<ngroup0; i++ ) reporterr( ">seq0[%d], len = %d\n%s\n", list0[i], strlen( seq[list0[i]] ), seq[list0[i]] );
		for( i=0; i<ngroup1; i++ ) reporterr( ">seq1[%d], len = %d\n%s\n", list1[i], strlen( seq[list1[i]] ), seq[list1[i]] );
		for( i=0; i<ngroup2; i++ ) reporterr( ">seq2[%d], len = %d\n%s\n", list2[i], strlen( seq[list2[i]] ), seq[list2[i]] );
		exit( 1 );
	}

#endif

//	for( i=0; i<njob; i++ ) if( mar[i] != 3 ) strcpy( seq[i], aseq[i] );
	for( i=0; i<ngroup0; i++ ) strcpy( seq[list0[i]], aseq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) strcpy( seq[list1[i]], aseq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) strcpy( seq[list2[i]], aseq[list2[i]] );

	free( mar );
	free( gaps );
	free( list0 );
	free( list1 );
	free( list2 );
	FreeCharMtx( mseq2 );
	FreeCharMtx( mseq1 ); // ? added 2012/02/12
	FreeCharMtx( mseq0 );
	FreeCharMtx( aseq ); // ? added 2012/02/12

}


static void reflectsmoothing( char *ref, int *mem, char **seq, int len )
{
	char *tmpseq;
	int i, j, k, p;

//	reporterr( "#### reflectsmoothing!!!!!\n" );

//	if( mem[1] != -1 ) reporterr( "original = %s\n", seq[mem[1]] );

	tmpseq = calloc( len+1, sizeof( char ) );


	for( j=1; (i=mem[j])!=-1; j++ )
	{
		eqpick( tmpseq, seq[i] );
		for( k=0, p=0; p<len; k++ )
		{
			while( ref[p] == '=' ) seq[i][p++] = '=';
			seq[i][p++] = tmpseq[k];
		}
	}
	free( tmpseq );

//	if( mem[1] != -1 ) reporterr( "output   = %s\n", seq[mem[1]] );

//	reporterr( "#### done!!!!!\n" );
}

static int smoothing1rightmulti( int len, char *ref ) // osoi!
{
	int i, j, k;
	int shiftfrom = -1;
	int shiftto = -1;
	int *hit;
	int val = 0, nhit = 0;

	hit = NULL;

//	reporterr( "ref (1rightmulti) = %s\n", ref );

	for( i=1, nhit=0; i<len-1; i++ ) // break nashi no baai, hidarihaji ha saigo
//	for( i=len-2; i>0; i-- ) // break ari no baai, migihajiha saigo
	{
		if( ref[i-1] == '+' && ( ref[i] != '+' && ref[i] != '=' ) && ref[i+1] == '=' )
		{
//			reporterr( "hit! i=%d, len=%d\n", i, len );
			hit = realloc( hit, (nhit+1) * sizeof( int ) );
			hit[nhit] = i;
			nhit += 1;
//			break;
		}
	}
	if( nhit == 0 ) return( 0 );


	for( k=0; k<nhit; k++ )
	{
		for( j=hit[k]+1; j<=len; j++ )
		{
			if( ref[j] != '=' )
			{
				shiftto = j-1;
				break;
			}
		}
		if( j == len && ref[len-1] == '=' )
		{
			reporterr( "hit[i].end = %d, j = len-1, skip!\n" );
			continue;
		}

		if( shiftto < len-1 && ref[shiftto+1] == '+' ) continue; // muda dakara

		val += 1;
		shiftfrom = hit[k];
		if( ref[shiftto] != '=' ) // atode sakujo 
		{
			reporterr( "Error in smoothing1left!\n" );
			exit( 1 );
		}
		ref[shiftto] = ref[shiftfrom];
		ref[shiftfrom] = '=';
	}
	free( hit );

//	reporterr( "ref (1rightmulti) = %s\n", ref );
	reporterr( " %d out of %d have been smoothed (right).\n", val, nhit );

//	if( nhit > 1 ) exit( 1 );
	return( val );
}

static int smoothing1leftmulti( int len, char *ref ) // osoi!
{
	int i, j, k;
	int shiftfrom = -1;
	int shiftto = -1;
	int *hit;
	int val = 0, nhit = 0;

	hit = NULL;

//	reporterr( "ref (1leftmulti) = %s\n", ref );

	for( i=1, nhit=0; i<len-1; i++ ) // break nashi no baai, hidarihaji ha saigo
//	for( i=len-2; i>0; i-- ) // break ari no baai, migihajiha saigo
	{
		if( ref[i-1] == '=' && ( ref[i] != '+' && ref[i] != '=' ) && ref[i+1] == '+' )
		{
//			reporterr( "hit! i=%d, len=%d\n", i, len );
			hit = realloc( hit, (nhit+1) * sizeof( int ) );
			hit[nhit] = i;
			nhit += 1;
//			break;
		}
	}
	if( nhit == 0 ) return( 0 );

	for( k=0; k<nhit; k++ )
	{
		for( j=hit[k]-1; j>-1; j-- )
		{
			if( ref[j] != '=' )
			{
				shiftto = j+1;
				break;
			}
		}
		if( j == -1 && ref[0] == '=' )
		{
			reporterr( "hit[i].end = %d, j = -1, skip!\n" );
			continue;
		}

		if( shiftto > 0 && ref[shiftto-1] == '+' ) continue; // muda dakara

		val += 1;
		shiftfrom = hit[k];
		if( ref[shiftto] != '=' ) // atode sakujo 
		{
			reporterr( "Error in smoothing1left!\n" );
			exit( 1 );
		}
		ref[shiftto] = ref[shiftfrom];
		ref[shiftfrom] = '=';
	
	}
	free( hit );

//	reporterr( "ref (1leftmulti) = %s\n", ref );
	reporterr( " %d out of %d have been smoothed (left).\n", val, nhit );


//	if( nhit > 1 ) exit( 1 );
	return( val );
}

void restorecommongapssmoothly( int njob, int n0, char **seq, int *ex1, int *ex2, int *gapmap, int alloclen, char gapchar )
{
	int *mem;
	char *tmpseq;
	char *cptr;
	int *iptr;
	int *tmpgapmap;
	int i, j, k, len, rep1, rep2, len1, klim, leninserted;
	int totalres;
	
	if( n0 == 0 ) return;


	mem = calloc( njob+1, sizeof( int ) ); // +1 ha iranai.
	intcpy( mem, ex1 );
	intcat( mem, ex2 );
//	tmpseq = calloc( alloclen+2, sizeof( char ) );
//	tmpgapmap = calloc( alloclen+2, sizeof( int ) );

#if 0 // iranai
	for( i=0; (k=mem[i])!=-1; i++ ) // iranai
		reporterr( "mem[%d] = %d\n", i, k ); // iranai
	if( i == njob ) // iranai
	{
		fprintf( stderr, "Error in restorecommongaps()\n" );
		free( mem );
		exit( 1 );
	}
#endif
	rep1 = ex1[0];
	rep2 = ex2[0];
	len = strlen( seq[rep1] );
	len1 = len+1;

	tmpseq = calloc( alloclen, sizeof( char ) );
	tmpgapmap = calloc( alloclen, sizeof( int ) );

#if 0
	reporterr( "\n" );
	reporterr( "seq[rep1] = %s\n", seq[rep1] );
	reporterr( "seq[rep2] = %s\n", seq[rep2] );
#endif

	for( k=0; (i=mem[k])!=-1; k++ )
	{
		cptr = tmpseq;
		for( j=0; j<len1; j++ )
		{
			klim = gapmap[j];
//			for( k=0; k<gapmap[j]; k++ )
			while( klim-- )
				*(cptr++) = '+'; // ???
			*(cptr++) = seq[i][j];
		}
		*cptr = 0;
		strcpy( seq[i], tmpseq );
	}
#if 0
	reporterr( "->\n" );
	reporterr( "seq[rep1] = \n%s\n", seq[rep1] );
	reporterr( "seq[rep2] = \n%s\n", seq[rep2] );
#endif

	leninserted = strlen( seq[rep1] );
#if 0
	reporterr( "gapmap =\n" );
	for(j=0; j<len1; j++) 
	{
		reporterr( "%d", gapmap[j] );
		for( i=gapmap[j]; i>0; i-- ) reporterr( "-" );
	}
	reporterr( "\n" );
#endif

#if 0
	resprev = 10000; // tekitou
	while( 1 )
	{
		res = 0;
//		reporterr( "\nsmoothing1right..\n" );
		res  = (0<smoothing1right( leninserted, seq[rep1], gapmap, seq, ex1 ));
//		reporterr( "done. res = %d\n", res );
//		reporterr( "smoothing1right..\n" );
		res += (0<smoothing1right( leninserted, seq[rep2], gapmap, seq, ex2 ));
//		reporterr( "done. res = %d\n", res );

//		reporterr( "smoothing1left..\n" );
		res += (0<smoothing1left( leninserted, seq[rep1], gapmap, seq, ex1 ));
//		reporterr( "done. res = %d\n", res );
//		reporterr( "smoothing1left..\n" );
		res += (0<smoothing1left( leninserted, seq[rep2], gapmap, seq, ex2 ));
//		reporterr( "done. res = %d\n", res );

		reporterr( " Smoothing .. %d \n", res );
		if( res >= resprev ) break;
//		if( res == 0 ) break;
		resprev = res;
	}
#else
	totalres = 0;
	totalres += smoothing1rightmulti( leninserted, seq[rep1] );
	totalres += smoothing1leftmulti( leninserted, seq[rep1] );
	if( totalres ) reflectsmoothing( seq[rep1], ex1, seq, leninserted );

	totalres = 0;
	totalres += smoothing1rightmulti( leninserted, seq[rep2] );
	totalres += smoothing1leftmulti( leninserted, seq[rep2] );
	if( totalres ) reflectsmoothing( seq[rep2], ex2, seq, leninserted );
#endif

	for( k=0; (i=mem[k])!=-1; k++ ) plus2gapchar( seq[i], gapchar );

#if 0
	reporterr( "->\n" );
	reporterr( "seq[rep1] = \n%s\n", seq[rep1] );
	reporterr( "seq[rep2] = \n%s\n", seq[rep2] );
	reporterr( "gapmap =\n" );
	for(j=0; j<len1; j++) 
	{
		reporterr( "%d", gapmap[j] );
		for( i=gapmap[j]; i>0; i-- ) reporterr( "-" );
	}
	reporterr( "\n" );
#endif

	iptr = tmpgapmap;
	for( j=0; j<len1; j++ )
	{
		*(iptr++) = gapmap[j];
		for( k=0; k<gapmap[j]; k++ )
			*(iptr++) = 0;
	}
	*iptr = -1;

	intcpy( gapmap, tmpgapmap );
//	iptr = tmpgapmap;
//	while( *iptr != -1 ) *gapmap++ = *iptr++;

	free( mem );
	free( tmpseq );
	free( tmpgapmap );
}

void restorecommongaps( int njob, int n0, char **seq, int *ex1, int *ex2, int *gapmap, int alloclen, char gapchar )
{
	int *mem;
	char *tmpseq;
	char *cptr;
	int *iptr;
	int *tmpgapmap;
	int i, j, k, len, rep, len1, klim;
	

	if( n0 == 0 ) return;


	mem = calloc( njob+1, sizeof( int ) ); // +1 ha iranai.
	intcpy( mem, ex1 );
	intcat( mem, ex2 );
//	tmpseq = calloc( alloclen+2, sizeof( char ) );
//	tmpgapmap = calloc( alloclen+2, sizeof( int ) );

#if 0 // iranai
	for( i=0; (k=mem[i])!=-1; i++ ) // iranai
		reporterr( "mem[%d] = %d\n", i, k ); // iranai
	if( i == njob ) // iranai
	{
		fprintf( stderr, "Error in restorecommongaps()\n" );
		free( mem );
		exit( 1 );
	}
#endif
	rep = mem[0];
	len = strlen( seq[rep] );
	len1 = len+1;

	tmpseq = calloc( alloclen, sizeof( char ) );
	tmpgapmap = calloc( alloclen, sizeof( int ) );



	for( k=0; (i=mem[k])!=-1; k++ )
	{
		cptr = tmpseq;
		for( j=0; j<len1; j++ )
		{
			klim = gapmap[j];
//			for( k=0; k<gapmap[j]; k++ )
			while( klim-- )
				*(cptr++) = gapchar; // ???
			*(cptr++) = seq[i][j];
		}
		*cptr = 0;
		strcpy( seq[i], tmpseq );
	}

	iptr = tmpgapmap;
	for( j=0; j<len1; j++ )
	{
		*(iptr++) = gapmap[j];
		for( k=0; k<gapmap[j]; k++ )
			*(iptr++) = 0;
	}
	*iptr = -1;

	iptr = tmpgapmap;
	while( *iptr != -1 ) *gapmap++ = *iptr++;

	free( mem );
	free( tmpseq );
	free( tmpgapmap );
}

#if 0 // mada
int deletenewinsertions_difflist( int on, int an, char **oseq, char **aseq, GapPos **difflist )
{
	int i, j, p, q, allgap, ndel, istart, pstart;
	int len = strlen( oseq[0] );
	char *eqseq, tmpc;

	reporterr( "In deletenewinsertions_difflist\n" );
	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );


	eqseq = calloc( len+1, sizeof( char ) );
	for( i=0; i<len; i++ )
	{
		allgap = 0;
		for( j=0; j<on; j++ )
		{
			tmpc = oseq[j][i];
			if( tmpc != '-' && tmpc != '=' ) break;
		}
		if( j == on ) 
			allgap = 1;

		if( allgap )
		{
			eqseq[i] = '=';
		}
		else
		{
			eqseq[i] = 'o';
		}
	}

	for( j=0; j<1; j++ ) reporterr( "\no                = %s\n", oseq[j] );
	reporterr( "\ne                = %s\n", eqseq );
	for( j=0; j<1; j++ ) reporterr( "a                = %s\n", aseq[j] );

	if( difflist )
	{
		for( j=0; j<an; j++ )
		{
			ndel = 0;
			for( i=0,q=0,p=0; i<len; i++ ) // 0 origin
			{
				tmpc = aseq[j][i];
				if( eqseq[i] != '=' ) p++;
				if( tmpc != '-' && tmpc != '=' ) q++;

#if 0
				if( eqseq[i] == '=' && ( tmpc != '-' && tmpc != '=' ) )
				{
					difflist[j] = realloc( difflist[j], sizeof( GapPos ) * (ndel+2) );
					difflist[j][ndel].pos = -p;
					difflist[j][ndel].len = 1;
					ndel++;
				}
				if( eqseq[i] != '=' && ( tmpc == '-' || tmpc == '=' ) ) 
				{
//					reporterr( "deleting %d-%d, %c\n", j, i, aseq[j][i] );
					difflist[j] = realloc( difflist[j], sizeof( GapPos ) * (ndel+2) );
					difflist[j][ndel].pos = p;
					difflist[j][ndel].len = 1;
					ndel++;
				}
#else
				istart = i;
				pstart = p;
				while( eqseq[i] == '=' && ( tmpc != '-' && tmpc != '=' ) )
				{
					i++;
					tmpc = aseq[j][i];
				}
				if( i != istart )
				{
					difflist[j] = realloc( difflist[j], sizeof( GapPos ) * (ndel+2) );
					difflist[j][ndel].pos = pstart;
					difflist[j][ndel].len = -(i-istart);
					ndel++;
				}
				istart = i;
				pstart = p;
				while( eqseq[i] != '=' && ( tmpc == '-' || tmpc == '=' ) ) 
				{
					i++;
					tmpc = aseq[j][i];
				}
				if( i != istart )
				{
					difflist[j] = realloc( difflist[j], sizeof( GapPos ) * (ndel+2) );
					difflist[j][ndel].pos = pstart;
					difflist[j][ndel].len = i-istart;
					ndel++;
				}
#endif
			}
			difflist[j][ndel].pos = -1;
			difflist[j][ndel].len = 0;
			for( i=0; ; i++ )
			{
				if( difflist[j][i].pos == -1 ) break;
				reporterr( "sequence %d,%d, pos=%d,len=%d\n", j, i, difflist[j][i].pos, difflist[j][i].len );
			}
		}

	}
exit( 1 );
	for( i=0,p=0; i<len; i++ )
	{

//		if( oseq[0][i] != '=' )
//		reporterr( "i=%d, p=%d, q=%d, originally, %c\n", i, p, q, originallygapped[p]);
//		if( eqseq[i] != '=' && originallygapped[p] != '-' ) // dame!!
		if( eqseq[i] != '=' )
		{
//			reporterr( "COPY! p=%d\n", p );
			if( p != i )
			{
				for( j=0; j<on; j++ ) oseq[j][p] = oseq[j][i];
				for( j=0; j<an; j++ ) aseq[j][p] = aseq[j][i];
			}
			p++;
		}
	}
//		reporterr( "deletemap        = %s\n", deletemap );
//		reporterr( "eqseq            = %s\n", eqseq );
//		reporterr( "originallygapped = %s\n", originallygapped );
	for( j=0; j<on; j++ ) oseq[j][p] = 0;
	for( j=0; j<an; j++ ) aseq[j][p] = 0;

	free( eqseq );

//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );

	return( i-p );
}
#endif

int deletenewinsertions_whole_eq( int on, int an, char **oseq, char **aseq, GapPos **deletelist )
{
	int i, j, p, q, allgap, ndel, qstart;
	int len = strlen( oseq[0] );
	char *eqseq, tmpc;

//	reporterr( "In deletenewinsertions_whole_eq\n" );
//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );

	eqseq = calloc( len+1, sizeof( char ) );
	for( i=0; i<len; i++ )
	{
		allgap = 0;
		for( j=0; j<on; j++ )
		{
			tmpc = oseq[j][i];
			if( tmpc != '-' && tmpc != '=' ) break;
		}
		if( j == on ) 
			allgap = 1;

		if( allgap )
		{
			eqseq[i] = '=';
		}
		else
		{
			eqseq[i] = 'o';
		}
	}
	eqseq[len] = 0; // hitsuyou

//	for( j=0; j<1; j++ ) reporterr( "\no                = %s\n", oseq[j] );
//	reporterr( "\ne                = %s\n", eqseq );
//	for( j=0; j<1; j++ ) reporterr( "a                = %s\n", aseq[j] );

	if( deletelist )
	{
		for( j=0; j<an; j++ )
		{
			ndel = 0;
			for( i=0,q=0; i<len; i++ )
			{
				tmpc = aseq[j][i];
#if 0
				if( tmpc != '-' && tmpc != '=' ) 
				{
					if( eqseq[i] == '=' )
					{
//						reporterr( "deleting %d-%d, %c\n", j, i, aseq[j][i] );
						deletelist[j] = realloc( deletelist[j], sizeof( GapPos ) * (ndel+2) );
						deletelist[j][ndel].pos = q;
						deletelist[j][ndel].len = 1;
						ndel++;
					}
					q++;
				}
#else
				qstart = q;
				while( (tmpc != '-' && tmpc != '=') && eqseq[i] == '=' )
				{
					i++;
					q++;
					tmpc = aseq[j][i];
				}
				if( qstart != q )
				{
//					reporterr( "deleting %d-%d, %c\n", j, i, aseq[j][i] );
					deletelist[j] = realloc( deletelist[j], sizeof( GapPos ) * (ndel+2) );
					deletelist[j][ndel].pos = qstart;
					deletelist[j][ndel].len = q-qstart;
					ndel++;
				}
				if( tmpc != '-' && tmpc != '=' ) q++;
#endif
			}
			deletelist[j][ndel].pos = -1;
			deletelist[j][ndel].len = 0;
		}
	}
	for( i=0,p=0; i<len; i++ )
	{

//		if( oseq[0][i] != '=' )
//		reporterr( "i=%d, p=%d, q=%d, originally, %c\n", i, p, q, originallygapped[p]);
//		if( eqseq[i] != '=' && originallygapped[p] != '-' ) // dame!!
		if( eqseq[i] != '=' )
		{
//			reporterr( "COPY! p=%d\n", p );
			if( p != i )
			{
				for( j=0; j<on; j++ ) oseq[j][p] = oseq[j][i];
				for( j=0; j<an; j++ ) aseq[j][p] = aseq[j][i];
			}
			p++;
		}
	}
//		reporterr( "deletemap        = %s\n", deletemap );
//		reporterr( "eqseq            = %s\n", eqseq );
//		reporterr( "originallygapped = %s\n", originallygapped );
	for( j=0; j<on; j++ ) oseq[j][p] = 0;
	for( j=0; j<an; j++ ) aseq[j][p] = 0;

	free( eqseq );

//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );

	return( i-p );
}

#if 1 // 2022/Jan. disttbfast to tbfast niha hitsuyou.
int deletenewinsertions_whole( int on, int an, char **oseq, char **aseq, GapPos **deletelist )
{
	int i, j, p, q, allgap, ndel, qstart;
	int len = strlen( oseq[0] );
	char *eqseq, tmpc;

//	reporterr( "In deletenewinsertions_whole\n" );
//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );

	eqseq = calloc( len+1, sizeof( char ) );
	for( i=0,p=0; i<len; i++ )
	{
		allgap = 0;
		for( j=0; j<on; j++ )
		{
			tmpc = oseq[j][i];
			if( tmpc != '-' ) break;
		}
		if( j == on ) 
			allgap = 1;

		if( allgap )
		{
			eqseq[i] = '=';
		}
		else
		{
			eqseq[i] = 'o';
		}
	}
	eqseq[len] = 0;

//	for( j=0; j<1; j++ ) reporterr( "\no                = %s\n", oseq[j] );
//	reporterr( "\ne                = %s\n", eqseq );
//	for( j=0; j<1; j++ ) reporterr( "a                = %s\n", aseq[j] );

	if( deletelist )
	{
		for( j=0; j<an; j++ )
		{
			ndel = 0;
			for( i=0,q=0; i<len; i++ )
			{
				tmpc = aseq[j][i];
#if 0
				if( tmpc != '-' ) 
				{
					if( eqseq[i] == '=' )
					{
//						reporterr( "deleting %d-%d, %c\n", j, i, aseq[j][i] );
						deletelist[j] = realloc( deletelist[j], sizeof( int ) * (ndel+2) );
						deletelist[j][ndel] = q;
						ndel++;
					}
					q++;
				}
#else
				qstart = q;
				while( tmpc != '-' && eqseq[i] == '=' )
				{
					i++;
					q++;
					tmpc = aseq[j][i];
				}
				if( qstart != q )
				{
//					reporterr( "deleting %d-%d, %c\n", j, i, aseq[j][i] );
					deletelist[j] = realloc( deletelist[j], sizeof( GapPos ) * (ndel+2) );
					deletelist[j][ndel].pos = qstart;
					deletelist[j][ndel].len = q-qstart;
					ndel++;
				}
				if( tmpc != '-' ) q++;
#endif
			}
			deletelist[j][ndel].pos = -1;
			deletelist[j][ndel].len = 0;
		}
	}
	for( i=0,p=0; i<len; i++ )
	{

//		if( oseq[0][i] != '=' )
//		reporterr( "i=%d, p=%d, q=%d, originally, %c\n", i, p, q, originallygapped[p]);
//		if( eqseq[i] != '=' && originallygapped[p] != '-' ) // dame!!
		if( eqseq[i] != '=' )
		{
//			reporterr( "COPY! p=%d\n", p );
			if( p != i )
			{
				for( j=0; j<on; j++ ) oseq[j][p] = oseq[j][i];
				for( j=0; j<an; j++ ) aseq[j][p] = aseq[j][i];
			}
			p++;
		}
	}
//		reporterr( "deletemap        = %s\n", deletemap );
//		reporterr( "eqseq            = %s\n", eqseq );
//		reporterr( "originallygapped = %s\n", originallygapped );
	for( j=0; j<on; j++ ) oseq[j][p] = 0;
	for( j=0; j<an; j++ ) aseq[j][p] = 0;

	free( eqseq );

//	for( j=0; j<on; j++ ) reporterr( "\no=%s\n", oseq[j] );
//	for( j=0; j<an; j++ ) reporterr( "a=%s\n", aseq[j] );
	return( i-p );

}


int maskoriginalgaps( char *repseq, char *originallygapped )
{
	int i, p;
	int len = strlen( repseq );
//	reporterr( "repseq = %s\n", repseq );
	for( i=0,p=0; i<len; i++ )
	{
		if( repseq[i] == '=' )
		{
			if( originallygapped[p] == '-' )
			{
				repseq[i] = '-';
				p++;
			}
		}
		else
		{
			p++;
		}
	}
	reporterr( "repseq = %s\n", repseq );
exit( 1 );
}

void restoregaponlysites( char *originallygapped, int n1, int n2, char **s1, char **s2, int rep )
{
	int i, j, p;
	char *tmpnew;
	int len;
	reporterr( "originallygapped = %s\n", originallygapped );
	reporterr( "s1[0]            = %s\n", s1[0] );
	reporterr( "s1[rep]          = %s\n", s1[rep] );
	reporterr( "s2[0]            = %s\n", s2[0] );
exit( 1 );

	tmpnew = calloc( strlen( originallygapped )+1, sizeof( char ) );
	len = strlen( s1[0] );

	for( i=0,p=0; i<len; i++ )
	{
		reporterr( "i=%d, p=%d, s[]=%c, o[]=%c\n", i, p, s1[0][i], originallygapped[p] );
		if( originallygapped[p] == 'o' )
		{
			tmpnew[p] = s1[0][i];
			p++;
		}
		while( originallygapped[p] == '-' )
		{
			tmpnew[p] = '-';
			p++;
		}
	}
	reporterr( "s1[0]            = %s\n", s1[0] );
	reporterr( "tmpnew           = %s\n", tmpnew );
	
}

#endif


int recordoriginalgaps( char *originallygapped, int n, char **s )
{
	int i, j;
	int len = strlen( s[0] );
	int v = 0;
	for( i=0; i<len; i++ )
	{
		for( j=0; j<n; j++ ) if( s[j][i] != '-' ) break;

		if( j == n ) 
			originallygapped[i] = '-';
		else
			originallygapped[i] = 'o';
	}
	originallygapped[i] = 0;
	return( v );
}

void restoreoriginalgaps( int n, char **seq, char *originalgaps )
{
	int i, j, p;
	int lenf = strlen( originalgaps );
	char *tmpseq = calloc( lenf+1, sizeof( char ) );

	for( i=0; i<n; i++ )
	{
		for( j=0,p=0; j<lenf; j++ )
		{
			if( originalgaps[j] == '-' )
				tmpseq[j] = '-';
			else
				tmpseq[j] = seq[i][p++];
		}
		strcpy( seq[i], tmpseq );
	}
	free( tmpseq );
}

void reconstructdeletemap( int nadd, char **addbk, GapPos **deletelist, char **realn, FILE *fp, char **name )
{
	int i, j, p, len, gaplen;
	char *gapped, *nameptr, *tmpptr;

	for( i=0; i<nadd; i++ )
	{
		len = strlen( addbk[i] );
		gapped = calloc( len+1, sizeof( char ) );
//		for( j=0; j<len; j++ ) gapped[j] = 'o'; // iranai
//		gapped[len] = 0; // iranai

		nameptr = name[i] + 1;
		if( outnumber )
			nameptr = strstr( nameptr, "_numo_e" ) + 8;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		fprintf( fp, ">%s\n", nameptr );
		fprintf( fp, "# letter, position in the original sequence, position in the reference alignment\n" );

#if 0
//		reporterr( "addbk[%d] = %s\n", i, addbk[i] );
		for( j=0; (p=deletelist[i][j])!=-1; j++ )
		{
//			reporterr( "deleting %d, %c\n", p, addbk[i][p] );
			gapped[p] = '-';
		}
#else
//		reporterr( "addbk[%d] = %s\n", i, addbk[i] );
		for( j=0; (p=deletelist[i][j].pos)!=-1; j++ )
		{
//			reporterr( "deleting %d, %c\n", p, addbk[i][p] );
			gaplen = deletelist[i][j].len;
			while( gaplen-- ) gapped[p++] = '-';
		}
#endif

//		reporterr( "addbk  = %s\n", addbk[i] );
//		reporterr( "gapped = %s\n", gapped );

		for( j=0,p=0; j<len; j++ )
		{
			while( realn[i][p] == '-' )
				p++;

			if( gapped[j] == '-' )
			{
				fprintf( fp, "%c, %d, -\n", addbk[i][j], j+1 ); // 1origin
			}
			else
			{
				fprintf( fp, "%c, %d, %d\n", addbk[i][j], j+1, p+1 ); // 1origin
				p++;
			}
		}
		free( gapped );
	}
}

#define MODIFYNAME 1 // MODIFYNAME mo --anysymbol to ryouritsu, 7.502-

void reconstructdeletemap_compact( int nadd, char **addbk, GapPos **deletelist, char **realn, FILE *fp, char **name )
{
	int i, j, p, len, nins, gaplen;
	char *gapped, *nameptr, *tmpptr;
	int status = 0;

#if MODIFYNAME
	char *newname, *insstr;
	newname = calloc( B+100, sizeof( char ) );
	insstr = calloc( 1000, sizeof( char ) );
#endif
	fprintf( fp, "# Insertion in added sequence > Position in reference\n" );
	for( i=0; i<nadd; i++ )
	{
		len = strlen( addbk[i] );
		gapped = calloc( len+1, sizeof( char ) );
//		reporterr( "len=%d", len );
//		for( j=0; j<len; j++ ) gapped[j] = 'o'; // iranai
//		gapped[len] = 0; // iranai

		nameptr = name[i] + 1;
		if( outnumber )
			nameptr = strstr( nameptr, "_numo_e" ) + 8;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		status = 0;


#if 0
//		reporterr( "addbk[%d] = %s\n", i, addbk[i] );
		for( j=0; (p=deletelist[i][j])!=-1; j++ )
		{
//			reporterr( "deleting %d, %c\n", p, addbk[i][p] );
			gapped[p] = '-';
			status = 1;
		}
#else
//		reporterr( "addbk[%d] = %s\n", i, addbk[i] );
		for( j=0; (p=deletelist[i][j].pos)!=-1; j++ )
		{
//			reporterr( "deleting %d-%d, %c\n", p, p+deletelist[i][j].len, addbk[i][p] );
			gaplen = deletelist[i][j].len;
			while( gaplen-- ) gapped[p++] = '-'; // origin??????????? 2022/Jan
			status = 1;
		}
#endif

//		reporterr( "addbk  = %s\n", addbk[i] );
//		reporterr( "gapped = %s\n", gapped );


		if( status == 0 ) 
		{
			free( gapped );
			continue;
		}

		fprintf( fp, ">%s\n", nameptr );

		status = -1;
#if MODIFYNAME
		insstr[0] = 0;
		nins = 0;
#endif
		for( j=0,p=0; j<len; j++ )
		{
			while( realn[i][p] == '-' )
				p++;

			if( gapped[j] == '-' )
			{
				if( status != 1 )	
				{
					status = 1;
					fprintf( fp, "%d%c - ", j+1, addbk[i][j] ); // 1origin
#if MODIFYNAME
					if( nins == 0 )
						sprintf( insstr, "%d%c-", j+1, addbk[i][j] );
					else if( nins == 1 )
						sprintf( insstr+strlen(insstr), "etc," );
					nins++;
#endif
				}
//				fprintf( fp, "%c, %d, -\n", addbk[i][j], j+1 ); // 1origin
			}
			else
			{
				if( status == 1 )	
				{
					fprintf( fp, "%d%c > %dv%d\n", j, addbk[i][j-1], p, p+1 ); // 1origin
#if MODIFYNAME
					if( nins == 1 ) sprintf( insstr+strlen(insstr), "%d%c,", j, addbk[i][j-1] );
#endif
				}
				status = 0;
//				fprintf( fp, "%c, %d, %d\n", addbk[i][j], j+1, p+1 ); // 1origin
				p++;
			}
		}
		if( status == 1 )	
		{
			fprintf( fp, "%d%c > %dv%d\n", j, addbk[i][j-1], p, p+1 ); // 1origin
#if MODIFYNAME
			if( nins == 1 ) sprintf( insstr+strlen(insstr), "%d%c,", j, addbk[i][j-1] );
#endif
		}
		free( gapped );

#if MODIFYNAME
		insstr[strlen(insstr)-1] = 0;
		strcpy( newname, name[i] );
		sprintf( newname+(nameptr-name[i]), "%dins:%s|%s", nins, insstr, nameptr );
		newname[B] = 0;
		strcpy( name[i], newname );
#endif
	}
#if MODIFYNAME
	free( newname );
	free( insstr );
#endif
}
