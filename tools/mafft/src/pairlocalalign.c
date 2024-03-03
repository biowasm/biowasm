#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SHISHAGONYU 0 // for debug

#define NODIST -9999

static char *whereispairalign;
static char *laraparams;
static char foldalignopt[1000];
static int stdout_align;
static int stdout_dist;
static int store_localhom;
static int store_dist;
static int laste;
static int lastm;
static int lastsubopt;
static int lastonce;

typedef struct _lastres
{
	int score;
	int start1;
	int start2;
	char *aln1;
	char *aln2;
} Lastres;

typedef struct _reg
{
	int start;
	int end;
} Reg;

typedef struct _aln
{
	int nreg;
	Reg *reg1;
	Reg *reg2;
} Aln;

typedef struct _lastresx
{
	int score;
	int naln;
	Aln *aln;
} Lastresx;

#ifdef enablemultithread
typedef struct _jobtable
{
	int i;
	int j;
} Jobtable;

typedef struct _thread_arg
{
	int thread_no;
	int njob;
	Jobtable *jobpospt;
	char **name;
	char **seq;
	char **dseq;
	int *thereisxineachseq;
	LocalHom **localhomtable;
	double **distancemtx;
	double *selfscore;
	char ***bpp;
	Lastresx **lastresx;
	int alloclen;
	int *targetmap;
	double **expdist;
	pthread_mutex_t *mutex_counter;
	pthread_mutex_t *mutex_stdout;
} thread_arg_t;
#endif

typedef struct _lastcallthread_arg
{
	int nq, nd;
	char **dseq;
	char **qseq;
	Lastresx **lastresx;
#ifdef enablemultithread
	int thread_no;
	int *kshare;
	pthread_mutex_t *mutex;
#endif
} lastcallthread_arg_t;

static void t2u( char *seq )
{
	while( *seq )
	{
		if     ( *seq == 'A' ) *seq = 'a';
		else if( *seq == 'a' ) *seq = 'a';
		else if( *seq == 'T' ) *seq = 'u';
		else if( *seq == 't' ) *seq = 'u';
		else if( *seq == 'U' ) *seq = 'u';
		else if( *seq == 'u' ) *seq = 'u';
		else if( *seq == 'G' ) *seq = 'g';
		else if( *seq == 'g' ) *seq = 'g';
		else if( *seq == 'C' ) *seq = 'c';
		else if( *seq == 'c' ) *seq = 'c';
		else *seq = 'n';
		seq++;
	}
}

static int removex( char *d, char *m )
{
	int val = 0;
	while( *m != 0 )
	{
		if( *m == 'X' || *m == 'x' ) 
		{
			m++;
			val++;
		}
		else 
		{
			*d++ = *m++;
		}
	}
	*d = 0;
	return( val );
}

static void putlocalhom_last( char *s1, char *s2, LocalHom *localhompt, Lastresx *lastresx, char korh )
{
	char *pt1, *pt2;
	int naln, nreg;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	LocalHom *tmppt2;
	LocalHom *localhompt0;
	Reg *rpt1, *rpt2;
	Aln *apt;
	int nlocalhom = 0;
	int len;

//	fprintf( stderr, "s1=%s\n", s1 );
//	fprintf( stderr, "s2=%s\n", s2 );


	naln = lastresx->naln;
	apt = lastresx->aln;

	if( naln == 0 ) return;
	while( naln-- )
	{
		rpt1 = apt->reg1;
		rpt2 = apt->reg2;
		nreg = apt->nreg;
		isumscore = 0;
		sumoverlap = 0;
		while( nreg-- )
		{
			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = rpt1->start;
			tmppt->start2 = rpt2->start;
			tmppt->end1   = rpt1->end;
			tmppt->end2   = rpt2->end;
			tmppt->korh   = 'h';
			if( rpt1 == apt->reg1 ) localhompt0 = tmppt; // ?
	
//			fprintf( stderr, "in putlocalhom, reg1: %d-%d (nreg=%d)\n", rpt1->start, rpt1->end, lastresx->nreg );
//			fprintf( stderr, "in putlocalhom, reg2: %d-%d (nreg=%d)\n", rpt2->start, rpt2->end, lastresx->nreg );
	
			len = tmppt->end1 - tmppt->start1 + 1;
	
//			fprintf( stderr, "tmppt->start1=%d\n", tmppt->start1 );
//			fprintf( stderr, "tmppt->start2=%d\n", tmppt->start2 );

//			fprintf( stderr, "s1+tmppt->start1=%*.*s\n", len, len, s1+tmppt->start1 );
//			fprintf( stderr, "s2+tmppt->start2=%*.*s\n", len, len, s2+tmppt->start2 );
	
			pt1 = s1 + tmppt->start1;
			pt2 = s2 + tmppt->start2;
			iscore = 0;
			while( len-- )
			{
				iscore += n_dis[(int)amino_n[(unsigned char)*pt1++]][(int)amino_n[(unsigned char)*pt2++]]; // - offset はいらないかも
//				fprintf( stderr, "len=%d, %c-%c, iscore(0) = %d\n", len, *(pt1-1), *(pt2-1), iscore );
			}
	
			if( divpairscore )
			{
				tmppt->overlapaa   = tmppt->end2-tmppt->start2+1;
				tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			}
			else
			{
				isumscore += iscore;
				sumoverlap += tmppt->end2-tmppt->start2+1;
			}
			rpt1++;
			rpt2++;
		}
#if 0
		fprintf( stderr, "iscore (1)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif

		if( !divpairscore )
		{
			for( tmppt2=localhompt0; tmppt2; tmppt2=tmppt2->next )
			{
				tmppt2->overlapaa = sumoverlap;
				tmppt2->opt = (double)isumscore * 5.8 / ( 600 * sumoverlap );
//				fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
			}
		}
		apt++;
	}
}

static int countcomma( char *s )
{
	int v = 0;
	while( *s ) if( *s++ == ',' ) v++;
	return( v );
}

static double recallpairfoldalign( char **mseq1, char **mseq2, int m1, int m2, int *of1pt, int *of2pt, int alloclen )
{
	static FILE *fp = NULL;
	double value;
	char *aln1;
	char *aln2;
	int of1tmp, of2tmp;

	if( fp == NULL )
	{
		fp = fopen( "_foldalignout", "r" );
		if( fp == NULL )
		{
			fprintf( stderr, "Cannot open _foldalignout\n" );
			exit( 1 );
		}
	}

	aln1 = calloc( alloclen, sizeof( char ) );
	aln2 = calloc( alloclen, sizeof( char ) );

	readpairfoldalign( fp, *mseq1, *mseq2, aln1, aln2, m1, m2, &of1tmp, &of2tmp, alloclen );

	if( strstr( foldalignopt, "-global") )
	{
		fprintf( stderr, "Calling G__align11\n" );
		value = G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );
		*of1pt = 0;
		*of2pt = 0;
	}
	else
	{
		fprintf( stderr, "Calling L__align11\n" );
		value = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, of1pt, of2pt );
	}

//	value = (double)naivepairscore11( *mseq1, *mseq2, penalty ); // nennnotame

	if( aln1[0] == 0 )
	{
		fprintf( stderr, "FOLDALIGN returned no alignment between %d and %d.  Sequence alignment is used instead.\n", m1+1, m2+1 );
	}
	else
	{
		strcpy( *mseq1, aln1 );
		strcpy( *mseq2, aln2 );
		*of1pt = of1tmp;
		*of2pt = of2tmp;
	}

//	value = naivepairscore11( *mseq1, *mseq2, penalty ); // v6.511 ha kore wo tsukau, global nomi dakara.

//	fclose( fp ); // saigo dake yatta houga yoi.

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );


	free( aln1 );
	free( aln2 );

	return( value );
}

static void block2reg( char *block, Reg *reg1, Reg *reg2, int start1, int start2 )
{
	Reg *rpt1, *rpt2;
	char *tpt, *npt;
	int pos1, pos2;
	int len, glen1, glen2;
	pos1 = start1;
	pos2 = start2;
	rpt1 = reg1;
	rpt2 = reg2;
	while( block )
	{
		block++;
//		fprintf( stderr, "block = %s\n", block );
		tpt = strchr( block, ':' );
		npt = strchr( block, ',' );
		if( !tpt || tpt > npt )
		{
			len = atoi( block );
			reg1->start = pos1;
			reg2->start = pos2;
			pos1 += len - 1;
			pos2 += len - 1;
			reg1->end = pos1;
			reg2->end = pos2;
//			fprintf( stderr, "in loop reg1: %d-%d\n", reg1->start, reg1->end );
//			fprintf( stderr, "in loop reg2: %d-%d\n", reg2->start, reg2->end );
			reg1++;
			reg2++;
		}
		else
		{
			sscanf( block, "%d:%d", &glen1, &glen2 );
			pos1 += glen1 + 1;
			pos2 += glen2 + 1;
		}
		block = npt;

	}
	reg1->start = reg1->end = reg2->start = reg2->end = -1;
	
	while( rpt1->start != -1 )
	{
//		fprintf( stderr, "reg1: %d-%d\n", rpt1->start, rpt1->end );
//		fprintf( stderr, "reg2: %d-%d\n", rpt2->start, rpt2->end );
		rpt1++;
		rpt2++;
	}
//	*apt1 = *apt2 = 0;
//	fprintf( stderr, "aln1 = %s\n", aln1 );
//	fprintf( stderr, "aln2 = %s\n", aln2 );
}


static void readlastresx_singleq( FILE *fp, int n1, int nameq, Lastresx **lastresx )
{
	char *gett;
	Aln *tmpaln;
	int prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
	int score, name1, start1, alnSize1, seqSize1;
	int        name2, start2, alnSize2, seqSize2;
	char strand1, strand2;
	int includeintoscore;
	gett = calloc( 10000, sizeof( char ) );

//	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
//	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

	while( 1 )
	{
		fgets( gett, 9999, fp );
		if( feof( fp ) ) break;
		if( gett[0] == '#' ) continue;
//		fprintf( stdout, "gett = %s\n", gett );
		if( gett[strlen(gett)-1] != '\n' )
		{
			fprintf( stderr, "Too long line?\n" );
			exit( 1 );
		}

		sscanf( gett, "%d %d %d %d %c %d %d %d %d %c %d", 
					&score, &name1, &start1, &alnSize1, &strand1, &seqSize1,
					        &name2, &start2, &alnSize2, &strand2, &seqSize2 );

		if( alg == 'R' && name2 <= name1 ) continue;
		if( name2 != nameq )
		{
			fprintf( stderr, "BUG!!!\n" );
			exit( 1 );
		}

//		if( lastresx[name1][name2].score ) continue; // dame!!!!


		prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 1 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 1 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 1 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 1 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
		if( prevnaln ) includeintoscore = 0;
		else includeintoscore = 1;
#endif
		if( !includeintoscore && !lastsubopt )
			continue;

		naln = prevnaln + 1;
		lastresx[name1][name2].naln = naln;
//		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

		if( ( tmpaln = (Aln *)realloc( lastresx[name1][name2].aln, (naln) * sizeof( Aln ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].aln\n" );
			exit( 1 );
		}
		else
			lastresx[name1][name2].aln = tmpaln;

		nreg = countcomma( gett )/2 + 1;
		lastresx[name1][name2].aln[prevnaln].nreg = nreg;
//		lastresx[name1][name2].aln[naln].nreg = -1;
//		lastresx[name1][name2].aln[naln].reg1 = NULL;
//		lastresx[name1][name2].aln[naln].reg2 = NULL;
//		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

		if( ( lastresx[name1][name2].aln[prevnaln].reg1 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

		if( ( lastresx[name1][name2].aln[prevnaln].reg2 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

//		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
//		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
		block2reg( strrchr( gett, '\t' ), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2 );

		if( includeintoscore )
		{
			if( lastresx[name1][name2].score ) score += penalty;
			lastresx[name1][name2].score += score;
		}

//		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
	}
	free( gett );
}

#ifdef enablemultithread
#if 0
static void readlastresx_group( FILE *fp, Lastresx **lastresx )
{
	char *gett;
	Aln *tmpaln;
	int prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
	int score, name1, start1, alnSize1, seqSize1;
	int        name2, start2, alnSize2, seqSize2;
	char strand1, strand2;
	int includeintoscore;
	gett = calloc( 10000, sizeof( char ) );

//	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
//	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

	while( 1 )
	{
		fgets( gett, 9999, fp );
		if( feof( fp ) ) break;
		if( gett[0] == '#' ) continue;
//		fprintf( stdout, "gett = %s\n", gett );
		if( gett[strlen(gett)-1] != '\n' )
		{
			fprintf( stderr, "Too long line?\n" );
			exit( 1 );
		}

		sscanf( gett, "%d %d %d %d %c %d %d %d %d %c %d", 
					&score, &name1, &start1, &alnSize1, &strand1, &seqSize1,
					        &name2, &start2, &alnSize2, &strand2, &seqSize2 );

		if( alg == 'R' && name2 <= name1 ) continue;

//		if( lastresx[name1][name2].score ) continue; // dame!!!!

		prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 3 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 3 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 3 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 3 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
		if( prevnaln ) includeintoscore = 0;
		else includeintoscore = 1;
#endif
		if( !includeintoscore && !lastsubopt )
			continue;

		naln = prevnaln + 1;
		lastresx[name1][name2].naln = naln;
//		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

		if( ( tmpaln = (Aln *)realloc( lastresx[name1][name2].aln, (naln) * sizeof( Aln ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].aln\n" );
			exit( 1 );
		}
		else
			lastresx[name1][name2].aln = tmpaln;



		nreg = countcomma( gett )/2 + 1;
		lastresx[name1][name2].aln[prevnaln].nreg = nreg;
//		lastresx[name1][name2].aln[naln].nreg = -1;
//		lastresx[name1][name2].aln[naln].reg1 = NULL;
//		lastresx[name1][name2].aln[naln].reg2 = NULL;
//		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

		if( ( lastresx[name1][name2].aln[prevnaln].reg1 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

		if( ( lastresx[name1][name2].aln[prevnaln].reg2 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

//		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
//		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
		block2reg( strrchr( gett, '\t' ), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2 );

		if( includeintoscore )
		{
			if( lastresx[name1][name2].score ) score += penalty;
			lastresx[name1][name2].score += score;
		}

//		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
	}
	free( gett );
}
#endif
#endif

static void readlastresx( FILE *fp, int n1, int n2, Lastresx **lastresx, char **seq1, char **seq2 )
{
	char *gett;
	Aln *tmpaln;
	int prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
	int score, name1, start1, alnSize1, seqSize1;
	int        name2, start2, alnSize2, seqSize2;
	char strand1, strand2;
	int includeintoscore;
	gett = calloc( 10000, sizeof( char ) );

//	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
//	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

	while( 1 )
	{
		fgets( gett, 9999, fp );
		if( feof( fp ) ) break;
		if( gett[0] == '#' ) continue;
//		fprintf( stdout, "gett = %s\n", gett );
		if( gett[strlen(gett)-1] != '\n' )
		{
			fprintf( stderr, "Too long line?\n" );
			exit( 1 );
		}

		sscanf( gett, "%d %d %d %d %c %d %d %d %d %c %d", 
					&score, &name1, &start1, &alnSize1, &strand1, &seqSize1,
					        &name2, &start2, &alnSize2, &strand2, &seqSize2 );

		if( alg == 'R' && name2 <= name1 ) continue;

//		if( lastresx[name1][name2].score ) continue; // dame!!!!

		prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 3 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 3 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 3 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 3 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
		if( prevnaln ) includeintoscore = 0;
		else includeintoscore = 1;
#endif
		if( !includeintoscore && !lastsubopt )
			continue;

		naln = prevnaln + 1;
		lastresx[name1][name2].naln = naln;
//		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

		if( ( tmpaln = (Aln *)realloc( lastresx[name1][name2].aln, (naln) * sizeof( Aln ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].aln\n" );
			exit( 1 );
		}
		else
			lastresx[name1][name2].aln = tmpaln;



		nreg = countcomma( gett )/2 + 1;
		lastresx[name1][name2].aln[prevnaln].nreg = nreg;
//		lastresx[name1][name2].aln[naln].nreg = -1;
//		lastresx[name1][name2].aln[naln].reg1 = NULL;
//		lastresx[name1][name2].aln[naln].reg2 = NULL;
//		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

		if( ( lastresx[name1][name2].aln[prevnaln].reg1 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

		if( ( lastresx[name1][name2].aln[prevnaln].reg2 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

//		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
//		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
		block2reg( strrchr( gett, '\t' ), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2 );

		if( includeintoscore )
		{
			if( lastresx[name1][name2].score ) score += penalty;
			lastresx[name1][name2].score += score;
		}

//		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
	}
	free( gett );
}

#ifdef enablemultithread
#if 0
static void *lastcallthread_group( void *arg )
{
	lastcallthread_arg_t *targ = (lastcallthread_arg_t *)arg;
	int k, i;
	int nq = targ->nq;
	int nd = targ->nd;
#ifdef enablemultithread
	int thread_no = targ->thread_no;
	int *kshare = targ->kshare; 
#endif
	Lastresx **lastresx = targ->lastresx;
	char **dseq = targ->dseq;
	char **qseq = targ->qseq;
	char command[5000];
	FILE *lfp;
	int msize;
	int klim;
	int qstart, qend, shou, amari;
	char kd[1000];

	if( nthread )
	{
		shou = nq / nthread;
		amari = nq - shou * nthread;
		fprintf( stderr, "shou: %d, amari: %d\n", shou, amari );

		qstart = thread_no * shou;
		if( thread_no - 1 < amari ) qstart += thread_no;
		else qstart += amari;

		qend = qstart + shou - 1;
		if( thread_no < amari ) qend += 1;
		fprintf( stderr, "%d: %d-%d\n", thread_no, qstart, qend );
	}
	k = -1;
	while( 1 )
	{
		if( nthread )
		{
			if( qstart > qend ) break;
			if( k == thread_no ) break;
			fprintf( stderr, "\n%d-%d / %d (thread %d)                    \n", qstart, qend, nq, thread_no );
			k = thread_no;
		}
		else
		{
			k++;
			if( k == nq ) break;
			fprintf( stderr, "\r%d / %d                    \r", k, nq );
		}

		if( alg == 'R' ) // if 'r' -> calllast_fast
		{
			fprintf( stderr, "Not supported\n" );
			exit( 1 );
		}
		else // 'r'
		{
			kd[0] = 0;
		}
		
		sprintf( command, "_q%d", k );
		lfp = fopen( command, "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open %s", command );
			exit( 1 );
		}
		for( i=qstart; i<=qend; i++ )
			fprintf( lfp, ">%d\n%s\n", i, qseq[i] );
		fclose( lfp );
	
//		if( alg == 'R' ) msize = MAX(10,k+nq);
//			else msize = MAX(10,nd+nq);
		if( alg == 'R' ) msize = MAX(10,k*lastm);
			else msize = MAX(10,nd*lastm);

//		fprintf( stderr, "Calling lastal from lastcallthread, msize = %d, k=%d\n", msize, k );
//		sprintf( command, "grep '>' _db%sd", kd );
//		system( command );
		sprintf( command, "%s/lastal -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db%sd _q%d > _lastres%d", whereispairalign, msize, laste, -penalty, -penalty_ex, kd, k, k );
		if( system( command ) ) exit( 1 );
	
		sprintf( command, "_lastres%d", k );
		lfp = fopen( command, "r" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot read _lastres%d", k );
			exit( 1 );
		}
//		readlastres( lfp, nd, nq, lastres, dseq, qseq );
//		fprintf( stderr, "Reading lastres\n" );
		readlastresx_group( lfp, lastresx );
		fclose( lfp );
	}
	return( NULL );
}
#endif
#endif

static void *lastcallthread( void *arg )
{
	lastcallthread_arg_t *targ = (lastcallthread_arg_t *)arg;
	int k, i;
	int nq = targ->nq;
	int nd = targ->nd;
#ifdef enablemultithread
	int thread_no = targ->thread_no;
	int *kshare = targ->kshare; 
#endif
	Lastresx **lastresx = targ->lastresx;
	char **dseq = targ->dseq;
	char **qseq = targ->qseq;
	char command[5000];
	FILE *lfp;
	int msize;
	int klim;
	char kd[1000];

	k = -1;
	while( 1 )
	{

#ifdef enablemultithread
		if( nthread )
		{
			pthread_mutex_lock( targ->mutex );
			k = *kshare;
			if( k == nq )
			{
				pthread_mutex_unlock( targ->mutex );
				break;
			}
			fprintf( stderr, "\r%d / %d (thread %d)                    \r", k, nq, thread_no );
			++(*kshare);
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			k++;
			if( k == nq ) break;
			fprintf( stderr, "\r%d / %d                    \r", k, nq );
		}

		if( alg == 'R' ) // if 'r' -> calllast_fast
		{
			klim = MIN( k, njob-nadd );
//			klim = k; // dochira demo yoi
			if( klim == k ) 
			{
				sprintf( command, "_db%dd", k );
				lfp = fopen( command, "w" );
				if( !lfp )
				{
					fprintf( stderr, "Cannot open _db." );
					exit( 1 );
				}
				for( i=0; i<klim; i++ ) fprintf( lfp, ">%d\n%s\n", i, dseq[i] );
				fclose( lfp );

//				sprintf( command, "md5sum _db%dd > /dev/tty", k );
//				system( command );

				if( dorp == 'd' ) 
					sprintf( command, "%s/lastdb _db%dd _db%dd", whereispairalign, k, k );
				else
					sprintf( command, "%s/lastdb -p _db%dd _db%dd", whereispairalign, k, k );
				system( command );
				sprintf( kd, "%d", k );
			}
			else // calllast_fast de tsukutta nowo riyou
			{
				kd[0] = 0;
//				fprintf( stderr, "klim=%d, njob=%d, nadd=%d, skip!\n", klim, njob, nadd );
			}
		}
		else // 'r'
		{
			kd[0] = 0;
		}
		
		sprintf( command, "_q%d", k );
		lfp = fopen( command, "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open %s", command );
			exit( 1 );
		}
		fprintf( lfp, ">%d\n%s\n", k, qseq[k] );
		fclose( lfp );
	
//		if( alg == 'R' ) msize = MAX(10,k+nq);
//			else msize = MAX(10,nd+nq);
		if( alg == 'R' ) msize = MAX(10,k*lastm);
			else msize = MAX(10,nd*lastm);

//		fprintf( stderr, "Calling lastal from lastcallthread, msize = %d, k=%d\n", msize, k );
//		sprintf( command, "grep '>' _db%sd", kd );
//		system( command );
		sprintf( command, "%s/lastal -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db%sd _q%d > _lastres%d", whereispairalign, msize, laste, -penalty, -penalty_ex, kd, k, k );
		if( system( command ) ) exit( 1 );
	
		sprintf( command, "_lastres%d", k );
		lfp = fopen( command, "r" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot read _lastres%d", k );
			exit( 1 );
		}
//		readlastres( lfp, nd, nq, lastres, dseq, qseq );
//		fprintf( stderr, "Reading lastres\n" );
		readlastresx_singleq( lfp, nd, k, lastresx );
		fclose( lfp );
	}
	return( NULL );
}


static void calllast_fast( int nd, char **dseq, int nq, char **qseq, Lastresx **lastresx )
{
	int i, j;
	FILE *lfp;
	char command[1000];

	lfp = fopen( "_scoringmatrixforlast", "w" );
	if( !lfp )
	{
		fprintf( stderr, "Cannot open _scoringmatrixforlast" );
		exit( 1 );
	}
	if( dorp == 'd' ) 
	{
		fprintf( lfp, "      " );
		for( j=0; j<4; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<4; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<4; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		}
	}
	else
	{
		fprintf( lfp, "      " );
		for( j=0; j<20; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<20; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<20; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		}
	}
	fclose( lfp );

//	if( alg == 'r' ) // if 'R' -> lastcallthread, kokonoha nadd>0 no toki nomi shiyou
	{
		sprintf( command, "_dbd" );
		lfp = fopen( command, "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open _dbd" );
			exit( 1 );
		}
		if( alg == 'R' )
			j = njob-nadd;
		else
			j = nd;
		for( i=0; i<j; i++ ) fprintf( lfp, ">%d\n%s\n", i, dseq[i] );

		fclose( lfp );
		if( dorp == 'd' ) 
			sprintf( command, "%s/lastdb _dbd _dbd", whereispairalign );
		else
			sprintf( command, "%s/lastdb -p _dbd _dbd", whereispairalign );
		system( command );
	}

#ifdef enablemultithread
	if( nthread )
	{
		pthread_t *handle;
		pthread_mutex_t mutex;
		lastcallthread_arg_t *targ;
		int *ksharept;
		targ = (lastcallthread_arg_t *)calloc( nthread, sizeof( lastcallthread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		ksharept = calloc( 1, sizeof(int) );
		*ksharept = 0;
		pthread_mutex_init( &mutex, NULL );
		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].kshare = ksharept;
			targ[i].nq = nq;
			targ[i].nd = nd;
			targ[i].dseq = dseq;
			targ[i].qseq = qseq;
			targ[i].lastresx = lastresx;
			targ[i].mutex = &mutex;
			pthread_create( handle+i, NULL, lastcallthread, (void *)(targ+i) );
		}

		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex );
		free( handle );
		free( targ );
		free( ksharept );
	}
	else
#endif
	{
		lastcallthread_arg_t *targ;
		targ = (lastcallthread_arg_t *)calloc( 1, sizeof( lastcallthread_arg_t ) );
		targ[0].nq = nq;
		targ[0].nd = nd;
		targ[0].dseq = dseq;
		targ[0].qseq = qseq;
		targ[0].lastresx = lastresx;
		lastcallthread( targ );
		free( targ );
	}

}

static void calllast_once( int nd, char **dseq, int nq, char **qseq, Lastresx **lastresx )
{
	int i, j;
	char command[5000];
	FILE *lfp;
	int msize;
	int res;

	fprintf( stderr, "nq=%d\n", nq );

	lfp = fopen( "_db", "w" );
	if( !lfp )
	{
		fprintf( stderr, "Cannot open _db" );
		exit( 1 );
	}
	for( i=0; i<nd; i++ ) fprintf( lfp, ">%d\n%s\n", i, dseq[i] );
	fclose( lfp );

	if( dorp == 'd' ) 
	{
		sprintf( command, "%s/lastdb _db _db", whereispairalign );
		system( command );
		lfp = fopen( "_scoringmatrixforlast", "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open _scoringmatrixforlast" );
			exit( 1 );
		}
		fprintf( lfp, "      " );
		for( j=0; j<4; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<4; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<4; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		}
		fclose( lfp );
#if 0
		sprintf( command, "lastex -s 2 -a %d -b %d -p _scoringmatrixforlast -E 10000 _db.prj _db.prj > _lastex", -penalty, -penalty_ex );
		system( command );
		lfp = fopen( "_lastex", "r" );
		fgets( command, 4999, lfp );
		fgets( command, 4999, lfp );
		fgets( command, 4999, lfp );
		fgets( command, 4999, lfp );
		laste = atoi( command );
		fclose( lfp );
		fprintf( stderr, "laste = %d\n", laste );
		sleep( 10 );
#else
//		laste = 5000;
#endif
	}
	else
	{
		sprintf( command, "%s/lastdb -p _db _db", whereispairalign );
		system( command );
		lfp = fopen( "_scoringmatrixforlast", "w" );
		if( !lfp )
		{
			fprintf( stderr, "Cannot open _scoringmatrixforlast" );
			exit( 1 );
		}
		fprintf( lfp, "      " );
		for( j=0; j<20; j++ ) fprintf( lfp, " %c ", amino[j] );
		fprintf( lfp, "\n" );
		for( i=0; i<20; i++ )
		{
			fprintf( lfp, "%c ", amino[i] );
			for( j=0; j<20; j++ ) fprintf( lfp, " %d ", n_dis[i][j] );
			fprintf( lfp, "\n" );
		}
		fclose( lfp );
//		fprintf( stderr, "Not written yet\n" );
	}

	lfp = fopen( "_q", "w" );
	if( !lfp )
	{
		fprintf( stderr, "Cannot open _q" );
		exit( 1 );
	}
	for( i=0; i<nq; i++ )
	{
		fprintf( lfp, ">%d\n%s\n", i, qseq[i] );
	}
	fclose( lfp );

	msize = MAX(10,nd*lastm);

//	fprintf( stderr, "Calling lastal from calllast_once, msize=%d\n", msize );
	sprintf( command, "%s/lastal -v -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db _q > _lastres", whereispairalign, msize, laste, -penalty, -penalty_ex );
//	sprintf( command, "lastal -v -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db _q > _lastres", 1, laste, -penalty, -penalty_ex );
//	sprintf( command, "lastal -v -e 40 -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db _q > _lastres", -penalty, -penalty_ex );
	res = system( command );
	if( res )
	{
		fprintf( stderr, "LAST aborted\n" );
		exit( 1 );
	}

	lfp = fopen( "_lastres", "r" );
	if( !lfp )
	{
		fprintf( stderr, "Cannot read _lastres" );
		exit( 1 );
	}
//	readlastres( lfp, nd, nq, lastres, dseq, qseq );
	fprintf( stderr, "Reading lastres\n" );
	readlastresx( lfp, nd, nq, lastresx, dseq, qseq );
	fclose( lfp );
}

static void callfoldalign( int nseq, char **mseq )
{
	FILE *fp;
	int i;
	int res;
	static char com[10000];

	for( i=0; i<nseq; i++ )
		t2u( mseq[i] );

	fp = fopen( "_foldalignin", "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _foldalignin\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ )
	{
		fprintf( fp, ">%d\n", i+1 );
		fprintf( fp, "%s\n", mseq[i] );
	}
	fclose( fp );

	sprintf( com, "env PATH=%s  foldalign210 %s _foldalignin > _foldalignout ", whereispairalign, foldalignopt );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in foldalign\n" );
		exit( 1 );
	}

}

static void calllara( int nseq, char **mseq, char *laraarg )
{
	FILE *fp;
	int i;
	int res;
	static char com[10000];

//	for( i=0; i<nseq; i++ )

	fp = fopen( "_larain", "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _larain\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ )
	{
		fprintf( fp, ">%d\n", i+1 );
		fprintf( fp, "%s\n", mseq[i] );
	}
	fclose( fp );


//	fprintf( stderr, "calling LaRA\n" );
	sprintf( com, "env PATH=%s:/bin:/usr/bin mafft_lara -i _larain -w _laraout -o _lara.params %s", whereispairalign, laraarg );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in lara\n" );
		exit( 1 );
	}
}

static double recalllara( char **mseq1, char **mseq2, int alloclen )
{
	static FILE *fp = NULL;
	static char *ungap1;
	static char *ungap2;
	static char *ori1;
	static char *ori2;
//	int res;
	static char com[10000];
	double value;


	if( fp == NULL )
	{
		fp = fopen( "_laraout", "r" );
		if( fp == NULL )
		{
			fprintf( stderr, "Cannot open _laraout\n" );
			exit( 1 );
		}
		ungap1 = AllocateCharVec( alloclen );
		ungap2 = AllocateCharVec( alloclen );
		ori1 = AllocateCharVec( alloclen );
		ori2 = AllocateCharVec( alloclen );
	}


	strcpy( ori1, *mseq1 );
	strcpy( ori2, *mseq2 );

	fgets( com, 999, fp );
	myfgets( com, 9999, fp );
	strcpy( *mseq1, com );
	myfgets( com, 9999, fp );
	strcpy( *mseq2, com );

	gappick0( ungap1, *mseq1 );
	gappick0( ungap2, *mseq2 );
	t2u( ungap1 );
	t2u( ungap2 );
	t2u( ori1 );
	t2u( ori2 );

	if( strcmp( ungap1, ori1 ) || strcmp( ungap2, ori2 ) )
	{
		fprintf( stderr, "SEQUENCE CHANGED!!\n" );
		fprintf( stderr, "*mseq1  = %s\n", *mseq1 );
		fprintf( stderr, "ungap1  = %s\n", ungap1 );
		fprintf( stderr, "ori1    = %s\n", ori1 );
		fprintf( stderr, "*mseq2  = %s\n", *mseq2 );
		fprintf( stderr, "ungap2  = %s\n", ungap2 );
		fprintf( stderr, "ori2    = %s\n", ori2 );
		exit( 1 );
	}

	value = (double)naivepairscore11( *mseq1, *mseq2, penalty );

//	fclose( fp ); // saigo dake yatta houga yoi.

	return( value );
}


static double calldafs_giving_bpp( char **mseq1, char **mseq2, char **bpp1, char **bpp2, int alloclen, int i, int j )
{
	FILE *fp;
	int res;
	char *com;
	double value;
	char *dirname;


	dirname = calloc( 100, sizeof( char ) );
	com = calloc( 1000, sizeof( char ) );
	sprintf( dirname, "_%d-%d", i, j );
	sprintf( com, "rm -rf %s", dirname );
	system( com );
	sprintf( com, "mkdir %s", dirname );
	system( com );


	sprintf( com, "%s/_bpporg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot write to %s/_bpporg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">a\n" );
	while( *bpp1 )
		fprintf( fp, "%s", *bpp1++ );

	fprintf( fp, ">b\n" );
	while( *bpp2 )
		fprintf( fp, "%s", *bpp2++ );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_bpporg > %s/_bpp", dirname, dirname );
	system( com ); // for cygwin, wakaran

	t2u( *mseq1 );
	t2u( *mseq2 );

	sprintf( com, "%s/_dafsinorg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_dafsinorg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">1\n" );
//	fprintf( fp, "%s\n", *mseq1 );
	write1seq( fp, *mseq1 );
	fprintf( fp, ">2\n" );
//	fprintf( fp, "%s\n", *mseq2 );
	write1seq( fp, *mseq2 );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_dafsinorg > %s/_dafsin", dirname, dirname );
	system( com ); // for cygwin, wakaran

	sprintf( com, "_dafssh%s", dirname );
	fp = fopen( com, "w" );
	fprintf( fp, "cd %s\n", dirname );
	fprintf( fp, "%s/dafs --mafft-in _bpp _dafsin > _dafsout 2>_dum\n", whereispairalign );
	fprintf( fp, "exit $tatus\n" );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < _dafssh%s > _dafssh%s.unix", dirname, dirname );
	system( com ); // for cygwin, wakaran

	sprintf( com, "sh _dafssh%s.unix 2>_dum%s", dirname, dirname );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in dafs\n" );
		exit( 1 );
	}

	sprintf( com, "%s/_dafsout", dirname );

	fp = fopen( com, "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_dafsout\n", dirname );
		exit( 1 );
	}

	myfgets( com, 999, fp ); // nagai kanousei ga arunode
	fgets( com, 999, fp );
	myfgets( com, 999, fp ); // nagai kanousei ga arunode
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq1 );
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq2 );

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

	value = (double)naivepairscore11( *mseq1, *mseq2, penalty );

#if 0
	sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
	if( system( com ) )
	{
		fprintf( stderr, "retrying to rmdir\n" );
		usleep( 2000 );
		system( com );
	}
#endif

	free( dirname );
	free( com );


	return( value );
}

static double callmxscarna_giving_bpp( char **mseq1, char **mseq2, char **bpp1, char **bpp2, int alloclen, int i, int j )
{
	FILE *fp;
	int res;
	char *com;
	double value;
	char *dirname;


	dirname = calloc( 100, sizeof( char ) );
	com = calloc( 1000, sizeof( char ) );
	sprintf( dirname, "_%d-%d", i, j );
	sprintf( com, "rm -rf %s", dirname );
	system( com );
	sprintf( com, "mkdir %s", dirname );
	system( com );


	sprintf( com, "%s/_bpporg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot write to %s/_bpporg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">a\n" );
	while( *bpp1 )
		fprintf( fp, "%s", *bpp1++ );

	fprintf( fp, ">b\n" );
	while( *bpp2 )
		fprintf( fp, "%s", *bpp2++ );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_bpporg > %s/_bpp", dirname, dirname );
	system( com ); // for cygwin, wakaran

	t2u( *mseq1 );
	t2u( *mseq2 );

	sprintf( com, "%s/_mxscarnainorg", dirname );
	fp = fopen( com, "w" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_mxscarnainorg\n", dirname );
		exit( 1 );
	}
	fprintf( fp, ">1\n" );
//	fprintf( fp, "%s\n", *mseq1 );
	write1seq( fp, *mseq1 );
	fprintf( fp, ">2\n" );
//	fprintf( fp, "%s\n", *mseq2 );
	write1seq( fp, *mseq2 );
	fclose( fp );

	sprintf( com, "tr -d '\\r' < %s/_mxscarnainorg > %s/_mxscarnain", dirname, dirname );
	system( com ); // for cygwin, wakaran

#if 0
	sprintf( com, "cd %s; %s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum", dirname, whereispairalign );
#else
	sprintf( com, "_mxscarnash%s", dirname );
	fp = fopen( com, "w" );
	fprintf( fp, "cd %s\n", dirname );
	fprintf( fp, "%s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum\n", whereispairalign );
	fprintf( fp, "exit $tatus\n" );
	fclose( fp );
//sleep( 10000 );

	sprintf( com, "tr -d '\\r' < _mxscarnash%s > _mxscarnash%s.unix", dirname, dirname );
	system( com ); // for cygwin, wakaran

	sprintf( com, "sh _mxscarnash%s.unix 2>_dum%s", dirname, dirname );
#endif
	res = system( com );
	if( res )
	{
		fprintf( stderr, "Error in mxscarna\n" );
		exit( 1 );
	}

	sprintf( com, "%s/_mxscarnaout", dirname );

	fp = fopen( com, "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s/_mxscarnaout\n", dirname );
		exit( 1 );
	}

	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq1 );
	fgets( com, 999, fp );
	load1SeqWithoutName_new( fp, *mseq2 );

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

	value = (double)naivepairscore11( *mseq1, *mseq2, penalty );

#if 0
	sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
	if( system( com ) )
	{
		fprintf( stderr, "retrying to rmdir\n" );
		usleep( 2000 );
		system( com );
	}
#endif

	free( dirname );
	free( com );


	return( value );
}

static void readhat4( FILE *fp, char ***bpp )
{
	char oneline[1000];
	int bppsize;
	int onechar;
//	double prob;
//	int posi, posj;

	bppsize = 0;
//	fprintf( stderr, "reading hat4\n" );
	onechar = getc(fp);
//	fprintf( stderr, "onechar = %c\n", onechar );
	if( onechar != '>' )
	{
		fprintf( stderr, "Format error\n" );
		exit( 1 );
	}
	ungetc( onechar, fp );
	fgets( oneline, 999, fp );
	while( 1 )
	{
		onechar = getc(fp);
		ungetc( onechar, fp );
		if( onechar == '>' || onechar == EOF )
		{
//			fprintf( stderr, "Next\n" );
			*bpp = realloc( *bpp, (bppsize+2) * sizeof( char * ) );
			(*bpp)[bppsize] = NULL;
			break;
		}
		fgets( oneline, 999, fp );
//		fprintf( stderr, "oneline=%s\n", oneline );
//		sscanf( oneline, "%d %d %lf", &posi, &posj, &prob );
//		fprintf( stderr, "%d %d -> %f\n", posi, posj, prob );
		*bpp = realloc( *bpp, (bppsize+2) * sizeof( char * ) );
		(*bpp)[bppsize] = calloc( 100, sizeof( char ) );
		strcpy( (*bpp)[bppsize], oneline );
		bppsize++;
	}
}

static void preparebpp( int nseq, char ***bpp )
{
	FILE *fp;
	int i;

	fp = fopen( "hat4", "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open hat4\n" );
		exit( 1 );
	}
	for( i=0; i<nseq; i++ )
		readhat4( fp, bpp+i );
	fclose( fp );
}

static void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	laste = 5000;
	lastm = 3;
	nadd = 0;
	lastsubopt = 0;
	lastonce = 0;
	foldalignopt[0] = 0;
	laraparams = NULL;
	inputfile = NULL;
	fftkeika = 0;
	pslocal = -1000.0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
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
	treemethod = 'x';
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	divpairscore = 0;
	stdout_align = 0;
	stdout_dist = 0;
	store_dist = 1;
	store_localhom = 1;
//	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_OP = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	ppenalty_EX = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	specificityconsideration = 0.0;
	usenaivescoreinsteadofalignmentscore = 0;
	specifictarget = 0;
	nwildcard = 0;

//	reporterr( "argc=%d\n", argc );
//	reporterr( "*argv=%s\n", *argv );
//	reporterr( "(*argv)[0]=%c\n", (*argv)[0] );
    while( --argc > 0 && (*++argv)[0] == '-' )
	{
//		reporterr( "(*argv)[0] in while loop = %s\n", (*argv) );
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
//					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'O':
					ppenalty_OP = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'E':
					ppenalty_EX = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
//					fprintf( stderr, "kimuraR = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					fprintf( stderr, "blosum %d\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
//					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
//					fprintf( stderr, "TM %d\n", pamN );
					--argc;
					goto nextoption;
#if 0
				case 'l':
					ppslocal = (int)( atof( *++argv ) * 1000 + 0.5 );
					pslocal = (int)( 600.0 / 1000.0 * ppslocal + 0.5);
//					fprintf( stderr, "ppslocal = %d\n", ppslocal );
//					fprintf( stderr, "pslocal = %d\n", pslocal );
					--argc;
					goto nextoption;
#else
				case 'l':
					if( atof( *++argv ) < 0.00001 ) store_localhom = 0;
					--argc;
					goto nextoption;
#endif
				case 'd':
					whereispairalign = *++argv;
					fprintf( stderr, "whereispairalign = %s\n", whereispairalign );
					--argc; 
					goto nextoption;
				case 'p':
					laraparams = *++argv;
					fprintf( stderr, "laraparams = %s\n", laraparams );
					--argc; 
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
//					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
#ifndef enablemultithread
					nthread = 0;
#endif
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
//					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'w':
					lastm = myatoi( *++argv );
					fprintf( stderr, "lastm = %d\n", lastm );
					--argc;
					goto nextoption;
				case 'e':
					laste = myatoi( *++argv );
					fprintf( stderr, "laste = %d\n", laste );
					--argc;
					goto nextoption;
				case 'u':
					specificityconsideration = (double)myatof( *++argv );
//					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
				case 'K': // Hontou ha iranai. disttbfast.c, tbfast.c to awaserutame.
					break;
				case 'c':
					stdout_dist = 1;
					break;
				case 'n':
					stdout_align = 1;
					break;
				case 'x':
					store_localhom = 0;
					store_dist = 0;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
#if 0
				case 'r':
					fmodel = -1;
					break;
#endif
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
#if 0
				case 'Q':
					calledByXced = 1;
					break;
				case 'x':
					disp = 1;
					break;
				case 'a':
					alg = 'a';
					break;
				case 'S':
					alg = 'S';
					break;
#endif
				case 'U':
					lastonce = 1;
					break;
				case 'S':
					lastsubopt = 1;
					break;
				case 't':
					alg = 't';
					store_localhom = 0;
					break;
				case 'L':
					alg = 'L';
					break;
				case 'Y':
					alg = 'Y'; // nadd>0 no toki nomi. moto no hairetsu to atarashii hairetsuno alignmnt -> L;
					break;
				case 'Z':
					usenaivescoreinsteadofalignmentscore = 1;
					break;
				case 's':
					alg = 's';
					break;
				case 'G':
					alg = 'G';
					break;
				case 'B':
					alg = 'B';
					break;
				case 'T':
					alg = 'T';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'R':
					alg = 'R';
					break;
				case 'r':
					alg = 'r'; // nadd>0 no toki nomi. moto no hairetsu to atarashii hairetsuno alignmnt -> R, last
					break;
				case 'N':
					alg = 'N';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'V':
					alg = 'V';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'y':
					divpairscore = 1;
					break;
				case '=':
					specifictarget = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
				case 'o':
//					foldalignopt = *++argv;
					strcat( foldalignopt, " " );
					strcat( foldalignopt, *++argv );
					fprintf( stderr, "foldalignopt = %s\n", foldalignopt );
					--argc; 
					goto nextoption;
#if 0
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
					break;
#endif
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
        fprintf( stderr, "pairlocalalign options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}

int countamino( char *s, int end )
{
	int val = 0;
	while( end-- )
		if( *s++ != '-' ) val++;
	return( val );
}

static double score2dist( double pscore, double selfscore1, double selfscore2)
{
	double val;
	double bunbo;
//	fprintf( stderr, "In score2dist\n" );

	if( (bunbo=MIN( selfscore1, selfscore2 )) == 0.0 )
		val = 2.0;
	else if( bunbo < pscore ) // mondai ari
		val = 0.0;
	else
		val = ( 1.0 - pscore / bunbo ) * 2.0;
	return( val );
}

#if enablemultithread
static void *athread( void *arg ) // alg='R', alg='r' -> tsukawarenai.
{
	thread_arg_t *targ = (thread_arg_t *)arg;
	int i, ilim, j, jst;
	int off1, off2, dum1, dum2, thereisx;
	int intdum;
	double pscore = 0.0; // by D.Mathog
	double *effarr1;
	double *effarr2;
	char **mseq1, **mseq2, **distseq1, **distseq2, **dumseq1, **dumseq2;
	char **aseq;
	double **dynamicmtx = NULL;
	double dist;
	double scoreoffset;

// thread_arg
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	Jobtable *jobpospt = targ->jobpospt;
	char **name = targ->name;
	char **seq = targ->seq;
	char **dseq = targ->dseq;
	int *thereisxineachseq = targ->thereisxineachseq;
	LocalHom **localhomtable = targ->localhomtable;
	double **distancemtx = targ->distancemtx;
	double *selfscore = targ->selfscore;
	char ***bpp = targ->bpp;
	Lastresx **lastresx = targ->lastresx;
	int alloclen = targ->alloclen;
	int *targetmap = targ->targetmap;
	double **expdist = targ->expdist;

//	fprintf( stderr, "thread %d start!\n", thread_no );

	effarr1 = AllocateDoubleVec( 1 );
	effarr2 = AllocateDoubleVec( 1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	if( alg == 'N' )
	{
		dumseq1 = AllocateCharMtx( 1, alloclen+10 );
		dumseq2 = AllocateCharMtx( 1, alloclen+10 );
	}
	distseq1 = AllocateCharMtx( 1, 0 );
	distseq2 = AllocateCharMtx( 1, 0 );
	aseq = AllocateCharMtx( 2, alloclen+10 );
	if( specificityconsideration > 0.0 ) dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );

	if( alg == 'Y' || alg == 'r' ) ilim = njob - nadd;
	else ilim = njob - 1;


	while( 1 )
	{
		pthread_mutex_lock( targ->mutex_counter );
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;

			if( alg == 'Y' || alg == 'r' ) jst = njob - nadd;
			else jst = i + 1;
			j = jst; 

			if( i == ilim )
			{
//				fprintf( stderr, "thread %d end!\n", thread_no );
				pthread_mutex_unlock( targ->mutex_counter );

				if( commonIP ) FreeIntMtx( commonIP );
				commonIP = NULL;
				if( commonJP ) FreeIntMtx( commonJP );
				commonJP = NULL;
				Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
				G__align11( NULL, NULL, NULL, 0, 0, 0 ); // 20130603
				G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
				L__align11( NULL, 0.0, NULL, NULL, 0, NULL, NULL );
				L__align11_noalign( NULL, NULL, NULL );
				genL__align11( NULL, NULL, NULL, 0, NULL, NULL );
				free( effarr1 );
				free( effarr2 );
				free( mseq1 );
				free( mseq2 );
				if( alg == 'N' )
				{
					FreeCharMtx( dumseq1 );
					FreeCharMtx( dumseq2 );
				}
				free( distseq1 );
				free( distseq2 );
				FreeCharMtx( aseq  );
				if( dynamicmtx ) FreeDoubleMtx( dynamicmtx  );
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex_counter );


//		if( j == i+1 || j % 100 == 0 ) 
		if( j == i+1 && i % 10 == 0 ) 
		{
			fprintf( stderr, "% 5d / %d (by thread %3d) \r", i, njob-nadd, thread_no );
//			fprintf( stderr, "% 5d - %5d / %d (thread %d)\n", i, j, njob, thread_no );
		}


		if( strlen( seq[i] ) == 0 || strlen( seq[j] ) == 0 )
		{
			if( store_dist )
			{
				if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = 3.0;
				else distancemtx[i][j-i] = 3.0;
			}
			if( stdout_dist) 
			{
				pthread_mutex_lock( targ->mutex_stdout );
				fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, 3.0 );
				pthread_mutex_unlock( targ->mutex_stdout );
			}
			continue;
		}

		strcpy( aseq[0], seq[i] );
		strcpy( aseq[1], seq[j] );
//		clus1 = conjuctionfortbfast( pair, i, aseq, mseq1, effarr1, effarr, indication1 );
//		clus2 = conjuctionfortbfast( pair, j, aseq, mseq2, effarr2, effarr, indication2 );
//		fprintf( stderr, "Skipping conjuction..\n" );

		effarr1[0] = 1.0;
		effarr2[0] = 1.0;
		mseq1[0] = aseq[0];
		mseq2[0] = aseq[1];

		thereisx = thereisxineachseq[i] + thereisxineachseq[j];
//		strcpy( distseq1[0], dseq[i] ); // nen no tame
//		strcpy( distseq2[0], dseq[j] ); // nen no tame
		distseq1[0] = dseq[i];
		distseq2[0] = dseq[j];

//		fprintf( stderr, "mseq1 = %s\n", mseq1[0] );
//		fprintf( stderr, "mseq2 = %s\n", mseq2[0] );
	
#if 0
		fprintf( stderr, "group1 = %.66s", indication1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		fprintf( stderr, "\n" );
#endif
//		for( l=0; l<clus1; l++ ) fprintf( stderr, "## STEP-eff for mseq1-%d %f\n", l, effarr1[l] );

		if( use_fft )
		{
			pscore = Falign( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, 1, 1, alloclen, &intdum, NULL, 0, NULL );
//			fprintf( stderr, "pscore (fft) = %f\n", pscore );
			off1 = off2 = 0;
		}
		else
		{
			switch( alg )
			{
				case( 'R' ):
					if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
						pscore = 0.0;
					else
						pscore = (double)lastresx[i][j].score; // all pair
					break;
				case( 'r' ):
					if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) )
						pscore = (double)lastresx[i][j-(njob-nadd)].score;
					else
						pscore = 0.0;
					break;
				case( 'L' ):
					if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
						pscore = 0.0;
					else
					{
						if( usenaivescoreinsteadofalignmentscore )
						{
							L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
							pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
						}
						else
						{
//							if( store_localhom )
							if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) )
							{
								pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
								if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // uwagaki
#if 1
								if( specificityconsideration > 0.0 )
								{
									if( expdist ) 
										dist = expdist[i][j];
									else
										dist = score2dist( pscore, selfscore[i], selfscore[j] );
									if( ( scoreoffset = dist2offset( dist ) ) < 0.0 )
									{
										makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
										strcpy( mseq1[0], seq[i] );
										strcpy( mseq2[0], seq[j] );
										L__align11( dynamicmtx, scoreoffset, mseq1, mseq2, alloclen, &off1, &off2 );
									}
								}
#endif
							}
							else
								pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 );
						}
					}
//					pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // CHUUI!!!!!!
					break;
				case( 'Y' ):
					if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) ) // new sequence vs exiting sequence nomi keisan
					{
						if( usenaivescoreinsteadofalignmentscore )
						{
							L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
							pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
						}
						else
						{
							if( store_localhom )
							{
								pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
								if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // uwagaki
							}
							else
								pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 );
						}
					}
					else
						pscore = 0.0;
					break;
				case( 'A' ):
					if( usenaivescoreinsteadofalignmentscore )
					{
						G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );
						pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
					}
					else
					{
//						if( store_localhom )
						if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) )
						{
							pscore = G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );
							if( thereisx ) pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
#if 1
							if( specificityconsideration > 0.0 )
							{
								if( expdist ) 
									dist = expdist[i][j];
								else
									dist = score2dist( pscore, selfscore[i], selfscore[j] );
//									dist = score2dist( L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ), selfscore[i], selfscore[j] ); // 2014/Feb/20
								if( dist2offset( dist ) < 0.0 )
								{
									makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
									strcpy( mseq1[0], seq[i] );
									strcpy( mseq2[0], seq[j] );
									G__align11( dynamicmtx, mseq1, mseq2, alloclen, outgap, outgap );
					
								}
//								pscore = (double)naivepairscore11( *mseq1, *mseq2, 0.0 );
							}
//
#endif
						}
						else
							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
					}
					off1 = off2 = 0;
					break;
				case( 'N' ):
					if( usenaivescoreinsteadofalignmentscore )
					{
						genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 );
						pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
					}
					else
					{
//						pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
						pscore = genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 );
						if( thereisx )
						{
							strcpy( dumseq1[0], distseq1[0] );
							strcpy( dumseq2[0], distseq2[0] );
							pscore = genL__align11( n_dis_consweight_multi, dumseq1, dumseq2, alloclen, &dum1, &dum2 ); // uwagaki
						}
#if 1
						if( specificityconsideration > 0.0 )
						{
							if( expdist ) 
								dist = expdist[i][j];
							else
								dist = score2dist( pscore, selfscore[i], selfscore[j] );
							if( dist2offset( dist ) < 0.0 )
							{
								makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
								strcpy( mseq1[0], seq[i] );
								strcpy( mseq2[0], seq[j] );
								genL__align11( dynamicmtx, mseq1, mseq2, alloclen, &off1, &off2 );
							}
						}
#endif
					}
					break;
				case( 't' ):
					off1 = off2 = 0;
//					pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
					pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // tsuneni distseq shiyou
					break;
				case( 's' ):
					pscore = callmxscarna_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
					off1 = off2 = 0;
					break;
				case( 'G' ):
					pscore = calldafs_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
					off1 = off2 = 0;
					break;
#if 0 
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'K' ):
					pscore = genG__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'H' ):
					pscore = recallpairfoldalign( mseq1, mseq2, i, j, &off1, &off2, alloclen );
					break;
				case( 'B' ):
				case( 'T' ):
					pscore = recalllara( mseq1, mseq2, alloclen );
					off1 = off2 = 0;
					break;
				case( 'M' ):
					pscore = MSalign11( mseq1, mseq2, alloclen );
					break;
#endif
				default:
					ErrorExit( "\n\nERROR IN SOURCE FILE\n\n" );
			}
		}

		if( alg == 't' || ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) ) // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
		{
#if SCOREOUT
			fprintf( stderr, "score = %10.2f (%d,%d)\n", pscore, i, j );
#endif
//			if( pscore > 0.0 && ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) ) x-ins-i de seido teika
			if( ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) )
			{
				if( !store_localhom )
					;
				else if( specifictarget && targetmap[i] == -1 && targetmap[j] == -1)
					;
				else if( alg == 'R' )
					putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j, lastresx[i]+j, 'h' );
				else if( alg == 'r' )
					putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), lastresx[i]+j-(njob-nadd), 'h' );// ?????
				else if( alg == 'H' )
					putlocalhom_ext( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
				else if( alg == 'Y' )
					putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
				else if( !specifictarget && alg != 'S' && alg != 'V' )
					putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-i, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
				else
//					putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ) );
				{
					if( targetmap[i] != -1 && targetmap[j] != -1 )
					{
						putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
						putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' ); // sukoshi muda.
					}
					else if( targetmap[j] != -1 )
						putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
					else if( targetmap[i] != -1 )
						putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
#if 0
					if( targetmap[i] != -1 )
						putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
					
					else if( targetmap[j] != -1 )
						putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
#endif
					else
					{
						reporterr( "okashii\n" );
						exit( 1 );
					}
				}
			}
			pscore = score2dist( pscore, selfscore[i], selfscore[j] );

//			pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 );
//			pscore = score2dist( pscore, selfscore[i], selfscore[j] );
//			reporterr( "->pscore = %f\n", pscore );

		}
		else
		{
			pscore = 2.0;
		}

#if 1 // mutex
		if( stdout_align )
		{
			pthread_mutex_lock( targ->mutex_stdout );
			if( alg != 't' )
			{
				fprintf( stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i+1, j+1, pscore );
				fprintf( stdout, ">%s\n", name[i] );
				write1seq( stdout, mseq1[0] );
				fprintf( stdout, ">%s\n", name[j] );
				write1seq( stdout, mseq2[0] );
				fprintf( stdout, "\n" );
			}
			pthread_mutex_unlock( targ->mutex_stdout );
		}
		if( stdout_dist )
		{
			pthread_mutex_lock( targ->mutex_stdout );
			if( j == i+1 ) fprintf( stdout, "%d %d d=%.3f\n", i+1, i+1, 0.0 );
			fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, pscore );
			pthread_mutex_unlock( targ->mutex_stdout );
		}
#endif // mutex
		if( store_dist )
		{
			if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = pscore;
			else distancemtx[i][j-i] = pscore;
		}
	}
}
#endif

static void pairalign( char **name, int *nlen, char **seq, char **aseq, char **dseq, int *thereisxineachseq, char **mseq1, char **mseq2, int alloclen, Lastresx **lastresx, double **distancemtx, LocalHom **localhomtable, double **expdist, int ngui )
{
	int i, j, ilim, jst, jj;
	int off1, off2, dum1, dum2, thereisx;
	double pscore = 0.0; // by D.Mathog
	FILE *hat2p, *hat3p;
//	double **distancemtx;
	double *selfscore;
	double *effarr1;
	double *effarr2;
	char *pt;
	char *hat2file = "hat2";
//	LocalHom **localhomtable = NULL, 
	LocalHom *tmpptr;
	int intdum;
	char ***bpp = NULL; // mxscarna no toki dake
	char **distseq1, **distseq2;
	char **dumseq1, **dumseq2;
	double dist;
	double scoreoffset;
	int ntarget;
	int *targetmap, *targetmapr;


	if( specifictarget )
	{
		targetmap = calloc( njob, sizeof( int ) );
		ntarget = 0;
		for( i=0; i<njob; i++ )
		{
			targetmap[i] = -1;
			if( !strncmp( name[i]+1, "_focus_", 7 ) )
				targetmap[i] = ntarget++;
		}
		targetmapr = calloc( ntarget, sizeof( int ) );
		for( i=0; i<njob; i++ )
			if( targetmap[i] != -1 ) targetmapr[targetmap[i]] = i;

		if( ntarget == 0 )
		{
			reporterr( "\n\nAdd '>_focus_' to the title lines of the sequences to be focused on.\n\n" );
			exit( 1 );
		}
		else
		{
			reporterr( "nfocus = %d \n", ntarget );
		}
	}
	else
	{
		ntarget = njob;
		targetmap = calloc( njob, sizeof( int ) );
		targetmapr = calloc( njob, sizeof( int ) );
		for( i=0; i<njob; i++ )
			targetmap[i] = targetmapr[i] = i;
	}

#if 0
	for( i=0; i<njob; i++ )
		reporterr( "targetmap[%d] = %d\n", i, targetmap[i] );
	for( i=0; i<ntarget; i++ )
		reporterr( "targetmapr[%d] = %d\n", i, targetmapr[i] );
#endif

	if( store_localhom && localhomtable == NULL )
	{
		if( alg == 'Y' || alg == 'r' )
		{
			ilim = njob - nadd;
			jst = nadd;
		}
		else
		{
			ilim = ntarget;
			jst = njob;
		}
		localhomtable = (LocalHom **)calloc( ilim, sizeof( LocalHom *) );
		for( i=0; i<ilim; i++)
		{
			localhomtable[i] = (LocalHom *)calloc( jst, sizeof( LocalHom ) );
			for( j=0; j<jst; j++)
			{
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1; 
				localhomtable[i][j].end2 = -1; 
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].next = NULL;
				localhomtable[i][j].nokori = 0;
				localhomtable[i][j].extended = -1;
				localhomtable[i][j].last = localhomtable[i]+j;
				localhomtable[i][j].korh = 'h';
			}
			if( !specifictarget && alg != 'Y' && alg != 'r' ) jst--;
		}
	}

	if( store_dist )
	{
		if( ngui == 0 )
		{
			if( alg == 'Y' || alg == 'r' )
//				distancemtx = AllocateDoubleMtx( njob, nadd );
				distancemtx = AllocateDoubleMtx( njob-nadd, nadd ); // 2020/Oct/23
			else
				distancemtx = AllocateDoubleHalfMtx( njob );
//				distancemtx = AllocateDoubleMtx( njob, njob );
		}
	}
	else distancemtx = NULL;

	if( alg == 'N' )
	{
		dumseq1 = AllocateCharMtx( 1, alloclen+10 );
		dumseq2 = AllocateCharMtx( 1, alloclen+10 );
	}
	distseq1 = AllocateCharMtx( 1, 0 ); // muda
	distseq2 = AllocateCharMtx( 1, 0 ); // muda

	selfscore = AllocateDoubleVec( njob );
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif


//	writePre( njob, name, nlen, aseq, 0 );

	reporterr( "All-to-all alignment.\n" );
	if( alg == 'R' )
	{
		fprintf( stderr, "Calling last (http://last.cbrc.jp/)\n" );
		if( lastonce )
			calllast_once( njob, seq, njob, seq, lastresx );
		else
			calllast_fast( njob, seq, njob, seq, lastresx );
		fprintf( stderr, "done.\n" );
//		nthread = 0; // igo multithread nashi
	}
	if( alg == 'r' )
	{
		fprintf( stderr, "Calling last (http://last.cbrc.jp/)\n" );
		fprintf( stderr, "nadd=%d\n", nadd );
#if 1 // last_fast ha, lastdb ga muda
		if( lastonce )
			calllast_once( njob-nadd, seq, nadd, seq+njob-nadd, lastresx );
		else
			calllast_fast( njob-nadd, seq, nadd, seq+njob-nadd, lastresx );
#else
		calllast_once( njob-nadd, seq, nadd, seq+njob-nadd, lastresx );
#endif

		fprintf( stderr, "nadd=%d\n", nadd );
		fprintf( stderr, "done.\n" );
//		nthread = 0; // igo multithread nashi
	}

	if( alg == 'H' )
	{
		fprintf( stderr, "Calling FOLDALIGN with option '%s'\n", foldalignopt );
		callfoldalign( njob, seq );
		fprintf( stderr, "done.\n" );
	}
	if( alg == 'B' )
	{
		fprintf( stderr, "Running LARA (Bauer et al. http://www.planet-lisa.net/)\n" );
		calllara( njob, seq, "" );
		fprintf( stderr, "done.\n" );
	}
	if( alg == 'T' )
	{
		fprintf( stderr, "Running SLARA (Bauer et al. http://www.planet-lisa.net/)\n" );
		calllara( njob, seq, "-s" );
		fprintf( stderr, "done.\n" );
	}
	if( alg == 's' )
	{
		fprintf( stderr, "Preparing bpp\n" );
//		bpp = AllocateCharCub( njob, nlenmax, 0 );
		bpp = calloc( njob, sizeof( char ** ) );
		preparebpp( njob, bpp );
		fprintf( stderr, "done.\n" );
		fprintf( stderr, "Running MXSCARNA (Tabei et al. http://www.ncrna.org/software/mxscarna)\n" );
	}
	if( alg == 'G' )
	{
		fprintf( stderr, "Preparing bpp\n" );
//		bpp = AllocateCharCub( njob, nlenmax, 0 );
		bpp = calloc( njob, sizeof( char ** ) );
		preparebpp( njob, bpp );
		fprintf( stderr, "done.\n" );
		fprintf( stderr, "Running DAFS (Sato et al. http://www.ncrna.org/)\n" );
	}

	for( i=0; i<njob; i++ )
	{
		pscore = 0.0;
		for( pt=seq[i]; *pt; pt++ )
			pscore += amino_dis[(unsigned char)*pt][(unsigned char)*pt];
		selfscore[i] = pscore;
//		fprintf( stderr, "selfscore[%d] = %f\n", i, selfscore[i] );
	}

#if enablemultithread
	if( nthread > 0 ) // alg=='r' || alg=='R' -> nthread:=0 (sukoshi ue)
	{
		Jobtable jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex_counter;
		pthread_mutex_t mutex_stdout;
		thread_arg_t *targ;

		if( alg == 'Y' || alg == 'r' ) jobpos.j = njob - nadd - 1;
		else jobpos.j = 0;
		jobpos.i = 0;

		targ = calloc( nthread, sizeof( thread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex_counter, NULL );
		pthread_mutex_init( &mutex_stdout, NULL );

		for( i=0; i<nthread; i++ )
		{
			targ[i].thread_no = i;
			targ[i].njob = njob;
			targ[i].jobpospt = &jobpos;
			targ[i].name = name;
			targ[i].seq = seq;
			targ[i].dseq = dseq;
			targ[i].thereisxineachseq = thereisxineachseq;
			targ[i].localhomtable = localhomtable;
			targ[i].distancemtx = distancemtx;
			targ[i].selfscore = selfscore;
			targ[i].bpp = bpp; 
			targ[i].lastresx = lastresx;
			targ[i].alloclen = alloclen;
			targ[i].expdist = expdist;
			targ[i].targetmap = targetmap;
			targ[i].mutex_counter = &mutex_counter;
			targ[i].mutex_stdout = &mutex_stdout;

//			athread( (void *)targ );
			pthread_create( handle+i, NULL, athread, (void *)(targ+i) );
//			pthread_create( handle+i, NULL, bthread, (void *)(targ+i) );
		}


		for( i=0; i<nthread; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex_counter );
		pthread_mutex_destroy( &mutex_stdout );
		free( handle );
		free( targ );
	}
	else
#endif
	{
		double **dynamicmtx = NULL;
		if( specificityconsideration > 0.0 ) dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );

		if( alg == 'Y' || alg == 'r' ) ilim = njob - nadd;
		else ilim = njob - 1;
		for( i=0; i<ilim; i++ ) 
		{
			if( stdout_dist) fprintf( stdout, "%d %d d=%.3f\n", i+1, i+1, 0.0 );
			fprintf( stderr, "% 5d / %d\r", i, njob-nadd );
			fflush( stderr );

			if( alg == 'Y' || alg == 'r' ) jst = njob - nadd;
			else jst = i + 1;
			for( j=jst; j<njob; j++ )
			{
	
				if( strlen( seq[i] ) == 0 || strlen( seq[j] ) == 0 )
				{
					if( store_dist ) 
					{
						if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = 3.0;
						else distancemtx[i][j-i] = 3.0;
					}
					if( stdout_dist) fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, 3.0 );
					continue;
				}
	
				strcpy( aseq[0], seq[i] );
				strcpy( aseq[1], seq[j] );
//				clus1 = conjuctionfortbfast( pair, i, aseq, mseq1, effarr1, effarr, indication1 );
//				clus2 = conjuctionfortbfast( pair, j, aseq, mseq2, effarr2, effarr, indication2 );
//				fprintf( stderr, "Skipping conjuction..\n" );

				effarr1[0] = 1.0;
				effarr2[0] = 1.0;
				mseq1[0] = aseq[0];
				mseq2[0] = aseq[1];

				thereisx = thereisxineachseq[i] + thereisxineachseq[j];
//				strcpy( distseq1[0], dseq[i] ); // nen no tame
//				strcpy( distseq2[0], dseq[j] ); // nen no tame
				distseq1[0] = dseq[i];
				distseq2[0] = dseq[j];

	//			fprintf( stderr, "mseq1 = %s\n", mseq1[0] );
	//			fprintf( stderr, "mseq2 = %s\n", mseq2[0] );
		
#if 0
				fprintf( stderr, "group1 = %.66s", indication1 );
				fprintf( stderr, "\n" );
				fprintf( stderr, "group2 = %.66s", indication2 );
				fprintf( stderr, "\n" );
#endif
	//			for( l=0; l<clus1; l++ ) fprintf( stderr, "## STEP-eff for mseq1-%d %f\n", l, effarr1[l] );
	
				if( use_fft )
				{
					pscore = Falign( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, 1, 1, alloclen, &intdum, NULL, 0, NULL );
//					fprintf( stderr, "pscore (fft) = %f\n", pscore );
					off1 = off2 = 0;
				}
				else
				{
					switch( alg )
					{
						case( 't' ):
//							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // tsuneni distseq shiyou
							off1 = off2 = 0;
							break;
						case( 'A' ):
							if( usenaivescoreinsteadofalignmentscore )
							{
								G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );
								pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
							}
							else
							{
//								if( store_localhom )
								if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) )
								{
									pscore = G__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, outgap, outgap );
									if( thereisx ) pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
#if 1
									if( specificityconsideration > 0.0 )
									{
										if( expdist ) 
											dist = expdist[i][j];
										else
											dist = score2dist( pscore, selfscore[i], selfscore[j] );
//										dist = score2dist( L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ), selfscore[i], selfscore[j] ); // 2014/Feb/20
//										reporterr( "dist(%d,%d)=%f\n", i, j, dist );
										if( dist2offset( dist ) < 0.0 )
										{
											makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
											strcpy( mseq1[0], seq[i] );
											strcpy( mseq2[0], seq[j] );
											G__align11( dynamicmtx, mseq1, mseq2, alloclen, outgap, outgap );
										}
//										pscore = (double)naivepairscore11( *mseq1, *mseq2, 0.0 );
									}
#endif
								}
								else
									pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // uwagaki
							}
							off1 = off2 = 0;
							break;
						case( 'N' ):
//							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, alloclen );
							if( usenaivescoreinsteadofalignmentscore )
							{
								genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 );
								pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
							}
							else
							{
								pscore = genL__align11( n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2 );
								if( thereisx )
								{
									strcpy( dumseq1[0], distseq1[0] );
									strcpy( dumseq2[0], distseq2[0] );
									pscore = genL__align11( n_dis_consweight_multi, dumseq1, dumseq2, alloclen, &dum1, &dum2 ); // uwagaki
								}
#if 1
								if( specificityconsideration > 0.0 )
								{
//									fprintf( stderr, "dist = %f\n", score2dist( pscore, selfscore[i], selfscore[j] ) );
									if( expdist ) 
										dist = expdist[i][j];
									else
										dist = score2dist( pscore, selfscore[i], selfscore[j] );
									if( dist2offset( dist ) < 0.0 )
									{
										makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
										strcpy( mseq1[0], seq[i] );
										strcpy( mseq2[0], seq[j] );
										genL__align11( dynamicmtx, mseq1, mseq2, alloclen, &off1, &off2 );
									}
								}
#endif
							}
							break;
						case( 'R' ):
							if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
								pscore = 0.0;
							else
								pscore = (double)lastresx[i][j].score; // all pair
							break;
						case( 'r' ):
							if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) )
								pscore = (double)lastresx[i][j-(njob-nadd)].score;
							else
								pscore = 0.0;
							break;
						case( 'L' ):
							if( nadd && njob-nadd <= j && njob-nadd <= i ) // new sequence doushi ha mushi
								pscore = 0.0;
							else
							{
								if( usenaivescoreinsteadofalignmentscore )
								{
									L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
									pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
								}
								else
								{
//									if( store_localhom )
									if( store_localhom && ( targetmap[i] != -1 || targetmap[j] != -1 ) )
									{
										pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 ); // all pair
										if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // all pair
#if 1
										if( specificityconsideration > 0.0 )
										{
											if( expdist ) 
												dist = expdist[i][j];
											else
												dist = score2dist( pscore, selfscore[i], selfscore[j] );
											if( ( scoreoffset = dist2offset( dist ) ) < 0.0 )
											{
												makedynamicmtx( dynamicmtx, n_dis_consweight_multi, 0.5 * dist ); // upgma ni awaseru.
												strcpy( mseq1[0], seq[i] );
												strcpy( mseq2[0], seq[j] );
												L__align11( dynamicmtx, scoreoffset, mseq1, mseq2, alloclen, &off1, &off2 );
											}
										}
#endif
									}
									else
										pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // all pair
								}
							}
//							pscore = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, distseq1, distseq2, alloclen ); // CHUUI!!!!!!
							break;
						case( 'Y' ):
							if( nadd == 0 || ( i < njob-nadd && njob-nadd <= j ) ) // new sequence vs exiting sequence nomi keisan
							{
								if( usenaivescoreinsteadofalignmentscore )
								{
									L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
									pscore = (double)naivepairscore11( mseq1[0], mseq2[0], 0.0 ); // uwagaki
								}
								else
								{
									if( store_localhom )
									{
										pscore = L__align11( n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2 );
										if( thereisx ) pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 ); // uwagaki
									}
									else
										pscore = L__align11_noalign( n_dis_consweight_multi, distseq1, distseq2 );
								}
							}
							else
								pscore = 0.0;
							break;
						case( 'a' ):
							pscore = Aalign( mseq1, mseq2, effarr1, effarr2, 1, 1, alloclen );
							off1 = off2 = 0;
							break;
#if 0
						case( 'K' ):
							pscore = genG__align11( mseq1, mseq2, alloclen );
							off1 = off2 = 0;
							break;
#endif
						case( 'H' ):
							pscore = recallpairfoldalign( mseq1, mseq2, i, j, &off1, &off2, alloclen );
							break;
						case( 'B' ):
						case( 'T' ):
							pscore = recalllara( mseq1, mseq2, alloclen );
							off1 = off2 = 0;
							break;
						case( 's' ):
							pscore = callmxscarna_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
							off1 = off2 = 0;
							break;
						case( 'G' ):
							pscore = calldafs_giving_bpp( mseq1, mseq2, bpp[i], bpp[j], alloclen, i, j );
							off1 = off2 = 0;
							break;
						case( 'M' ):
							pscore = MSalign11( mseq1, mseq2, alloclen );
							break;
						default:
							ErrorExit( "ERROR IN SOURCE FILE" );
					}
				}
	
				if( alg == 't' || ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) ) // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
				{
#if SCOREOUT
					fprintf( stderr, "score = %10.2f (%d,%d)\n", pscore, i, j );
#endif
//					if( pscore > 0.0 && ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) ) // x-ins-i de seido teika
					if( ( nadd == 0 || ( alg != 'Y' && alg != 'r' ) || ( i < njob-nadd && njob-nadd <= j ) ) )
					{
						if( !store_localhom )
							;
						else if( specifictarget && targetmap[i] == -1 && targetmap[j] == -1)
							;
						else if( alg == 'R' )
							putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j, lastresx[i]+j, 'h' );
						else if( alg == 'r' )
							putlocalhom_last( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), lastresx[i]+j-(njob-nadd), 'h' );// ?????
						else if( alg == 'H' )
							putlocalhom_ext( mseq1[0], mseq2[0], localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
						else if( alg == 'Y' )
							putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-(njob-nadd), off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
						else if( !specifictarget && alg != 'S' && alg != 'V' )
							putlocalhom2( mseq1[0], mseq2[0], localhomtable[i]+j-i, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
						else
						{
							if( targetmap[i] != -1 && targetmap[j] != -1 )
							{
								putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
								putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' ); // sukoshi muda.
							}
							else if( targetmap[j] != -1 )
								putlocalhom2( mseq2[0], mseq1[0], localhomtable[targetmap[j]]+i, off2, off1, (int)pscore, strlen( mseq2[0] ), 'h' );
							else if( targetmap[i] != -1 )
								putlocalhom2( mseq1[0], mseq2[0], localhomtable[targetmap[i]]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'h' );
							else
							{
								reporterr( "okashii\n" );
								exit( 1 );
							}
						}
					}

					pscore = score2dist( pscore, selfscore[i], selfscore[j] );
				}
				else
				{
					pscore = 2.0;
				}
	
				if( stdout_align )
				{
					if( alg != 't' )
					{
						fprintf( stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i+1, j+1, pscore );
						fprintf( stdout, ">%s\n", name[i] );
						write1seq( stdout, mseq1[0] );
						fprintf( stdout, ">%s\n", name[j] );
						write1seq( stdout, mseq2[0] );
						fprintf( stdout, "\n" );
					}
				}
				if( stdout_dist ) fprintf( stdout, "%d %d d=%.3f\n", i+1, j+1, pscore );
				if( store_dist) 
				{
					if( alg == 'Y' || alg == 'r' ) distancemtx[i][j-(njob-nadd)] = pscore;
					else distancemtx[i][j-i] = pscore;
				}
			}
		}
		if( dynamicmtx ) FreeDoubleMtx( dynamicmtx );
	}


	if( store_dist && ngui == 0 )
	{
		hat2p = fopen( hat2file, "w" );
		if( !hat2p ) ErrorExit( "Cannot open hat2." );
		if( alg == 'Y' || alg == 'r' )
			WriteHat2_part_pointer( hat2p, njob, nadd, name, distancemtx );
		else
//			WriteHat2_pointer( hat2p, njob, name, distancemtx );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, distancemtx ); // jissiha double
		fclose( hat2p );
	}

	hat3p = fopen( "hat3", "w" );
	if( !hat3p ) ErrorExit( "Cannot open hat3." );
	if( store_localhom && ngui == 0 )
	{

		fprintf( stderr, "\n\n##### writing hat3\n" );
		if( alg == 'Y' || alg == 'r' )
			ilim = njob-nadd;	
		else if( specifictarget )
			ilim = ntarget;
		else
			ilim = njob-1;	
		for( i=0; i<ilim; i++ ) 
		{
			if( alg == 'Y' || alg == 'r' )
			{
				jst = njob-nadd;
				jj = 0;
			}
			else if( specifictarget )
			{
				jst = 0;
				jj = 0;
			}
			else
			{
				jst = i;
				jj = 0;
			}
			for( j=jst; j<njob; j++, jj++ )
			{
				for( tmpptr=localhomtable[i]+jj; tmpptr; tmpptr=tmpptr->next )
				{
//					fprintf( stderr, "j=%d, jj=%d\n", j, jj );
					if( tmpptr->opt == -1.0 ) continue;
// tmptmptmptmptmp
//					if( alg == 'B' || alg == 'T' )
//						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, 1.0, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, (void *)tmpptr->next ); 
//					else
					if( targetmap[j] == -1 || targetmap[i] < targetmap[j] )
						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d h\n", targetmapr[i], j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2 );
//						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d h\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2+1, tmpptr->end2+1 ); // zettai dame!!!!
				}
			}
		}
//		if( ngui == 0 )
//		{
#if DEBUG
			fprintf( stderr, "calling FreeLocalHomTable\n" );
#endif
			if( alg == 'Y' || alg == 'r' )
				FreeLocalHomTable_part( localhomtable, (njob-nadd), nadd );
			else if( specifictarget )
				FreeLocalHomTable_part( localhomtable, ntarget, njob );
			else
				FreeLocalHomTable_half( localhomtable, njob );
#if DEBUG
			fprintf( stderr, "done. FreeLocalHomTable\n" );
#endif
//		}
	}
	fclose( hat3p );

	if( alg == 's' )
	{
		char **ptpt;
		for( i=0; i<njob; i++ )
		{
			ptpt = bpp[i];
			while( 1 )
			{
				if( *ptpt ) free( *ptpt );
				else break;
				ptpt++;
			}
			free( bpp[i] );
		}
		free( bpp );
	}
	free( selfscore );
	free( effarr1 );
	free( effarr2 );
	if( alg == 'N' )
	{
		FreeCharMtx( dumseq1 );
		FreeCharMtx( dumseq2 );
	}
	free( distseq1 );
	free( distseq2 );
	if( store_dist && ngui == 0 ) 
	{
		if( alg == 'Y' || alg == 'r' )
			FreeDoubleMtx( distancemtx ); // 2020/Oct/23
		else
			FreeDoubleHalfMtx( distancemtx, njob );
	}

	free( targetmap );
	free( targetmapr );
}


int pairlocalalign( int ngui, int lgui, char **namegui, char **seqgui, double **distancemtx, LocalHom **localhomtable, int argc, char **argv, double **expdist )
{
	int  *nlen, *thereisxineachseq;
	char **name, **seq;
	char **mseq1, **mseq2;
	char **aseq;
	char **bseq;
	char **dseq;
	int i, j, k;
	FILE *infp;
	char c;
	int alloclen;
	Lastresx **lastresx;

//	reporterr( "argc=%d, argv[0]=%s\n", argc, argv[0] );

	arguments( argc, argv );


	if( !ngui )
	{
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
	
		if( njob < 2 )
		{
			fprintf( stderr, "At least 2 sequences should be input!\n"
							 "Only %d sequence found.\n", njob ); 
			exit( 1 );
		}
		if( njob > M )
		{
			fprintf( stderr, "The number of sequences must be < %d\n", M );
			fprintf( stderr, "Please try --6merpair --addfragments for such large data.\n" );
			exit( 1 );
		}
	}

	if( ( alg == 'r' || alg == 'R' ) && dorp == 'p' )
	{
		fprintf( stderr, "Not yet supported\n" );
		exit( 1 );
	}

	alloclen = nlenmax*2;
	if( ngui ) 
	{
		seq = seqgui;
		name = namegui;
	}
	else
	{
		seq = AllocateCharMtx( njob, alloclen+10 );
		name = AllocateCharMtx( njob, B );
	}

	aseq = AllocateCharMtx( 2, alloclen+10 );
	bseq = AllocateCharMtx( njob, alloclen+10 );
	dseq = AllocateCharMtx( njob, alloclen+10 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	nlen = AllocateIntVec( njob );
	thereisxineachseq = AllocateIntVec( njob );


	if( alg == 'R' )
	{
		lastresx = calloc( njob+1, sizeof( Lastresx * ) );
		for( i=0; i<njob; i++ ) 
		{
			lastresx[i] = calloc( njob+1, sizeof( Lastresx ) ); // muda
			for( j=0; j<njob; j++ ) 
			{
				lastresx[i][j].score = 0;
				lastresx[i][j].naln = 0;
				lastresx[i][j].aln = NULL;
			}
			lastresx[i][njob].naln = -1;
		}
		lastresx[njob] = NULL;
	}
	else if( alg == 'r' )
	{
//		fprintf( stderr, "Allocating lastresx (%d), njob=%d, nadd=%d\n", njob-nadd+1, njob, nadd );
		lastresx = calloc( njob-nadd+1, sizeof( Lastresx * ) );
		for( i=0; i<njob-nadd; i++ )
		{
//			fprintf( stderr, "Allocating lastresx[%d]\n", i );
			lastresx[i] = calloc( nadd+1, sizeof( Lastresx ) );
			for( j=0; j<nadd; j++ ) 
			{
//				fprintf( stderr, "Initializing lastresx[%d][%d]\n", i, j );
				lastresx[i][j].score = 0;
				lastresx[i][j].naln = 0;
				lastresx[i][j].aln = NULL;
			}
			lastresx[i][nadd].naln = -1;
		}
		lastresx[njob-nadd] = NULL;
	}
	else
		lastresx = NULL;

#if 0
	Read( name, nlen, seq );
#else
	if( !ngui ) 
	{
		readData_pointer( infp, name, nlen, seq );
		fclose( infp );
	}
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

//	writePre( njob, name, nlen, seq, 0 );

	//reporterr( "expdist=%p\n", expdist );

	if( dorp == 'p' && scoremtx == 1 && nblosum > 0 ) // protein, not text.  hitsuyou?
	{
		for( i=0; i<njob; i++ ) 
		{
			gappick0( bseq[i], seq[i] );
			thereisxineachseq[i] = removex( dseq[i], bseq[i] );
		}
	}
	else // text, dna
	{
		for( i=0; i<njob; i++ ) 
		{
			gappick0( bseq[i], seq[i] );
			strcpy( dseq[i], bseq[i] );
			thereisxineachseq[i] = 0;
		}
	}

	pairalign( name, nlen, bseq, aseq, dseq, thereisxineachseq, mseq1, mseq2, alloclen, lastresx, distancemtx, localhomtable, expdist, ngui );

	fprintf( trap_g, "done.\n" );
#if DEBUG
	fprintf( stderr, "closing trap_g\n" );
#endif
	fclose( trap_g );
	fclose( prep_g );

//	writePre( njob, name, nlen, aseq, !contin );
#if 0
	writeData( stdout, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif
	SHOWVERSION;

	if( stdout_dist && nthread > 1 )
	{
		fprintf( stderr, "\nThe order of distances is not identical to that in the input file, because of the parallel calculation.  Reorder them by yourself, using sort -n -k 2 | sort -n -k 1 -s\n" );
	}
	if( stdout_align && nthread > 1 )
	{
		fprintf( stderr, "\nThe order of pairwise alignments is not identical to that in the input file, because of the parallel calculation.  Reorder them by yourself.\n" );
	}

#if 1
	if( lastresx ) 
	{
		for( i=0; lastresx[i]; i++ ) 
		{
			for( j=0; lastresx[i][j].naln!=-1; j++ ) 
			{
				for( k=0; k<lastresx[i][j].naln; k++ )
				{
					free( lastresx[i][j].aln[k].reg1 );
					free( lastresx[i][j].aln[k].reg2 );
				}
				free( lastresx[i][j].aln );
			}
			free( lastresx[i] );
		}
		free( lastresx );
	}
#endif
	if( ngui == 0 ) 
	{
		FreeCharMtx( seq );
		FreeCharMtx( name );
	}
	FreeCharMtx( aseq );
	FreeCharMtx( bseq );
	FreeCharMtx( dseq );
	free( mseq1 );
	free( mseq2 );
	free( nlen );
	free( thereisxineachseq );
	freeconstants();

	if( !ngui )
	{
		FreeCommonIP();
	}
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // 20130603
	G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
	L__align11( NULL, 0.0, NULL, NULL, 0, NULL, NULL );
	L__align11_noalign( NULL, NULL, NULL );
	genL__align11( NULL, NULL, NULL, 0, NULL, NULL );

#if SHISHAGONYU
	if( ngui )
	{
		char buf[100];
		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
		{
			sprintf( buf, "%5.3f", distancemtx[i][j-i] );
			distancemtx[i][j-i] = 0.0;
			sscanf( buf, "%lf", distancemtx[i]+j-i );
//			distancemtx[i][j-i] = 0.001 * (int)(distancemtx[i][j-i] * 1000 + 0.5);
		}

	}
#endif


	return( 0 );
}

void pairalign_node( int njob, int nlenmax, char **name, char **seq, int ***topol, double **len, Treedep *dep, int treegiven, int treeout ) // used by nodepair
{
/* test */
	double *selfscore;
//	int ***topol;
//	double **len;
//	Treedep *dep;
	char *fn;
	int i, ntarget;
	int *uselh;
	int *nfilesfornode;
	char **bseq;
	char **dseq;
	int alloclen = nlenmax * 2;
	int alignmentlength;
	FILE *fp;

	bseq = AllocateCharMtx( njob, alloclen+10 );
	dseq = AllocateCharMtx( njob, alloclen+10 );
	uselh = calloc( njob, sizeof( int ) );
	nfilesfornode = calloc( njob-1, sizeof( int ) );

	for( i=0; i<njob; i++ ) 
	{
		gappick0( bseq[i], seq[i] );
		removex( dseq[i], bseq[i] );
	}

	if( nadd )
	{
		if( addprofile )
		{
			reporterr( "--addprofile is not yet supported\n" );
			exit( 1 );
		}
		alignmentlength = strlen( seq[0] );
		for( i=njob-nadd-1; i>0; i-- )
		{
			if( alignmentlength != strlen( seq[i] ) )
			{
				fprintf( stderr, "#################################################################################\n" );
				fprintf( stderr, "# ERROR!                                                                       \n" );
				fprintf( stderr, "# For the --add option, the original%4d sequences must be aligned              \n", njob-nadd );
				fprintf( stderr, "#################################################################################\n" );
				exit( 1 );
			}
		}
	}

	if( specifictarget )
	{
		reporterr( "specifictarget\n" );

		ntarget = 0;
		for( i=0; i<njob; i++ )
		{
			uselh[i] = 0;
			if( !strncmp( name[i]+1, "_focus_", 7 ) )
			{
				uselh[i] = 1;
				ntarget++;
			}
		}

		if( ntarget == 0 )
		{
			reporterr( "\n\nAdd '>_focus_' to the title lines of the sequences to be focused on.\n\n" );
			exit( 1 );
		}
		else
		{
			reporterr( "nfocus = %d \n", ntarget );
		}
	}
	else
	{
		ntarget = njob;
//		targetmap = calloc( njob, sizeof( int ) );
//		targetmapr = calloc( njob, sizeof( int ) );

//		for( i=0; i<njob; i++ ) targetmap[i] = i;
//		if( alg != 'Y' && alg != 'r' )
//			stringshuffle( targetmap, njob );

#if 0
		for( i=0; i<njob; i++ ) uselh[i] = 1;
#else

		char *tmpseq;
		Lennum *tmpstr;
		tmpseq = (char *)calloc( nlenmax+1, sizeof( char ) );
		tmpstr = (Lennum *)calloc( njob, sizeof( Lennum ) );
		for( i=0; i<njob; i++ )
		{
			gappick0( tmpseq, seq[i] );
			tmpstr[i].len = strlen( tmpseq );
			tmpstr[i].num = i;
//			reporterr( "nlen[%d] = %d\n", i, nlen[i] );
//			reporterr( "strlen = %d\n", strlen( seq[i] ) );
//			reporterr( "seq = %s\n", seq[i] );
		}

		limitlh( uselh, tmpstr, njob, lhlimit );

		free( tmpseq );
		free( tmpstr );
#endif

//		for( i=0; i<njob; i++ )
//			reporterr( "targetmap[%d] = %d\n", i, targetmap[i] );
	}

//	for( i=0; i<njob; i++ ) reporterr( "uselh[%d] = %d\n", i, uselh[i] );


	selfscore = AllocateDoubleVec( njob );
	fn = calloc( 100, sizeof( char ) );
//	topol = AllocateIntCub( njob, 2, 0 );
//	len = AllocateFloatMtx( njob, 2 );
//	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );
	system( "rm -rf hat3dir" ); // toriaezu
	system( "mkdir hat3dir" ); // toriaezu

	for( i=0; i<njob-1; i+=HAT3NODEBLOCK ) 
	{
		sprintf( fn, "mkdir \"hat3dir/%d-\"", i ); // windows de slash wo tsukau tame niju in-youfu hitsuyou
		system( fn );
	}
	free( fn );

	for( i=0; i<njob; i++ ) selfscore[i] = (double)naivepairscore11( seq[i], seq[i], 0.0 ); // dseq??



#if 0
	remove( "hat1node" );
	if( (fd1 = open("hat1node",O_RDWR|O_CREAT,0600)) == -1)
	{
		reporterr( "failed to open hat1node\n" );
		exit( 1 );
	}

	remove( "hat0node" );
	if( (fd0 = open("hat0node",O_RDWR|O_CREAT,0600)) == -1)
	{
		reporterr( "failed to open hat0node\n" );
		exit( 1 );
	}

	remove( "hat2node" ); // test only. ato de kesu.
	if( (fd2 = open("hat2node",O_RDWR|O_CREAT,0600)) == -1)
	{
		reporterr( "failed to open hat2node\n" );
		exit( 1 );
	}
#endif

//	reporterr( "Call compacttreedpdist here.\n" );
//	compacttreedpdist( njob, seq, dseq, selfscore, topol, len, name, dep, 1, 1, 2*nlenmax+1, uselh ); // dame
//	compacttreedpdist( njob, bseq, dseq, selfscore, topol, len, name, dep, 1, 2*nlenmax+1, uselh, nfilesfornode, treegiven );
	compacttreedpdist( njob, bseq, dseq, selfscore, topol, len, name, dep, treeout, 2*nlenmax+1, uselh, nfilesfornode, treegiven );


	fp = fopen( "hat3dir/tree", "wb" ); // window no tame wb
	treeout_bin( fp, njob, topol, len, dep, nfilesfornode );
	fclose( fp );

	fp = fopen( "hat3dir/uselh", "wb" ); // nenno tame 
	uselhout( fp, njob, uselh );
	fclose( fp );

#if 0
	if( close( fd0 ) != 0 )
	{
		reporterr( "error in close( hat0node )\n" );
		exit( 1 );
	}

	if( close( fd1 ) != 0 )
	{
		reporterr( "error in close( hat1node )\n" );
		exit( 1 );
	}

	if( close( fd2 ) != 0 )
	{
		reporterr( "error in close( hat2node )\n" );
		exit( 1 );
	}
#endif

//	FreeFloatMtx( len );
//	free( dep );
	free( selfscore );
//	free( targetmap );
//	free( targetmapr );
	free( uselh );
	free( nfilesfornode );
	FreeCharMtx( bseq );
	FreeCharMtx( dseq );
//	FreeIntCub( topol ); topol = NULL; // koko?

/* test owari */
}
