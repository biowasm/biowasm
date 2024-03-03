#include "mltaln.h"
#include "dp.h"

#define DEBUG 0
#define XXXXXXX    0
#define USE_PENALTY_EX  1

#define TERMGAPFAC 0.0
#define TERMGAPFAC_EX 0.0

#if 1

#if 0
static void match_calc_mtx_codon( double **codonmtx, double **mtx, double *match, char **s1, char **s2, int **codonseq1, int **codonseq2, int i1, int lgth2 )  
// gstart, gend wo tsukatte coding ka kakuninn suru beki
{
	char *seq2 = s2[0];
	double *doubleptr = mtx[(unsigned char)s1[0][i1]];
	int codonid1, codonid2;
//	int codonid1p, codonid2p;
//	int codonid1pp, codonid2pp;
	double codonmatch;
//	int *pt0, *pt1, *pt2;
	int *pt0;

	codonid1 = codonseq1[0][i1];
//	codonid1p = codonseq1[1][i1];
//	codonid1pp = codonseq1[2][i1];
	pt0 = codonseq2[0];
//	pt1 = codonseq2[1];
//	pt2 = codonseq2[2];
//	posin2 = 0;
	while( lgth2-- )
	{
//		reporterr( "%c-%c -> %f\n", s1[0][i1], *seq2, doubleptr[(unsigned char)*seq2] );
		
// track ha at de kokoromiru
		codonid2 = *pt0++;
//		codonid2p = *pt1++;
//		codonid2pp = *pt2++;
//		posin2++;
//
//		reporterr( "%c%c%c(%d)-%c%c%c(%d)\n", s1[0][i1], s1[0][i1+1], s1[0][i1+2], codonid1, *seq2, *(seq2+1), *(seq2+2), codonid2 );
		codonmatch = 0.0;
		if( codonid1 > -1 && codonid2 > -1 ) codonmatch = codonmtx[codonid1][codonid2] * 600.0;
//		if( codonid1p > -1 && codonid2p > -1 ) codonmatch -= codonmtx[codonid1p][codonid2p] * 600;
//		if( codonid1pp > -1 && codonid2pp > -1 ) codonmatch -= codonmtx[codonid1pp][codonid2pp] * 600;
		*match++ = doubleptr[(unsigned char)*seq2++] + codonmatch;
	}
}
#else
static void match_calc_mtx_codon( double **codonmtx, double **mtx, double *match, char **s1, char **s2, int **codonseq1, int **codonseq2, int i1, int lgth2 )  
// gstart, gend wo tsukatte coding ka kakuninn suru beki
{
	char *seq2 = s2[0];
	double *doubleptr = mtx[(unsigned char)s1[0][i1]];
	int codonid1, codonid2;
	double codonmatch;
	int *pt0;

	codonid1 = codonseq1[0][i1];
	pt0 = codonseq2[0];
	while( lgth2-- )
	{
//		reporterr( "%c-%c -> %f\n", s1[0][i1], *seq2, doubleptr[(unsigned char)*seq2] );
		
		codonid2 = *pt0++;
//		reporterr( "%c%c%c(%d)-%c%c%c(%d)\n", s1[0][i1], s1[0][i1+1], s1[0][i1+2], codonid1, *seq2, *(seq2+1), *(seq2+2), codonid2 );
		codonmatch = 0.0;
		if( codonid1 > -1 && codonid2 > -1 ) codonmatch = codonmtx[codonid1][codonid2] * 600.0; // codonid1>-1 ha iranai.
		*match++ = doubleptr[(unsigned char)*seq2++] + codonmatch;
	}
}
#endif

static void match_calc_mtx( double **mtx, double *match, char **s1, char **s2, int i1, int lgth2 ) 
{
	char *seq2 = s2[0];
	double *doubleptr = mtx[(unsigned char)s1[0][i1]];

	while( lgth2-- )
		*match++ = doubleptr[(unsigned char)*seq2++];
}
#else
static void match_calc( double *match, char **s1, char **s2, int i1, int lgth2 )
{
	int j;

	for( j=0; j<lgth2; j++ )
		match[j] = amino_dis[(*s1)[i1]][(*s2)[j]];
}
#endif

static double Atracking( double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp,
						int tailgp,
						int *warpis, int *warpjs, int warpbase )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
//	char gap[] = "-";
	char *gap;
	gap = newgapstr;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	double wm, g;
	double fpenalty = (double)penalty;
	double fpenalty_ex = (double)penalty_ex;


#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }


//	if( tailgp == 1 || ijp[lgth1][lgth2] >= warpbase )
	if( tailgp == 1 )
		;
	else
	{
#if 1
//		reporterr( "lastverticalw[lgth1-1] = %f\n", lastverticalw[lgth1-1] );
//		reporterr( "lasthorizontalw[lgth2-1] = %f\n", lasthorizontalw[lgth2-1] );
		wm = lasthorizontalw[lgth2-1] - 1.0; // lasthorizontalw[lgth2-1] yori kanarazu chiisai.
		for( j=lgth2-2; j>=0; j-- )
		{
			if( (g=lasthorizontalw[j]+ ( fpenalty * TERMGAPFAC + fpenalty_ex * (lgth2-1-j) * TERMGAPFAC_EX ) ) > wm )
			{
				wm = g;
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
		for( i=lgth1-2; i>=0; i-- )
		{
			if( ( g=lastverticalw[i]+ ( fpenalty * TERMGAPFAC + fpenalty_ex * (lgth1-1-i) * TERMGAPFAC_EX ) ) > wm )
			{
				wm = g;
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		if( lasthorizontalw[lgth2-1] > wm )  // score ga onaji baai erabarenai
		{
			wm = lasthorizontalw[lgth2-1];
			iin = lgth1-1; jin = lgth2-1;
			ijp[lgth1][lgth2] = 0;
		}
#else
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
#endif
	}



	mseq1[0] += lgth1+lgth2;
	*mseq1[0] = 0;
	mseq2[0] += lgth1+lgth2;
	*mseq2[0] = 0;



	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] >= warpbase )
		{
//			fprintf( stderr, "WARP!\n" );
			ifi = warpis[ijp[iin][jin]-warpbase]; 
			jfi = warpjs[ijp[iin][jin]-warpbase];
		}
		else if( ijp[iin][jin] < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijp[iin][jin];
		}
		else if( ijp[iin][jin] > 0 )
		{
			ifi = iin-ijp[iin][jin]; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}

		if( ifi == -warpbase && jfi == -warpbase )
		{
			l = iin;
			while( --l >= 0 ) 
			{
				*--mseq1[0] = seq1[0][l];
				*--mseq2[0] = *gap;
				k++;
			}
			l= jin;
			while( --l >= 0 )
			{
				*--mseq1[0] = *gap;
				*--mseq2[0] = seq2[0][l];
				k++;
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l > 0 ) 
			{
				*--mseq1[0] = seq1[0][ifi+l];
				*--mseq2[0] = *gap;
				k++;
			}
			l= jin - jfi;
			while( --l > 0 )
			{
				*--mseq1[0] = *gap;
				*--mseq2[0] = seq2[0][jfi+l];
				k++;
			}
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--mseq1[0] = seq1[0][ifi];
		*--mseq2[0] = seq2[0][jfi];
		k++;
		iin = ifi; jin = jfi;
	}

//	fprintf( stderr, "%s\n", mseq1[0] );
//	fprintf( stderr, "%s\n", mseq2[0] );
	return( wm );
}


double G__align11psg( double **codonmtx, double **n_dynamicmtx, char **seq1, char **seq2, int alloclen, int headgp, int tailgp, double *gstart, double *gend )
{



//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lastj;
	int lgth1, lgth2;
	int resultlen;
	double wm, wmo;   /* int ?????? */
	double g;
	double *currentw, *previousw;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	double fpenalty_tmp;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
	double fpenalty_ex_i;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi = 0.0;
	static TLS double *m = NULL;
	static TLS int **ijp = NULL;
	static TLS int mpi = 0;
	static TLS int *mp = NULL;
	static TLS double *w1 = NULL;
	static TLS double *w2 = NULL;
	static TLS double *match = NULL;
	static TLS double *initverticalw = NULL;    /* kufuu sureba iranai */
	static TLS double *lastverticalw = NULL;    /* kufuu sureba iranai */
	static TLS char **mseq1 = NULL;
	static TLS char **mseq2 = NULL;
	static TLS char **mseq = NULL;
	static TLS int **intwork = NULL;
	static TLS double **doublework = NULL;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx = NULL; // ??
	static TLS int **codonseq1 = NULL;
	static TLS int **codonseq2 = NULL;

	int *warpis = NULL;
	int *warpjs = NULL;
	int *warpi = NULL;
	int *warpj = NULL;
	int *prevwarpi = NULL;
	int *prevwarpj = NULL;
	double *wmrecords = NULL;
	double *prevwmrecords = NULL;
	int warpn = 0;
	int warpbase;
	double curm = 0.0;
	double *wmrecordspt, *wmrecords1pt, *prevwmrecordspt;
	int *warpipt, *warpjpt;

	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			if( mseq1 ) free( mseq1 ); mseq1 = NULL;
			if( mseq2 ) free( mseq2 ); mseq2 = NULL;
			if( w1 ) FreeFloatVec( w1 ); w1 = NULL;
			if( w2 ) FreeFloatVec( w2 ); w2 = NULL;
			if( match ) FreeFloatVec( match ); match = NULL;
			if( initverticalw ) FreeFloatVec( initverticalw ); initverticalw = NULL;
			if( lastverticalw ) FreeFloatVec( lastverticalw ); lastverticalw = NULL;

			if( m ) FreeFloatVec( m ); m = NULL;
			if( mp ) FreeIntVec( mp ); mp = NULL;

			if( mseq ) FreeCharMtx( mseq ); mseq = NULL;



			if( doublework ) FreeFloatMtx( doublework ); doublework = NULL;
			if( intwork ) FreeIntMtx( intwork ); intwork = NULL;

			if( amino_dynamicmtx ) FreeDoubleMtx( amino_dynamicmtx ); amino_dynamicmtx = NULL;

			if( codonseq1 ) FreeIntMtx( codonseq1 ); codonseq1 = NULL;
			if( codonseq2 ) FreeIntMtx( codonseq2 ); codonseq2 = NULL;
		}
		orlgth1 = 0;
		orlgth2 = 0;
		return( 0.0 );
	}



	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );



	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;
	if( trywarp )
	{
//		fprintf( stderr, "IN G__align11\n" );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1.\n" );
			exit( 1 );
		}

		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}

#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

#if 1
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		seq1[0][lgth2] = 0;
		while( lgth2 ) seq1[0][--lgth2] = *newgapstr;
//		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		seq2[0][lgth1] = 0;
		while( lgth1 ) seq2[0][--lgth1] = *newgapstr;
//		fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
		return( 0.0 );
	}
#endif


	wm = 0.0;

	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( 2, 0 ); // 2020/Apr
		mseq2 = AllocateCharMtx( 2, 0 ); // 2020/Apr
	}



	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );



			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );

			if( codonseq1 ) FreeIntMtx( codonseq1 );
			if( codonseq2 ) FreeIntMtx( codonseq2 );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( 2, ll1+ll2 ); // 2020/Apr


		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		amino_dynamicmtx = AllocateDoubleMtx( 0x100, 0x100 );

		if( codonscore )
		{
			codonseq1 = AllocateIntMtx( 3, ll1 ); // Only codonseq1[0] is used at this point
			codonseq2 = AllocateIntMtx( 3, ll2 ); // Only codonseq2[0] is used at this point
		}

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}
    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];


	mseq1[0] = mseq[0];
	mseq2[0] = mseq[1];


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeIntMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;

	if( codonscore )
	{
// codonseq[1..2] are not used yet.
		for( i=0; i<3&&i<lgth1; i++ )
		{
				codonseq1[0][i] = -1;
				codonseq1[1][i] = -1;
				codonseq1[2][i] = -1;
		}
		for( i=0; i<3&&i<lgth2; i++ )
		{
				codonseq2[0][i] = -1;
				codonseq2[1][i] = -1;
				codonseq2[2][i] = -1;
		}
		for( i=3; i<lgth1; i++ )
		{
			codonseq1[2][i] = codon2id( seq1[0]+i );
			codonseq1[1][i] = codon2id( seq1[0]+i-1 );
			codonseq1[0][i] = codon2id( seq1[0]+i-2 );
		}
		for( i=3; i<lgth2; i++ )
		{
			codonseq2[2][i] = codon2id( seq2[0]+i );
			codonseq2[1][i] = codon2id( seq2[0]+i-1 );
			codonseq2[0][i] = codon2id( seq2[0]+i-2 );
		}
	}
	else
	{
		codonseq1 = NULL;
		codonseq2 = NULL;
	}

#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	if( codonscore )
	{

//		for( i=0; i<5; i++ ) for( j=0; j<5; j++ ) reporterr( "score[%c][%c]=%f\n", amino[i], amino[j], amino_dynamicmtx[amino[i]][amino[j]] );
//		exit( 1 );
		match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );
//		if( !codonpos || ! (gstart[0] == 0.5 && gend[0] == 0.5) ) // consider codon when unknown or coding
//			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, 0, lgth2 );
#if 1
		if( gstart[0] == 0.5 && gend[0] > 0.5 ) // only at 3rd position
			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, 0, lgth2 );
		else
			match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );
#else
		match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, 0, lgth2, gstart[0] == 0.5 && gend[0] > 0.5 );
#endif
	}
	else
	{
		match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );
		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );
	}

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
#if USE_PENALTY_EX
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
#if USE_PENALTY_EX
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
	}
	else
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
	}


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1]; mp[j] = 0;
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif


	for( i=1; i<lasti; i++ )
	{
//		reporterr( "seq1[0][%d]=%c, gstart[i]=%f, gend[i]=%f\n", i, seq1[0][i], gstart[i], gend[i] );
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];


//		if( codonscore && ( !codonpos || !( gstart[i] == 0.5 && gend[i] == 0.5 ) ) ) // unknown or coding
//			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, i, lgth2 );
#if 1
		if( codonscore && ( gstart[i] == 0.5 && gend[i] > 0.5 ) ) // only at 3rd position in seq1
			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, i, lgth2 );
		else
			match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#else
		if( codonscore ) // only at 3rd position in seq1
			match_calc_mtx_codon( codonmtx, amino_dynamicmtx, currentw, seq1, seq2, codonseq1, codonseq2, i, lgth2, gstart[i] == 0.5 && gend[i] > 0.5 );
		else
			match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0]; mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}

		if( i < lgth1 ) fpenalty_ex_i = fpenalty_ex;
		else fpenalty_ex_i = 0.0; // 2018/May/11

		for( j=1; j<lastj; j++ )
		{

			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty*gend[i]) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
//				reporterr( "hit! jump from %d to %d: %c->%c\n", j, mpi, seq2[0][j], seq2[0][mpi] );
			}
			if( (g=*prept+fpenalty*gstart[i-1]) >= mi )
//			if( (g=*prept) > mi )
			{
				mi = g;
				mpi = j-1;
//				reporterr( "hit! jump to %d: ->%c:%c\n", mpi, seq1[0][i-1], seq2[0][mpi] );
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex_i;
//			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty*gend[i]) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			if( (g=*prept+fpenalty*gstart[i-1]) >= *mjpt )
//			if( (g=*prept) > *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			if( j < lgth2 ) // 2018/May/11
				m[j] += fpenalty_ex;
#endif
#if 1
			if( trywarp )
			{
				fpenalty_tmp = fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] );
//				fprintf( stderr, "fpenalty_shift = %f\n", fpenalty_tmp );

//				fprintf( stderr, "\n\n\nwarp to %c-%c (%d-%d) from %c-%c (%d-%d) ? prevwmrecords[%d] = %f + %f <- wm = %f\n", seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], prevwarpi[j-1], prevwarpj[j-1], seq1[0][i], seq2[0][j], i, j, j, prevwmrecords[j-1], fpenalty_tmp, wm );
//				if( (g=prevwmrecords[j-1] + fpenalty_shift )> wm )
				if( ( g=*prevwmrecordspt++ + fpenalty_tmp )> wm ) // naka ha osokute kamawanai
				{
//					fprintf( stderr, "Yes! Warp!! from %d-%d (%c-%c) to  %d-%d (%c-%c) fpenalty_tmp = %f! warpn = %d\n", i, j, seq1[0][i], seq2[0][j-1], prevwarpi[j-1], prevwarpj[j-1],seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], fpenalty_tmp, warpn );
					if( warpn && prevwarpi[j-1] == warpis[warpn-1] && prevwarpj[j-1] == warpjs[warpn-1] )
					{
						*ijppt = warpbase + warpn - 1;
					}
					else
					{
						*ijppt = warpbase + warpn;
						warpis = realloc( warpis, sizeof(int) * ( warpn+1 ) );
						warpjs = realloc( warpjs, sizeof(int) * ( warpn+1 ) );
						warpis[warpn] = prevwarpi[j-1];
						warpjs[warpn] = prevwarpj[j-1];
						warpn++;
					}
					wm = g;
				}
				else
				{
				}

				curm = *curpt + wm;
	
//				fprintf( stderr, "###### curm = %f at %c-%c, i=%d, j=%d\n", curm, seq1[0][i], seq2[0][j], i, j ); 
	
//				fprintf( stderr, "copy from i, j-1? %f > %f?\n", wmrecords[j-1], curm );
//				if( wmrecords[j-1] > wmrecords[j] )
				if( *wmrecords1pt > *wmrecordspt )
				{
//					fprintf( stderr, "yes\n" );
//					wmrecords[j] = wmrecords[j-1];
					*wmrecordspt = *wmrecords1pt;
//					warpi[j] = warpi[j-1];
//					warpj[j] = warpj[j-1];
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
//					fprintf( stderr, "warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
				}
//				else
//				{
//					fprintf( stderr, "no\n" );
//				}
	
//				fprintf( stderr, " curm = %f at %c-%c\n", curm, seq1[0][i], seq2[0][j] ); 
//				fprintf( stderr, " wmrecords[%d] = %f\n", j, wmrecords[j] ); 
//				fprintf( stderr, "replace?\n" );
	
//				if( curm > wmrecords[j] )
				if( curm > *wmrecordspt )
				{
//					fprintf( stderr, "yes at %d-%d (%c-%c), replaced warp: warpi[j]=%d, warpj[j]=%d warpn=%d, wmrecords[j] = %f -> %f\n", i, j, seq1[0][i], seq2[0][j], i, j, warpn, wmrecords[j], curm );
//					wmrecords[j] = curm;
					*wmrecordspt = curm;
//					warpi[j] = i;
//					warpj[j] = j;
					*warpipt = i;
					*warpjpt = j;
				}
//				else
//				{
//					fprintf( stderr, "No! warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
//				}
//				fprintf( stderr, "%d-%d (%c-%c) curm = %5.0f, wmrecords[j]=%f\n", i, j, seq1[0][i], seq2[0][j], curm, wmrecords[j] );
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}
#endif

			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}


	if( trywarp )
	{
//		fprintf( stderr, "\nwm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}

	wmo = Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, tailgp, warpis, warpjs, warpbase );
	if( !tailgp ) wm = wmo;

//	reporterr( "wm (after tracking) = %f\n", wm );
	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );


	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	strcpy( seq1[0], mseq1[0] );
	strcpy( seq2[0], mseq2[0] );
#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm = %f\n", wm );
#endif
	return( wm );
}
double G__align11( double **n_dynamicmtx, char **seq1, char **seq2, int alloclen, int headgp, int tailgp )
{

//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lastj;
	int lgth1, lgth2;
	int resultlen;
	double wm, wmo;   /* int ?????? */
	double g;
	double *currentw, *previousw;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	double fpenalty_tmp;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
	double fpenalty_ex_i;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi = 0.0;
	static TLS double *m = NULL;
	static TLS int **ijp = NULL;
	static TLS int mpi = 0;
	static TLS int *mp = NULL;
	static TLS double *w1 = NULL;
	static TLS double *w2 = NULL;
	static TLS double *match = NULL;
	static TLS double *initverticalw = NULL;    /* kufuu sureba iranai */
	static TLS double *lastverticalw = NULL;    /* kufuu sureba iranai */
	static TLS char **mseq1 = NULL;
	static TLS char **mseq2 = NULL;
	static TLS char **mseq = NULL;
	static TLS int **intwork = NULL;
	static TLS double **doublework = NULL;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx = NULL; // ??

	int *warpis = NULL;
	int *warpjs = NULL;
	int *warpi = NULL;
	int *warpj = NULL;
	int *prevwarpi = NULL;
	int *prevwarpj = NULL;
	double *wmrecords = NULL;
	double *prevwmrecords = NULL;
	int warpn = 0;
	int warpbase;
	double curm = 0.0;
	double *wmrecordspt, *wmrecords1pt, *prevwmrecordspt;
	int *warpipt, *warpjpt;


	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			if( mseq1 ) free( mseq1 ); mseq1 = NULL;
			if( mseq2 ) free( mseq2 ); mseq2 = NULL;
			if( w1 ) FreeFloatVec( w1 ); w1 = NULL;
			if( w2 ) FreeFloatVec( w2 ); w2 = NULL;
			if( match ) FreeFloatVec( match ); match = NULL;
			if( initverticalw ) FreeFloatVec( initverticalw ); initverticalw = NULL;
			if( lastverticalw ) FreeFloatVec( lastverticalw ); lastverticalw = NULL;

			if( m ) FreeFloatVec( m ); m = NULL;
			if( mp ) FreeIntVec( mp ); mp = NULL;

			if( mseq ) FreeCharMtx( mseq ); mseq = NULL;



			if( doublework ) FreeFloatMtx( doublework ); doublework = NULL;
			if( intwork ) FreeIntMtx( intwork ); intwork = NULL;

			if( amino_dynamicmtx ) FreeDoubleMtx( amino_dynamicmtx ); amino_dynamicmtx = NULL;
		}
		orlgth1 = 0;
		orlgth2 = 0;
		return( 0.0 );
	}



	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;
	if( trywarp )
	{
//		fprintf( stderr, "IN G__align11\n" );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1.\n" );
			exit( 1 );
		}

		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}

#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

#if 1
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		seq1[0][lgth2] = 0;
		while( lgth2 ) seq1[0][--lgth2] = *newgapstr;
//		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		seq2[0][lgth1] = 0;
		while( lgth1 ) seq2[0][--lgth1] = *newgapstr;
//		fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
		return( 0.0 );
	}
#endif


	wm = 0.0;

	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( 2, 0 ); // 2020/Apr
		mseq2 = AllocateCharMtx( 2, 0 ); // 2020/Apr
	}



	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );



			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( 2, ll1+ll2 ); // 2020/Apr


		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		amino_dynamicmtx = AllocateDoubleMtx( 0x100, 0x100 );

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}
    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(unsigned char)amino[i]][(unsigned char)amino[j]] = (double)n_dynamicmtx[i][j];


	mseq1[0] = mseq[0];
	mseq2[0] = mseq[1];


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeIntMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;


#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );
	match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
#if USE_PENALTY_EX
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
#if USE_PENALTY_EX
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
	}
	else
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX // 2018/Apr/22
			currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
		}
	}


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1]; mp[j] = 0;
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif

	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0]; mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}

		if( i < lgth1 ) fpenalty_ex_i = fpenalty_ex;
		else fpenalty_ex_i = 0.0; // 2018/May/11

		for( j=1; j<lastj; j++ )
		{

			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			if( (g=*prept) >= mi )
//			if( (g=*prept) > mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex_i;
//			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			if( (g=*prept) >= *mjpt )
//			if( (g=*prept) > *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			if( j < lgth2 ) // 2018/May/11
				m[j] += fpenalty_ex;
#endif
#if 1
			if( trywarp )
			{
				fpenalty_tmp = fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] );
//				fprintf( stderr, "fpenalty_shift = %f\n", fpenalty_tmp );

//				fprintf( stderr, "\n\n\nwarp to %c-%c (%d-%d) from %c-%c (%d-%d) ? prevwmrecords[%d] = %f + %f <- wm = %f\n", seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], prevwarpi[j-1], prevwarpj[j-1], seq1[0][i], seq2[0][j], i, j, j, prevwmrecords[j-1], fpenalty_tmp, wm );
//				if( (g=prevwmrecords[j-1] + fpenalty_shift )> wm )
				if( ( g=*prevwmrecordspt++ + fpenalty_tmp )> wm ) // naka ha osokute kamawanai
				{
//					fprintf( stderr, "Yes! Warp!! from %d-%d (%c-%c) to  %d-%d (%c-%c) fpenalty_tmp = %f! warpn = %d\n", i, j, seq1[0][i], seq2[0][j-1], prevwarpi[j-1], prevwarpj[j-1],seq1[0][prevwarpi[j-1]], seq2[0][prevwarpj[j-1]], fpenalty_tmp, warpn );
					if( warpn && prevwarpi[j-1] == warpis[warpn-1] && prevwarpj[j-1] == warpjs[warpn-1] )
					{
						*ijppt = warpbase + warpn - 1;
					}
					else
					{
						*ijppt = warpbase + warpn;
						warpis = realloc( warpis, sizeof(int) * ( warpn+1 ) );
						warpjs = realloc( warpjs, sizeof(int) * ( warpn+1 ) );
						warpis[warpn] = prevwarpi[j-1];
						warpjs[warpn] = prevwarpj[j-1];
						warpn++;
					}
					wm = g;
				}
				else
				{
				}

				curm = *curpt + wm;
	
//				fprintf( stderr, "###### curm = %f at %c-%c, i=%d, j=%d\n", curm, seq1[0][i], seq2[0][j], i, j ); 
	
//				fprintf( stderr, "copy from i, j-1? %f > %f?\n", wmrecords[j-1], curm );
//				if( wmrecords[j-1] > wmrecords[j] )
				if( *wmrecords1pt > *wmrecordspt )
				{
//					fprintf( stderr, "yes\n" );
//					wmrecords[j] = wmrecords[j-1];
					*wmrecordspt = *wmrecords1pt;
//					warpi[j] = warpi[j-1];
//					warpj[j] = warpj[j-1];
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
//					fprintf( stderr, "warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
				}
//				else
//				{
//					fprintf( stderr, "no\n" );
//				}
	
//				fprintf( stderr, " curm = %f at %c-%c\n", curm, seq1[0][i], seq2[0][j] ); 
//				fprintf( stderr, " wmrecords[%d] = %f\n", j, wmrecords[j] ); 
//				fprintf( stderr, "replace?\n" );
	
//				if( curm > wmrecords[j] )
				if( curm > *wmrecordspt )
				{
//					fprintf( stderr, "yes at %d-%d (%c-%c), replaced warp: warpi[j]=%d, warpj[j]=%d warpn=%d, wmrecords[j] = %f -> %f\n", i, j, seq1[0][i], seq2[0][j], i, j, warpn, wmrecords[j], curm );
//					wmrecords[j] = curm;
					*wmrecordspt = curm;
//					warpi[j] = i;
//					warpj[j] = j;
					*warpipt = i;
					*warpjpt = j;
				}
//				else
//				{
//					fprintf( stderr, "No! warpi[j]=%d, warpj[j]=%d wmrecords[j] = %f\n", warpi[j], warpj[j], wmrecords[j] );
//				}
//				fprintf( stderr, "%d-%d (%c-%c) curm = %5.0f, wmrecords[j]=%f\n", i, j, seq1[0][i], seq2[0][j], curm, wmrecords[j] );
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}
#endif

			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}


	if( trywarp )
	{
//		fprintf( stderr, "\nwm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}

	wmo = Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, tailgp, warpis, warpjs, warpbase );
	if( !tailgp ) wm = wmo;

//	reporterr( "wm (after tracking) = %f\n", wm );
	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );


	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	strcpy( seq1[0], mseq1[0] );
	strcpy( seq2[0], mseq2[0] );
#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm = %f\n", wm );
#endif

	return( wm );
}

double G__align11_noalign( double **n_dynamicmtx, int penal, int penal_ex, char **seq1, char **seq2, int alloclen )
/* warp mitaiou */
{
//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
//	int resultlen;
	double wm;   /* int ?????? */
	double g;
	double *currentw, *previousw;
	double fpenalty = (double)penal;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penal_ex;
	double fpenalty_ex_i;
#endif
#if 1
	double *wtmp;
	double *mjpt, *prept, *curpt;
//	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS int **intwork;
	static TLS double **doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx;

	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );
			free( m );


			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );
		}
		return( 0.0 );
	}


	wm = 0.0;



	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );



#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );




			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );

			FreeDoubleMtx( amino_dynamicmtx );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );



		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 


//		amino_dynamicmtx = AllocateDoubleMtx( 0x80, 0x80 ); 
		amino_dynamicmtx = AllocateDoubleMtx( 0x100, 0x100 );  // 2017/Nov.  constants.c no 'charsize' wo global hensuu nishita houga yoi?
#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}


    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(int)amino[i]][(int)amino[j]] = (double)n_dynamicmtx[i][j];


#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );


	match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );

	if( 1 ) // tsuneni outgap-1
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
#if USE_PENALTY_EX // 2018/Apr/23
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
#if USE_PENALTY_EX // 2018/Apr/23
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1];
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( 1 ) lasti = lgth1+1; else lasti = lgth1; // tsuneni outgap-1

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif

	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0];

		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;

		if( i < lgth1 ) fpenalty_ex_i = fpenalty_ex;
		else fpenalty_ex_i = 0.0; // 2018/May/11

		for( j=1; j<lgth2+1; j++ )
		{
			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty) > wm )
			{
				wm = g;
			}
//			if( (g=*prept) >= mi )
			if( (g=*prept) > mi ) // onaji hazu
			{
				mi = g;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex_i;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty) > wm )
			{
				wm = g;
			}
//			if( (g=*prept) >= *mjpt )
			if( (g=*prept) > *mjpt ) // onaji hazu
			{
				*mjpt = g;
			}
#if USE_PENALTY_EX
			if( j < lgth2 ) // 2018/May/11
				m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt++ += wm;
			mjpt++;
			prept++;
		}
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error
	}

#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm (noalign) = %f\n", wm );
#endif

	return( wm );
}
