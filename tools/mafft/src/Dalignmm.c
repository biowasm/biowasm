#include "mltaln.h"
#include "dp.h"

#define DEBUG 0

#define WMCHECK 1
#define ALGZSTRAIGHT 0
#define ALGZGAP 0
#define USEGAPLENMTX 0
#define USEGAPLENHALF 0
#define FREEFREQUENTLY 1

#define IDATEND 0


#define MACHIGAI 0
#define OUTGAP0TRY 0
#define XXXXXXX    0
#define USE_PENALTY_EX  0
#define FASTMATCHCALC 1
#define SLOW 0

#define zero 0
#define one 1

#if USEGAPLENHALF
#define USEGAPLENHALFORMTX 1
#endif
#if USEGAPLENMTX
#define USEGAPLENHALFORMTX 1
#endif


#if WMCHECK
static int PFACERROR = 0;
#endif


static TLS double **impmtx = NULL;
static TLS int impalloclen = 0;

double imp_match_out_scD( int i1, int j1 )
{
//	fprintf( stderr, "imp+match = %f\n", impmtx[i1][j1] * fastathreshold );
//	fprintf( stderr, "val = %f\n", impmtx[i1][j1] );
	return( impmtx[i1][j1] );
}

typedef struct _gaplenvec
{
	int relend;
#if USEGAPLENHALFORMTX
	int idatend;
#endif
	int idatnext;
	int idatprev;
	int npat;
	int len; // sukoshi muda.
	double freq;
} Gaplen;

#if 0
static void imp_match_out_vead_gapmap( double *imp, int i1, int lgth2, int *gapmap2 )
{
#if FASTMATCHCALC
	double *pt = impmtx[i1];
	int *gapmappt = gapmap2;
	while( lgth2-- )
		*imp++ += pt[*gapmappt++];
#else
	int j;
	double *pt = impmtx[i1];
	for( j=0; j<lgth2; j++ )
		*imp++ += pt[gapmap2[j]];
#endif
}


static void imp_match_out_vead_tate_gapmap( double *imp, int j1, int lgth1, int *gapmap1 )
{
#if FASTMATCHCALC
	int *gapmappt = gapmap1;
	while( lgth1-- )
		*imp++ += impmtx[*gapmappt++][j1];
#else
	int i;
	for( i=0; i<lgth1; i++ )
		*imp++ += impmtx[gapmap1[i]][j1];
#endif
}
#endif

static void imp_match_out_vead( double *imp, int i1, int lgth2 )
{
#if FASTMATCHCALC 
	double *pt = impmtx[i1];
	while( lgth2-- )
		*imp++ += *pt++;
#else
	int j;
	double *pt = impmtx[i1];
	for( j=0; j<lgth2; j++ )
		*imp++ += pt[j];
#endif
}
static void imp_match_out_vead_tate( double *imp, int j1, int lgth1 )
{
	int i;
	for( i=0; i<lgth1; i++ )
		*imp++ += impmtx[i][j1];
}

void imp_rnaD( int nseq1, int nseq2, char **seq1, char **seq2, double *eff1, double *eff2, RNApair ***grouprna1, RNApair ***grouprna2, int *gapmap1, int *gapmap2, RNApair *pair )
{
	foldrna( nseq1, nseq2, seq1, seq2, eff1, eff2, grouprna1, grouprna2, impmtx, gapmap1, gapmap2, pair );
}


void imp_match_init_strictD( double *imp, int clus1, int clus2, int lgth1, int lgth2, char **seq1, char **seq2, double *eff1, double *eff2, double *eff1_kozo, double *eff2_kozo, LocalHom ***localhom, char *swaplist, int forscore, int *orinum1, int *orinum2, int *uselh, int *seedinlh1, int *seedinlh2, int nodeid, int nfiles )
{
//	int i, j, k1, k2, tmpint, start1, start2, end1, end2;
//	double effij;
//	double effij_kozo;
//	double effijx;
//	char *pt, *pt1, *pt2;
//	static TLS char *nocount1 = NULL;
//	static TLS char *nocount2 = NULL;
//	LocalHom *tmpptr;

	if( seq1 == NULL )
	{
		if( impmtx ) FreeFloatMtx( impmtx );
		impmtx = NULL;
//		if( nocount1 ) free( nocount1 );
//		nocount1 = NULL;
//		if( nocount2 ) free( nocount2 );
//		nocount2 = NULL;
		
		return;
	}

	if( impalloclen < lgth1 + 2 || impalloclen < lgth2 + 2 )
	{
		if( impmtx ) FreeFloatMtx( impmtx );
//		if( nocount1 ) free( nocount1 );
//		if( nocount2 ) free( nocount2 );
		impalloclen = MAX( lgth1, lgth2 ) + 2;
		impmtx = AllocateFloatMtx( impalloclen, impalloclen );
//		nocount1 = AllocateCharVec( impalloclen );
//		nocount2 = AllocateCharVec( impalloclen );
	}

	if( nodeid == -1 )
		fillimp( impmtx, imp, clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, eff1_kozo, eff2_kozo, localhom, swaplist, forscore, orinum1, orinum2 );
	else
		fillimp_file( impmtx, imp, clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, eff1_kozo, eff2_kozo, localhom, swaplist, forscore, orinum1, orinum2, uselh, seedinlh1, seedinlh2, nodeid, nfiles );
}





static void match_calc_del( int **which, double ***matrices, double *match, int n1, char **seq1, double *eff1, int n2, char **seq2, double *eff2, int i1, int lgth2, int mid, int nmask, int *mask1, int *mask2 ) 
{
// osoi!
	int i, j, k, m;
	int c1, c2;
//	fprintf( stderr, "\nmatch_calc_dynamicmtx... %d", i1 );
//	fprintf( stderr, "\nseq1[0]=%s\n", seq1[0] );
//	fprintf( stderr, "\nseq2[0]=%s\n", seq2[0] );
//	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
//	{
//		if( flip ) reporterr( "in match_calc_slow, which[%d][%d] = %d\n", j, i, which[j][i] );
//		else       reporterr( "in match_calc_slow, which[%d][%d] = %d\n", i, j, which[i][j] );
//	}
	for( k=0; k<lgth2; k++ )
	{
		for( m=0; m<nmask; m++ )
		{
			i = mask1[m];
			j = mask2[m];
//			reporterr( "Deleting %d-%d (c=%d)\n", i, j, mid );
//			if( k==0 ) fprintf( stderr, "pairoffset[%d][%d] = %f\n", i, j, po );
			c1 = amino_n[(unsigned char)seq1[i][i1]];
			c2 = amino_n[(unsigned char)seq2[j][k]];
//			reporterr( "k=%d, c1=%d, c2=%d, seq1[i][i1]=%c, seq2[%d][%d]=%c\n", k, c1, c2, seq1[i][i1], j, k, seq2[j][k] );
			if( seq1[i][i1] == '-' || seq2[j][k] == '-' ) continue;
			if( c1 < 0 || c2 < 0 ) continue;
//			fprintf( stderr, "c1=%d, c2=%d\n", c1, c2 );
//			fprintf( stderr, "match[k] = %f -> ", match[k], mid );
			match[k] -= matrices[mid][c1][c2] * eff1[i] * eff2[j];
//			fprintf( stderr, "match[k] = %f (mid=%d)\n", match[k], mid );
		}
	}
//	fprintf( stderr, "done\n" );
	return;
}

#if SLOW
static void match_calc_slow( int **which, double ***matrices, double *match, int n1, char **seq1, double *eff1, int n2, char **seq2, double *eff2, int i1, int lgth2, double **doublework, int **intwork, int initialize, int flip ) 
{
// osoi!
	int i, j, k;
	int c1, c2;
	int mid;
//	fprintf( stderr, "\nmatch_calc_dynamicmtx... %d", i1 );
//	fprintf( stderr, "\nseq1[0]=%s\n", seq1[0] );
//	fprintf( stderr, "\nseq2[0]=%s\n", seq2[0] );
//	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
//	{
//		if( flip ) reporterr( "in match_calc_slow, which[%d][%d] = %d\n", j, i, which[j][i] );
//		else       reporterr( "in match_calc_slow, which[%d][%d] = %d\n", i, j, which[i][j] );
//	}
	for( k=0; k<lgth2; k++ )
	{
		match[k] = 0.0;
		for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
		{
			if( flip ) mid = which[j][i];
			else       mid = which[i][j];
//			if( k==0 ) fprintf( stderr, "pairoffset[%d][%d] = %f\n", i, j, po );
			c1 = amino_n[(unsigned char)seq1[i][i1]];
			c2 = amino_n[(unsigned char)seq2[j][k]];
			if( seq1[i][i1] == '-' || seq2[j][k] == '-' ) continue;
			if( c1 < 0 || c2 < 0 ) continue;
//			fprintf( stderr, "c1=%d, c2=%d\n", c1, c2 );
			if( flip ) 
				match[k] += matrices[mid][c1][c2] * eff1[i] * eff2[j];
			else
				match[k] += matrices[mid][c1][c2] * eff1[i] * eff2[j];
//			fprintf( stderr, "match[k] = %f (which=%d)\n", match[k], mid );
		}
	}
//	fprintf( stderr, "done\n" );
	return;
}
#endif

static void fillzero( double *s, int l )
{
	while( l-- ) *s++ = 0.0;
}


static void match_calc_add( double **scoreingmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
#if FASTMATCHCALC
//	fprintf( stderr, "\nmatch_calc... %d", i1 );
	int j, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *matchpt, *cpmxpdpt, **cpmxpdptpt;
	int *cpmxpdnpt, **cpmxpdnptpt;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[j][count] = cpmx2[l][j];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	{
		for( l=0; l<nalphabets; l++ )
		{
			scarr[l] = 0.0;
			for( j=0; j<nalphabets; j++ )
//				scarr[l] += n_dis[j][l] * cpmx1[j][i1];
//				scarr[l] += n_dis_consweight_multi[j][l] * cpmx1[j][i1];
				scarr[l] += scoreingmtx[j][l] * cpmx1[j][i1];
		}
		matchpt = match;
		cpmxpdnptpt = cpmxpdn;
		cpmxpdptpt = cpmxpd;
		while( lgth2-- )
		{
//			*matchpt = 0.0;
			cpmxpdnpt = *cpmxpdnptpt++;
			cpmxpdpt = *cpmxpdptpt++;
			while( *cpmxpdnpt>-1 )
				*matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
			matchpt++;
		} 
	}
	free( scarr );
//	fprintf( stderr, "done\n" );
#else
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
// simple
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}
	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
//			scarr[l] += n_dis[k][l] * cpmx1[k][i1];
//			scarr[l] += n_dis_consweight_multi[k][l] * cpmx1[k][i1];
			scarr[l] += scoreingmtx[k][l] * cpmx1[k][i1];
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[k][j]>-1; k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
	free( scarr );
#endif
}

static void match_calc( double **n_dynamicmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
#if FASTMATCHCALC
//	fprintf( stderr, "\nmatch_calc... %d", i1 );
	int j, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *matchpt, *cpmxpdpt, **cpmxpdptpt;
	int *cpmxpdnpt, **cpmxpdnptpt;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[j][count] = cpmx2[l][j];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	{
		for( l=0; l<nalphabets; l++ )
		{
			scarr[l] = 0.0;
			for( j=0; j<nalphabets; j++ )
//				scarr[l] += n_dis[j][l] * cpmx1[j][i1];
//				scarr[l] += n_dis_consweight_multi[j][l] * cpmx1[j][i1];
				scarr[l] += n_dynamicmtx[j][l] * cpmx1[j][i1];
		}
		matchpt = match;
		cpmxpdnptpt = cpmxpdn;
		cpmxpdptpt = cpmxpd;
		while( lgth2-- )
		{
			*matchpt = 0.0;
			cpmxpdnpt = *cpmxpdnptpt++;
			cpmxpdpt = *cpmxpdptpt++;
			while( *cpmxpdnpt>-1 )
				*matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
			matchpt++;
		} 
	}
	free( scarr );
//	fprintf( stderr, "done\n" );
#else
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
// simple
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}
	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
//			scarr[l] += n_dis[k][l] * cpmx1[k][i1];
//			scarr[l] += n_dis_consweight_multi[k][l] * cpmx1[k][i1];
			scarr[l] += n_dynamicmtx[k][l] * cpmx1[k][i1];
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[k][j]>-1; k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
	free( scarr );
#endif
}

static void Atracking_localhom( double *impwmpt, double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int *warpis, int *warpjs, int warpbase )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
	double wm;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	gt1bk = AllocateCharVec( lgth1+lgth2+1 );
	gt2bk = AllocateCharVec( lgth1+lgth2+1 );

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( outgap == 1 )
		;
	else
	{
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
	}

    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	gaptable1 = gt1bk + lgth1+lgth2;
	*gaptable1 = 0;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable2 = 0;

	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	*impwmpt = 0.0;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] >= warpbase )
		{
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
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
			}
			l= jin;
			while( --l >= 0 )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l ) 
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
			}
			l= jin - jfi;
			while( --l )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
			}
		}
		if( iin == lgth1 || jin == lgth2 )
			;
		else
		{
			*impwmpt += (double)imp_match_out_scD( iin, jin );

//		fprintf( stderr, "impwm = %f (iin=%d, jin=%d) seq1=%c, seq2=%c\n", *impwmpt, iin, jin, seq1[0][iin], seq2[0][jin] );
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;
	}

	for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );
	for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );

	free( gt1bk );
	free( gt2bk );
}

static double Atracking( double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int tailgp,
						int *warpis, int *warpjs, int warpbase )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
	double wm;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	gt1bk = AllocateCharVec( lgth1+lgth2+1 );
	gt2bk = AllocateCharVec( lgth1+lgth2+1 );

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( tailgp == 1 )
		;
	else
	{
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
	}

    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	gaptable1 = gt1bk + lgth1+lgth2;
	*gaptable1 = 0;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable2 = 0;

	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] >= warpbase )
		{
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
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
			}
			l= jin;
			while( --l >= 0 )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l ) 
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
			}
			l= jin - jfi;
			while( --l )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
			}
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;
	}

	for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );
	for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );

	free( gt1bk );
	free( gt2bk );

	return( 0.0 );
}

#if 1
static void cleargaplens( Gaplen ****gaplens )
{
	int i;
//	int j, k;
	Gaplen ****ptptptpt, ***ptptpt, ***ptptptmv, **ptptmv, **ptpt;
	for( i=0; i<12; i++ )
	{
		ptptptpt = gaplens+i;

		ptptptmv = *ptptptpt;
		while( *(ptptpt=ptptptmv++) )
//		for( j=0; *(ptptpt=*ptptptpt+j); j++ )
		{
			ptptmv = *ptptpt;
			while( *(ptpt=ptptmv++)!=(Gaplen *)1 )
//			for( k=0; gaplens[i][j][k]!=(Gaplen *)1; k++ )
//			for( k=0; *(ptpt=*ptptpt+k)!=(Gaplen *)1; k++ )
			{
//				if( gaplens[i][j][k] ) free( gaplens[i][j][k] );
				if( *ptpt ) free( *ptpt );
				*ptpt = NULL;
			}
//			free( gaplens[i][j] );
			free( *ptptpt );
//			gaplens[i][j] = NULL;
			*ptptpt = NULL;
		}
	}
}
#else
static void cleargaplens( Gaplen ****gaplens )
{
	int i, j, k;
	for( i=0; i<12; i++ )
	{
		for( j=0; gaplens[i][j]; j++ )
		{
			for( k=0; gaplens[i][j][k]!=(Gaplen *)1; k++ )
			{
				if( gaplens[i][j][k] ) free( gaplens[i][j][k] );
				gaplens[i][j][k] = NULL;
			}
			free( gaplens[i][j] );
			gaplens[i][j] = NULL;
		}
	}
}
#endif

#if USEGAPLENHALFORMTX
static void FreeGaplenMtxReport( Gaplen **mtx )
{
	int i;
	if( mtx == NULL ) return;

	for( i=0; ; i++ )
	{
		reporterr( "i=%d, mtx[i] = %p\n", i, mtx[i] );
		if( mtx[i] )
		{
			if( mtx[i] == (Gaplen *)1 ) break;
			free( mtx[i] ); mtx[i] = NULL;
		}
		mtx[i] = NULL;
	}
	free( mtx );
	mtx = NULL;
}
#endif

static void FreeGaplenMtx( Gaplen **mtx, int inclfreq )
{
	int i;
	if( mtx == NULL ) return;

	for( i=0; ; i++ )
	{
		if( mtx[i] )
		{
			if( mtx[i] == (Gaplen *)1 ) break;

#if 0
			if( inclfreq ) 
			{
//				reporterr( "inclfreq=%d\n", inclfreq );
				for( j=0; mtx[i][j].relend==0; j++ ) 
				{
//					reporterr( "j=%d\n", j );
//					reporterr( "Free! freq\n" );
					if( mtx[i][j].freq ) 
					{
						free( mtx[i][j].freq );
					}
					mtx[i][j].freq = NULL;
				}
			}
#endif

			free( mtx[i] ); mtx[i] = NULL;
		}
	}
	free( mtx );
	mtx = NULL;
}

#if USEGAPLENHALFORMTX
static void FreeGaplenCubgaplenReport( Gaplen ***cub )
{
	int i;
	if( cub == NULL ) return;

	for( i=0; cub[i]; i++ )
	{
		reporterr( "i=%d, cub[i]=%p\n", i, cub[i] );
		FreeGaplenMtx( cub[i], 0 );
		cub[i] = NULL;
	}
	free( cub );
	cub = NULL;
}
#endif

static void FreeGaplenCub( Gaplen ***cub )
{
	int i;
	if( cub == NULL ) return;

	for( i=0; cub[i]; i++ )
	{
		FreeGaplenMtx( cub[i], 0 );
		cub[i] = NULL;
	}
	free( cub );
	cub = NULL;
}

static int strralpha( const char *s, const char *first )
{
	int v = 0;
	s--;
	while( s >= first )
	{
		if( *s-- != '-' ) return( v );
		v++;
	}
	if( s == first-1 ) return( v );
	return( -1 );
}

static void fillgaplen( Gaplen **mtx, int l )
{
	int i, j, n, k, len, pos, idatnext;
	double freq;
	for( i=0; i<=l; i++ )
	{
//		reporterr( "i=%d\n", i );
		if( mtx[i] == NULL ) continue;
		for( n=0; (len=mtx[i][n].len)!=-1; n++ )
		{
			freq = mtx[i][n].freq;
			idatnext = n;
			for( j=0; j<len; j++ )
			{
//				reporterr( "n=%d, j=%d, i=%d, len=%d\n", n, j, i, len );
				pos = i-1-j;
//				reporterr( "pos = %d\n", pos );
				if( mtx[pos] == NULL ) 
				{
					mtx[pos] = calloc( 2, sizeof( Gaplen ) );
					mtx[pos][0].len = -1;
#if USEGAPLENHALFORMTX
					mtx[pos][0].idatend = -1;
#endif
					mtx[pos][0].idatnext = -1;
					mtx[pos][0].idatprev = -1;
					mtx[pos][0].relend = -1;
					mtx[pos][0].freq = 0.0;
					mtx[pos][0].npat = -1;
					k = 0;
				}
				else
				{
					k = mtx[pos][0].npat;
				}
				mtx[pos] = realloc( mtx[pos], sizeof( Gaplen ) * ( k + 2 ) );
//				mtx[pos][k].len = -100000; // tsukawanai!
				mtx[pos][k].len = len; // compact no toki tsukau!
				mtx[pos][k].relend = j+1;
				mtx[pos][k].freq = freq;// tsukawanai! mtx[i][n].freq;
				mtx[pos][k].idatnext = idatnext;
				mtx[pos][k].idatprev =-1;
				mtx[pos][k].npat = -1;
				mtx[pos][k+1].len = -1;
				mtx[pos][k+1].idatnext = -1;
				mtx[pos][k+1].relend = -1;
				mtx[pos][k+1].freq = 0.0;
				mtx[pos][k+1].npat = -1;
				mtx[pos][0].npat = k+1;
#if USEGAPLENHALFORMTX
				mtx[pos][k].idatend = n;
				mtx[pos][k+1].idatend = -1;
#endif

				mtx[pos+1][idatnext].idatprev = k; // kanarazu aru?

				idatnext = k;
			}
		}
	}
}

static int gapvariety( int n, int l, char **seq )
{
	int i, j, gl, *known, nknown, val;
	known = calloc( l+1, sizeof( int ) );
//	for( i=0; i<n; i++ ) reporterr( "seq[%d] = %s\n", i, seq[i] );

	val = 0;
	for( j=0; j<=l; j++ )
	{
		for( i=0; i<j; i++ ) known[i] = 0;
		nknown = 0;
		for( i=0; i<n; i++ )
		{
			if( seq[i][j] == '-' ) continue;

			gl = strralpha( seq[i]+j, seq[i] );
//			reporterr( "gl = %d\n", gl );
			if( gl > 0 )
			{
				if( known[gl] )
				{
					;
				}
				else
				{
					nknown++;
				}
			}
		}
		val += nknown;
	}
	free( known );

	return( val );
}


static void gaplencount( int n, int l, Gaplen **mtx, char **seq, double *eff )
{
	int i, j, k, gl, *known, nknown;
	known = calloc( l+1, sizeof( int ) );
//	for( i=0; i<n; i++ ) reporterr( "seq[%d] = %s\n", i, seq[i] );

	for( j=0; j<=l; j++ )
	{
		if( mtx[j] ) 
		{
			if( mtx[j] == (Gaplen *)1 ) break;
			for( k=0; mtx[j][k].relend==0; k++ ) 
			{
#if 0
//				reporterr( "j=%d\n", j );
//				reporterr( "Free! freq\n" );
				if( mtx[j][k].freq ) 
				{
					free( mtx[j][k].freq );
				}
				mtx[j][k].freq = NULL;
#endif
			}
			free( mtx[j] );
		}
		mtx[j] = NULL;
	}

	for( j=0; j<=l; j++ )
	{
		for( i=0; i<j; i++ ) known[i] = 0;
		nknown = 0;
		for( i=0; i<n; i++ )
		{
			if( seq[i][j] == '-' ) continue;

			gl = strralpha( seq[i]+j, seq[i] );
//			reporterr( "gl = %d\n", gl );
			if( gl > 0 )
			{
				if( known[gl] )
				{
//					reporterr( "gl=%d, Known!\n", gl );
					for( k=0; mtx[j][k].len!=-1; k++ ) if( mtx[j][k].len == gl ) break;
					if( mtx[j][k].len == -1 )
					{
						reporterr( "Unexpected error!\n" );
						exit( 1 );
					}
					mtx[j][k].freq += eff[i];
				}
				else
				{
//					reporterr( "gl=%d, First!\n", gl );
					mtx[j] = realloc( mtx[j], sizeof( Gaplen ) * (nknown+2) );
					mtx[j][nknown].len = gl;
					mtx[j][nknown].relend = 0;
					mtx[j][nknown].freq = eff[i];
					mtx[j][nknown].idatnext = -2;
					mtx[j][nknown+1].len = -1;
					mtx[j][nknown+1].idatnext = -1;
					mtx[j][nknown+1].relend = -1;
					mtx[j][nknown+1].freq = 0.0;
					mtx[j][nknown+1].npat = -1;
#if USEGAPLENHALFORMTX
					mtx[j][nknown].idatend = nknown;
					mtx[j][nknown+1].idatend = -1;
#endif
					known[gl]++;
					nknown++;
					mtx[j][0].npat = nknown;
				}
			}
		}
	}
	fillgaplen( mtx, l );
#if 0
	reporterr( "Gaplen:\n" );
	for( i=0; i<=l; i++ )
	{
//		reporterr( "i=%d, gaplen[i] = %p\n", i, mtx[i] );
		if( mtx[i] ) 
		{
			for( j=0; mtx[i][j].len!=-1; j++ )
				reporterr( "i=%d, len = %d, relend = %d, freq = %f\n", i, mtx[i][j].len, mtx[i][j].relend, mtx[i][j].freq );
		}
	}

#endif

	free( known );
}


#if DEBUG
static void showgaplen( Gaplen **mtx, int seqlen )
{
	int i, l;
#if USEGAPLENHALFORMTX
	int id, pos;
#endif
//	for( i=0; i<=seqlen; i++ )
	for( i=0; ; i++ )
	{
//		reporterr( "chain[%d] = %d\n", i, chain[i] );
		if( mtx[i] == NULL ) continue;
		if( mtx[i] == (Gaplen *)1 ) break;
		for( l=0; mtx[i][l].idatnext!=-1; l++ )
		{
#if USEGAPLENHALFORMTX
			reporterr( "i=%d, l=%d, len=%d, relend=%d, idatend=%d, idatnext=%d, idatprev=%d, freq=%f\n", i, l, mtx[i][l].len, mtx[i][l].relend, mtx[i][l].idatend, mtx[i][l].idatnext, mtx[i][l].idatprev, mtx[i][l].freq );
			pos = mtx[i][l].relend;
			id = mtx[i][l].idatend;
			if( mtx[i+pos] == NULL )
			{
//				reporterr( "Error in SOURCE\n" );
				reporterr( ".len and .freq were lost when i=%d!\n", i );
//				exit( 1 );
			}
#else
			reporterr( "i=%d, l=%d, len=%d, relend=%d, idatnext=%d, idatprev=%d, freq=%f\n", i, l, mtx[i][l].len, mtx[i][l].relend, mtx[i][l].idatnext, mtx[i][l].idatprev, mtx[i][l].freq );
#endif
		}
	}
}
#endif

#if WMCHECK
static int pairgapcount( char *s1, char *s2 )
{
	char **tmpseq;
	int i, len, st, k;
	int v = 0;
	
	len = strlen( s1 );
	tmpseq = calloc( sizeof( char * ), 2 );
	tmpseq[0] = malloc( sizeof( char ) * ( len + 1 ) );
	tmpseq[1] = malloc( sizeof( char ) * ( len + 1 ) );

	strcpy( tmpseq[0], s1 );
	strcpy( tmpseq[1], s2 );

	commongappick( 2, tmpseq );
	len = strlen( tmpseq[0] );


	for( k=0; k<2; k++ )
	{
		st = 0;
		for( i=0; i<len; i++ )
		{
			if( tmpseq[k][i] == '-' )
			{
				if( st == 0 )
				{
					v++;
					st = 1;
				}
			}
			else
			{
				st = 0;
			}
		}
	}
	free( tmpseq[0] );
	free( tmpseq[1] );
	free( tmpseq );

	return( v );
}
#endif

static double calcpfac_gap_noidatend( Gaplen **gaplen1, Gaplen **gaplen2, int newgaplen, int i, int j, char *seq1, char *seq2, int disp ) // seq1 to seq2 ha debug you
{
#if 1
	double pfac, pfac1, pfac10, pfac2;
	int k, l, pos1, pos2;
	Gaplen *gaplen1i, *gaplen2j, *g1, *g2;
	
#if 0 // .len no shouryaku ni taiou shiteinai
	int gl; 
	if( disp )
	{
		reporterr( "calcpfac_gap_noidatend, %c (%d) - %c (%d)\n", seq1[i], i, seq2[j], j );
		reporterr( "newgaplen = %d\n", newgaplen );
		reporterr( "In calcpfac_gap, gaplen1[%d(%c)] = \n", i, seq1[i] );
		for( k=0; gaplen1[i]&&(id1=gaplen1[i][k].idatnext)!=-1; k++ )
		{
			pos1 = gaplen1[i][k].relend;
			reporterr( ".len=%d, .relend=%d, .freq=%f\n", gaplen1[i+pos1][id1].len, gaplen1[i][k].relend, gaplen1[i+pos1][id1].freq[0] );
		}
		reporterr( "In calcpfac_gap, gaplen2[%d(%c)] = \n", j, seq2[j] );
		showgaplen( gaplen2, strlen(seq2) );
		for( k=0; gaplen2[j]&&(id2=gaplen2[j][k].idatnext)!=-1; k++ )
		{
			pos2 = gaplen2[j][k].relend;
			reporterr( ".len=%d, .relend=%d, .freq=%f\n", gaplen2[j+pos2][id2].len, gaplen2[j][k].relend, gaplen2[j+pos2][id2].freq[0] );
		}
	}
#endif
	gaplen2j = gaplen2[j];
	gaplen1i = gaplen1[i];

	pfac = 0.0;
	pfac1 = 0.0;
	pfac10 = 0.0;
	if( gaplen1i ) for( k=0; (g1=gaplen1i+k)->idatnext!=-1; k++ )
	{
		pos1 = g1->relend;
		if( pos1 != 0 )
		{
			pfac2 = 0.0;
			if( gaplen2j ) for( l=0; (g2=gaplen2j+l)->idatnext!=-1; l++ ) 
			{
				pos2 = g2->relend;
				if( pos2 == 0 && g2->len >g1->len - (pos1) + newgaplen ) 
				{
					pfac2 += g2->freq;
//					reporterr( "hit! pfac2=%f, .freq=%f\n", pfac2, gaplen2[j][l].freq );
				}
//				else
//					reporterr( "does not hit! pfac2=%f, gaplen1[i][k].len=%d, gaplen[i][k].relend=%d, newgaplen=%d\n", pfac2, gaplen1[i][k].len, gaplen1[i][k].relend, newgaplen );
			}
			pfac += pfac2 * g1->freq;
			pfac1 += g1->freq;
		}
		else if( pos1 == 0 )
		{
			pfac2 = 1.0;
			if( gaplen2j ) for( l=0; (g2=gaplen2j+l)->idatnext!=-1; l++ ) 
			{
				pos2 = g2->relend;
				if( pos2 == 0 && g2->len == g1->len+newgaplen ) pfac2 -= g2->freq;// kokode shuryou suru gap, gaplen1 ha kangaenai.
				if( pos2 != 0 && g2->len - (pos2-1) > g1->len+newgaplen ) pfac2 -= g2->freq;// keizoku suru gap, gaplen1 ha kangaenai.
			}
//			reporterr( "pfac2 in line 1056 = %f\n", pfac2 );
			pfac += pfac2 * g1->freq;
			pfac10 += g1->freq;
		}
	}
#if DEBUG
	reporterr( "pfac1 (step2) = %f\n", pfac1 );
	reporterr( "pfac10 (step2) = %f\n", pfac10 );
	reporterr( "pfac (step2) = %f\n", pfac );
#endif

	pfac1 = 1.0 - pfac1 - pfac10;
	pfac2 = 1.0;
	if( gaplen2j ) for( l=0; (g2=gaplen2j+l)->idatnext!=-1; l++ ) 
	{
		pos2 = g2->relend;
		if( pos2 == 0 && g2->len == newgaplen ) pfac2 -= g2->freq;// kokode shuryou suru gap, gaplen1 ha kangaenai.
		if( pos2 != 0 && g2->len - (pos2-1) > newgaplen ) pfac2 -= g2->freq;// keizoku suru gap, gaplen1 ha kangaenai.
	}
#if DEBUG
	reporterr( "pfac1 (type3) = %f\n", pfac1 );
	reporterr( "pfac2 (type3) = %f\n", pfac2 );
	reporterr( "pfac (step3) = %f\n", pfac );
#endif
	pfac += pfac1 * pfac2;
#if DEBUG
	reporterr( "incomplete pfac = %f, pfac1,pfac2 (%c%d,%c%d) = %f, %f\n", pfac, seq1[i], i, seq2[j], j, pfac1, pfac2 );
#endif


	return( pfac );

#else

	double pfac, pfac1, pfac10, pfac2;
	int k, l, pos1, pos2, id1, id2;
	Gaplen *gaplen1i, *gaplen2j;
	
#if 0 // .len no shouryaku ni taiou shiteinai
	int gl; 
	if( disp )
	{
		reporterr( "calcpfac_gap_noidatend, %c (%d) - %c (%d)\n", seq1[i], i, seq2[j], j );
		reporterr( "newgaplen = %d\n", newgaplen );
		reporterr( "In calcpfac_gap, gaplen1[%d(%c)] = \n", i, seq1[i] );
		for( k=0; gaplen1[i]&&(id1=gaplen1[i][k].idatend)!=-1; k++ )
		{
			pos1 = gaplen1[i][k].relend;
			reporterr( ".len=%d, .relend=%d, .freq=%f\n", gaplen1[i+pos1][id1].len, gaplen1[i][k].relend, gaplen1[i+pos1][id1].freq[0] );
		}
		reporterr( "In calcpfac_gap, gaplen2[%d(%c)] = \n", j, seq2[j] );
		showgaplen( gaplen2, strlen(seq2) );
		for( k=0; gaplen2[j]&&(id2=gaplen2[j][k].idatend)!=-1; k++ )
		{
			pos2 = gaplen2[j][k].relend;
			reporterr( ".len=%d, .relend=%d, .freq=%f\n", gaplen2[j+pos2][id2].len, gaplen2[j][k].relend, gaplen2[j+pos2][id2].freq[0] );
		}
	}
#endif
	gaplen2j = gaplen2[j];
	gaplen1i = gaplen1[i];

	pfac = 0.0;
	pfac1 = 0.0;
	pfac10 = 0.0;
	if( gaplen1i ) for( k=0; (gaplen1i[k].idatnext)!=-1; k++ )
	{
		pos1 = gaplen1i[k].relend;
		if( pos1 != 0 )
		{
			pfac2 = 0.0;
			if( gaplen2j ) for( l=0; (gaplen2j[l].idatnext)!=-1; l++ ) 
			{
				pos2 = gaplen2j[l].relend;
				if( pos2 == 0 && gaplen2j[l].len > gaplen1i[k].len - (pos1) + newgaplen ) 
				{
					pfac2 += gaplen2j[l].freq;
//					reporterr( "hit! pfac2=%f, .freq=%f\n", pfac2, gaplen2[j][l].freq );
				}
//				else
//					reporterr( "does not hit! pfac2=%f, gaplen1[i][k].len=%d, gaplen[i][k].relend=%d, newgaplen=%d\n", pfac2, gaplen1[i][k].len, gaplen1[i][k].relend, newgaplen );
			}
			pfac += pfac2 * gaplen1i[k].freq;
			pfac1 += gaplen1i[k].freq;
		}
		else if( pos1 == 0 )
		{
			pfac2 = 1.0;
			if( gaplen2j ) for( l=0; (gaplen2j[l].idatnext)!=-1; l++ ) 
			{
				pos2 = gaplen2j[l].relend;
				if( pos2 == 0 && gaplen2j[l].len == gaplen1i[k].len+newgaplen ) pfac2 -= gaplen2j[l].freq;// kokode shuryou suru gap, gaplen1 ha kangaenai.
				if( pos2 != 0 && gaplen2j[l].len - (pos2-1) > gaplen1i[k].len+newgaplen ) pfac2 -= gaplen2j[l].freq;// keizoku suru gap, gaplen1 ha kangaenai.
			}
//			reporterr( "pfac2 in line 1056 = %f\n", pfac2 );
			pfac += pfac2 * gaplen1i[k].freq;
			pfac10 += gaplen1i[k].freq;
		}
	}
#if DEBUG
	reporterr( "pfac1 (step2) = %f\n", pfac1 );
	reporterr( "pfac10 (step2) = %f\n", pfac10 );
	reporterr( "pfac (step2) = %f\n", pfac );
#endif

	pfac1 = 1.0 - pfac1 - pfac10;
	pfac2 = 1.0;
	if( gaplen2j ) for( l=0; (gaplen2j[l].idatnext)!=-1; l++ ) 
	{
		pos2 = gaplen2j[l].relend;
		if( pos2 == 0 && gaplen2j[l].len == newgaplen ) pfac2 -= gaplen2j[l].freq;// kokode shuryou suru gap, gaplen1 ha kangaenai.
		if( pos2 != 0 && gaplen2j[l].len - (pos2-1) > newgaplen ) pfac2 -= gaplen2j[l].freq;// keizoku suru gap, gaplen1 ha kangaenai.
	}
#if DEBUG
	reporterr( "pfac1 (type3) = %f\n", pfac1 );
	reporterr( "pfac2 (type3) = %f\n", pfac2 );
	reporterr( "pfac (step3) = %f\n", pfac );
#endif
	pfac += pfac1 * pfac2;
#if DEBUG
	reporterr( "incomplete pfac = %f, pfac1,pfac2 (%c%d,%c%d) = %f, %f\n", pfac, seq1[i], i, seq2[j], j, pfac1, pfac2 );
#endif


	return( pfac );

#endif
}

#if USEGAPLENHALFORMTX

static double calcpfac_gap_incomplete( Gaplen **gaplen1, Gaplen **gaplen2, int newgaplen, int i, int j, char *seq1, char *seq2, int disp ) // seq1 to seq2 ha debug you
{
	double pfac, pfac1, pfac10, pfac2;
	int k, l, pos1, pos2, id1, id2;
	Gaplen *gapend1, *gapend2;
	Gaplen *gaplen1i, *gaplen2j;
	
#if 0 // .len no shouryaku ni taiou shiteinai
	int gl; 
	if( disp )
	{
		reporterr( "calcpfac_gap_incomplete, %c (%d) - %c (%d)\n", seq1[i], i, seq2[j], j );
		reporterr( "newgaplen = %d\n", newgaplen );
		reporterr( "In calcpfac_gap, gaplen1[%d(%c)] = \n", i, seq1[i] );
		for( k=0; gaplen1[i]&&(id1=gaplen1[i][k].idatend)!=-1; k++ )
		{
			pos1 = gaplen1[i][k].relend;
			reporterr( ".len=%d, .relend=%d, .freq=%f\n", gaplen1[i+pos1][id1].len, gaplen1[i][k].relend, gaplen1[i+pos1][id1].freq[0] );
		}
		reporterr( "In calcpfac_gap, gaplen2[%d(%c)] = \n", j, seq2[j] );
		showgaplen( gaplen2, strlen(seq2) );
		for( k=0; gaplen2[j]&&(id2=gaplen2[j][k].idatend)!=-1; k++ )
		{
			pos2 = gaplen2[j][k].relend;
			reporterr( ".len=%d, .relend=%d, .freq=%f\n", gaplen2[j+pos2][id2].len, gaplen2[j][k].relend, gaplen2[j+pos2][id2].freq[0] );
		}
	}
#endif
	gaplen2j = gaplen2[j];
	gaplen1i = gaplen1[i];

	pfac = 0.0;
	pfac1 = 0.0;
	pfac10 = 0.0;
	if( gaplen1i ) for( k=0; (id1=gaplen1i[k].idatend)!=-1; k++ )
	{
		pos1 = gaplen1i[k].relend;
		gapend1 = gaplen1[i+pos1]+id1;
		if( pos1 != 0 )
		{
			pfac2 = 0.0;
			if( gaplen2j ) for( l=0; (id2=gaplen2j[l].idatend)!=-1; l++ ) 
			{
				pos2 = gaplen2j[l].relend;
				gapend2 = gaplen2[j+pos2]+id2;
//				if( pos2 == 0 && gapend2->len + 1 > gapend1->len - (pos1-1) + newgaplen ) 
				if( pos2 == 0 && gapend2->len > gapend1->len - (pos1) + newgaplen ) 
				{
					pfac2 += gapend2->freq;
//					reporterr( "hit! pfac2=%f, .freq=%f\n", pfac2, gaplen2[j][l].freq );
				}
//				else
//					reporterr( "does not hit! pfac2=%f, gaplen1[i][k].len=%d, gaplen[i][k].relend=%d, newgaplen=%d\n", pfac2, gaplen1[i][k].len, gaplen1[i][k].relend, newgaplen );
			}
			pfac += pfac2 * gapend1->freq;
			pfac1 += gapend1->freq;
		}
		else if( pos1 == 0 )
		{
			pfac2 = 1.0;
			if( gaplen2j ) for( l=0; (id2=gaplen2j[l].idatend)!=-1; l++ ) 
			{
				pos2 = gaplen2j[l].relend;
				gapend2 = gaplen2[j+pos2]+id2;
				if( pos2 == 0 && gapend2->len == gapend1->len+newgaplen ) pfac2 -= gapend2->freq;// kokode shuryou suru gap, gaplen1 ha kangaenai.
				if( pos2 != 0 && gapend2->len - (pos2-1) > gapend1->len+newgaplen ) pfac2 -= gapend2->freq;// keizoku suru gap, gaplen1 ha kangaenai.
			}
//			reporterr( "pfac2 in line 1056 = %f\n", pfac2 );
			pfac += pfac2 * gapend1->freq;
			pfac10 += gapend1->freq;
		}
	}
#if DEBUG
	reporterr( "pfac1 (step2) = %f\n", pfac1 );
	reporterr( "pfac10 (step2) = %f\n", pfac10 );
	reporterr( "pfac (step2) = %f\n", pfac );
#endif

	pfac1 = 1.0 - pfac1 - pfac10;
	pfac2 = 1.0;
	if( gaplen2j ) for( l=0; (id2=gaplen2j[l].idatend)!=-1; l++ ) 
	{
		pos2 = gaplen2j[l].relend;
		gapend2 = gaplen2[j+pos2]+id2;
		if( pos2 == 0 && gapend2->len == newgaplen ) pfac2 -= gapend2->freq;// kokode shuryou suru gap, gaplen1 ha kangaenai.
		if( pos2 != 0 && gapend2->len - (pos2-1) > newgaplen ) pfac2 -= gapend2->freq;// keizoku suru gap, gaplen1 ha kangaenai.
	}
#if DEBUG
	reporterr( "pfac1 (type3) = %f\n", pfac1 );
	reporterr( "pfac2 (type3) = %f\n", pfac2 );
	reporterr( "pfac (step3) = %f\n", pfac );
#endif
	pfac += pfac1 * pfac2;
#if DEBUG
	reporterr( "incomplete pfac = %f, pfac1,pfac2 (%c%d,%c%d) = %f, %f\n", pfac, seq1[i], i, seq2[j], j, pfac1, pfac2 );
#endif


	return( pfac );

}

static double calcpfac_gapex( Gaplen **gaplen1, Gaplen **gaplen2, int i, int j, int newgaplen, char *seq1, char *seq2, int disp )
{
	double pfac, pfac1, pfac2, pfac10;
	int k, l, id1, id2, pos1, pos2;
	Gaplen *gapend1, *gapend2;
	Gaplen *gaplen1i, *gaplen2j;

	gaplen1i = gaplen1[i];
	gaplen2j = gaplen2[j];

	pfac = 0.0;
	pfac2 = 0.0;
//	for( k=0; gaplen2[j]&&(gl=gaplen2[j][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
	if( gaplen2j ) for( k=0; (id2=gaplen2j[k].idatend)!=-1; k++ ) // ososugi!  hash ni atode henkou
	{
#if DEBUG
		int gl;
		pos2 = gaplen2j[k].relend;
		id2 = gaplen2j[k].idatend;
		gl = gaplen2[j+pos2][id2].len;
		if( disp ) reporterr( "gaplen2[][].len=%d, .relend=%d, .freq=%f\n", gaplen2[j+pos2][id2].len, gaplen2[j][k].relend, gaplen2[j+pos2][id2].freq );
		if( disp ) reporterr( "gl = %d, newgaplen=%d\n", gl, newgaplen );
#endif
		if( (pos2=gaplen2[j][k].relend) != 0 ) continue;

		gapend2 = gaplen2[j+pos2]+id2;
		pfac1 = 1.0;
		pfac10 = 1.0;
		if( gaplen1i ) for( l=0; (id1=gaplen1i[l].idatend)!=-1; l++ ) // ososugi!  hash ni atode henkou
		{
			pos1 = gaplen1i[l].relend;
			gapend1 = gaplen1[i+pos1]+id1;
			pfac10 -= gapend1->freq;
#if DEBUG
			if( disp ) reporterr( "gaplen1[][].len=%d, .relend=%d, .freq=%f\n", gaplen1[i+pos1][id1].len, gaplen1[i][l].relend, gaplen1[i+pos1][id1].freq );
#endif
			if( newgaplen + gapend1->len - (pos1) > gapend2->len - (pos2) ) pfac1 -= gapend1->freq;
//			reporterr( "pfac1 = %f\n", pfac1 );
		}
		pfac += pfac1 * gapend2->freq;


/* ???? */
		if( newgaplen >= gapend2->len - (pos2-1) )  // >= or >??
		{
			pfac -= pfac10 * gapend2->freq;
//			reporterr( "Hit! pfac1 = %f\n", pfac1 );
		}
/* ???? */


//		if( gaplen2[j][k].relend == -1 ) pfac += gaplen2[j][k].freq;
	}

	return( pfac );
}

static double calcpfac( Gaplen **gaplen1, Gaplen **gaplen2, int i, int j, char *seq1, char *seq2, int disp ) // seq1 to seq2 ha debug you
{
	double pfac, pfac1, pfac2;
	int k, l, pos1, pos2, id1, id2;
	Gaplen *gapend1, *gapend2;
	Gaplen *gaplen1i, *gaplen2j;

	gaplen1i = gaplen1[i];
	gaplen2j = gaplen2[j];

#if DEBUG
	if( disp )
	{
		reporterr( "seq1[0] = %s\n", seq1 );
		reporterr( "seq2[0] = %s\n", seq2 );
		reporterr( "i,j=%d,%d\n", i, j );
	
		reporterr( "In calcpfac(), gaplen1[%d(%c)] = \n", i, seq1[i] );
//		showgaplen( gaplen1, seqlen( seq1 ) );
		for( k=0; gaplen1[i]&&(id1=gaplen1[i][k].idatend)!=-1; k++ )
		{
			pos1 = gaplen1[i][k].relend;
			reporterr( "pos1=%d, id1=%d\n", pos1, id1 );
			reporterr( ".len=%d, .relend=%d, .freq=%f\n", gaplen1[i+pos1][id1].len, gaplen1[i][k].relend, gaplen1[i+pos1][id1].freq );
		}
	
		reporterr( "In calcpfac(), gaplen2[%d(%c)] = \n", j, seq2[j] );
//		showgaplen( gaplen2, seqlen( seq2 ) );
		for( k=0; gaplen2[j]&&(id2=gaplen2[j][k].idatend)!=-1; k++ )
		{
			pos2 = gaplen2[j][k].relend;
			reporterr( "j=%d, k=%d, id2=%d, pos2=%d\n", j, k, id2, pos2 );
			reporterr( ".len=%d, .relend=%d\n", gaplen2[j+pos2][id2].len, gaplen2[j][k].relend );
			reporterr( ".freq=%f\n", gaplen2[j+pos2][id2].freq );
		}
	}
#endif 

	pfac1 = pfac2 = 0.0;
	if( gaplen1i ) for( k=0; (id1=gaplen1i[k].idatend)!=-1; k++ )
	{
		if( (pos1=gaplen1i[k].relend) == 0 ) pfac1 += gaplen1[i+pos1][id1].freq;
	}

	if( gaplen2j ) for( l=0; (id2=gaplen2j[l].idatend)!=-1; l++ ) // ososugi!  hash ni atode henkou
	{
		if( (pos2=gaplen2j[l].relend) == 0 ) pfac2 += gaplen2[j+pos2][id2].freq;
	}
#if DEBUG
	reporterr( "\n\nInitial pfac1,pfac2 (%c%d,%c%d) = %f, %f\n", seq1[i], i, seq2[j], j, pfac1, pfac2 );
#endif
	pfac = pfac1 * pfac2 + pfac1 * (1-pfac2) + pfac2 * (1-pfac1);
#if DEBUG
	reporterr( "\n\nInitial pfac (%d,%d) = %f\n", i, j, pfac );
#endif

#if 1
//	if( pfac ) reporterr( "i,j=%d,%d, Cancel (eq len)? pfac = %f -> ", i, j, pfac );
	if( gaplen1i ) for( k=0; (id1=gaplen1i[k].idatend)!=-1; k++ ) // ososugi!  hash ni atode henkou
	{
		pos1=gaplen1i[k].relend;
		gapend1 = gaplen1[i+pos1]+id1;
		if( gaplen2j ) for( l=0; (id2=gaplen2j[l].idatend)!=-1; l++ ) // ososugi!  hash ni atode henkou
		{
			pos2 = gaplen2j[l].relend;
			gapend2 = gaplen2[j+pos2]+id2;
			if     ( pos1 == 0 && pos2 == 0 && gapend1->len == gapend2->len ) pfac -= gapend1->freq * gapend2->freq;
			else if( pos1 == 0 && pos2 != 0 && gapend2->len - (pos2-1) > gapend1->len ) pfac -= gapend1->freq * gapend2->freq;
			else if( pos1 != 0 && pos2 == 0 && gapend1->len - (pos1-1) > gapend2->len ) pfac -= gapend1->freq * gapend2->freq;
		}
	}

#if DEBUG
	reporterr( "\n\nFinal pfac1,pfac2 (%c%d,%c%d, straight) = %f\n\n", seq1[i], i, seq2[j], j, pfac );
#endif
#else
#endif
	return( pfac );
}
#endif

static double calcpfac_gapex_noidatend( Gaplen **gaplen1, Gaplen **gaplen2, int i, int j, int newgaplen, char *seq1, char *seq2, int disp )
{
#if 1
	double pfac, pfac1, pfac2, pfac10;
	int k, l,  pos1, pos2;
	Gaplen *gaplen1i, *gaplen2j, *g1, *g2;

	gaplen1i = gaplen1[i];
	gaplen2j = gaplen2[j];

	pfac = 0.0;
	pfac2 = 0.0;
	if( gaplen2j ) for( k=0; (g2=gaplen2j+k)->idatnext!=-1; k++ )
	{
#if DEBUG
		int gl;
		pos2 = gaplen2j[k].relend;
		gl = gaplen2j[k].len;
		if( disp ) reporterr( "gaplen2[][].len=%d, .relend=%d, .freq=%f\n", gaplen2[j][k].len, gaplen2[j][k].relend, gaplen2[j][k].freq );
		if( disp ) reporterr( "gl = %d, newgaplen=%d\n", gl, newgaplen );
#endif
		if( (pos2=g2->relend) != 0 ) continue;

		pfac1 = 1.0;
		pfac10 = 1.0;
		if( gaplen1i ) for( l=0; (g1=gaplen1i+l)->idatnext!=-1; l++ )
		{
			pos1 = g1->relend;
			pfac10 -= g1->freq;
#if DEBUG
			if( disp ) reporterr( "gaplen1[][].len=%d, .relend=%d, .freq=%f\n", gaplen1[i][l].len, gaplen1[i][l].relend, gaplen1[i][l].freq );
#endif
			if( newgaplen + g1->len - (pos1) > g2->len - (pos2) ) pfac1 -= g1->freq;
//			reporterr( "pfac1 = %f\n", pfac1 );
		}
		pfac += pfac1 * g2->freq;


/* ???? */
		if( newgaplen >= g2->len - (pos2-1) )  // >= or >??
		{
			pfac -= pfac10 * g2->freq;
//			reporterr( "Hit! pfac1 = %f\n", pfac1 );
		}
/* ???? */


//		if( gaplen2[j][k].relend == -1 ) pfac += gaplen2[j][k].freq;
	}

	return( pfac );
#else
	double pfac, pfac1, pfac2, pfac10;
	int k, l, id1, id2, pos1, pos2;
	Gaplen *gaplen1i, *gaplen2j;

	gaplen1i = gaplen1[i];
	gaplen2j = gaplen2[j];

	pfac = 0.0;
	pfac2 = 0.0;
//	for( k=0; gaplen2[j]&&(gl=gaplen2[j][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
	if( gaplen2j ) for( k=0; (gaplen2j[k].idatnext)!=-1; k++ ) // ososugi!  hash ni atode henkou
	{
#if DEBUG
		int gl;
		pos2 = gaplen2j[k].relend;
		gl = gaplen2j[k].len;
		if( disp ) reporterr( "gaplen2[][].len=%d, .relend=%d, .freq=%f\n", gaplen2[j][k].len, gaplen2[j][k].relend, gaplen2[j][k].freq );
		if( disp ) reporterr( "gl = %d, newgaplen=%d\n", gl, newgaplen );
#endif
		if( (pos2=gaplen2[j][k].relend) != 0 ) continue;

		pfac1 = 1.0;
		pfac10 = 1.0;
		if( gaplen1i ) for( l=0; (gaplen1i[l].idatnext)!=-1; l++ ) // ososugi!  hash ni atode henkou
		{
			pos1 = gaplen1i[l].relend;
			pfac10 -= gaplen1i[l].freq;
#if DEBUG
			if( disp ) reporterr( "gaplen1[][].len=%d, .relend=%d, .freq=%f\n", gaplen1[i][l].len, gaplen1[i][l].relend, gaplen1[i][l].freq );
#endif
			if( newgaplen + gaplen1i[l].len - (pos1) > gaplen2j[k].len - (pos2) ) pfac1 -= gaplen1i[l].freq;
//			reporterr( "pfac1 = %f\n", pfac1 );
		}
		pfac += pfac1 * gaplen2j[k].freq;


/* ???? */
		if( newgaplen >= gaplen2j[k].len - (pos2-1) )  // >= or >??
		{
			pfac -= pfac10 * gaplen2j[k].freq;
//			reporterr( "Hit! pfac1 = %f\n", pfac1 );
		}
/* ???? */


//		if( gaplen2[j][k].relend == -1 ) pfac += gaplen2[j][k].freq;
	}

	return( pfac );

#endif
}


static double calcpfacnoidatend( Gaplen **gaplen1, Gaplen **gaplen2, int i, int j, char *seq1, char *seq2, int disp ) // seq1 to seq2 ha debug you
{
	double pfac, pfac1, pfac2;
	int k, l, pos1, pos2;
	Gaplen *gaplen1i, *gaplen2j, *g1, *g2;

	gaplen1i = gaplen1[i];
	gaplen2j = gaplen2[j];

#if DEBUG
	if( disp )
	{
		reporterr( "seq1[0] = %s\n", seq1 );
		reporterr( "seq2[0] = %s\n", seq2 );
		reporterr( "i,j=%d,%d\n", i, j );
	
		reporterr( "In calcpfacnoidatend(), gaplen1[%d(%c)] = \n", i, seq1[i] );
		showgaplen( gaplen1, seqlen( seq1 ) );
		for( k=0; gaplen1[i]&&gaplen1[i][k].idatnext!=-1; k++ )
		{
			pos1 = gaplen1[i][k].relend;
			reporterr( ".len=%d, .relend=%d, .freq=%f (i=%d)\n", gaplen1[i][k].len, gaplen1[i][k].relend, gaplen1[i][k].freq, i );
		}
	
		reporterr( "In calcpfacnoidatend(), gaplen2[%d(%c)] = \n", j, seq2[j] );
		showgaplen( gaplen2, seqlen( seq2 ) );
		for( k=0; gaplen2[j]&&gaplen2[j][k].idatnext!=-1; k++ )
		{
			pos2 = gaplen2[j][k].relend;
			reporterr( ".len=%d, .relend=%d (j=%d)\n", gaplen2[j][k].len, gaplen2[j][k].relend, j );
			reporterr( ".freq=%f\n", gaplen2[j][k].freq );
		}
	}
#endif 

#if 1
	pfac1 = pfac2 = 0.0;
	if( gaplen1i ) for( k=0; (g1=gaplen1i+k)->idatnext!=-1; k++ )
	{
		if( (pos1=g1->relend) == 0 ) pfac1 += g1->freq;
	}

	if( gaplen2j ) for( l=0; (g2=gaplen2j+l)->idatnext!=-1; l++ ) // ososugi!  hash ni atode henkou
	{
		if( (pos2=g2->relend) == 0 ) pfac2 += g2->freq;
	}
#if DEBUG
	reporterr( "\n\nInitial pfac1,pfac2 (%c%d,%c%d) = %f, %f\n", seq1[i], i, seq2[j], j, pfac1, pfac2 );
#endif
	pfac = pfac1 * pfac2 + pfac1 * (1-pfac2) + pfac2 * (1-pfac1);
#if DEBUG
	reporterr( "\n\nInitial pfac (%d,%d) = %f\n", i, j, pfac );
#endif

//	if( pfac ) reporterr( "i,j=%d,%d, Cancel (eq len)? pfac = %f -> ", i, j, pfac );
	if( gaplen1i ) for( k=0; (g1=gaplen1i+k)->idatnext!=-1; k++ ) // ososugi!  hash ni atode henkou
	{
		pos1=g1->relend;
		if( gaplen2j ) for( l=0; (g2=gaplen2j+l)->idatnext!=-1; l++ ) // ososugi!  hash ni atode henkou
		{
			pos2 = gaplen2j[l].relend;
			if     ( pos1 == 0 && pos2 == 0 && g1->len == g2->len ) pfac -= g1->freq * g2->freq;
			else if( pos1 == 0 && pos2 != 0 && g2->len - (pos2-1) > g1->len ) pfac -= g1->freq * g2->freq;
			else if( pos1 != 0 && pos2 == 0 && g1->len - (pos1-1) > g2->len ) pfac -= g1->freq * g2->freq;
		}
	}

#else

	pfac1 = pfac2 = 0.0;
	if( gaplen1i ) for( k=0; (gaplen1i[k].idatnext)!=-1; k++ )
	{
		if( gaplen1i[k].relend == 0 ) pfac1 += gaplen1[i][k].freq;
	}

	if( gaplen2j ) for( l=0; (gaplen2j[l].idatnext)!=-1; l++ ) // ososugi!  hash ni atode henkou
	{
		if( gaplen2j[l].relend == 0 ) pfac2 += gaplen2[j][l].freq;
	}
#if DEBUG
	reporterr( "\n\nInitial pfac1,pfac2 (%c%d,%c%d) = %f, %f\n", seq1[i], i, seq2[j], j, pfac1, pfac2 );
#endif
	pfac = pfac1 * pfac2 + pfac1 * (1-pfac2) + pfac2 * (1-pfac1);
#if DEBUG
	reporterr( "\n\nInitial pfac (%d,%d) = %f\n", i, j, pfac );
#endif

#if 1
//	if( pfac ) reporterr( "i,j=%d,%d, Cancel (eq len)? pfac = %f -> ", i, j, pfac );
	if( gaplen1i ) for( k=0; (gaplen1i[k].idatnext)!=-1; k++ ) // ososugi!  hash ni atode henkou
	{
		pos1=gaplen1i[k].relend;
		if( gaplen2j ) for( l=0; (gaplen2j[l].idatnext)!=-1; l++ ) // ososugi!  hash ni atode henkou
		{
			pos2 = gaplen2j[l].relend;
			if     ( pos1 == 0 && pos2 == 0 && gaplen1i[k].len == gaplen2j[l].len ) pfac -= gaplen1i[k].freq * gaplen2j[l].freq;
			else if( pos1 == 0 && pos2 != 0 && gaplen2j[l].len - (pos2-1) > gaplen1i[k].len ) pfac -= gaplen1i[k].freq * gaplen2j[l].freq;
			else if( pos1 != 0 && pos2 == 0 && gaplen1i[k].len - (pos1-1) > gaplen2j[l].len ) pfac -= gaplen1i[k].freq * gaplen2j[l].freq;
		}
	}
#endif
#endif

#if DEBUG
	reporterr( "\n\nFinal pfac1,pfac2 (%c%d,%c%d, straight) = %f\n\n", seq1[i], i, seq2[j], j, pfac );
#endif
	return( pfac );
}


static void extendgaplencompactx( Gaplen **cpy, Gaplen **orig, int start )
{
	Gaplen *opt, *cpt;
	int l, id;
#if DEBUG
	Gaplen cpybk;
#endif

//	if( start < 0 ) start = 0;

	if( orig[start] == NULL )
	{
		if( cpy[start] ) 
		{
			free( cpy[start] );
			cpy[start] = NULL;
		}
		return;
	}


#if DEBUG
	reporterr( "At first, cpy -> \n" );
	showgaplen( cpy, 100 );
	reporterr( "Look at %d \n", start );
#endif

	if( cpy[start] == NULL )
	{
		l = orig[start][0].npat;
	
		cpy[start] = realloc( cpy[start], (l+2) * sizeof( Gaplen ) );
	
#if 0
		for( l=0; (gl=orig[start][l].idatend)!=-1; l++ )
			cpy[start][l] = orig[start][l]; // freq ha pointer de copy
		cpy[start][l] = orig[start][l]; // dekiru?
#else
		for( opt = orig[start],cpt = cpy[start]; opt->idatnext!=-1; )
			*cpt++ = *opt++;
		*cpt = *opt;
#endif
	}

#if DEBUG
	cpybk = cpy[start][0];
#endif

#if 0
	for( l=0; (opt=orig[start]+l)->idatend!=-1; l++ )
	{
		if( (pos=opt->relend) == 0 ) continue;

		if( cpy[posplus=start+pos] != NULL )
		{
			id = opt->idatend;
//			reporterr( "cpy[%d][%d].len: %d -> %d (relend=%d)\n", start, l, cpy[start][l].len, cpy[posplus][id].len, pos );
			cpy[start][l].len = cpy[posplus][id].len; // Ato de posplus wo tsukawanaiyouni henkou.
			continue; // HITSUYOU!!!
		}
		else
		{
//			reporterr( "cpy[%d][%d].len: %d (relend=%d)\n", start, l, cpy[start][l].len, pos );
		}

#if 0
		for( k=0; orig[start+pos][k].idatend!=-1; k++ )
			;
#else
		optplus = orig[posplus];
		k = optplus->npat;
#endif
	

		cptplus = cpy[posplus] = realloc( cpy[posplus], (k+2) * sizeof( Gaplen ) );
//		cptplus = realloc( cptplus, (k+2) * sizeof( Gaplen ) );

#if 0
		for( k=0; optplus[k].idatend!=-1; k++ )
		{
			cptplus[k] = optplus[k]; // dekiru?
		}
		cptplus[k] = optplus[k]; // dekiru?
#else
		while( optplus->idatend!=-1 ) *cptplus++ = *optplus++;
		*cptplus = *optplus;
#endif
	}
#endif


	if( start == 0 ) return;
	if( cpy[start-1] == NULL ) return;

#if DEBUG
	reporterr( "cpy -> \n" );
	showgaplen( cpy, 100 );
	reporterr( "Look at %d \n", start );
#endif

	for( l=0; orig[start][l].idatnext!=-1; l++ )
	{
		if( (id=orig[start][l].idatprev) == -1 ) continue;

//		if( cpy[start][l].relend != 0 ) cpy[start][l].len = cpy[start-1][id].len; // Shinchou ni
		cpy[start][l].len = cpy[start-1][id].len; // Shinchou ni

//		if( cpy[start][l].len != cpy[start-1][id].len )
#if DEBUG
		if( 1 || cpy[start][l].len != cpy[start-1][id].len )
		{
			reporterr( "Check!! cpy[%d][%d].len=%d, but [start-1][].len=%d, relend=%d\n", start, l, cpy[start][l].len, cpy[start-1][id].len, cpy[start][l].relend );
			reporterr( "orig[%d][%d].len=%d, relend=%d\n", start, l, orig[start][l].len, orig[start][l].relend );
			reporterr( "cpybk.len=%d, relend=%d\n", cpybk.len, cpybk.relend );

		}
		else
		{
//			reporterr( "OK, cpy[%d][%d].len=%d, relend=%d\n", start, l, cpy[start][l].len, cpy[start][l].relend );
		}
#endif
	}

}


#if USEGAPLENHALFORMTX
static void extendgaplenpartly( Gaplen **cpy, Gaplen **orig, int start, int end )
{
	int i, l, gl, extrascope;
	Gaplen *pt;

	if( start < 0 ) start = 0;
//	for( i=start; i<=end; i++ )
//	{
//		if( cpy[i] == (Gaplen *)1 )
//		{
//			end = i-1;
//			break;
////			reporterr( "Okashii! i=%d\n", i );
////			exit( 1 );
//		}
//		if( cpy[i] ) free( cpy[i] );
//		cpy[i] = NULL;
//	}


	extrascope = 0;
#if 0
	for( i=start; i<=end; i++ ) if( orig[i] )
	{ 
		for( pt=orig[i]; (pt->idatend)!=-1; )
		{
			if( (gl=pt++->relend) > extrascope ) extrascope = i+gl-end+1; 
		}
//		extrascope = 10; // Kinji

	}
#else
	if( orig[end] )
	{ 
		for( pt=orig[end]; (pt->idatend)!=-1; )
		{
			if( (gl=pt++->relend) > extrascope ) extrascope = gl; 
		}
//		extrascope = 10; // Kinji

	}
#endif
	end += extrascope;

	for( i=start; i<=end; i++ )
	{
		if( cpy[i] != NULL ) continue;

		if( orig[i] == NULL )
		{
			if( cpy[i] ) free( cpy[i] ); // muda dakedo 
			cpy[i] = NULL;
			continue;
		}

		for( l=0; (gl=orig[i][l].idatend)!=-1; l++ )
			;

		cpy[i] = realloc( cpy[i], (l+2) * sizeof( Gaplen ) );
//		cpy[i] = calloc( sizeof( Gaplen ), l+2 );

		for( l=0; (gl=orig[i][l].idatend)!=-1; l++ )
		{
#if 1
			cpy[i][l] = orig[i][l]; // freq ha pointer de copy
#else
			cpy[i][l].len = gl;
			cpy[i][l].relend = orig[i][l].relend;
			cpy[i][l].freq = orig[i][l].freq;
			cpy[i][l].gapidatend = orig[i][l].gapidatend;
#endif

//			reporterr( "i=%d, l=%d, len=%d, freq=%f, relend=%d\n", i, l, cpy[i][l].len, cpy[i][l].freq, cpy[i][l].relend );
		}
		cpy[i][l] = orig[i][l]; // dekiru?
//		cpy[i][l].relend = -1;
//		cpy[i][l].len = -1;
	}

}
#endif

static void duplicategaplencompactx( Gaplen **cpy, Gaplen **orig, int maxlen, int start, int end )
{
	int i, l;


	if( start < 0 ) start = 0;
	for( i=start; i<=end; i++ )
	{
//		reporterr( "i=%d / %d\n", i, maxlen );
		if( cpy[i] == (Gaplen *)1 )
		{
			end = i-1;
			break;
//			reporterr( "Okashii! i=%d\n", i );
//			exit( 1 );
		}
		if( cpy[i] ) free( cpy[i] );
		cpy[i] = NULL;
	}

	for( i=start; i<=end; i++ )
	{
		if( orig[i] == NULL )
		{
			if( cpy[i] ) free( cpy[i] ); // muda dakedo 
			cpy[i] = NULL;
			continue;
		}

#if 0
		for( l=0; (gl=orig[i][l].idatend)!=-1; l++ )
			;
#else
		l = orig[i][0].npat;
#endif

		cpy[i] = realloc( cpy[i], (l+2) * sizeof( Gaplen ) );
//		cpy[i] = calloc( sizeof( Gaplen ), l+2 );

		for( l=0; orig[i][l].idatnext!=-1; l++ )
		{
			cpy[i][l] = orig[i][l]; // freq ha pointer de copy
//			reporterr( "i=%d, l=%d, len=%d, freq=%f, relend=%d\n", i, l, cpy[i][l].len, cpy[i][l].freq, cpy[i][l].relend );
		}
		cpy[i][l] = orig[i][l]; // dekiru?
//		cpy[i][l].relend = -1;
//		cpy[i][l].len = -1;
	}

	return;
}



#if USEGAPLENHALFORMTX
static void duplicategaplenpartly( Gaplen **cpy, Gaplen **orig, int start, int end )
{
	int i, l, gl, extrascope;
	Gaplen *pt;

	if( start < 0 ) start = 0;
	for( i=start; i<=end; i++ )
	{
		if( cpy[i] == (Gaplen *)1 )
		{
			end = i-1;
			break;
//			reporterr( "Okashii! i=%d\n", i );
//			exit( 1 );
		}
		if( cpy[i] ) free( cpy[i] );
		cpy[i] = NULL;
	}


	extrascope = 0;
#if 0
	for( i=start; i<=end; i++ ) if( orig[i] )
	{ 
		for( pt=orig[i]; (pt->idatend)!=-1; )
		{
			if( (gl=pt++->relend) > extrascope ) extrascope = i+gl-end+1; 
		}
//		extrascope = 10; // Kinji

	}
#else
	if( orig[end] )
	{ 
		for( pt=orig[end]; (pt->idatend)!=-1; )
		{
			if( (gl=pt++->relend) > extrascope ) extrascope = gl; 
		}
//		extrascope = 10; // Kinji

	}
#endif
	end += extrascope;

	for( i=start; i<=end; i++ )
	{
		if( orig[i] == NULL )
		{
			if( cpy[i] ) free( cpy[i] ); // muda dakedo 
			cpy[i] = NULL;
			continue;
		}

		for( l=0; (gl=orig[i][l].idatend)!=-1; l++ )
			;

		cpy[i] = realloc( cpy[i], (l+2) * sizeof( Gaplen ) );
//		cpy[i] = calloc( sizeof( Gaplen ), l+2 );

		for( l=0; (gl=orig[i][l].idatend)!=-1; l++ )
		{
#if 1
			cpy[i][l] = orig[i][l]; // freq ha pointer de copy
#else
			cpy[i][l].len = gl;
			cpy[i][l].relend = orig[i][l].relend;
			cpy[i][l].freq = orig[i][l].freq;
			cpy[i][l].gapidatend = orig[i][l].gapidatend;
#endif

//			reporterr( "i=%d, l=%d, len=%d, freq=%f, relend=%d\n", i, l, cpy[i][l].len, cpy[i][l].freq, cpy[i][l].relend );
		}
		cpy[i][l] = orig[i][l]; // dekiru?
//		cpy[i][l].relend = -1;
//		cpy[i][l].len = -1;
	}

}
#endif


static void gaplenextendnoidatend( Gaplen **cpy, int gapstartpos, int insertionlen )
{
	int l, id, idn, pos, len;

#if 0
//	reporterr( "inserting %d gaps at position %d\n", insertionlen, gapstartpos );
	for( l=0; cpy[gapstartpos] && (id=cpy[gapstartpos][l].idatend) !=-1; l++ )
	{
		pos = cpy[gapstartpos][l].relend;
		cpy[gapstartpos+pos][id].len += insertionlen;
	}
#endif

#if 1
	for( l=0; cpy[gapstartpos] && (id=cpy[gapstartpos][l].idatnext) !=-1; l++ )
	{
		len = cpy[gapstartpos][l].len + insertionlen;
//		reporterr( "ext\n" );
		for( pos=gapstartpos, idn=l; cpy[pos] != NULL && cpy[pos][idn].relend != 0; pos++ )
		{
//			reporterr( "%d, plus %d %d->%d\n", pos, insertionlen, cpy[pos][idn].len, cpy[pos][idn].len+insertionlen );
			cpy[pos][idn].len = len;
			idn = cpy[pos][idn].idatnext;
//			if( pos == gapstartpos + 1 ) break;
			break;
		}
//		reporterr( "end\n" );

		idn = cpy[gapstartpos][l].idatprev;
		if( gapstartpos != 0 && idn != -1 && cpy[gapstartpos-1] ) cpy[gapstartpos-1][idn].len = len;
	}
#endif
}


#if USEGAPLENHALFORMTX

static void gaplenextend( Gaplen **cpy, int gapstartpos, int insertionlen )
{
	int l, id, idn, pos, len;

#if 1
//	reporterr( "inserting %d gaps at position %d\n", insertionlen, gapstartpos );
	for( l=0; cpy[gapstartpos] && (id=cpy[gapstartpos][l].idatend) !=-1; l++ )
	{
		pos = cpy[gapstartpos][l].relend;
		cpy[gapstartpos+pos][id].len += insertionlen;
	}
#endif

#if 1
	for( l=0; cpy[gapstartpos] && (id=cpy[gapstartpos][l].idatend) !=-1; l++ )
	{
		len = cpy[gapstartpos][l].len + insertionlen;
//		reporterr( "ext\n" );
		for( pos=gapstartpos, idn=l; cpy[pos] != NULL && cpy[pos][idn].relend != 0; pos++ )
		{
//			reporterr( "%d, plus %d %d->%d\n", pos, insertionlen, cpy[pos][idn].len, cpy[pos][idn].len+insertionlen );
			cpy[pos][idn].len = len;
			idn = cpy[pos][idn].idatnext;
//			if( pos == gapstartpos + 1 ) break;
//			break;
		}
//		reporterr( "end\n" );

		idn = cpy[gapstartpos][l].idatprev;
		if( gapstartpos != 0 && idn != -1 && cpy[gapstartpos-1] ) cpy[gapstartpos-1][idn].len = len;
	}
#endif
}
#endif

static void copygaplencompactx( Gaplen **cpy, Gaplen **orig, int seqlen, int gapstartpos, int insertionlen, int posincopy, int posinori )
{
	Gaplen *pt, *cpt;



#if DEBUG
	reporterr( "At the head of copygaplencompactx, cpy=\n" );
	showgaplen( cpy+posincopy, 100 );
	reporterr( "At the head of copygaplencompactx, orig=\n" );
	showgaplen( orig+posinori, 100 );
	reporterr( "posinori=%d\n", posinori );
#endif

	if( orig[posinori] == NULL ) return;

//	for( pt=orig[posinori],cpt=cpy[posincopy]; pt->relend==0; ) // zenhan ni relend=0 ga matomatteirukara.
	for( pt=orig[posinori],cpt=cpy[posincopy]; pt->idatnext!=-1; ) // kouhan mo copy
	{
		cpt++->len = pt++->len;
	}		


#if 0
	for( l=0; (id=orig[posinori][l].idatend)!=-1; l++ )
	{
		pos = orig[posinori][l].relend;
		if( pos == 0 ) continue;
		if( orig[posinori+pos] == NULL )
		{
			reporterr( "Okashii\n" );
			PFACERROR = 1;
			continue;
		}

#if 0
		for( k=0; orig[posinori+pos][k].relend==0; k++ ) // zenhan dake
		{
			cpy[posincopy+pos][k].len = orig[posinori+pos][k].len; // dekiru?
		}
#else
		cpy[posincopy+pos][id].len = orig[posinori+pos][id].len; // dekiru?
#endif
	}
#endif


	if( gapstartpos == -1 ) gapstartpos = posincopy;
	gaplenextendnoidatend( cpy, gapstartpos, insertionlen );

#if DEBUG
	reporterr( "At the end of copygaplencompactx, cpy=\n" );
	showgaplen( cpy+posincopy, 100 );
#endif
}


#if USEGAPLENHALF
static void copygaplenrestricted_zurasu( Gaplen **cpy, Gaplen **orig, int seqlen, int gapstartpos, int insertionlen, int startincopy, int endincopy, int startinori, int endinori )
{
	int i, extrascope, gl, j;
	int zure, newend;
	Gaplen *pt, *cpt;
//	int ncopied = 0;

#if 0
// mae houkou nimo renzoku gap de enchou suru hitsuyou ga aru to omou.
	for( i=startinori-1; 0<=i&&i<=seqlen; i-- )
	{
//		reporterr( "i=%d\n", i );
		if( orig[i] == NULL ) break;
		for( pt=orig[i],cpt=cpy[i]; (gl=pt++->len)!=-1; ) cpt++->len = gl;
	}
#endif

	zure = startincopy - startinori; // end ha check shinai

//	int ncopied = 0;
	if( orig[endinori] )
	{ 
		extrascope = 0;
		for( pt=orig[endinori]; (pt->idatend)!=-1; )
		{
			if( (gl=pt++->relend) > extrascope ) extrascope = gl; 
		}
//		extrascope = 10; // Kinji

		newend = endinori + extrascope;
	}
	else newend = endinori;

//	reporterr( "ncopy = %d\n", newend - startinori );
//
#if 0 // extra end wo shizen ni kimereba iranai
	if( newend > seqlen ) newend = seqlen; 
//	if( startinori < 0 ) startinori = 0; 
#endif

	for( i=startinori, j=startincopy; i<=newend; i++, j++ )
	{
		if( orig[i] == NULL ) continue;

//		ncopied += 1;

#if 0
		for( pt=orig[i],cpt=cpy[i]; (gl=pt++->len)!=-1; )
			cpt++->len = gl;
#else
		for( pt=orig[i],cpt=cpy[j]; pt->relend==0; ) // zenhan ni relend=0 ga matomatteirukara.
//		int k;
//		for( k=0; orig[i][k].relend==0; k++ ) // zenhan ni relend=0 ga matomatteirukara.
		{
			cpt++->len = pt++->len;
//			reporterr( "i=%d, k=%d\n", i, k );
//			cpy[i][k].len = orig[i][k].len;
		}		
#endif
	}



#if 0
	for( i=0; i<=seqlen; i++ )
	{
		for( l=0; cpy[i]&&(gl=cpy[i][l].len)!=-1; l++ )
			reporterr( "after copy, i=%d, l=%d, len=%d, freq=%f, relend=%d\n", i, l, cpy[i][l].len, cpy[i][l].freq, cpy[i][l].relend );
	}
#endif

	if( gapstartpos < 0 ) return;

	gaplenextend( cpy, gapstartpos, insertionlen );



//	return;


// TEST
//	for( i=endinori+1; i<=newend; i++ )
	for( i=endincopy+1; i<=newend+zure; i++ )
	{
		if( cpy[i] == NULL ) continue;
		for( j=0; cpy[i][j].idatend!=-1; j++ )
		{
			if( cpy[i][j].relend == 0 ) 
			{
				break;
			}
		}
		if( cpy[i][j].idatend == -1 )
		{
			free( cpy[i] );
			cpy[i] = NULL;
		}
	}






#if 0
	reporterr( "\n" );
	for( i=0; i<=seqlen; i++ )
	{
		for( l=0; cpy[i]&&(gl=cpy[i][l].len)!=-1; l++ )
			reporterr( "after add, i=%d, l=%d, len=%d, freq=%f, relend=%d\n", i, l, cpy[i][l].len, cpy[i][l].freq, cpy[i][l].relend );
	}
#endif
}
#endif

#if USEGAPLENHALFORMTX
static void copygaplenrestricted( Gaplen **cpy, Gaplen **orig, int seqlen, int gapstartpos, int insertionlen, int scopestart, int scopeend )
{
	int i, j, extrascope, gl, endinori, newend;
	Gaplen *pt, *cpt;
//	int ncopied = 0;

#if 0
// mae houkou nimo renzoku gap de enchou suru hitsuyou ga aru to omou.
	for( i=scopestart-1; 0<=i&&i<=seqlen; i-- )
	{
//		reporterr( "i=%d\n", i );
		if( orig[i] == NULL ) break;
		for( pt=orig[i],cpt=cpy[i]; (gl=pt++->len)!=-1; ) cpt++->len = gl;
	}
#endif

//	int ncopied = 0;
	endinori = scopeend;
	if( orig[scopeend] )
	{ 
		extrascope = 0;
		for( pt=orig[scopeend]; (pt->idatend)!=-1; )
		{
			if( (gl=pt++->relend) > extrascope ) extrascope = gl; 
		}
//		extrascope = 10; // Kinji

		scopeend += extrascope;
	}
	newend = scopeend;

//	reporterr( "ncopy = %d\n", scopeend - scopestart );
//
#if 0 // extra end wo shizen ni kimereba iranai
	if( scopeend > seqlen ) scopeend = seqlen; 
//	if( scopestart < 0 ) scopestart = 0; 
#endif

	if( scopestart < 0 ) scopestart = 0; 
	for( i=scopestart; i<=scopeend; i++ )
	{
		if( orig[i] == NULL ) continue;

//		ncopied += 1;

#if 0
		for( pt=orig[i],cpt=cpy[i]; (gl=pt++->len)!=-1; )
			cpt++->len = gl;
#else
		for( pt=orig[i],cpt=cpy[i]; pt->relend==0; ) // zenhan ni relend=0 ga matomatteirukara.
//		int k;
//		for( k=0; orig[i][k].relend==0; k++ ) // zenhan ni relend=0 ga matomatteirukara.
		{
			cpt++->len = pt++->len;
//			reporterr( "i=%d, k=%d\n", i, k );
//			cpy[i][k].len = orig[i][k].len;
		}		
#endif
	}



#if 0
	for( i=0; i<=seqlen; i++ )
	{
		for( l=0; cpy[i]&&(gl=cpy[i][l].len)!=-1; l++ )
			reporterr( "after copy, i=%d, l=%d, len=%d, freq=%f, relend=%d\n", i, l, cpy[i][l].len, cpy[i][l].freq, cpy[i][l].relend );
	}
#endif

	if( gapstartpos < 0 ) return;

	gaplenextend( cpy, gapstartpos, insertionlen );

	return;

// TEST extra scope de tsukaunoha end dake?
	for( i=endinori+1; i<=newend; i++ )
	{
		if( cpy[i] == NULL ) continue;
		for( j=0; cpy[i][j].idatend!=-1; j++ )
		{
			if( cpy[i][j].relend == 0 ) break;
		}
		if( cpy[i][j].idatend == -1 )
		{
			free( cpy[i] );
			cpy[i] = NULL;
		}
	}



#if 0
	reporterr( "\n" );
	for( i=0; i<=seqlen; i++ )
	{
		for( l=0; cpy[i]&&(gl=cpy[i][l].len)!=-1; l++ )
			reporterr( "after add, i=%d, l=%d, len=%d, freq=%f, relend=%d\n", i, l, cpy[i][l].len, cpy[i][l].freq, cpy[i][l].relend );
	}
#endif
}
#endif

#if 1
static void freegaplenpartly( Gaplen **mtx, int startpos, int endpos )
{
	int i;
	Gaplen **pt;
	if( startpos < 0 ) startpos = 0;

	for( i=startpos; i<=endpos; i++ )
	{
		if( *(pt=mtx+i) == (Gaplen *)1 ) break;
		if( *pt ) free( *pt );
		*pt = NULL;
	}
}
#else
static void freegaplenpartly( Gaplen **mtx, int startpos, int endpos )
{
	int i;
	if( startpos < 0 ) startpos = 0;

	for( i=startpos; i<=endpos; i++ )
	{
		if( mtx[i] == (Gaplen *)1 ) break;
		if( mtx[i] ) free( mtx[i] );
		mtx[i] = NULL;
	}
}
#endif


double D__align( double **n_dynamicmtx, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, int constraint, double *impmatch, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{

//	int k;
	register int i, j;




	int lasti, lastj;      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
//	double fpenalty = (double)penalty;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS int **ijp;
	static TLS int mpi, *mp;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS char **mseq1;
	static TLS char **mseq2;
	static TLS char **mseq;
	static TLS double **cpmx1;
	static TLS double **cpmx2;
	static TLS int **intwork;
	static TLS double **doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;
#if USEGAPLENHALF
	Gaplen ****gaplen1half = NULL; // NULL ga iru to omou.
	Gaplen ****gaplen2half = NULL; // NULL ga iru to omou.
#endif
#if USEGAPLENMTX
	Gaplen ****gaplen1mtx = NULL; // NULL ga iru to omou.
	Gaplen ****gaplen2mtx = NULL; // NULL ga iru to omou.
#endif
	static TLS Gaplen **gaplen1 = NULL; // NULL ga iru to omou.
	static TLS Gaplen **gaplen2 = NULL; // NULL ga iru to omou.
	static TLS Gaplen ***gaplen1jprev = NULL;
	static TLS Gaplen ***gaplen2jprev = NULL;
	static TLS Gaplen ***gaplen1jcurr = NULL;
	static TLS Gaplen ***gaplen2jcurr = NULL;
	static TLS Gaplen ***gaplen1icurr = NULL;
	static TLS Gaplen ***gaplen2icurr = NULL;
	static TLS Gaplen ***gaplen1jbestkamo = NULL;
	static TLS Gaplen ***gaplen2jbestkamo = NULL;
	static TLS Gaplen ***gaplen1ibestkamo = NULL;
	static TLS Gaplen ***gaplen2ibestkamo = NULL;
	static TLS Gaplen ***gaplen1jbest = NULL;
	static TLS Gaplen ***gaplen2jbest = NULL;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	static TLS Gaplen ****gaplens = NULL;

	Gaplen ***gaplentmp = NULL;
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
	int k;
	double pfac, pfactmp;
	int newgaplen;

//	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s, %f\n", seq1[i], eff1[i] );
//	for( i=0; i<jcyc; i++ ) fprintf( stderr, "%s, %f\n", seq2[i], eff2[i] );

//	reporterr( "\nsizeof(gaplen) = %d\n", sizeof( Gaplen ) );
//	reporterr( "\nsizeof(int) = %d\n", sizeof( int ) );
//	reporterr( "\nsizeof(double) = %d\n", sizeof( double ) );
//	reporterr( "\nsizeof(double*) = %d\n", sizeof( double * ) );


	if( seq1 == NULL )
	{
		if( orlgth1 )
		{
//			fprintf( stderr, "## Freeing local arrays in D__align\n" );
			orlgth1 = 0;
			orlgth2 = 0;

			imp_match_init_strictD( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL, 0, 0 );

			free( mseq1 );
			free( mseq2 );
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );



			free( gaplens );
			if( gaplen1ibestkamo ) FreeGaplenCub( gaplen1ibestkamo ); gaplen1ibestkamo = NULL;
			if( gaplen2ibestkamo ) FreeGaplenCub( gaplen2ibestkamo ); gaplen2ibestkamo = NULL;
			if( gaplen1icurr ) FreeGaplenCub( gaplen1icurr ); gaplen1icurr = NULL;
			if( gaplen2icurr ) FreeGaplenCub( gaplen2icurr ); gaplen2icurr = NULL;
   	                   
			if( gaplen1jprev ) FreeGaplenCub( gaplen1jprev ); gaplen1jprev = NULL;
			if( gaplen2jprev ) FreeGaplenCub( gaplen2jprev ); gaplen2jprev = NULL;
			if( gaplen1jcurr ) FreeGaplenCub( gaplen1jcurr ); gaplen1jcurr = NULL;
			if( gaplen2jcurr ) FreeGaplenCub( gaplen2jcurr ); gaplen2jcurr = NULL;
			if( gaplen1jbestkamo ) FreeGaplenCub( gaplen1jbestkamo ); gaplen1jbestkamo = NULL;
			if( gaplen2jbestkamo ) FreeGaplenCub( gaplen2jbestkamo ); gaplen2jbestkamo = NULL;
			if( gaplen1jbest ) FreeGaplenCub( gaplen1jbest ); gaplen1jbest = NULL;
			if( gaplen2jbest ) FreeGaplenCub( gaplen2jbest ); gaplen2jbest = NULL;
			if( gaplen1 ) FreeGaplenMtx( gaplen1, 1 ); gaplen1 = NULL;
			if( gaplen2 ) FreeGaplenMtx( gaplen2, 1 ); gaplen2 = NULL;
		}
		else
		{
//			fprintf( stderr, "## Not allocated\n" );
		}
		return( 0.0 );
	}


	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );


	reporterr( "%d x %d sequences, len=%d, %d\n", icyc, jcyc, lgth1, lgth2 );


#if 0
	if( lgth1 == 0 || lgth2 == 0 )
	{
		fprintf( stderr, "WARNING (Aalignmm): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		for( i=0; i<icyc; i++ )
		{
			j = lgth2;
			seq1[i][j] = 0;
			while( j ) seq1[i][--j] = *newgapstr;
//			fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
		}
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		for( i=0; i<jcyc; i++ )
		{
			j = lgth1;
			seq2[i][j] = 0;
			while( j ) seq2[i][--j] = *newgapstr;
//			fprintf( stderr, "seq2[i] = %s\n", seq2[i] );
		}
		return( 0.0 );
	}

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;



	if( trywarp )
	{
//		reporterr( "Not supported yet!\n" );
//		exit( 1 );
//		fprintf( stderr, "IN D__align, penalty_shift = %d\n", penalty_shift );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1 to allow shift.\n" );
			exit( 1 );
		}
		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}


#if 0
	fprintf( stderr, "####  eff in SA+++align\n" );
	fprintf( stderr, "####  seq1[0] = %s\n", seq1[0] );
	fprintf( stderr, "####  strlen( seq1[0] ) = %d\n", strlen( seq1[0] ) );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
	fprintf( stderr, "####  seq2[0] = %s\n", seq2[0] );
	fprintf( stderr, "####  strlen( seq2[0] ) = %d\n", strlen( seq2[0] ) );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "eff2[%d] = %f\n", i, eff2[i] );
#endif
	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 0 );
		mseq2 = AllocateCharMtx( njob, 0 );
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

			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );


			free( gaplens );

			if( gaplen1ibestkamo ) FreeGaplenCub( gaplen1ibestkamo ); gaplen1ibestkamo = NULL;
			if( gaplen2ibestkamo ) FreeGaplenCub( gaplen2ibestkamo ); gaplen2ibestkamo = NULL;
			if( gaplen1icurr ) FreeGaplenCub( gaplen1icurr ); gaplen1icurr = NULL;
			if( gaplen2icurr ) FreeGaplenCub( gaplen2icurr ); gaplen2icurr = NULL;
   	                   
			if( gaplen1jcurr ) FreeGaplenCub( gaplen1jcurr ); gaplen1jcurr = NULL;
			if( gaplen1jprev ) FreeGaplenCub( gaplen1jprev ); gaplen1jprev = NULL;
			if( gaplen2jcurr ) FreeGaplenCub( gaplen2jcurr ); gaplen2jcurr = NULL;
			if( gaplen2jprev ) FreeGaplenCub( gaplen2jprev ); gaplen2jprev = NULL;
			if( gaplen1jbestkamo ) FreeGaplenCub( gaplen1jbestkamo ); gaplen1jbestkamo = NULL;
			if( gaplen2jbestkamo ) FreeGaplenCub( gaplen2jbestkamo ); gaplen2jbestkamo = NULL;
			if( gaplen1jbest ) FreeGaplenCub( gaplen1jbest ); gaplen1jbest = NULL;
			if( gaplen2jbest ) FreeGaplenCub( gaplen2jbest ); gaplen2jbest = NULL;
			if( gaplen1 ) FreeGaplenMtx( gaplen1, 1 ); gaplen1 = NULL;
			if( gaplen2 ) FreeGaplenMtx( gaplen2, 1 ); gaplen2 = NULL;


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

		mseq = AllocateCharMtx( njob, ll1+ll2 );

		cpmx1 = AllocateFloatMtx( nalphabets, ll1+2 );
		cpmx2 = AllocateFloatMtx( nalphabets, ll2+2 );

#if FASTMATCHCALC
		doublework = AllocateFloatMtx( MAX( ll1, ll2 )+2, nalphabets ); 
		intwork = AllocateIntMtx( MAX( ll1, ll2 )+2, nalphabets+1 ); 
#else
		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
#endif

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;


//		reporterr( "Allocating gaplen1 and gaplen2\n" );
		gaplen1 = (Gaplen ** )calloc( ll1+2, sizeof( Gaplen * ) );
		gaplen1[ll1+1] = (Gaplen *)1;
		gaplen2 = (Gaplen ** )calloc( ll2+2, sizeof( Gaplen * ) );
		gaplen2[ll2+1] = (Gaplen *)1;


//		reporterr( "Allocating gaplen*\n" );
		gaplen1ibestkamo = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen2ibestkamo = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen1icurr = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen2icurr = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen1jbestkamo = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jbestkamo = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen1jbest = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jbest = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen1jcurr = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jcurr = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen1jprev = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jprev = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );

		gaplens = calloc( sizeof( Gaplen ***), 12 );
		gaplens[0] = gaplen1ibestkamo;
		gaplens[1] = gaplen2ibestkamo;
		gaplens[2] = gaplen1icurr;
		gaplens[3] = gaplen2icurr;
		gaplens[4] = gaplen1jbestkamo;
		gaplens[5] = gaplen2jbestkamo;
		gaplens[6] = gaplen1jbest;
		gaplens[7] = gaplen2jbest;
		gaplens[8] = gaplen1jcurr;
		gaplens[9] = gaplen2jcurr;
		gaplens[10] = gaplen1jprev;
		gaplens[11] = gaplen2jprev;
//		reporterr( "Allocation end\n" );
	}

	{
		int ll1 = lgth1;
		int ll2 = lgth2;

//		reporterr( "Allocating gaplen*i\n" );
		for(i=0; i<ll1+1; i++ ) 
		{
			gaplen1ibestkamo[i] = (Gaplen **)calloc( ll1+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1-i; k++ ) gaplen1ibestkamo[i][k] = NULL;
			gaplen1ibestkamo[i][k] = (Gaplen *)1;

			gaplen2ibestkamo[i] = (Gaplen **)calloc( ll2+2, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1; k++ ) gaplen2ibestkamo[i][k] = NULL;
			gaplen2ibestkamo[i][k] = (Gaplen *)1;

			gaplen1icurr[i] = (Gaplen **)calloc( ll1+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1-i; k++ ) gaplen1icurr[i][k] = NULL;
			gaplen1icurr[i][k] = (Gaplen *)1;

			gaplen2icurr[i] = (Gaplen **)calloc( ll2+2, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1; k++ ) gaplen2icurr[i][k] = NULL;
			gaplen2icurr[i][k] = (Gaplen *)1;
		}
		gaplen1ibestkamo[ll1+1] = NULL;
		gaplen2ibestkamo[ll1+1] = NULL;
		gaplen1icurr[ll1+1] = NULL;
		gaplen2icurr[ll1+1] = NULL;
	
//		reporterr( "Allocating gaplen*j\n" );
		for(i=0; i<ll2+1; i++ ) 
		{
			gaplen1jbestkamo[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jbestkamo[i][k] = NULL;
			gaplen1jbestkamo[i][k] = (Gaplen *)1;

			gaplen2jbestkamo[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jbestkamo[i][k] = NULL;
			gaplen2jbestkamo[i][k] = (Gaplen *)1;

			gaplen1jbest[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jbest[i][k] = NULL;
			gaplen1jbest[i][k] = (Gaplen *)1;

			gaplen2jbest[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jbest[i][k] = NULL;
			gaplen2jbest[i][k] = (Gaplen *)1;

			gaplen1jcurr[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jcurr[i][k] = NULL;
			gaplen1jcurr[i][k] = (Gaplen *)1;

			gaplen2jcurr[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jcurr[i][k] = NULL;
			gaplen2jcurr[i][k] = (Gaplen *)1;

			gaplen1jprev[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jprev[i][k] = NULL;
			gaplen1jprev[i][k] = (Gaplen *)1;

			gaplen2jprev[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jprev[i][k] = NULL;
			gaplen2jprev[i][k] = (Gaplen *)1;

		}
		gaplen1jbestkamo[ll2+1] = NULL;
		gaplen2jbestkamo[ll2+1] = NULL;
		gaplen1jbest[ll2+1] = NULL;
		gaplen2jbest[ll2+1] = NULL;
		gaplen1jcurr[ll2+1] = NULL;
		gaplen2jcurr[ll2+1] = NULL;
		gaplen1jprev[ll2+1] = NULL;
		gaplen2jprev[ll2+1] = NULL;
	}


#if USEGAPLENMTX
/* maikai allocate */

	reporterr( "Allocating gaplenmtx1\n" );
	gaplen1mtx = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen1mtx[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen1mtx[i][j] = (Gaplen **)calloc( lgth1+2, sizeof( Gaplen * ) );
			for( k=0; k<lgth1+1; k++ ) gaplen1mtx[i][j][k] = NULL;
			gaplen1mtx[i][j][k] = (Gaplen *)1;
		}
		gaplen1mtx[i][j] = NULL;
	}
	gaplen1mtx[i] = NULL;

	reporterr( "Allocating gaplenmtx2\n" );
	gaplen2mtx = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen2mtx[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen2mtx[i][j] = (Gaplen **)calloc( lgth2+2, sizeof( Gaplen * ) );
			for( k=0; k<lgth2+1; k++ ) gaplen2mtx[i][j][k] = NULL;
			gaplen2mtx[i][j][k] = (Gaplen *)1;
		}
		gaplen2mtx[i][j] = NULL;
	}
	gaplen2mtx[i] = NULL;

#endif

#if USEGAPLENHALF
	reporterr( "Allocating gaplenhalf1\n" );
	gaplen1half = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen1half[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen1half[i][j] = (Gaplen **)calloc( lgth1+2 - i, sizeof( Gaplen * ) );
			for( k=0; k<lgth1+1-i; k++ ) gaplen1half[i][j][k] = NULL;
			gaplen1half[i][j][k] = (Gaplen *)1;
		}
		gaplen1half[i][j] = NULL;
	}
	gaplen1half[i] = NULL;

	reporterr( "Allocating gaplenhalf2\n" );
	gaplen2half = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen2half[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen2half[i][j] = (Gaplen **)calloc( lgth2+2 - j, sizeof( Gaplen * ) );
			for( k=0; k<lgth2+1-j; k++ ) gaplen2half[i][j][k] = NULL;
			gaplen2half[i][j][k] = (Gaplen *)1;
		}
		gaplen2half[i][j] = NULL;
	}
	gaplen2half[i] = NULL;
#endif


/* maikai allocate */


	for( i=0; i<icyc; i++ )
	{
		mseq1[i] = mseq[i];
		seq1[i][lgth1] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		mseq2[j] = mseq[icyc+j];
		seq2[j][lgth2] = 0;
	}


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
	{
		double t = 0.0;
		for( i=0; i<icyc; i++ )
			t += eff1[i];
	fprintf( stderr, "## totaleff = %f\n", t );
	}
#endif

	cpmx_calc_new( seq1, cpmx1, eff1, lgth1, icyc );
	cpmx_calc_new( seq2, cpmx2, eff2, lgth2, jcyc );


//	reporterr( "Counting gaplen\n" );
	gaplencount( icyc, lgth1, gaplen1, seq1, eff1 );
	gaplencount( jcyc, lgth2, gaplen2, seq2, eff2 );
#if DEBUG
	reporterr( "group1 = \n" );
	showgaplen( gaplen1, lgth1 );
	reporterr( "group2 = \n" );
	showgaplen( gaplen2, lgth2 );
#endif
//	reporterr( "done.\n" );


	for( i=0; i<lgth1+1; i++ ) for( j=0; j<lgth2+1; j++ )
	{
#if USEGAPLENMTX
//		duplicategaplen( gaplen1mtx[i][j], gaplen1, lgth1 );
//		duplicategaplen( gaplen2mtx[i][j], gaplen2, lgth2 );

//		duplicategaplenpartly( gaplen2mtx[i][j], gaplen2, j-0, lgth2 ); // anzen
//		duplicategaplenpartly( gaplen1mtx[i][j], gaplen1, i-0, lgth1 ); // anzen
		duplicategaplenpartly( gaplen1mtx[i][j], gaplen1, i-0, i ); // iranaikamo
		duplicategaplenpartly( gaplen2mtx[i][j], gaplen2, j-0, j ); // iranaikamo
#endif
#if USEGAPLENHALF

//		duplicategaplenpartly( gaplen1half[i][j], gaplen1+i, 0, lgth1-i ); // KOKO de setsuyaku dekiru to omou.
//		duplicategaplenpartly( gaplen2half[i][j], gaplen2+j, 0, lgth2-j ); //  originally, j-1, lgth2
		duplicategaplenpartly( gaplen1half[i][j], gaplen1+i, 0, 0 ); // test
		duplicategaplenpartly( gaplen2half[i][j], gaplen2+j, 0, 0 ); // test
#endif
	}


#if USEGAPLENMTX
	reporterr( "Duplicating gaplen*mtx\n" );
	for( i=0; i<lgth1+1; i++ )
	{
//		addnewgaplen( gaplen1mtx[i][0], gaplen1, gaplen1, lgth1, -1, 0 );
//		addnewgaplen( gaplen2mtx[i][0], gaplen2, gaplen2, lgth2, 0, i );
//		duplicategaplenpartly( gaplen1mtx[i][0], gaplen1, 0, lgth1 );
//		duplicategaplenpartly( gaplen2mtx[i][0], gaplen2, 0, lgth2 );
		copygaplenrestricted( gaplen2mtx[i][0], gaplen2, lgth2, 0, i, 0, 0 );
	}
#endif
#if USEGAPLENHALF
	reporterr( "Duplicating gaplen*mtx\n" );
	for( i=0; i<lgth1+1; i++ )
	{
		copygaplenrestricted( gaplen2half[i][0], gaplen2, lgth2, 0, i, 0, 0 );
	}
#endif



	for( i=0; i<1; i++ )
	{
//		duplicategaplencompactx( gaplen1icurr[i], gaplen1, lgth1, i-0, lgth1 ); //originally, 0, lgth1
//
//		duplicategaplencompactx( gaplen1icurr[i], gaplen1+i, lgth1-i, 0, lgth1-i ); // half
		duplicategaplencompactx( gaplen1icurr[i], gaplen1+i, lgth1-i, 0, 1 ); //  0, 1  hitsuyou


//		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, lgth2 ); // ichiou zenbu
		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, 0 );


		copygaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, i, 0, 0 ); // -> zurasu -> error?


//		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1, lgth1, 0, 1 );
//		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1+i, lgth1-i, 0, 1 ); // half
		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1+i, lgth1-i, 0, 0 ); // half
//		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, lgth2 );
		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, 0 );
//		copygaplenrestricted( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, 0 ); // -> zurasu -> error?
//		copygaplenrestricted_zurasu( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, lgth2, 0, lgth2 ); // -> zurasu -> error?
		copygaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, 0 ); // -> zurasu -> error?
	}

//	reporterr( "Duplicating gaplen*j*curr \n" );
//	int nduplicated = 0;
	for( j=0; j<lgth2+1; j++ )
//	for( j=0; j<1; j++ )
	{
#if USEGAPLENMTX
//		addnewgaplen( gaplen1mtx[0][j], gaplen1, gaplen1, lgth1, 0, j );
//		addnewgaplen( gaplen2mtx[0][j], gaplen2, gaplen2, lgth2, -1, 0 );
//		duplicategaplenpartly( gaplen1mtx[0][j], gaplen1, 0, lgth1 );
//		duplicategaplenpartly( gaplen2mtx[0][j], gaplen2, 0, lgth2 );
		copygaplenrestricted( gaplen1mtx[0][j], gaplen1, lgth1, 0, j, 0, 0 );
#endif

#if USEGAPLENHALF
		copygaplenrestricted( gaplen1half[0][j], gaplen1, lgth1, 0, j, 0, 0 );
#endif
//		reporterr( "1jcurr?\n" );
//		duplicategaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, lgth1 ); // test
		duplicategaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, 0 ); // dame?
//		reporterr( "done\n" );
//		duplicategaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, 0 ); //test

//		duplicategaplencompactx( gaplen2jcurr[j], gaplen2, lgth2, j-0, lgth2 ); // full
//		duplicategaplencompactx( gaplen2jcurr[j], gaplen2+j, lgth2-j, 0, lgth2-j ); //half! KOKO????
//reporterr( "starting suspicious duplication\n" );
		duplicategaplencompactx( gaplen2jcurr[j], gaplen2+j, lgth2-j, 0, 0 ); //half!
//reporterr( "starting suspicious copy\n" );
		copygaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, j, 0, 0 ); // TEST
//reporterr( "finished\n" );

//		reporterr( "Checking gaplen1jcurr[%d]\n", j );
//		checkgaplen( gaplen1jcurr[j], 100 );
//		reporterr( "Checking gaplen2jcurr[%d]\n", j );
//		checkgaplen( gaplen2jcurr[j], 100 );
	}

//	reporterr( "nduplicated (corrected) = %d\n", nduplicated );

//	reporterr( "Duplicating gaplen*j*prev \n\n" );
	for( j=0; j<lgth2+1; j++ ) // allocate nominotame, atode uwagaki
	{
//		duplicategaplencompactx( gaplen1jprev[j], gaplen1, lgth1, 0, lgth1 );
		duplicategaplencompactx( gaplen1jprev[j], gaplen1, lgth1, 0, 0 ); // TEST
//		duplicategaplencompactx( gaplen2jprev[j], gaplen2, lgth2, j-0, lgth2 ); // originally, 0,lgth2
//		duplicategaplencompactx( gaplen2jprev[j], gaplen2+j, lgth2-j, 0, lgth2-j ); // half
		duplicategaplencompactx( gaplen2jprev[j], gaplen2+j, lgth2-j, 0, 0 ); // half


		copygaplencompactx( gaplen1jprev[j], gaplen1, lgth1, 0, j, 0, 0 ); // wasuretetakamo


	}


//	reporterr( "Duplicating gaplen*j*best \n\n" );

	for( j=0; j<lgth2+1; j++ )
//	for( j=0; j<1; j++ )
	{
//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, lgth1 ); // KOKO
//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 0 ); // test
//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 1 );
		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 0 );


//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 1 );
//		duplicategaplencompactx( gaplen2jbestkamo[j], gaplen2, lgth2, j-0, j+1 ); // originally, 0, j+1
		duplicategaplencompactx( gaplen2jbestkamo[j], gaplen2+j, lgth2-j, 0, 1 ); // half!
		copygaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, j, 0, 0 ); // TEST

//		duplicategaplencompactx( gaplen1jbest[j], gaplen1, lgth1, 0, lgth1 ); // KOKO
		duplicategaplencompactx( gaplen1jbest[j], gaplen1, lgth1, 0, 0 ); // test
//		duplicategaplencompactx( gaplen2jbest[j], gaplen2, lgth2,j-0, j+1 ); // originally, 0,j+1
		duplicategaplencompactx( gaplen2jbest[j], gaplen2+j, lgth2-j, 0, 1 ); // half!
		copygaplencompactx( gaplen1jbest[j], gaplen1, lgth1, 0, j, 0, 0 ); // TEST


	}

//	reporterr( "Duplication end\n" );


#if 0
	reporterr( "Checking gaplen1icurr\n" );
	checkgaplen( gaplen1icurr[0], 100 );
	reporterr( "Checking gaplen2icurr\n" );
	checkgaplen( gaplen2icurr[0], 100 );
#endif




//	showgaplen( gaplen1jcurr[50], lgth2 );

	currentw = w1;
	previousw = w2;

	match_calc( n_dynamicmtx, initverticalw, cpmx2, cpmx1, 0, lgth1, doublework, intwork, 1 );
	if( constraint )
		imp_match_out_vead_tate( initverticalw, 0, lgth1 ); // 060306

	match_calc( n_dynamicmtx, currentw, cpmx1, cpmx2, 0, lgth2, doublework, intwork, 1 );
	if( constraint )
		imp_match_out_vead( currentw, 0, lgth2 ); // 060306
#if 0 // -> tbfast.c // impossible
	if( localhom )
		imp_match_calc( n_dynamicmtx, currentw, icyc, jcyc, lgth1, lgth2, seq1, seq2, eff1, eff2, localhom, 1, 0 );

#endif

	for( j=1; j<lgth2+1; j++ )
	{
		pfac = calcpfac_gap_noidatend( gaplen1, gaplen2, j, 0, j, seq1[0], seq2[0], 0 ); 
//		reporterr( "computing initial end gap penalty for %c-%c, i=0, j=%d, pfac=%f\n", seq1[0][0], seq2[0][j], j, pfac );
//		reporterr( "%c-%c, i=0, j=%d, currentw[j]=%f -> ", seq1[0][0], seq2[0][j], j, currentw[j] );
		currentw[j] += fpenalty * pfac; // tekitou
//		reporterr( " %f\n", currentw[j] );
	}
	for( i=1; i<lgth1+1; i++ )
	{
		pfac = calcpfac_gap_noidatend( gaplen2, gaplen1, i, 0, i, seq2[0], seq1[0], 0 );
//		reporterr( "computing initial end gap penalty for %c-%c, i=%d, j=0, pfac=%f\n", seq1[0][i], seq2[0][0], i, pfac );
		initverticalw[i] += fpenalty * pfac; // tekitou
	}



	for( j=1; j<lgth2+1; ++j ) 
	{


#if ALGZGAP
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2[j-1]; mp[j] = 0;;
#else
		pfac = calcpfac_gapex_noidatend( gaplen2, gaplen1, j, 1, j, seq2[0], seq1[0], 1 );
#if DEBUG
		reporterr( "%c-%c, INITIAL jgap extension check, pfac = %f\n\n", seq1[0][j], '=', pfac );
#endif
		m[j] = currentw[j-1] + fpenalty * pfac; 
		mp[j] = 0;
#endif
	}
	if( lgth2 == 0 )
		lastverticalw[0] = 0.0; // Falign kara yobaretatoki kounarukanousei ari
	else
		lastverticalw[0] = currentw[lgth2-1];

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;


	for( i=1; i<lasti; i++ )
	{
//		reporterr( "i = %d\n", i );

//		reporterr( "err1? i=%d/%d\n", i, lgth1 );
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
			cleargaplens( gaplens );
//			fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
			*chudanres = 1;
			return( -1.0 );
		}
#endif


		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if 1
		gaplentmp = gaplen1jprev;
		gaplen1jprev = gaplen1jcurr;
		gaplen1jcurr = gaplentmp;

		gaplentmp = gaplen2jprev;
		gaplen2jprev = gaplen2jcurr;
		gaplen2jcurr = gaplentmp;

#if DEBUG
		reporterr( "Entering a small j loop, i=%d\n", i );
		for( j=1; j<lgth2+1; j++ )
		{
			reporterr( "before j loop, i=%d, gaplen2jcurr[%d] = \n", i, j );
			showgaplen( gaplen2jcurr[j], 100 );
			reporterr( "\n" );
			reporterr( "before j loop, i=%d, gaplen2prev[%d] = \n", i, j );
			showgaplen( gaplen2jprev[j], 100 );
			reporterr( "\n" );
		}
#endif
#else
		
		reporterr( "Entering a small j loop, i=%\n", i );
		for( j=1; j<lgth2+1; j++ )
		{
//			addnewgaplen( gaplen1jprev[j], gaplen1jcurr[j], gaplen1, lgth1, -1, 0 );
//			addnewgaplen( gaplen2jprev[j], gaplen2jcurr[j], gaplen2, lgth2, -1, 0 );
			reporterr( "err1? j=%d/%d\n", j, lgth2 );
			copygaplencompactx( gaplen1jprev[j-1], gaplen1jcurr[j-1], lgth1, -1, 0, i-1, i-1 ); // TEST
			reporterr( "err1? j=%d/%d\n", j, lgth2 );
			copygaplencompactx( gaplen2jprev[j-1], gaplen2jcurr[j-1], lgth2, -1, 0, j-1, j-1 ); // TETS
#if DEBUG
			reporterr( "before j loop, i=%d, gaplen2jcurr[%d] = \n", i, j );
			showgaplen( gaplen2jcurr[j], 100 );
			reporterr( "\n" );
			reporterr( "before j loop, i=%d, gaplen2prev[%d] = \n", i, j );
			showgaplen( gaplen2jprev[j], 100 );
			reporterr( "\n" );
#endif
		}
#endif

//		reporterr( "err2? i=%d/%d\n", i, lgth1 );

//		duplicategaplencompactx( gaplen1icurr[i], gaplen1, lgth1, i, i+1 ); //originally 0, i+1
//		reporterr( "gaplen+0=\n");
//		showgaplen( gaplen1, 10 );
//		reporterr( "i=%d, lgth1=%d, lgth1-i=%d, gaplen+i-1=\n", i, lgth1, lgth1-i );
//		showgaplen( gaplen1+i-1, 100 );
		duplicategaplencompactx( gaplen1icurr[i], gaplen1+i, lgth1-i, 0, 1 ); // half!!
//		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, lgth2 ); // KOKO
		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, 0 ); // test
		copygaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, i, 0, 0 ); // IRU? TEST



//		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1, lgth1, i, i+1 ); //originally 0, i+1
		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1+i, lgth1-i, 0, 1 ); //half
//		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, lgth2 ); // ORIGINALLY, 0, lgth2
		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, 0 ); // ORIGINALLY, 0, lgth2
//		copygaplenrestricted( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, lgth2, 0, 0 ); // IRU? // TEST
		copygaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, 0 ); // IRU? // TEST

		extendgaplencompactx( gaplen1jprev[0], gaplen1, i ); // ???


//		addnewgaplen( gaplen1jprev[0], gaplen1icurr[i-1], gaplen1, lgth1, -1, 0 );
//		addnewgaplen( gaplen2jprev[0], gaplen2icurr[i-1], gaplen2, lgth2, -1, 0 );
//		copygaplenrestricted( gaplen1jprev[0], gaplen1icurr[i-1], lgth1, -1, 0, i, i ); // i-1, i da to omou.
		copygaplencompactx( gaplen1jprev[0], gaplen1icurr[i-1], lgth1-i, -1, 0, i, 1 ); // half? lgth1-i?
//		copygaplenrestricted( gaplen2jprev[0], gaplen2icurr[i-1], lgth2, -1, 0, 0, 0 );
		copygaplencompactx( gaplen2jprev[0], gaplen2icurr[i-1], lgth2-j, -1, 0, 0, 0 ); // half?? lgth2-j?


		match_calc( n_dynamicmtx, currentw, cpmx1, cpmx2, i, lgth2, doublework, intwork, 0 );
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
		if( constraint )
		{
//			fprintf( stderr, "Calling imp_match_calc (o) lgth = %d, i = %d\n", lgth1, i );
#if  0
			imp_match_out_vead( currentw, i, lgth2 );
#else
			imp_match_out_vead( currentw, i, lgth2 );
#endif
		}
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

#if 0
		fprintf( stderr, "%c ", seq1[0][i] );
		for( j=0; j<lgth2+1; j++ )
		{
			fprintf( stderr, "%5.0f ", currentw[j] );
		}
		fprintf( stderr, "\n"  );
#endif
	
//		mi = previousw[0] + ogcp2[1]; mpi = 0;





#if ALGZGAP
		mi = previousw[0] + ogcp2[1] * gapfreq1[i-1]; mpi=0;
#else
		pfac = calcpfac_gapex_noidatend( gaplen1, gaplen2, i, 1, i, seq1[0], seq2[0], 1 );
#if DEBUG
		reporterr( "%c-%c, INITIAL igap extension check, pfac = %f\n\n", '=', seq2[0][j], pfac );
#endif
		mi = previousw[0] + fpenalty * pfac; 
		mpi=0;
#endif
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




//		reporterr( "\n\ni=%d, %c\n", i, seq1[0][i] );
		for( j=1; j<lastj; j++ )
		{

#if DEBUG
			reporterr( "*****   %c%d-%c%d   ******* \n", seq1[0][i], i, seq2[0][j], j );
			reporterr( "mpi=%d (%c), *mpjpt=%d (%c)\n", mpi, seq2[0][mpi], *mpjpt, seq1[0][*mpjpt] );
#endif


// Hitsuyou na bun dake tsuika
#if USEGAPLENMTX
			extendgaplenpartly( gaplen1mtx[i-1][mpi], gaplen1, i, i );
			extendgaplenpartly( gaplen2mtx[i-1][mpi], gaplen2, j, j );
			extendgaplenpartly( gaplen1mtx[*mpjpt][j-1], gaplen1, i, i );
			extendgaplenpartly( gaplen2mtx[*mpjpt][j-1], gaplen2, j, j );
			extendgaplenpartly( gaplen1mtx[i-1][j-1], gaplen1, i, i );
			extendgaplenpartly( gaplen2mtx[i-1][j-1], gaplen2, j, j );
#endif
#if USEGAPLENHALF
			extendgaplenpartly( gaplen1half[i-1][mpi], gaplen1+i-1, 1, 1 );
			extendgaplenpartly( gaplen2half[i-1][mpi], gaplen2+mpi, j-mpi, j-mpi );
			extendgaplenpartly( gaplen1half[*mpjpt][j-1], gaplen1+*mpjpt, i-*mpjpt, i-*mpjpt );
			extendgaplenpartly( gaplen2half[*mpjpt][j-1], gaplen2+j-1, 1, 1 );
			extendgaplenpartly( gaplen1half[i-1][j-1], gaplen1+i-1, 1, 1 );
			extendgaplenpartly( gaplen2half[i-1][j-1], gaplen2+j-1, 1, 1 );

#endif

//			reporterr( "extending gaplen1icurr\n" );
			extendgaplencompactx( gaplen1icurr[i-1], gaplen1+i-1, 1 ); // iruhazu
//			reporterr( "extending gaplen2icurr\n" );
			extendgaplencompactx( gaplen2icurr[i-1], gaplen2, j ); // iruhazu
//			reporterr( "extending gaplen1jprev[j-1], j-1=%d\n", j-1 );
			extendgaplencompactx( gaplen1jprev[j-1], gaplen1, i );
//			reporterr( "extending gaplen1jcurr, j-1=%d\n", j-1 );
			extendgaplencompactx( gaplen1jcurr[j-1], gaplen1, i );
//			reporterr( "extending gaplen2jprev\n" );
			extendgaplencompactx( gaplen2jprev[j-1], gaplen2+j-1, 1 );
//			reporterr( "extending gaplen2jcurr\n" );
			extendgaplencompactx( gaplen2jcurr[j-1], gaplen2+j-1, 1 );
//			reporterr( "extending gaplen1jbest[j-1]\n" );
			extendgaplencompactx( gaplen1jbest[j-1], gaplen1, i );
//			reporterr( "extending gaplen1jbestkamo[j-1]\n" );
			extendgaplencompactx( gaplen1jbestkamo[j-1], gaplen1, i );
//			reporterr( "extending gaplen1jprev[mpi], j-1=%d\n", j-1 );
			extendgaplencompactx( gaplen1jprev[mpi], gaplen1, i ); // full
//			reporterr( "extending gaplen2jprev[mpi]\n" );
			extendgaplencompactx( gaplen2jprev[mpi], gaplen2+mpi, j-mpi ); // half
//			reporterr( "extending gaplen1ibestkamo[i-1]\n" );
			extendgaplencompactx( gaplen1ibestkamo[i-1], gaplen1+i-1, 1 );
//			reporterr( "extending gaplen2ibestkamo[i-1]\n" );
			extendgaplencompactx( gaplen2ibestkamo[i-1], gaplen2, j );


#if DEBUG
			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jcurr[j]=\n", i, j );
			showgaplen( gaplen2jcurr[j], 100 );

			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jcurr[j-1]=\n", i, j );
			showgaplen( gaplen2jcurr[j-1], 100 );


			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jprev[j]=\n", i, j );
			showgaplen( gaplen2jprev[j], 100 );

			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jprev[j-1]=\n", i, j );
			showgaplen( gaplen2jprev[j-1], 100 );
#endif


#ifdef xxxenablemultithread
//			fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
			if( chudanpt && *chudanpt != chudanref ) 
			{
				cleargaplens( gaplens );
//				fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
				*chudanres = 1;
				return( -1.0 );
			}
#endif
#if USEGAPLENHALF
// i==248, j==80 wo check
#if DEBUG80
			if( j==80 )
			{
				reporterr( "When i==%d, j==%d,\n", i, j );
				reporterr( "gaplen2jprev[j-1=%d]=\n", j-1 );
				showgaplen( gaplen2jprev[j-1], lgth1 );
				reporterr( "gaplen2half[i-1=%d][j-1=%d]=\n", i-1, j-1 );
				showgaplen( gaplen2half[i-1][j-1], lgth1 );
			}
			if( j==79 )
			{
				reporterr( "When i==%d, j==%d,\n", i, j );
				reporterr( "gaplen2jprev[j-1=%d]=\n", j-1 );
				showgaplen( gaplen2jprev[j-1], lgth1 );
				reporterr( "gaplen2half[i-1=%d][j-1=%d]=\n", i-1, j-1 );
				showgaplen( gaplen2half[i-1][j-1], lgth1 );
			}
#endif
#endif


//			pfac = calcpfac( gaplen1jprev[j-1], gaplen2jprev[j-1], i, j, seq1[0], seq2[0] );
//reporterr( "#### COMPACT, i,j=%d,%d\n", i, j );
			pfac = calcpfacnoidatend( gaplen1jprev[j-1], gaplen2jprev[j-1], i, 1, seq1[0], seq2[0]+j, one ); // 1j->full, 2j->half
#if USEGAPLENMTX
//reporterr( "#### FULL, i,j=%d,%d\n", i, j );
			pfactmp = calcpfac( gaplen1mtx[i-1][j-1], gaplen2mtx[i-1][j-1], i, j, seq1[0], seq2[0], one );
#endif
#if USEGAPLENHALF
//reporterr( "#### HALF, i,j=%d/%d,%d/%d\n", i, lgth1, j, lgth2 );
//			showgaplen( gaplen2half[i-1][j-1], lgth2-j );
			pfactmp = calcpfac( gaplen1half[i-1][j-1], gaplen2half[i-1][j-1], 1, 1, seq1[0]+i, seq2[0]+j, zero );
#endif
#if USEGAPLENMTX + USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(straight) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
				PFACERROR = 1;
				exit( 1 );
			}
#endif
//if( i==50 && j==135 ) exit( 1 );


//			reporterr( "i,j=%d,%d, *prept = %f\n", i, j, *prept );

#if ALGZSTRAIGHT
			wm = *prept;  // Machigai!!
#else
			wm = *prept + fpenalty * pfac;
#endif
			*ijppt = 0;


#if DEBUG
			if( i == j )
			{
				fprintf( stderr, "\n i=%d, j=%d %c, %c ", i, j, seq1[0][i], seq2[0][j] );
				fprintf( stderr, "%5.0f, pfac for straight =%f\n", wm, pfac );
			}
#endif
			newgaplen = j-mpi-1;


//if( i == 53 && j == 93 ) exit( 1 );




//			pfac = calcpfac_gap_incomplete( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], newgaplen, i, j, seq1[0], seq2[0], 0 ); // i-1
			pfac = calcpfac_gap_noidatend( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], newgaplen, 1, j, seq1[0]+i-1, seq2[0], 0 ); // i-1
#if USEGAPLENMTX
			pfactmp = calcpfac_gap_incomplete( gaplen1mtx[i-1][mpi], gaplen2mtx[i-1][mpi], newgaplen, i, j, seq1[0], seq2[0], 1 );
#endif
#if USEGAPLENHALF
			pfactmp = calcpfac_gap_incomplete( gaplen1half[i-1][mpi], gaplen2half[i-1][mpi], newgaplen, 1, j-mpi, seq1[0]+i-1, seq2[0]+mpi, 1 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(igap) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
				PFACERROR = 1;
			}
#endif


#if DEBUG
			reporterr( "%c-%c pfac for igap end incomplete = %f\n", seq1[0][i], seq2[0][j], pfac );
			reporterr( "mi when igap end checking = %f\n", mi );
			reporterr( "wm = %f, mi+fpenalty*pfac=%f\n", wm, mi+fpenalty*pfac );
#endif


#if ALGZGAP
			if( (g=mi+*fgcp2pt*gf1va) > wm )
#else
			if( (g=mi+fpenalty*pfac) > wm )
#endif
			{
				wm = g;
				*ijppt = -( j - mpi );
#if DEBUG80
				reporterr( "Selected as igap end! wm = %f, mi = %f\n", wm, mi );
				fprintf( stderr, "Jump from %d-%d (%c-%c) to %d (%c-%c)!\n", i, j, seq1[0][i], seq2[0][j], mpi, seq1[0][i-1], seq2[0][mpi] );
#endif
			}


#if 0 
			fprintf( stderr, "%5.0f->", wm );
#endif
//			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )

#if 0
//			reporterr( "Checking %c, (current pos = %c)\n", seq2[0][j+1], seq2[0][j] );
			sfac = 1.0;
			for( k=0; gaplen2[j+1]&&(gl=gaplen2[j+1][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
			{
//				reporterr( ".len = %d, .relend = %d\n", gaplen2[j+1][k].len, gaplen2[j+1][k].relend );
				if( gl - 1 == gaplen2[j+1][k].relend ) 
				{
					sfac -= gaplen2[j+1][k].freq;
//					reporterr( "Hit! sfac = %f\n", sfac );
				}
			}
			sfac2 = 1.0;
			for( k=0; gaplen1[i+1]&&(gl=gaplen1[i+1][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
				if( gaplen1[i+1][k].relend != -1 ) sfac2 -= gaplen1[i+1][k].freq;
			sfac *= sfac2;
#else
//			sfac = 0.0;
#endif



#if ALGZGAP
			if( (g=*prept+*ogcp2pt*gf1vapre) >= mi )
#else
//			if( (g=*prept + fpenalty * sfac ) >= mi )
			if( (g=*prept ) >= mi )
#endif
			{
//				mpibk = mpi;
//				mi = g - fpenalty * sfac;
				mi = g;
				mpi = j-1;
#if DEBUG80
				reporterr( "Selected as igap start! %c%d-%c%d, mi=%f, g=%f\n", seq1[0][i-1], i-1, seq2[0][mpi], mpi, mi, g );
#endif

#if FREEFREQUENTLY
//				freegaplenpartly( gaplen1ibestkamo[i-1], 0, i-1 );
				freegaplenpartly( gaplen2ibestkamo[i-1], j-3, j-2 );
#endif
//				freegaplenpartly( gaplen1jprev[mpibk], 0, lgth2 ); // full
//				freegaplenpartly( gaplen2jprev[mpibk], 0, lgth2-mpibk ); // half
//				if( gaplen1jprev[mpibk] ) FreeGaplenMtx( gaplen1jprev[mpibk], 0 );
//				gaplen1jprev[mpibk] = NULL;
//				if( gaplen2jprev[mpibk] ) FreeGaplenMtx( gaplen2jprev[mpibk], 0 );
//				gaplen2jprev[mpibk] = NULL;


//				addnewgaplen( gaplen1ibestkamo[i-1], gaplen1jprev[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2ibestkamo[i-1], gaplen2jprev[j-1], gaplen2, lgth2, -1, 0 );
//				copygaplenrestricted( gaplen1ibestkamo[i-1], gaplen1jprev[j-1], lgth1, -1, 0, i, i ); // i-1, i
				copygaplencompactx( gaplen1ibestkamo[i-1], gaplen1jprev[j-1], lgth1, -1, 0, 1, i ); // half
//				copygaplenrestricted( gaplen2ibestkamo[i-1], gaplen2jprev[j-1], lgth2, -1, 0, j, j ); // mpi, j
				copygaplencompactx( gaplen2ibestkamo[i-1], gaplen2jprev[j-1], lgth2, -1, 0, j, 1 ); //half


			}






//			reporterr( "g=%f, *prept=%f, mi=%f\n", g, *prept, mi );


#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if ALGZGAP
			pfac = 0.0; // CHUUI!
#else

//			pfac = calcpfac_gapex( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], i, j, j-mpi, seq1[0], seq2[0], 1 ); // i-1 
			pfac = calcpfac_gapex_noidatend( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], 1, j, j-mpi, seq1[0]+i, seq2[0], 1 ); // 1ibest->half, 2ibest->full
#if USEGAPLENMTX
			pfactmp = calcpfac_gapex( gaplen1mtx[i-1][mpi], gaplen2mtx[i-1][mpi], i, j, j-mpi, seq1[0], seq2[0], 1 );
#endif
#if USEGAPLENHALF
			pfactmp = calcpfac_gapex( gaplen1half[i-1][mpi], gaplen2half[i-1][mpi], 1, j-mpi, j-mpi, seq1[0]+i, seq2[0]+mpi, 1 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(igapex) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
				PFACERROR = 1;
			}
#endif







#if DEBUG
			reporterr( "%c-%c, igap extension check, pfac = %f\n\n", '=', seq2[0][j], pfac );
#endif
#endif
//			reporterr( "mi = %f -> ", mi );
			mi += fpenalty * pfac;
//			reporterr( "mi = %f\n", mi );


//			reporterr( "using %d-%d, %d, %d\n", *mpjpt, j-1, i, j );
			newgaplen = i-*mpjpt-1;
//			pfac = calcpfac_gap_incomplete( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], newgaplen, j, i, seq2[0], seq1[0], 0 ); // j-1 deha???


			pfac = calcpfac_gap_noidatend( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], newgaplen, 1, i, seq2[0]+j-1, seq1[0], 1 ); // 2jbestkamo->half, 1jbestkamo->full
#if USEGAPLENMTX
			pfactmp = calcpfac_gap_incomplete( gaplen2mtx[*mpjpt][j-1], gaplen1mtx[*mpjpt][j-1], newgaplen, j, i, seq2[0], seq1[0], 1 );
#endif
#if USEGAPLENHALF
			pfactmp = calcpfac_gap_incomplete( gaplen2half[*mpjpt][j-1], gaplen1half[*mpjpt][j-1], newgaplen, 1, i-*mpjpt, seq2[0]+j-1, seq1[0]+*mpjpt, 1 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(jgap) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
//				exit( 1 );
				PFACERROR = 1;
			}
#endif

#if ALGZGAP
			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )
#else
			if( (g=*mjpt + fpenalty*pfac) > wm )
#endif
			{
				wm = g;
				*ijppt = +( i - *mpjpt );


#if FREEFREQUENTLY
				freegaplenpartly( gaplen1jbest[j-1], i-3, i-2 );
//				freegaplenpartly( gaplen2jbest[j-1], j-3, j-2 );
#endif


#if DEBUG
				reporterr( "Selected as jgap end!, pfac = %f\n", pfac );
				fprintf( stderr, "Jump from %d (%c) to %d (%c)!\n", j, seq1[0][j], *mpjpt, seq1[0][*mpjpt] );
#endif
//				addnewgaplen( gaplen1jbest[j-1], gaplen1jbestkamo[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jbest[j-1], gaplen2jbestkamo[j-1], gaplen2, lgth2, -1, 0 );
				copygaplencompactx( gaplen1jbest[j-1], gaplen1jbestkamo[j-1], lgth1, -1, 0, i, i );// *mpjpt, i
//				copygaplenrestricted( gaplen2jbest[j-1], gaplen2jbestkamo[j-1], lgth2, -1, 0, j, j ); // j-1, j
				copygaplencompactx( gaplen2jbest[j-1], gaplen2jbestkamo[j-1], lgth2, -1, 0, 1, 1 ); // half!




			}


//			extendgaplenpartly( gaplen1jbest[j-1], gaplen1, i, i ); // tmptmptmp
//			extendgaplenpartly( gaplen2jbest[j-1], gaplen2, 0, 0 ); // tmptmptmp

#if 0
			sfac = 1.0;
			for( l=0; gaplen1[i+1]&&(gl=gaplen1[i+1][l].len)!=-1; l++ ) // ososugi!  hash ni atode henkou
				if( gl - 1 == gaplen1[i+1][l].relend ) sfac -= gaplen1[i+1][l].freq;
			sfac2 = 1.0;
			for( k=0; gaplen2[j+1]&&(gl=gaplen2[j+1][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
				if( gaplen2[j+1][k].relend != -1 ) sfac2 -= gaplen2[j+1][k].freq;
			sfac *= sfac2;
#else
//			sfac = 0.0;
#endif

#if DEBUG
			reporterr( " (jgap start check i=%d) -> *prept=%f, *mjpt=%f\n", i, seq1[0][i], seq2[0][j], *prept, *mjpt );
#endif

#if ALGZGAP
			if( (g=*prept+ ogcp1va* *gf2ptpre) >= *mjpt )
#else
//			if( (g=*prept + fpenalty * sfac ) >= *mjpt )
			if( (g=*prept ) >= *mjpt )
#endif
			{
//				*mjpt = g - fpenalty * sfac;
				*mjpt = g;
				*mpjpt = i-1;
#if DEBUG
				reporterr( "Selected as jgap start!\n" );
#endif


#if FREEFREQUENTLY
				freegaplenpartly( gaplen1jbestkamo[j-1], i-3, i-2 );
//				freegaplenpartly( gaplen2jbestkamo[j-1], j-3, j-2 );
#endif


//				addnewgaplen( gaplen1jbestkamo[j-1], gaplen1jprev[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], gaplen2, lgth2, -1, 0 );
//				reporterr( "copying gaplen1jbestkamo[%d-1] from galpen1jprev, j=%d, i=%d\n", j, j, i );
				copygaplencompactx( gaplen1jbestkamo[j-1], gaplen1jprev[j-1], lgth1, -1, 0, i, i ); // *mpjpt, i
//				copygaplenrestricted( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], lgth2, -1, 0, j, j ); // j-1, j
//				copygaplencompactx( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], lgth2, -1, 0, j, 1 ); // half!
//				reporterr( "copying gaplen2jbestkamo[%d-1] from galpen2jprev\n", j );
				copygaplencompactx( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], lgth2-j, -1, 0, 1, 1 ); // ryouhou half!


//				if( j==2 && i==1 ) exit( 1 );



			}

//			extendgaplenpartly( gaplen1ibestkamo[i-1], gaplen1, 0, 0 ); // tmptmptmp
//			extendgaplenpartly( gaplen2ibestkamo[i-1], gaplen2, j, j ); // tmptmptmp


//			extendgaplenpartly( gaplen1jbestkamo[j-1], gaplen1, i, i ); // tmptmptmp
//			extendgaplenpartly( gaplen2jbestkamo[j-1], gaplen2, 0, 0 ); // tmptmptmp


#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if ALGZGAP
			pfac = 0.0;
#else

//			pfactmp = calcpfac_gapex( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], j, i, i-*mpjpt, seq2[0], seq1[0], 0 ); // j-1 
			pfactmp = calcpfac_gapex_noidatend( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], 1, i, i-*mpjpt, seq2[0]+j, seq1[0], 0 ); // 2jbestkamo->half, 1jbestkamo->full
#if USEGAPLENMTX
			pfac = calcpfac_gapex( gaplen2mtx[*mpjpt][j-1], gaplen1mtx[*mpjpt][j-1], j, i, i-*mpjpt, seq2[0], seq1[0], 0 );
#endif
#if USEGAPLENHALF
			pfac = calcpfac_gapex( gaplen2half[*mpjpt][j-1], gaplen1half[*mpjpt][j-1], 1, i-*mpjpt, i-*mpjpt, seq2[0]+j, seq1[0]+*mpjpt, 0 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(jgapex) pfac=%f, but pfactmp=%f (i,j=%d,%d) diff=%f\n", pfac, pfactmp, i, j, pfac-pfactmp );
//				exit( 1 );
				PFACERROR = 1;
			}
#endif
			pfac = pfactmp;
#if DEBUG
			reporterr( "%c-%c, jgap extension check (j=%d), pfac = %f\n", seq1[0][i], '=', j, pfac );
#endif
#endif
			m[j] += fpenalty * pfac;



			if( trywarp )
			{
#if USE_PENALTY_EX
				if( ( g=*prevwmrecordspt++ + fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] ) ) > wm ) // naka ha osokute kamawanai
#else
				if( ( g=*prevwmrecordspt++ + fpenalty_shift ) > wm ) // naka ha osokute kamawanai
#endif
				{
//					fprintf( stderr, "WARP!!\n" );
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
	
#if 0
				fprintf( stderr, "%5.0f ", wm );
#endif
				curm = *curpt + wm;
	
				if( *wmrecords1pt > *wmrecordspt )
				{
					*wmrecordspt = *wmrecords1pt;
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
				}
				if( curm > *wmrecordspt )
				{
					*wmrecordspt = curm;
					*warpipt = i;
					*warpjpt = j;
				}
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}

#if DEBUG
			reporterr( "extention-x 1j???, before extention-x, j=%d\n", j );
			showgaplen( gaplen1jcurr[j], 100 );
#endif
			extendgaplencompactx( gaplen1jcurr[j], gaplen1, i );

#if DEBUG
			reporterr( "after extention-x\n" );
			showgaplen( gaplen1jcurr[j], 100 );
			reporterr( "extention-x 2j???\n" );
#endif
			extendgaplencompactx( gaplen2jcurr[j], gaplen2+j, 0 );


#if 1
			if( *ijppt < 0 )
			{
#if DEBUG
				reporterr( "Path: %d-%d->%d-%d, i=%d (%c), j=%d (%c), *ijppt=%d\n", i, j, i-1, j+*ijppt, i, seq1[0][i], j, seq2[0][j], *ijppt );
				reporterr( "Inserting %d gaps to gaplen1 and copying gaplen2 (%c%d-%c%d)\n", -*ijppt-1, seq1[0][i], i, seq2[0][j], j );
#endif
#if USEGAPLENMTX
//				addnewgaplen( gaplen1mtx[i][j], gaplen1mtx[i-1][j+*ijppt], gaplen1, lgth1, i, -*ijppt-1 );
//				addnewgaplen( gaplen2mtx[i][j], gaplen2mtx[i-1][j+*ijppt], gaplen2, lgth2, -1, 0 );
				copygaplenrestricted( gaplen1mtx[i][j], gaplen1mtx[i-1][j+*ijppt], lgth1, i, -*ijppt-1, i, i ); // 0, lgth1
				copygaplenrestricted( gaplen2mtx[i][j], gaplen2mtx[i-1][j+*ijppt], lgth2, -1, 0, j, j ); // 0, lgth2
#endif
#if USEGAPLENHALF
				copygaplenrestricted_zurasu( gaplen1half[i][j], gaplen1half[i-1][j+*ijppt], lgth1, 0, -*ijppt-1, 0, 0, 1, 1 ); // 0, lgth1
				copygaplenrestricted_zurasu( gaplen2half[i][j], gaplen2half[i-1][j+*ijppt], lgth2, -1, 0, 0, 0, -*ijppt, -*ijppt ); // 0, lgth2
#endif
//				addnewgaplen( gaplen1jcurr[j], gaplen1jprev[j+*ijppt], gaplen1, lgth1, i, -*ijppt-1 );
//				addnewgaplen( gaplen2jcurr[j], gaplen2jprev[j+*ijppt], gaplen2, lgth2, -1, 0 );
//				reporterr( "copying gaplen1jcurr from gaplen1jbest, with a %d insertion\n", -*ijppt-1 ); 
				copygaplencompactx( gaplen1jcurr[j], gaplen1jprev[j+*ijppt], lgth1, i, -*ijppt-1, i, i ); // scope: i+*ijppt+1, i ?
//				reporterr( "copy end\n" ); 
//				copygaplenrestricted( gaplen2jcurr[j], gaplen2jprev[j+*ijppt], lgth2, -1, 0, j, j );
				copygaplencompactx( gaplen2jcurr[j], gaplen2jprev[j+*ijppt], lgth2, -1, 0, 0, -*ijppt ); // half! ryouho zureteru
			}
			else if( *ijppt > 0 )
			{
#if DEBUG
				reporterr( "Path: %d-%d->%d-%d, i=%d (%c), j=%d (%c), *ijppt=%d\n", i, j, i-*ijppt, j-1, i, seq1[0][i], j, seq2[0][j], *ijppt );
				reporterr( "Copying gaplen1 inserting %d gaps to gaplen2 (%c%d-%c%d)\n", *ijppt-1, seq1[0][i], i, seq2[0][j], j );
#endif
#if USEGAPLENMTX
//				addnewgaplen( gaplen1mtx[i][j], gaplen1mtx[i-*ijppt][j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2mtx[i][j], gaplen2mtx[i-*ijppt][j-1], gaplen2, lgth2, j, *ijppt-1 );
				copygaplenrestricted( gaplen1mtx[i][j], gaplen1mtx[i-*ijppt][j-1], lgth1, -1, 0, i, i ); // 0, lgth1
				copygaplenrestricted( gaplen2mtx[i][j], gaplen2mtx[i-*ijppt][j-1], lgth2, j, *ijppt-1, j, j ); // 0, lgth2
#endif
#if USEGAPLENHALF
				copygaplenrestricted_zurasu( gaplen1half[i][j], gaplen1half[i-*ijppt][j-1], lgth1, -1, 0, 0, 0, *ijppt, *ijppt ); // 0, lgth1
				copygaplenrestricted_zurasu( gaplen2half[i][j], gaplen2half[i-*ijppt][j-1], lgth2, 0, *ijppt-1, 0, 0, 1, 1 ); // 0, lgth2
#endif
//				addnewgaplen( gaplen1jcurr[j], gaplen1jbest[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jcurr[j], gaplen2jbest[j-1], gaplen2, lgth2, j, *ijppt-1 );
				copygaplencompactx( gaplen1jcurr[j], gaplen1jbest[j-1], lgth1, -1, 0, i, i );
//				copygaplenrestricted( gaplen2jcurr[j], gaplen2jbest[j-1], lgth2, j, *ijppt-1, j, j ); // j-*ijppt+1?
//				copygaplenrestricted_zurasu( gaplen2jcurr[j], gaplen2jbest[j-1], lgth2, 0, *ijppt-1, 0, 0, j, j ); // 2jcurr->half, but 2jbest->full, imanotokoro
				copygaplencompactx( gaplen2jcurr[j], gaplen2jbest[j-1], lgth2, 0, *ijppt-1, 0, 1 ); //ryouhou half

			}
			else
#endif
			{
#if DEBUG
				reporterr( "Path: %d-%d->%d-%d, i=%d (%c), j=%d (%c), *ijppt=%d\n", i, j, i-1, j-1, i, seq1[0][i], j, seq2[0][j], *ijppt );
				reporterr( "Copying gaplen1 and gaplen2 (%c%d-%c%d)\n", seq1[0][i], i, seq2[0][j], j );
#endif	
#if USEGAPLENMTX
//				addnewgaplen( gaplen1mtx[i][j], gaplen1mtx[i-1][j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2mtx[i][j], gaplen2mtx[i-1][j-1], gaplen2, lgth2, -1, 0 );
				copygaplenrestricted( gaplen1mtx[i][j], gaplen1mtx[i-1][j-1], lgth1, -1, 0, i, i ); // 0, lgth1
				copygaplenrestricted( gaplen2mtx[i][j], gaplen2mtx[i-1][j-1], lgth2, -1, 0, j, j ); // 0, lgth2
#endif
#if USEGAPLENHALF
				copygaplenrestricted_zurasu( gaplen1half[i][j], gaplen1half[i-1][j-1], lgth1, -1, 0, 0, 0, 1, 1 ); // 0, lgth1
				copygaplenrestricted_zurasu( gaplen2half[i][j], gaplen2half[i-1][j-1], lgth2, -1, 0, 0, 0, 1, 1 ); // 0, lgth2
#endif
//				addnewgaplen( gaplen1jcurr[j], gaplen1jprev[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jcurr[j], gaplen2jprev[j-1], gaplen2, lgth2, -1, 0 );
				copygaplencompactx( gaplen1jcurr[j], gaplen1jprev[j-1], lgth1, -1, 0, i, i );
//				copygaplenrestricted( gaplen2jcurr[j], gaplen2jprev[j-1], lgth2, -1, 0, j, j );
				copygaplencompactx( gaplen2jcurr[j], gaplen2jprev[j-1], lgth2, -1, 0, 0, 1 ); // half
			}

#if DEBUG
			reporterr( "at the end of j loop, gaplen1jcurr[%d] = \n", j );
			showgaplen( gaplen1jcurr[j], 100 );
			reporterr( "at the end of j loop, gaplen1prev[%d] = \n", j );
			showgaplen( gaplen1jprev[j], 100 );
#endif



#if 1
			freegaplenpartly( gaplen1jcurr[j-1], i-3, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // half!
			freegaplenpartly( gaplen1jbestkamo[j-1], i-3, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbestkamo[j-1], j-3, j-2 ); // -1 dehanaku -2??
			freegaplenpartly( gaplen1jbest[j-1], i-3, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbest[j-1], j-3, j-2 ); // -1 dehanaku -2??
#else
			freegaplenpartly( gaplen1jprev[j-1], 0, i-2 ); // -1 dehanaku -2??
			freegaplenpartly( gaplen1jcurr[j-1], 0, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // half!
			freegaplenpartly( gaplen1jbestkamo[j-1], 0, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbestkamo[j-1], j-3, j-2 ); // -1 dehanaku -2??
			freegaplenpartly( gaplen1jbest[j-1], 0, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbest[j-1], j-3, j-2 ); // -1 dehanaku -2??
#endif


#if USEGAPLENMTX
//			freegaplenpartly( gaplen1mtx[i-1][j-1], 0, i-2 );
//			freegaplenpartly( gaplen2mtx[i-1][j-1], 0, j-2 );
#endif


			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
		lastverticalw[i] = currentw[lgth2-1];

#if 1
//		freegaplenpartly( gaplen1icurr[i-1], i-1, i-1 );
		freegaplenpartly( gaplen1icurr[i-1], 0, lgth1-i );
		freegaplenpartly( gaplen2icurr[i-1], 0, lgth2 );
//		freegaplenpartly( gaplen1ibestkamo[i-1], i-1, i-1 );
		freegaplenpartly( gaplen1ibestkamo[i-1], 0, lgth1-i );
		freegaplenpartly( gaplen2ibestkamo[i-1], 0, lgth2 );
#endif

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
#if 0
		fprintf( stderr, "i=%d, %15.5f \n", i, wm );
#endif
//if( i == 2 ) exit( 1 );
	}

	if( trywarp )
	{
//		fprintf( stderr, "wm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}


#if OUTGAP0TRY
	if( !outgap )
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * ( lgth2 - j ) / 2.0;
		for( i=1; i<lgth1+1; i++ )
			lastverticalw[i] -= offset * ( lgth1 - i  / 2.0);
	}
#endif
		
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/
	if( constraint )
	{
		Atracking_localhom( impmatch, currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, warpis, warpjs, warpbase );
	}
	else
		Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase );

	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );

//	fprintf( stderr, "### impmatch = %f\n", *impmatch );

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );
#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
#endif

//	reporterr( "clearing\n" );
	cleargaplens( gaplens );

#if USEGAPLENMTX
/* maikai free */
	reporterr( "Freeing!\n" );
	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen1mtx[i][j] ) FreeGaplenMtx( gaplen1mtx[i][j], 0 );
			gaplen1mtx[i][j] = NULL;
		}
		free( gaplen1mtx[i] );
		gaplen1mtx[i] = NULL;
	}
	free( gaplen1mtx );
	gaplen1mtx = NULL;

	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen2mtx[i][j] ) FreeGaplenMtx( gaplen2mtx[i][j], 0 );
			gaplen2mtx[i][j] = NULL;
		}
		free( gaplen2mtx[i] );
		gaplen2mtx[i] = NULL;
	}
	free( gaplen2mtx );
	gaplen2mtx = NULL;
#endif


#if USEGAPLENHALF
	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen1half[i][j] ) FreeGaplenMtx( gaplen1half[i][j], 0 );
			gaplen1half[i][j] = NULL;
		}
		free( gaplen1half[i] );
		gaplen1half[i] = NULL;
	}
	free( gaplen1half );
	gaplen1half = NULL;

	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen2half[i][j] ) FreeGaplenMtx( gaplen2half[i][j], 0 );
			gaplen2half[i][j] = NULL;
		}
		free( gaplen2half[i] );
		gaplen2half[i] = NULL;
	}
	free( gaplen2half );
	gaplen2half = NULL;
#endif
/* maikai free */


#if WMCHECK
	fprintf( stderr, "wm = %f\n", wm - *impmatch);
	fprintf( stderr, "*impmatch = %f\n", *impmatch);

	int kenzan = 0;
	for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
	{
		kenzan += pairgapcount( mseq1[i], mseq2[j] );
	}


	reporterr( "kenzan = %d -> %f\n", kenzan, (double)kenzan /( icyc*jcyc ) );

	double pairscore, nogappairscore;
	char **pseq;
	pseq = AllocateCharMtx( 2, strlen( seq1[0] ) + 1 );
	pairscore = nogappairscore = 0.0;
	for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
	{
		strcpy( pseq[0], seq1[i] );
		strcpy( pseq[1], seq2[j] );
		commongappick( 2, pseq );
		pairscore += eff1[i] * eff2[j] * naivepairscore11_dynmtx( n_dynamicmtx, pseq[0], pseq[1], penalty );
		nogappairscore += eff1[i] * eff2[j] * naivepairscore11_dynmtx( n_dynamicmtx, pseq[0], pseq[1], 0 );
	}

	FreeCharMtx( pseq );
	reporterr( "pairscore = %f\n", (double)pairscore );
	reporterr( "pairscore-nogappairscore = %f\n", (double)(pairscore-nogappairscore) );
	reporterr( "pairscore-nogappairscore / penalty = %f\n", (double)(pairscore-nogappairscore)/(double)(fpenalty) );
	reporterr( "diff = %f\n\n", (pairscore - wm + *impmatch ) / fpenalty );

#if 1
	if( ( !trywarp && fabs( pairscore - wm +*impmatch ) > 0.01 ) || PFACERROR ) // abs() -> fabs(), 2019/Jan/25
//	if( abs( pairscore - wm +*impmatch ) > 0.01 )
#else
	if( abs( pairscore - wm +*impmatch ) > 0.01 )
#endif
//	if( abs( pairscore - wm +*impmatch ) > 0.01 )
	{
		for( i=0; i<icyc; i++ )
			printf( ">group1\n%s\n", seq1[i] );
		for( j=0; j<jcyc; j++ )
			printf( ">group2\n%s\n", seq2[j] );
		exit( 1 );
	}
#else
	reporterr( "\n" );
#endif

#if 0
//	if( strlen( seq1[0] ) - lgth1 > 100 && icyc > 1 || strlen( seq2[0] ) - lgth2 > 100 & jcyc > 1 )
	if( strstr( seq1[0], "LNDDAT" ) && icyc == 1 || strstr( seq2[0], "LNDDAT" ) && jcyc==1)
	{
		for( i=0; i<icyc; i++ )
			printf( ">group1\n%s\n", seq1[i] );
		for( j=0; j<jcyc; j++ )
			printf( ">group2\n%s\n", seq2[j] );
		exit( 1 );
	}
#endif


	return( wm );
}


double D__align_ls( double **n_dynamicmtx, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, int constraint, double *impmatch, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
{
	int v1, v2;
	double val;

#if 1
	v1 = gapvariety( icyc, strlen( seq1[0] ), seq1 );
	v2 = gapvariety( jcyc, strlen( seq2[0] ), seq2 );
#else
	v1 = icyc;
	v2 = jcyc;
#endif

//	reporterr( "\nicyc,jcyc = %d,%d\n", icyc, jcyc );
	reporterr( " v1,v2 = %d,%d\n", v1, v2 );

	if( v1 >= v2 )
	{
		val = D__align( n_dynamicmtx, seq1, seq2, eff1, eff2, icyc, jcyc, alloclen, constraint, impmatch, sgap1, sgap2, egap1, egap2, chudanpt, chudanref, chudanres, headgp, tailgp );
	}
	else
	{
		val = D__align( n_dynamicmtx, seq2, seq1, eff2, eff1, jcyc, icyc, alloclen, constraint, impmatch, sgap2, sgap1, egap2, egap1, chudanpt, chudanref, chudanres, headgp, tailgp );
	}
	return val;
}



double D__align_gapmap( char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, int constraint, double *impmatch, int *gapmap1, int *gapmap2 )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
	fprintf( stderr, "Unexpected error.  Please contact katoh@ifrec.osaka-u.ac.jp\n" );
	exit( 1 );
}


double D__align_variousdist( int **which, double ***matrices, double **n_dynamicmtx, char **seq1, char **seq2, double *eff1, double *eff2, double **eff1s, double **eff2s, int icyc, int jcyc, int alloclen, int constraint, double *impmatch, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{

//	int k;
	register int i, j, c;
	int lasti, lastj;      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
//	double fpenalty = (double)penalty;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS int **ijp;
	static TLS int mpi, *mp;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS char **mseq1;
	static TLS char **mseq2;
	static TLS char **mseq;
	static TLS double ***cpmx1s;
	static TLS double ***cpmx2s;
	static TLS int ***intwork;
	static TLS double ***doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;
#if USEGAPLENHALF
	Gaplen ****gaplen1half = NULL; // NULL ga iru to omou.
	Gaplen ****gaplen2half = NULL; // NULL ga iru to omou.
#endif
#if USEGAPLENMTX
	Gaplen ****gaplen1mtx = NULL; // NULL ga iru to omou.
	Gaplen ****gaplen2mtx = NULL; // NULL ga iru to omou.
#endif
	static TLS Gaplen **gaplen1 = NULL; // NULL ga iru to omou.
	static TLS Gaplen **gaplen2 = NULL; // NULL ga iru to omou.
	static TLS Gaplen ***gaplen1jprev = NULL;
	static TLS Gaplen ***gaplen2jprev = NULL;
	static TLS Gaplen ***gaplen1jcurr = NULL;
	static TLS Gaplen ***gaplen2jcurr = NULL;
	static TLS Gaplen ***gaplen1icurr = NULL;
	static TLS Gaplen ***gaplen2icurr = NULL;
	static TLS Gaplen ***gaplen1jbestkamo = NULL;
	static TLS Gaplen ***gaplen2jbestkamo = NULL;
	static TLS Gaplen ***gaplen1ibestkamo = NULL;
	static TLS Gaplen ***gaplen2ibestkamo = NULL;
	static TLS Gaplen ***gaplen1jbest = NULL;
	static TLS Gaplen ***gaplen2jbest = NULL;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	static TLS Gaplen ****gaplens = NULL;

	Gaplen ***gaplentmp = NULL;
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
	int k;
	double pfac, pfactmp;
	int newgaplen;
	int **masklist1 = NULL, **masklist2 = NULL;
	int *nmask;

//	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s, %f\n", seq1[i], eff1[i] );
//	for( i=0; i<jcyc; i++ ) fprintf( stderr, "%s, %f\n", seq2[i], eff2[i] );

//	reporterr( "\nsizeof(gaplen) = %d\n", sizeof( Gaplen ) );
//	reporterr( "\nsizeof(int) = %d\n", sizeof( int ) );
//	reporterr( "\nsizeof(double) = %d\n", sizeof( double ) );
//	reporterr( "\nsizeof(double*) = %d\n", sizeof( double * ) );


	if( seq1 == NULL )
	{
		if( orlgth1 )
		{
//			fprintf( stderr, "## Freeing local arrays in D__align\n" );
			orlgth1 = 0;
			orlgth2 = 0;

			imp_match_init_strictD( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL, 0, 0 );

			free( mseq1 );
			free( mseq2 );
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatCub( cpmx1s );
			FreeFloatCub( cpmx2s );

			FreeFloatCub( doublework );
			FreeIntCub( intwork );



			free( gaplens );
			if( gaplen1ibestkamo ) FreeGaplenCub( gaplen1ibestkamo ); gaplen1ibestkamo = NULL;
			if( gaplen2ibestkamo ) FreeGaplenCub( gaplen2ibestkamo ); gaplen2ibestkamo = NULL;
			if( gaplen1icurr ) FreeGaplenCub( gaplen1icurr ); gaplen1icurr = NULL;
			if( gaplen2icurr ) FreeGaplenCub( gaplen2icurr ); gaplen2icurr = NULL;
   	                   
			if( gaplen1jprev ) FreeGaplenCub( gaplen1jprev ); gaplen1jprev = NULL;
			if( gaplen2jprev ) FreeGaplenCub( gaplen2jprev ); gaplen2jprev = NULL;
			if( gaplen1jcurr ) FreeGaplenCub( gaplen1jcurr ); gaplen1jcurr = NULL;
			if( gaplen2jcurr ) FreeGaplenCub( gaplen2jcurr ); gaplen2jcurr = NULL;
			if( gaplen1jbestkamo ) FreeGaplenCub( gaplen1jbestkamo ); gaplen1jbestkamo = NULL;
			if( gaplen2jbestkamo ) FreeGaplenCub( gaplen2jbestkamo ); gaplen2jbestkamo = NULL;
			if( gaplen1jbest ) FreeGaplenCub( gaplen1jbest ); gaplen1jbest = NULL;
			if( gaplen2jbest ) FreeGaplenCub( gaplen2jbest ); gaplen2jbest = NULL;
			if( gaplen1 ) FreeGaplenMtx( gaplen1, 1 ); gaplen1 = NULL;
			if( gaplen2 ) FreeGaplenMtx( gaplen2, 1 ); gaplen2 = NULL;

		}
		else
		{
//			fprintf( stderr, "## Not allocated\n" );
		}
		return( 0.0 );
	}


#if SLOW
	nmask = calloc( maxdistclass, sizeof( int ) );
#else
	masklist1 = AllocateIntMtx( maxdistclass, 0 );
	masklist2 = AllocateIntMtx( maxdistclass, 0 );
	nmask = calloc( maxdistclass, sizeof( int ) );

	for( c=0; c<maxdistclass; c++ )
	{
		for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
		{
			if( eff1s[c][i] * eff2s[c][j] != 0.0 )
			{
//				reporterr( "%d-%d, c=%d\n", i, j, c );
				if( c != which[i][j] )
				{
					masklist1[c] = realloc( masklist1[c], sizeof( int ) * nmask[c]+1 );
					masklist2[c] = realloc( masklist2[c], sizeof( int ) * nmask[c]+1 );

					masklist1[c][nmask[c]] = i;
					masklist2[c][nmask[c]] = j;
					nmask[c]++;
				}
			}
		}
	}
	for( c=0; c<maxdistclass; c++ ) if( nmask[c] ) break;
	if( c<maxdistclass ) reporterr( "Found a complex grouping. This step may be a bit slow.\n" );
#endif


	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );


	reporterr( "%d x %d sequences, len=%d, %d\n", icyc, jcyc, lgth1, lgth2 );


#if 0
	if( lgth1 == 0 || lgth2 == 0 )
	{
		fprintf( stderr, "WARNING (Aalignmm): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		for( i=0; i<icyc; i++ )
		{
			j = lgth2;
			seq1[i][j] = 0;
			while( j ) seq1[i][--j] = *newgapstr;
//			fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
		}
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		for( i=0; i<jcyc; i++ )
		{
			j = lgth1;
			seq2[i][j] = 0;
			while( j ) seq2[i][--j] = *newgapstr;
//			fprintf( stderr, "seq2[i] = %s\n", seq2[i] );
		}
		return( 0.0 );
	}

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;



	if( trywarp )
	{
//		reporterr( "Not supported yet!\n" );
//		exit( 1 );
//		fprintf( stderr, "IN D__align, penalty_shift = %d\n", penalty_shift );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1 to allow shift.\n" );
			exit( 1 );
		}
		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}


#if 0
	fprintf( stderr, "####  eff in SA+++align\n" );
	fprintf( stderr, "####  seq1[0] = %s\n", seq1[0] );
	fprintf( stderr, "####  strlen( seq1[0] ) = %d\n", strlen( seq1[0] ) );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
	fprintf( stderr, "####  seq2[0] = %s\n", seq2[0] );
	fprintf( stderr, "####  strlen( seq2[0] ) = %d\n", strlen( seq2[0] ) );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "eff2[%d] = %f\n", i, eff2[i] );
#endif
	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 0 );
		mseq2 = AllocateCharMtx( njob, 0 );
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

			FreeFloatCub( cpmx1s );
			FreeFloatCub( cpmx2s );

			FreeFloatCub( doublework );
			FreeIntCub( intwork );


			free( gaplens );

			if( gaplen1ibestkamo ) FreeGaplenCub( gaplen1ibestkamo ); gaplen1ibestkamo = NULL;
			if( gaplen2ibestkamo ) FreeGaplenCub( gaplen2ibestkamo ); gaplen2ibestkamo = NULL;
			if( gaplen1icurr ) FreeGaplenCub( gaplen1icurr ); gaplen1icurr = NULL;
			if( gaplen2icurr ) FreeGaplenCub( gaplen2icurr ); gaplen2icurr = NULL;
   	                   
			if( gaplen1jcurr ) FreeGaplenCub( gaplen1jcurr ); gaplen1jcurr = NULL;
			if( gaplen1jprev ) FreeGaplenCub( gaplen1jprev ); gaplen1jprev = NULL;
			if( gaplen2jcurr ) FreeGaplenCub( gaplen2jcurr ); gaplen2jcurr = NULL;
			if( gaplen2jprev ) FreeGaplenCub( gaplen2jprev ); gaplen2jprev = NULL;
			if( gaplen1jbestkamo ) FreeGaplenCub( gaplen1jbestkamo ); gaplen1jbestkamo = NULL;
			if( gaplen2jbestkamo ) FreeGaplenCub( gaplen2jbestkamo ); gaplen2jbestkamo = NULL;
			if( gaplen1jbest ) FreeGaplenCub( gaplen1jbest ); gaplen1jbest = NULL;
			if( gaplen2jbest ) FreeGaplenCub( gaplen2jbest ); gaplen2jbest = NULL;
			if( gaplen1 ) FreeGaplenMtx( gaplen1, 1 ); gaplen1 = NULL;
			if( gaplen2 ) FreeGaplenMtx( gaplen2, 1 ); gaplen2 = NULL;

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

		mseq = AllocateCharMtx( njob, ll1+ll2 );

		cpmx1s = AllocateFloatCub( maxdistclass, nalphabets, ll1+2 );
		cpmx2s = AllocateFloatCub( maxdistclass, nalphabets, ll2+2 );

		doublework = AllocateFloatCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets ); 
		intwork = AllocateIntCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets+1 ); 

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;


//		reporterr( "Allocating gaplen1 and gaplen2\n" );
		gaplen1 = (Gaplen ** )calloc( ll1+2, sizeof( Gaplen * ) );
		gaplen1[ll1+1] = (Gaplen *)1;
		gaplen2 = (Gaplen ** )calloc( ll2+2, sizeof( Gaplen * ) );
		gaplen2[ll2+1] = (Gaplen *)1;


//		reporterr( "Allocating gaplen*\n" );
		gaplen1ibestkamo = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen2ibestkamo = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen1icurr = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen2icurr = (Gaplen ***)calloc( (ll1+2), sizeof( Gaplen **) );
		gaplen1jbestkamo = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jbestkamo = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen1jbest = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jbest = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen1jcurr = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jcurr = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen1jprev = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );
		gaplen2jprev = (Gaplen ***)calloc( (ll2+2), sizeof( Gaplen **) );

		gaplens = calloc( sizeof( Gaplen ***), 12 );
		gaplens[0] = gaplen1ibestkamo;
		gaplens[1] = gaplen2ibestkamo;
		gaplens[2] = gaplen1icurr;
		gaplens[3] = gaplen2icurr;
		gaplens[4] = gaplen1jbestkamo;
		gaplens[5] = gaplen2jbestkamo;
		gaplens[6] = gaplen1jbest;
		gaplens[7] = gaplen2jbest;
		gaplens[8] = gaplen1jcurr;
		gaplens[9] = gaplen2jcurr;
		gaplens[10] = gaplen1jprev;
		gaplens[11] = gaplen2jprev;
//		reporterr( "Allocation end\n" );
	}

	{
		int ll1 = lgth1;
		int ll2 = lgth2;

//		reporterr( "Allocating gaplen*i\n" );
		for(i=0; i<ll1+1; i++ ) 
		{
			gaplen1ibestkamo[i] = (Gaplen **)calloc( ll1+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1-i; k++ ) gaplen1ibestkamo[i][k] = NULL;
			gaplen1ibestkamo[i][k] = (Gaplen *)1;

			gaplen2ibestkamo[i] = (Gaplen **)calloc( ll2+2, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1; k++ ) gaplen2ibestkamo[i][k] = NULL;
			gaplen2ibestkamo[i][k] = (Gaplen *)1;

			gaplen1icurr[i] = (Gaplen **)calloc( ll1+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1-i; k++ ) gaplen1icurr[i][k] = NULL;
			gaplen1icurr[i][k] = (Gaplen *)1;

			gaplen2icurr[i] = (Gaplen **)calloc( ll2+2, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1; k++ ) gaplen2icurr[i][k] = NULL;
			gaplen2icurr[i][k] = (Gaplen *)1;
		}
		gaplen1ibestkamo[ll1+1] = NULL;
		gaplen2ibestkamo[ll1+1] = NULL;
		gaplen1icurr[ll1+1] = NULL;
		gaplen2icurr[ll1+1] = NULL;
	
//		reporterr( "Allocating gaplen*j\n" );
		for(i=0; i<ll2+1; i++ ) 
		{
			gaplen1jbestkamo[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jbestkamo[i][k] = NULL;
			gaplen1jbestkamo[i][k] = (Gaplen *)1;

			gaplen2jbestkamo[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jbestkamo[i][k] = NULL;
			gaplen2jbestkamo[i][k] = (Gaplen *)1;

			gaplen1jbest[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jbest[i][k] = NULL;
			gaplen1jbest[i][k] = (Gaplen *)1;

			gaplen2jbest[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jbest[i][k] = NULL;
			gaplen2jbest[i][k] = (Gaplen *)1;

			gaplen1jcurr[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jcurr[i][k] = NULL;
			gaplen1jcurr[i][k] = (Gaplen *)1;

			gaplen2jcurr[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jcurr[i][k] = NULL;
			gaplen2jcurr[i][k] = (Gaplen *)1;

			gaplen1jprev[i] = (Gaplen **)calloc( ll1+2, sizeof( Gaplen * ) );
			for( k=0; k<ll1+1; k++ ) gaplen1jprev[i][k] = NULL;
			gaplen1jprev[i][k] = (Gaplen *)1;

			gaplen2jprev[i] = (Gaplen **)calloc( ll2+2-i, sizeof( Gaplen * ) );
			for( k=0; k<ll2+1-i; k++ ) gaplen2jprev[i][k] = NULL;
			gaplen2jprev[i][k] = (Gaplen *)1;

		}
		gaplen1jbestkamo[ll2+1] = NULL;
		gaplen2jbestkamo[ll2+1] = NULL;
		gaplen1jbest[ll2+1] = NULL;
		gaplen2jbest[ll2+1] = NULL;
		gaplen1jcurr[ll2+1] = NULL;
		gaplen2jcurr[ll2+1] = NULL;
		gaplen1jprev[ll2+1] = NULL;
		gaplen2jprev[ll2+1] = NULL;
	}


#if USEGAPLENMTX
/* maikai allocate */

	reporterr( "Allocating gaplenmtx1\n" );
	gaplen1mtx = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen1mtx[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen1mtx[i][j] = (Gaplen **)calloc( lgth1+2, sizeof( Gaplen * ) );
			for( k=0; k<lgth1+1; k++ ) gaplen1mtx[i][j][k] = NULL;
			gaplen1mtx[i][j][k] = (Gaplen *)1;
		}
		gaplen1mtx[i][j] = NULL;
	}
	gaplen1mtx[i] = NULL;

	reporterr( "Allocating gaplenmtx2\n" );
	gaplen2mtx = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen2mtx[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen2mtx[i][j] = (Gaplen **)calloc( lgth2+2, sizeof( Gaplen * ) );
			for( k=0; k<lgth2+1; k++ ) gaplen2mtx[i][j][k] = NULL;
			gaplen2mtx[i][j][k] = (Gaplen *)1;
		}
		gaplen2mtx[i][j] = NULL;
	}
	gaplen2mtx[i] = NULL;

#endif

#if USEGAPLENHALF
	reporterr( "Allocating gaplenhalf1\n" );
	gaplen1half = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen1half[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen1half[i][j] = (Gaplen **)calloc( lgth1+2 - i, sizeof( Gaplen * ) );
			for( k=0; k<lgth1+1-i; k++ ) gaplen1half[i][j][k] = NULL;
			gaplen1half[i][j][k] = (Gaplen *)1;
		}
		gaplen1half[i][j] = NULL;
	}
	gaplen1half[i] = NULL;

	reporterr( "Allocating gaplenhalf2\n" );
	gaplen2half = (Gaplen ****)calloc( (lgth1+2), sizeof( Gaplen ***) );
	for(i=0; i<lgth1+1; i++ ) gaplen2half[i] = (Gaplen ***)calloc( lgth2+2, sizeof( Gaplen ** ) );
	for(i=0; i<lgth1+1; i++ ) 
	{
		for(j=0; j<lgth2+1; j++ ) 
		{
			gaplen2half[i][j] = (Gaplen **)calloc( lgth2+2 - j, sizeof( Gaplen * ) );
			for( k=0; k<lgth2+1-j; k++ ) gaplen2half[i][j][k] = NULL;
			gaplen2half[i][j][k] = (Gaplen *)1;
		}
		gaplen2half[i][j] = NULL;
	}
	gaplen2half[i] = NULL;
#endif


/* maikai allocate */


	for( i=0; i<icyc; i++ )
	{
		mseq1[i] = mseq[i];
		seq1[i][lgth1] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		mseq2[j] = mseq[icyc+j];
		seq2[j][lgth2] = 0;
	}


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
	{
		double t = 0.0;
		for( i=0; i<icyc; i++ )
			t += eff1[i];
	fprintf( stderr, "## totaleff = %f\n", t );
	}
#endif


#if SLOW
#else
//	cpmx_calc_new( seq1, cpmx1, eff1, lgth1, icyc );
//	cpmx_calc_new( seq2, cpmx2, eff2, lgth2, jcyc );
	for( c=0; c<maxdistclass; c++ )
	{
		cpmx_calc_new( seq1, cpmx1s[c], eff1s[c], lgth1, icyc );
		cpmx_calc_new( seq2, cpmx2s[c], eff2s[c], lgth2, jcyc );
	}
#endif


//	reporterr( "Counting gaplen\n" );
	gaplencount( icyc, lgth1, gaplen1, seq1, eff1 );
	gaplencount( jcyc, lgth2, gaplen2, seq2, eff2 );
#if DEBUG
	reporterr( "group1 = \n" );
	showgaplen( gaplen1, lgth1 );
	reporterr( "group2 = \n" );
	showgaplen( gaplen2, lgth2 );
#endif
//	reporterr( "done.\n" );


	for( i=0; i<lgth1+1; i++ ) for( j=0; j<lgth2+1; j++ )
	{
#if USEGAPLENMTX
//		duplicategaplen( gaplen1mtx[i][j], gaplen1, lgth1 );
//		duplicategaplen( gaplen2mtx[i][j], gaplen2, lgth2 );

//		duplicategaplenpartly( gaplen2mtx[i][j], gaplen2, j-0, lgth2 ); // anzen
//		duplicategaplenpartly( gaplen1mtx[i][j], gaplen1, i-0, lgth1 ); // anzen
		duplicategaplenpartly( gaplen1mtx[i][j], gaplen1, i-0, i ); // iranaikamo
		duplicategaplenpartly( gaplen2mtx[i][j], gaplen2, j-0, j ); // iranaikamo
#endif
#if USEGAPLENHALF

//		duplicategaplenpartly( gaplen1half[i][j], gaplen1+i, 0, lgth1-i ); // KOKO de setsuyaku dekiru to omou.
//		duplicategaplenpartly( gaplen2half[i][j], gaplen2+j, 0, lgth2-j ); //  originally, j-1, lgth2
		duplicategaplenpartly( gaplen1half[i][j], gaplen1+i, 0, 0 ); // test
		duplicategaplenpartly( gaplen2half[i][j], gaplen2+j, 0, 0 ); // test
#endif
	}


#if USEGAPLENMTX
	reporterr( "Duplicating gaplen*mtx\n" );
	for( i=0; i<lgth1+1; i++ )
	{
//		addnewgaplen( gaplen1mtx[i][0], gaplen1, gaplen1, lgth1, -1, 0 );
//		addnewgaplen( gaplen2mtx[i][0], gaplen2, gaplen2, lgth2, 0, i );
//		duplicategaplenpartly( gaplen1mtx[i][0], gaplen1, 0, lgth1 );
//		duplicategaplenpartly( gaplen2mtx[i][0], gaplen2, 0, lgth2 );
		copygaplenrestricted( gaplen2mtx[i][0], gaplen2, lgth2, 0, i, 0, 0 );
	}
#endif
#if USEGAPLENHALF
	reporterr( "Duplicating gaplen*mtx\n" );
	for( i=0; i<lgth1+1; i++ )
	{
		copygaplenrestricted( gaplen2half[i][0], gaplen2, lgth2, 0, i, 0, 0 );
	}
#endif



	for( i=0; i<1; i++ )
	{
//		duplicategaplencompactx( gaplen1icurr[i], gaplen1, lgth1, i-0, lgth1 ); //originally, 0, lgth1
//
//		duplicategaplencompactx( gaplen1icurr[i], gaplen1+i, lgth1-i, 0, lgth1-i ); // half
		duplicategaplencompactx( gaplen1icurr[i], gaplen1+i, lgth1-i, 0, 1 ); //  0, 1  hitsuyou


//		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, lgth2 ); // ichiou zenbu
		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, 0 );


		copygaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, i, 0, 0 ); // -> zurasu -> error?


//		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1, lgth1, 0, 1 );
//		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1+i, lgth1-i, 0, 1 ); // half
		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1+i, lgth1-i, 0, 0 ); // half
//		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, lgth2 );
		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, 0 );
//		copygaplenrestricted( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, 0 ); // -> zurasu -> error?
//		copygaplenrestricted_zurasu( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, lgth2, 0, lgth2 ); // -> zurasu -> error?
		copygaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, 0 ); // -> zurasu -> error?
	}

//	reporterr( "Duplicating gaplen*j*curr \n" );
//	int nduplicated = 0;
	for( j=0; j<lgth2+1; j++ )
//	for( j=0; j<1; j++ )
	{
#if USEGAPLENMTX
//		addnewgaplen( gaplen1mtx[0][j], gaplen1, gaplen1, lgth1, 0, j );
//		addnewgaplen( gaplen2mtx[0][j], gaplen2, gaplen2, lgth2, -1, 0 );
//		duplicategaplenpartly( gaplen1mtx[0][j], gaplen1, 0, lgth1 );
//		duplicategaplenpartly( gaplen2mtx[0][j], gaplen2, 0, lgth2 );
		copygaplenrestricted( gaplen1mtx[0][j], gaplen1, lgth1, 0, j, 0, 0 );
#endif

#if USEGAPLENHALF
		copygaplenrestricted( gaplen1half[0][j], gaplen1, lgth1, 0, j, 0, 0 );
#endif
//		reporterr( "1jcurr?\n" );
//		duplicategaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, lgth1 ); // test
		duplicategaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, 0 ); // dame?
//		reporterr( "done\n" );
//		duplicategaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, 0 ); //test

//		duplicategaplencompactx( gaplen2jcurr[j], gaplen2, lgth2, j-0, lgth2 ); // full
//		duplicategaplencompactx( gaplen2jcurr[j], gaplen2+j, lgth2-j, 0, lgth2-j ); //half! KOKO????
//reporterr( "starting suspicious duplication\n" );
		duplicategaplencompactx( gaplen2jcurr[j], gaplen2+j, lgth2-j, 0, 0 ); //half!
//reporterr( "starting suspicious copy\n" );
		copygaplencompactx( gaplen1jcurr[j], gaplen1, lgth1, 0, j, 0, 0 ); // TEST
//reporterr( "finished\n" );

//		reporterr( "Checking gaplen1jcurr[%d]\n", j );
//		checkgaplen( gaplen1jcurr[j], 100 );
//		reporterr( "Checking gaplen2jcurr[%d]\n", j );
//		checkgaplen( gaplen2jcurr[j], 100 );
	}

//	reporterr( "nduplicated (corrected) = %d\n", nduplicated );

//	reporterr( "Duplicating gaplen*j*prev \n\n" );
	for( j=0; j<lgth2+1; j++ ) // allocate nominotame, atode uwagaki
	{
//		duplicategaplencompactx( gaplen1jprev[j], gaplen1, lgth1, 0, lgth1 );
		duplicategaplencompactx( gaplen1jprev[j], gaplen1, lgth1, 0, 0 ); // TEST
//		duplicategaplencompactx( gaplen2jprev[j], gaplen2, lgth2, j-0, lgth2 ); // originally, 0,lgth2
//		duplicategaplencompactx( gaplen2jprev[j], gaplen2+j, lgth2-j, 0, lgth2-j ); // half
		duplicategaplencompactx( gaplen2jprev[j], gaplen2+j, lgth2-j, 0, 0 ); // half


		copygaplencompactx( gaplen1jprev[j], gaplen1, lgth1, 0, j, 0, 0 ); // wasuretetakamo


	}


//	reporterr( "Duplicating gaplen*j*best \n\n" );

	for( j=0; j<lgth2+1; j++ )
//	for( j=0; j<1; j++ )
	{
//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, lgth1 ); // KOKO
//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 0 ); // test
//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 1 );
		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 0 );


//		duplicategaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, 1 );
//		duplicategaplencompactx( gaplen2jbestkamo[j], gaplen2, lgth2, j-0, j+1 ); // originally, 0, j+1
		duplicategaplencompactx( gaplen2jbestkamo[j], gaplen2+j, lgth2-j, 0, 1 ); // half!
		copygaplencompactx( gaplen1jbestkamo[j], gaplen1, lgth1, 0, j, 0, 0 ); // TEST

//		duplicategaplencompactx( gaplen1jbest[j], gaplen1, lgth1, 0, lgth1 ); // KOKO
		duplicategaplencompactx( gaplen1jbest[j], gaplen1, lgth1, 0, 0 ); // test
//		duplicategaplencompactx( gaplen2jbest[j], gaplen2, lgth2,j-0, j+1 ); // originally, 0,j+1
		duplicategaplencompactx( gaplen2jbest[j], gaplen2+j, lgth2-j, 0, 1 ); // half!
		copygaplencompactx( gaplen1jbest[j], gaplen1, lgth1, 0, j, 0, 0 ); // TEST


	}

//	reporterr( "Duplication end\n" );


#if 0
	reporterr( "Checking gaplen1icurr\n" );
	checkgaplen( gaplen1icurr[0], 100 );
	reporterr( "Checking gaplen2icurr\n" );
	checkgaplen( gaplen2icurr[0], 100 );
#endif




//	showgaplen( gaplen1jcurr[50], lgth2 );

	currentw = w1;
	previousw = w2;

//	match_calc( n_dynamicmtx, initverticalw, cpmx2, cpmx1, 0, lgth1, doublework, intwork, 1 );
#if SLOW
	match_calc_slow( which, matrices, initverticalw, jcyc, seq2, eff2, icyc, seq1, eff1, 0, lgth1, *doublework, *intwork, 1, 1 );
//	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%d - %f\n", i, initverticalw[i] );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
	{
//		fprintf( stderr, "c=%d matrices[c][W][W] = %f\n", c, matrices[c][amino_n['W']][amino_n['W']] );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "seq1[i] = %c, cpmx1s[c][3][%d] = %f\n", seq1[0][i], i, cpmx1s[c][3][i] );
//		for( i=0; i<lgth2; i++ ) fprintf( stderr, "seq2[i] = %c, cpmx2s[c][3][%d] = %f\n", seq2[0][i], i, cpmx2s[c][3][i] );
		match_calc_add( matrices[c], initverticalw, cpmx2s[c], cpmx1s[c], 0, lgth1, doublework[c], intwork[c], 1 );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "c=%d, %d - %f\n", c, i, initverticalw[i] );

		if( nmask[c] ) match_calc_del( which, matrices, initverticalw, jcyc, seq2, eff2, icyc, seq1, eff1, 0, lgth1, c, nmask[c], masklist2[c], masklist1[c] );
	}
#endif
//	reporterr( "initverticalw = \n" );
//	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%d - %f\n", i, initverticalw[i] );


	if( constraint )
		imp_match_out_vead_tate( initverticalw, 0, lgth1 ); // 060306

//	match_calc( n_dynamicmtx, currentw, cpmx1, cpmx2, 0, lgth2, doublework, intwork, 1 );
#if SLOW
	match_calc_slow( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, 0, lgth2, *doublework, *intwork, 1, 0 );
//	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%d - %f\n", i, currentw[i] );
//	exit( 1 );
#else
	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
	{
		match_calc_add( matrices[c], currentw, cpmx1s[c], cpmx2s[c], 0, lgth2, doublework[c], intwork[c], 1 );
		if( nmask[c] ) match_calc_del( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, 0, lgth2, c, nmask[c], masklist1[c], masklist2[c] );
	}
#endif
//	reporterr( "currentw = \n" );
//	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%d - %f\n", i, currentw[i] );





//	exit( 1 );


	if( constraint )
		imp_match_out_vead( currentw, 0, lgth2 ); // 060306
#if 0 // -> tbfast.c // impossible
	if( localhom )
		imp_match_calc( n_dynamicmtx, currentw, icyc, jcyc, lgth1, lgth2, seq1, seq2, eff1, eff2, localhom, 1, 0 );

#endif

	for( j=1; j<lgth2+1; j++ )
	{
		pfac = calcpfac_gap_noidatend( gaplen1, gaplen2, j, 0, j, seq1[0], seq2[0], 0 ); 
//		reporterr( "computing initial end gap penalty for %c-%c, i=0, j=%d, pfac=%f\n", seq1[0][0], seq2[0][j], j, pfac );
//		reporterr( "%c-%c, i=0, j=%d, currentw[j]=%f -> ", seq1[0][0], seq2[0][j], j, currentw[j] );
		currentw[j] += fpenalty * pfac; // tekitou
//		reporterr( " %f\n", currentw[j] );
	}
	for( i=1; i<lgth1+1; i++ )
	{
		pfac = calcpfac_gap_noidatend( gaplen2, gaplen1, i, 0, i, seq2[0], seq1[0], 0 );
//		reporterr( "computing initial end gap penalty for %c-%c, i=%d, j=0, pfac=%f\n", seq1[0][i], seq2[0][0], i, pfac );
		initverticalw[i] += fpenalty * pfac; // tekitou
	}



	for( j=1; j<lgth2+1; ++j ) 
	{


#if ALGZGAP
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2[j-1]; mp[j] = 0;;
#else
		pfac = calcpfac_gapex_noidatend( gaplen2, gaplen1, j, 1, j, seq2[0], seq1[0], 1 );
#if DEBUG
		reporterr( "%c-%c, INITIAL jgap extension check, pfac = %f\n\n", seq1[0][j], '=', pfac );
#endif
		m[j] = currentw[j-1] + fpenalty * pfac; 
		mp[j] = 0;
#endif
	}
	if( lgth2 == 0 )
		lastverticalw[0] = 0.0; // Falign kara yobaretatoki kounarukanousei ari
	else
		lastverticalw[0] = currentw[lgth2-1];

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;


	for( i=1; i<lasti; i++ )
	{
//		reporterr( "i = %d\n", i );

//		reporterr( "err1? i=%d/%d\n", i, lgth1 );
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
			cleargaplens( gaplens );
			if( masklist1 ) FreeIntMtx( masklist1 ); masklist1 = NULL;
			if( masklist2 ) FreeIntMtx( masklist2 ); masklist2 = NULL;
			if( nmask ) free( nmask ); nmask = NULL;
//			fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
			*chudanres = 1;
			return( -1.0 );
		}
#endif


		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if 1
		gaplentmp = gaplen1jprev;
		gaplen1jprev = gaplen1jcurr;
		gaplen1jcurr = gaplentmp;

		gaplentmp = gaplen2jprev;
		gaplen2jprev = gaplen2jcurr;
		gaplen2jcurr = gaplentmp;

#if DEBUG
		reporterr( "Entering a small j loop, i=%d\n", i );
		for( j=1; j<lgth2+1; j++ )
		{
			reporterr( "before j loop, i=%d, gaplen2jcurr[%d] = \n", i, j );
			showgaplen( gaplen2jcurr[j], 100 );
			reporterr( "\n" );
			reporterr( "before j loop, i=%d, gaplen2prev[%d] = \n", i, j );
			showgaplen( gaplen2jprev[j], 100 );
			reporterr( "\n" );
		}
#endif
#else
		
		reporterr( "Entering a small j loop, i=%\n", i );
		for( j=1; j<lgth2+1; j++ )
		{
//			addnewgaplen( gaplen1jprev[j], gaplen1jcurr[j], gaplen1, lgth1, -1, 0 );
//			addnewgaplen( gaplen2jprev[j], gaplen2jcurr[j], gaplen2, lgth2, -1, 0 );
			reporterr( "err1? j=%d/%d\n", j, lgth2 );
			copygaplencompactx( gaplen1jprev[j-1], gaplen1jcurr[j-1], lgth1, -1, 0, i-1, i-1 ); // TEST
			reporterr( "err1? j=%d/%d\n", j, lgth2 );
			copygaplencompactx( gaplen2jprev[j-1], gaplen2jcurr[j-1], lgth2, -1, 0, j-1, j-1 ); // TETS
#if DEBUG
			reporterr( "before j loop, i=%d, gaplen2jcurr[%d] = \n", i, j );
			showgaplen( gaplen2jcurr[j], 100 );
			reporterr( "\n" );
			reporterr( "before j loop, i=%d, gaplen2prev[%d] = \n", i, j );
			showgaplen( gaplen2jprev[j], 100 );
			reporterr( "\n" );
#endif
		}
#endif

//		reporterr( "err2? i=%d/%d\n", i, lgth1 );

//		duplicategaplencompactx( gaplen1icurr[i], gaplen1, lgth1, i, i+1 ); //originally 0, i+1
//		reporterr( "gaplen+0=\n");
//		showgaplen( gaplen1, 10 );
//		reporterr( "i=%d, lgth1=%d, lgth1-i=%d, gaplen+i-1=\n", i, lgth1, lgth1-i );
//		showgaplen( gaplen1+i-1, 100 );
		duplicategaplencompactx( gaplen1icurr[i], gaplen1+i, lgth1-i, 0, 1 ); // half!!
//		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, lgth2 ); // KOKO
		duplicategaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, 0 ); // test
		copygaplencompactx( gaplen2icurr[i], gaplen2, lgth2, 0, i, 0, 0 ); // IRU? TEST



//		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1, lgth1, i, i+1 ); //originally 0, i+1
		duplicategaplencompactx( gaplen1ibestkamo[i], gaplen1+i, lgth1-i, 0, 1 ); //half
//		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, lgth2 ); // ORIGINALLY, 0, lgth2
		duplicategaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, 0 ); // ORIGINALLY, 0, lgth2
//		copygaplenrestricted( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, lgth2, 0, 0 ); // IRU? // TEST
		copygaplencompactx( gaplen2ibestkamo[i], gaplen2, lgth2, 0, i, 0, 0 ); // IRU? // TEST

		extendgaplencompactx( gaplen1jprev[0], gaplen1, i ); // ???


//		addnewgaplen( gaplen1jprev[0], gaplen1icurr[i-1], gaplen1, lgth1, -1, 0 );
//		addnewgaplen( gaplen2jprev[0], gaplen2icurr[i-1], gaplen2, lgth2, -1, 0 );
//		copygaplenrestricted( gaplen1jprev[0], gaplen1icurr[i-1], lgth1, -1, 0, i, i ); // i-1, i da to omou.
		copygaplencompactx( gaplen1jprev[0], gaplen1icurr[i-1], lgth1-i, -1, 0, i, 1 ); // half? lgth1-i?
//		copygaplenrestricted( gaplen2jprev[0], gaplen2icurr[i-1], lgth2, -1, 0, 0, 0 );
		copygaplencompactx( gaplen2jprev[0], gaplen2icurr[i-1], lgth2-j, -1, 0, 0, 0 ); // half?? lgth2-j?


//		match_calc( n_dynamicmtx, currentw, cpmx1, cpmx2, i, lgth2, doublework, intwork, 0 );
#if SLOW
		match_calc_slow( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, i, lgth2, *doublework, *intwork, 0, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
		{
			match_calc_add( matrices[c], currentw, cpmx1s[c], cpmx2s[c], i, lgth2, doublework[c], intwork[c], 0 );
			if( nmask[c] ) match_calc_del( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, i, lgth2, c, nmask[c], masklist1[c], masklist2[c] );
		}
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
		if( constraint )
		{
//			fprintf( stderr, "Calling imp_match_calc (o) lgth = %d, i = %d\n", lgth1, i );
#if  0
			imp_match_out_vead( currentw, i, lgth2 );
#else
			imp_match_out_vead( currentw, i, lgth2 );
#endif
		}
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

#if 0
		fprintf( stderr, "%c ", seq1[0][i] );
		for( j=0; j<lgth2+1; j++ )
		{
			fprintf( stderr, "%5.0f ", currentw[j] );
		}
		fprintf( stderr, "\n"  );
#endif
	
//		mi = previousw[0] + ogcp2[1]; mpi = 0;





#if ALGZGAP
		mi = previousw[0] + ogcp2[1] * gapfreq1[i-1]; mpi=0;
#else
		pfac = calcpfac_gapex_noidatend( gaplen1, gaplen2, i, 1, i, seq1[0], seq2[0], 1 );
#if DEBUG
		reporterr( "%c-%c, INITIAL igap extension check, pfac = %f\n\n", '=', seq2[0][j], pfac );
#endif
		mi = previousw[0] + fpenalty * pfac; 
		mpi=0;
#endif
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




//		reporterr( "\n\ni=%d, %c\n", i, seq1[0][i] );
		for( j=1; j<lastj; j++ )
		{

#if DEBUG
			reporterr( "*****   %c%d-%c%d   ******* \n", seq1[0][i], i, seq2[0][j], j );
			reporterr( "mpi=%d (%c), *mpjpt=%d (%c)\n", mpi, seq2[0][mpi], *mpjpt, seq1[0][*mpjpt] );
#endif


// Hitsuyou na bun dake tsuika
#if USEGAPLENMTX
			extendgaplenpartly( gaplen1mtx[i-1][mpi], gaplen1, i, i );
			extendgaplenpartly( gaplen2mtx[i-1][mpi], gaplen2, j, j );
			extendgaplenpartly( gaplen1mtx[*mpjpt][j-1], gaplen1, i, i );
			extendgaplenpartly( gaplen2mtx[*mpjpt][j-1], gaplen2, j, j );
			extendgaplenpartly( gaplen1mtx[i-1][j-1], gaplen1, i, i );
			extendgaplenpartly( gaplen2mtx[i-1][j-1], gaplen2, j, j );
#endif
#if USEGAPLENHALF
			extendgaplenpartly( gaplen1half[i-1][mpi], gaplen1+i-1, 1, 1 );
			extendgaplenpartly( gaplen2half[i-1][mpi], gaplen2+mpi, j-mpi, j-mpi );
			extendgaplenpartly( gaplen1half[*mpjpt][j-1], gaplen1+*mpjpt, i-*mpjpt, i-*mpjpt );
			extendgaplenpartly( gaplen2half[*mpjpt][j-1], gaplen2+j-1, 1, 1 );
			extendgaplenpartly( gaplen1half[i-1][j-1], gaplen1+i-1, 1, 1 );
			extendgaplenpartly( gaplen2half[i-1][j-1], gaplen2+j-1, 1, 1 );

#endif

//			reporterr( "extending gaplen1icurr\n" );
			extendgaplencompactx( gaplen1icurr[i-1], gaplen1+i-1, 1 ); // iruhazu
//			reporterr( "extending gaplen2icurr\n" );
			extendgaplencompactx( gaplen2icurr[i-1], gaplen2, j ); // iruhazu
//			reporterr( "extending gaplen1jprev[j-1], j-1=%d\n", j-1 );
			extendgaplencompactx( gaplen1jprev[j-1], gaplen1, i );
//			reporterr( "extending gaplen1jcurr, j-1=%d\n", j-1 );
			extendgaplencompactx( gaplen1jcurr[j-1], gaplen1, i );
//			reporterr( "extending gaplen2jprev\n" );
			extendgaplencompactx( gaplen2jprev[j-1], gaplen2+j-1, 1 );
//			reporterr( "extending gaplen2jcurr\n" );
			extendgaplencompactx( gaplen2jcurr[j-1], gaplen2+j-1, 1 );
//			reporterr( "extending gaplen1jbest[j-1]\n" );
			extendgaplencompactx( gaplen1jbest[j-1], gaplen1, i );
//			reporterr( "extending gaplen1jbestkamo[j-1]\n" );
			extendgaplencompactx( gaplen1jbestkamo[j-1], gaplen1, i );
//			reporterr( "extending gaplen1jprev[mpi], j-1=%d\n", j-1 );
			extendgaplencompactx( gaplen1jprev[mpi], gaplen1, i ); // full
//			reporterr( "extending gaplen2jprev[mpi]\n" );
			extendgaplencompactx( gaplen2jprev[mpi], gaplen2+mpi, j-mpi ); // half
//			reporterr( "extending gaplen1ibestkamo[i-1]\n" );
			extendgaplencompactx( gaplen1ibestkamo[i-1], gaplen1+i-1, 1 );
//			reporterr( "extending gaplen2ibestkamo[i-1]\n" );
			extendgaplencompactx( gaplen2ibestkamo[i-1], gaplen2, j );


#if DEBUG
			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jcurr[j]=\n", i, j );
			showgaplen( gaplen2jcurr[j], 100 );

			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jcurr[j-1]=\n", i, j );
			showgaplen( gaplen2jcurr[j-1], 100 );


			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jprev[j]=\n", i, j );
			showgaplen( gaplen2jprev[j], 100 );

			reporterr( "at the head of j loop, i,j=%d,%d, gaplen2jprev[j-1]=\n", i, j );
			showgaplen( gaplen2jprev[j-1], 100 );
#endif


#ifdef xxxenablemultithread
//			fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
			if( chudanpt && *chudanpt != chudanref ) 
			{
				cleargaplens( gaplens );
				if( masklist1 ) FreeIntMtx( masklist1 ); masklist1 = NULL;
				if( masklist2 ) FreeIntMtx( masklist2 ); masklist2 = NULL;
				if( nmask ) free( nmask ); nmask = NULL;
//				fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
				*chudanres = 1;
				return( -1.0 );
			}
#endif
#if USEGAPLENHALF
// i==248, j==80 wo check
#if DEBUG80
			if( j==80 )
			{
				reporterr( "When i==%d, j==%d,\n", i, j );
				reporterr( "gaplen2jprev[j-1=%d]=\n", j-1 );
				showgaplen( gaplen2jprev[j-1], lgth1 );
				reporterr( "gaplen2half[i-1=%d][j-1=%d]=\n", i-1, j-1 );
				showgaplen( gaplen2half[i-1][j-1], lgth1 );
			}
			if( j==79 )
			{
				reporterr( "When i==%d, j==%d,\n", i, j );
				reporterr( "gaplen2jprev[j-1=%d]=\n", j-1 );
				showgaplen( gaplen2jprev[j-1], lgth1 );
				reporterr( "gaplen2half[i-1=%d][j-1=%d]=\n", i-1, j-1 );
				showgaplen( gaplen2half[i-1][j-1], lgth1 );
			}
#endif
#endif


//			pfac = calcpfac( gaplen1jprev[j-1], gaplen2jprev[j-1], i, j, seq1[0], seq2[0] );
//reporterr( "#### COMPACT, i,j=%d,%d\n", i, j );
			pfac = calcpfacnoidatend( gaplen1jprev[j-1], gaplen2jprev[j-1], i, 1, seq1[0], seq2[0]+j, one ); // 1j->full, 2j->half
#if USEGAPLENMTX
//reporterr( "#### FULL, i,j=%d,%d\n", i, j );
			pfactmp = calcpfac( gaplen1mtx[i-1][j-1], gaplen2mtx[i-1][j-1], i, j, seq1[0], seq2[0], one );
#endif
#if USEGAPLENHALF
//reporterr( "#### HALF, i,j=%d/%d,%d/%d\n", i, lgth1, j, lgth2 );
//			showgaplen( gaplen2half[i-1][j-1], lgth2-j );
			pfactmp = calcpfac( gaplen1half[i-1][j-1], gaplen2half[i-1][j-1], 1, 1, seq1[0]+i, seq2[0]+j, zero );
#endif
#if USEGAPLENMTX + USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(straight) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
				PFACERROR = 1;
				exit( 1 );
			}
#endif
//if( i==50 && j==135 ) exit( 1 );


//			reporterr( "i,j=%d,%d, *prept = %f\n", i, j, *prept );

#if ALGZSTRAIGHT
			wm = *prept;  // Machigai!!
#else
			wm = *prept + fpenalty * pfac;
#endif
			*ijppt = 0;


#if DEBUG
			if( i == j )
			{
				fprintf( stderr, "\n i=%d, j=%d %c, %c ", i, j, seq1[0][i], seq2[0][j] );
				fprintf( stderr, "%5.0f, pfac for straight =%f\n", wm, pfac );
			}
#endif
			newgaplen = j-mpi-1;


//if( i == 53 && j == 93 ) exit( 1 );




//			pfac = calcpfac_gap_incomplete( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], newgaplen, i, j, seq1[0], seq2[0], 0 ); // i-1
			pfac = calcpfac_gap_noidatend( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], newgaplen, 1, j, seq1[0]+i-1, seq2[0], 0 ); // i-1
#if USEGAPLENMTX
			pfactmp = calcpfac_gap_incomplete( gaplen1mtx[i-1][mpi], gaplen2mtx[i-1][mpi], newgaplen, i, j, seq1[0], seq2[0], 1 );
#endif
#if USEGAPLENHALF
			pfactmp = calcpfac_gap_incomplete( gaplen1half[i-1][mpi], gaplen2half[i-1][mpi], newgaplen, 1, j-mpi, seq1[0]+i-1, seq2[0]+mpi, 1 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(igap) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
				PFACERROR = 1;
			}
#endif


#if DEBUG
			reporterr( "%c-%c pfac for igap end incomplete = %f\n", seq1[0][i], seq2[0][j], pfac );
			reporterr( "mi when igap end checking = %f\n", mi );
			reporterr( "wm = %f, mi+fpenalty*pfac=%f\n", wm, mi+fpenalty*pfac );
#endif


#if ALGZGAP
			if( (g=mi+*fgcp2pt*gf1va) > wm )
#else
			if( (g=mi+fpenalty*pfac) > wm )
#endif
			{
				wm = g;
				*ijppt = -( j - mpi );
#if DEBUG80
				reporterr( "Selected as igap end! wm = %f, mi = %f\n", wm, mi );
				fprintf( stderr, "Jump from %d-%d (%c-%c) to %d (%c-%c)!\n", i, j, seq1[0][i], seq2[0][j], mpi, seq1[0][i-1], seq2[0][mpi] );
#endif
			}


#if 0 
			fprintf( stderr, "%5.0f->", wm );
#endif
//			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )

#if 0
//			reporterr( "Checking %c, (current pos = %c)\n", seq2[0][j+1], seq2[0][j] );
			sfac = 1.0;
			for( k=0; gaplen2[j+1]&&(gl=gaplen2[j+1][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
			{
//				reporterr( ".len = %d, .relend = %d\n", gaplen2[j+1][k].len, gaplen2[j+1][k].relend );
				if( gl - 1 == gaplen2[j+1][k].relend ) 
				{
					sfac -= gaplen2[j+1][k].freq;
//					reporterr( "Hit! sfac = %f\n", sfac );
				}
			}
			sfac2 = 1.0;
			for( k=0; gaplen1[i+1]&&(gl=gaplen1[i+1][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
				if( gaplen1[i+1][k].relend != -1 ) sfac2 -= gaplen1[i+1][k].freq;
			sfac *= sfac2;
#else
//			sfac = 0.0;
#endif



#if ALGZGAP
			if( (g=*prept+*ogcp2pt*gf1vapre) >= mi )
#else
//			if( (g=*prept + fpenalty * sfac ) >= mi )
			if( (g=*prept ) >= mi )
#endif
			{
//				mpibk = mpi;
//				mi = g - fpenalty * sfac;
				mi = g;
				mpi = j-1;
#if DEBUG80
				reporterr( "Selected as igap start! %c%d-%c%d, mi=%f, g=%f\n", seq1[0][i-1], i-1, seq2[0][mpi], mpi, mi, g );
#endif

#if FREEFREQUENTLY
//				freegaplenpartly( gaplen1ibestkamo[i-1], 0, i-1 );
				freegaplenpartly( gaplen2ibestkamo[i-1], j-3, j-2 );
#endif
//				freegaplenpartly( gaplen1jprev[mpibk], 0, lgth2 ); // full
//				freegaplenpartly( gaplen2jprev[mpibk], 0, lgth2-mpibk ); // half
//				if( gaplen1jprev[mpibk] ) FreeGaplenMtx( gaplen1jprev[mpibk], 0 );
//				gaplen1jprev[mpibk] = NULL;
//				if( gaplen2jprev[mpibk] ) FreeGaplenMtx( gaplen2jprev[mpibk], 0 );
//				gaplen2jprev[mpibk] = NULL;


//				addnewgaplen( gaplen1ibestkamo[i-1], gaplen1jprev[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2ibestkamo[i-1], gaplen2jprev[j-1], gaplen2, lgth2, -1, 0 );
//				copygaplenrestricted( gaplen1ibestkamo[i-1], gaplen1jprev[j-1], lgth1, -1, 0, i, i ); // i-1, i
				copygaplencompactx( gaplen1ibestkamo[i-1], gaplen1jprev[j-1], lgth1, -1, 0, 1, i ); // half
//				copygaplenrestricted( gaplen2ibestkamo[i-1], gaplen2jprev[j-1], lgth2, -1, 0, j, j ); // mpi, j
				copygaplencompactx( gaplen2ibestkamo[i-1], gaplen2jprev[j-1], lgth2, -1, 0, j, 1 ); //half


			}






//			reporterr( "g=%f, *prept=%f, mi=%f\n", g, *prept, mi );


#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if ALGZGAP
			pfac = 0.0; // CHUUI!
#else

//			pfac = calcpfac_gapex( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], i, j, j-mpi, seq1[0], seq2[0], 1 ); // i-1 
			pfac = calcpfac_gapex_noidatend( gaplen1ibestkamo[i-1], gaplen2ibestkamo[i-1], 1, j, j-mpi, seq1[0]+i, seq2[0], 1 ); // 1ibest->half, 2ibest->full
#if USEGAPLENMTX
			pfactmp = calcpfac_gapex( gaplen1mtx[i-1][mpi], gaplen2mtx[i-1][mpi], i, j, j-mpi, seq1[0], seq2[0], 1 );
#endif
#if USEGAPLENHALF
			pfactmp = calcpfac_gapex( gaplen1half[i-1][mpi], gaplen2half[i-1][mpi], 1, j-mpi, j-mpi, seq1[0]+i, seq2[0]+mpi, 1 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(igapex) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
				PFACERROR = 1;
			}
#endif







#if DEBUG
			reporterr( "%c-%c, igap extension check, pfac = %f\n\n", '=', seq2[0][j], pfac );
#endif
#endif
//			reporterr( "mi = %f -> ", mi );
			mi += fpenalty * pfac;
//			reporterr( "mi = %f\n", mi );


//			reporterr( "using %d-%d, %d, %d\n", *mpjpt, j-1, i, j );
			newgaplen = i-*mpjpt-1;
//			pfac = calcpfac_gap_incomplete( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], newgaplen, j, i, seq2[0], seq1[0], 0 ); // j-1 deha???


			pfac = calcpfac_gap_noidatend( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], newgaplen, 1, i, seq2[0]+j-1, seq1[0], 1 ); // 2jbestkamo->half, 1jbestkamo->full
#if USEGAPLENMTX
			pfactmp = calcpfac_gap_incomplete( gaplen2mtx[*mpjpt][j-1], gaplen1mtx[*mpjpt][j-1], newgaplen, j, i, seq2[0], seq1[0], 1 );
#endif
#if USEGAPLENHALF
			pfactmp = calcpfac_gap_incomplete( gaplen2half[*mpjpt][j-1], gaplen1half[*mpjpt][j-1], newgaplen, 1, i-*mpjpt, seq2[0]+j-1, seq1[0]+*mpjpt, 1 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(jgap) pfac=%f, but pfactmp=%f (i,j=%d,%d)\n", pfac, pfactmp, i, j );
//				exit( 1 );
				PFACERROR = 1;
			}
#endif

#if ALGZGAP
			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )
#else
			if( (g=*mjpt + fpenalty*pfac) > wm )
#endif
			{
				wm = g;
				*ijppt = +( i - *mpjpt );


#if FREEFREQUENTLY
				freegaplenpartly( gaplen1jbest[j-1], i-3, i-2 );
//				freegaplenpartly( gaplen2jbest[j-1], j-3, j-2 );
#endif


#if DEBUG
				reporterr( "Selected as jgap end!, pfac = %f\n", pfac );
				fprintf( stderr, "Jump from %d (%c) to %d (%c)!\n", j, seq1[0][j], *mpjpt, seq1[0][*mpjpt] );
#endif
//				addnewgaplen( gaplen1jbest[j-1], gaplen1jbestkamo[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jbest[j-1], gaplen2jbestkamo[j-1], gaplen2, lgth2, -1, 0 );
				copygaplencompactx( gaplen1jbest[j-1], gaplen1jbestkamo[j-1], lgth1, -1, 0, i, i );// *mpjpt, i
//				copygaplenrestricted( gaplen2jbest[j-1], gaplen2jbestkamo[j-1], lgth2, -1, 0, j, j ); // j-1, j
				copygaplencompactx( gaplen2jbest[j-1], gaplen2jbestkamo[j-1], lgth2, -1, 0, 1, 1 ); // half!




			}


//			extendgaplenpartly( gaplen1jbest[j-1], gaplen1, i, i ); // tmptmptmp
//			extendgaplenpartly( gaplen2jbest[j-1], gaplen2, 0, 0 ); // tmptmptmp

#if 0
			sfac = 1.0;
			for( l=0; gaplen1[i+1]&&(gl=gaplen1[i+1][l].len)!=-1; l++ ) // ososugi!  hash ni atode henkou
				if( gl - 1 == gaplen1[i+1][l].relend ) sfac -= gaplen1[i+1][l].freq;
			sfac2 = 1.0;
			for( k=0; gaplen2[j+1]&&(gl=gaplen2[j+1][k].len)!=-1; k++ ) // ososugi!  hash ni atode henkou
				if( gaplen2[j+1][k].relend != -1 ) sfac2 -= gaplen2[j+1][k].freq;
			sfac *= sfac2;
#else
//			sfac = 0.0;
#endif

#if DEBUG
			reporterr( " (jgap start check i=%d) -> *prept=%f, *mjpt=%f\n", i, seq1[0][i], seq2[0][j], *prept, *mjpt );
#endif

#if ALGZGAP
			if( (g=*prept+ ogcp1va* *gf2ptpre) >= *mjpt )
#else
//			if( (g=*prept + fpenalty * sfac ) >= *mjpt )
			if( (g=*prept ) >= *mjpt )
#endif
			{
//				*mjpt = g - fpenalty * sfac;
				*mjpt = g;
				*mpjpt = i-1;
#if DEBUG
				reporterr( "Selected as jgap start!\n" );
#endif


#if FREEFREQUENTLY
				freegaplenpartly( gaplen1jbestkamo[j-1], i-3, i-2 );
//				freegaplenpartly( gaplen2jbestkamo[j-1], j-3, j-2 );
#endif


//				addnewgaplen( gaplen1jbestkamo[j-1], gaplen1jprev[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], gaplen2, lgth2, -1, 0 );
//				reporterr( "copying gaplen1jbestkamo[%d-1] from galpen1jprev, j=%d, i=%d\n", j, j, i );
				copygaplencompactx( gaplen1jbestkamo[j-1], gaplen1jprev[j-1], lgth1, -1, 0, i, i ); // *mpjpt, i
//				copygaplenrestricted( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], lgth2, -1, 0, j, j ); // j-1, j
//				copygaplencompactx( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], lgth2, -1, 0, j, 1 ); // half!
//				reporterr( "copying gaplen2jbestkamo[%d-1] from galpen2jprev\n", j );
				copygaplencompactx( gaplen2jbestkamo[j-1], gaplen2jprev[j-1], lgth2-j, -1, 0, 1, 1 ); // ryouhou half!


//				if( j==2 && i==1 ) exit( 1 );



			}

//			extendgaplenpartly( gaplen1ibestkamo[i-1], gaplen1, 0, 0 ); // tmptmptmp
//			extendgaplenpartly( gaplen2ibestkamo[i-1], gaplen2, j, j ); // tmptmptmp


//			extendgaplenpartly( gaplen1jbestkamo[j-1], gaplen1, i, i ); // tmptmptmp
//			extendgaplenpartly( gaplen2jbestkamo[j-1], gaplen2, 0, 0 ); // tmptmptmp


#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if ALGZGAP
			pfac = 0.0;
#else

//			pfactmp = calcpfac_gapex( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], j, i, i-*mpjpt, seq2[0], seq1[0], 0 ); // j-1 
			pfactmp = calcpfac_gapex_noidatend( gaplen2jbestkamo[j-1], gaplen1jbestkamo[j-1], 1, i, i-*mpjpt, seq2[0]+j, seq1[0], 0 ); // 2jbestkamo->half, 1jbestkamo->full
#if USEGAPLENMTX
			pfac = calcpfac_gapex( gaplen2mtx[*mpjpt][j-1], gaplen1mtx[*mpjpt][j-1], j, i, i-*mpjpt, seq2[0], seq1[0], 0 );
#endif
#if USEGAPLENHALF
			pfac = calcpfac_gapex( gaplen2half[*mpjpt][j-1], gaplen1half[*mpjpt][j-1], 1, i-*mpjpt, i-*mpjpt, seq2[0]+j, seq1[0]+*mpjpt, 0 );
#endif
#if USEGAPLENMTX || USEGAPLENHALF
			if( pfac != pfactmp )
			{
				reporterr( "(jgapex) pfac=%f, but pfactmp=%f (i,j=%d,%d) diff=%f\n", pfac, pfactmp, i, j, pfac-pfactmp );
//				exit( 1 );
				PFACERROR = 1;
			}
#endif
			pfac = pfactmp;
#if DEBUG
			reporterr( "%c-%c, jgap extension check (j=%d), pfac = %f\n", seq1[0][i], '=', j, pfac );
#endif
#endif
			m[j] += fpenalty * pfac;



			if( trywarp )
			{
#if USE_PENALTY_EX
				if( ( g=*prevwmrecordspt++ + fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] ) ) > wm ) // naka ha osokute kamawanai
#else
				if( ( g=*prevwmrecordspt++ + fpenalty_shift ) > wm ) // naka ha osokute kamawanai
#endif
				{
//					fprintf( stderr, "WARP!!\n" );
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
	
#if 0
				fprintf( stderr, "%5.0f ", wm );
#endif
				curm = *curpt + wm;
	
				if( *wmrecords1pt > *wmrecordspt )
				{
					*wmrecordspt = *wmrecords1pt;
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
				}
				if( curm > *wmrecordspt )
				{
					*wmrecordspt = curm;
					*warpipt = i;
					*warpjpt = j;
				}
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}

#if DEBUG
			reporterr( "extention-x 1j???, before extention-x, j=%d\n", j );
			showgaplen( gaplen1jcurr[j], 100 );
#endif
			extendgaplencompactx( gaplen1jcurr[j], gaplen1, i );

#if DEBUG
			reporterr( "after extention-x\n" );
			showgaplen( gaplen1jcurr[j], 100 );
			reporterr( "extention-x 2j???\n" );
#endif
			extendgaplencompactx( gaplen2jcurr[j], gaplen2+j, 0 );


#if 1
			if( *ijppt < 0 )
			{
#if DEBUG
				reporterr( "Path: %d-%d->%d-%d, i=%d (%c), j=%d (%c), *ijppt=%d\n", i, j, i-1, j+*ijppt, i, seq1[0][i], j, seq2[0][j], *ijppt );
				reporterr( "Inserting %d gaps to gaplen1 and copying gaplen2 (%c%d-%c%d)\n", -*ijppt-1, seq1[0][i], i, seq2[0][j], j );
#endif
#if USEGAPLENMTX
//				addnewgaplen( gaplen1mtx[i][j], gaplen1mtx[i-1][j+*ijppt], gaplen1, lgth1, i, -*ijppt-1 );
//				addnewgaplen( gaplen2mtx[i][j], gaplen2mtx[i-1][j+*ijppt], gaplen2, lgth2, -1, 0 );
				copygaplenrestricted( gaplen1mtx[i][j], gaplen1mtx[i-1][j+*ijppt], lgth1, i, -*ijppt-1, i, i ); // 0, lgth1
				copygaplenrestricted( gaplen2mtx[i][j], gaplen2mtx[i-1][j+*ijppt], lgth2, -1, 0, j, j ); // 0, lgth2
#endif
#if USEGAPLENHALF
				copygaplenrestricted_zurasu( gaplen1half[i][j], gaplen1half[i-1][j+*ijppt], lgth1, 0, -*ijppt-1, 0, 0, 1, 1 ); // 0, lgth1
				copygaplenrestricted_zurasu( gaplen2half[i][j], gaplen2half[i-1][j+*ijppt], lgth2, -1, 0, 0, 0, -*ijppt, -*ijppt ); // 0, lgth2
#endif
//				addnewgaplen( gaplen1jcurr[j], gaplen1jprev[j+*ijppt], gaplen1, lgth1, i, -*ijppt-1 );
//				addnewgaplen( gaplen2jcurr[j], gaplen2jprev[j+*ijppt], gaplen2, lgth2, -1, 0 );
//				reporterr( "copying gaplen1jcurr from gaplen1jbest, with a %d insertion\n", -*ijppt-1 ); 
				copygaplencompactx( gaplen1jcurr[j], gaplen1jprev[j+*ijppt], lgth1, i, -*ijppt-1, i, i ); // scope: i+*ijppt+1, i ?
//				reporterr( "copy end\n" ); 
//				copygaplenrestricted( gaplen2jcurr[j], gaplen2jprev[j+*ijppt], lgth2, -1, 0, j, j );
				copygaplencompactx( gaplen2jcurr[j], gaplen2jprev[j+*ijppt], lgth2, -1, 0, 0, -*ijppt ); // half! ryouho zureteru
			}
			else if( *ijppt > 0 )
			{
#if DEBUG
				reporterr( "Path: %d-%d->%d-%d, i=%d (%c), j=%d (%c), *ijppt=%d\n", i, j, i-*ijppt, j-1, i, seq1[0][i], j, seq2[0][j], *ijppt );
				reporterr( "Copying gaplen1 inserting %d gaps to gaplen2 (%c%d-%c%d)\n", *ijppt-1, seq1[0][i], i, seq2[0][j], j );
#endif
#if USEGAPLENMTX
//				addnewgaplen( gaplen1mtx[i][j], gaplen1mtx[i-*ijppt][j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2mtx[i][j], gaplen2mtx[i-*ijppt][j-1], gaplen2, lgth2, j, *ijppt-1 );
				copygaplenrestricted( gaplen1mtx[i][j], gaplen1mtx[i-*ijppt][j-1], lgth1, -1, 0, i, i ); // 0, lgth1
				copygaplenrestricted( gaplen2mtx[i][j], gaplen2mtx[i-*ijppt][j-1], lgth2, j, *ijppt-1, j, j ); // 0, lgth2
#endif
#if USEGAPLENHALF
				copygaplenrestricted_zurasu( gaplen1half[i][j], gaplen1half[i-*ijppt][j-1], lgth1, -1, 0, 0, 0, *ijppt, *ijppt ); // 0, lgth1
				copygaplenrestricted_zurasu( gaplen2half[i][j], gaplen2half[i-*ijppt][j-1], lgth2, 0, *ijppt-1, 0, 0, 1, 1 ); // 0, lgth2
#endif
//				addnewgaplen( gaplen1jcurr[j], gaplen1jbest[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jcurr[j], gaplen2jbest[j-1], gaplen2, lgth2, j, *ijppt-1 );
				copygaplencompactx( gaplen1jcurr[j], gaplen1jbest[j-1], lgth1, -1, 0, i, i );
//				copygaplenrestricted( gaplen2jcurr[j], gaplen2jbest[j-1], lgth2, j, *ijppt-1, j, j ); // j-*ijppt+1?
//				copygaplenrestricted_zurasu( gaplen2jcurr[j], gaplen2jbest[j-1], lgth2, 0, *ijppt-1, 0, 0, j, j ); // 2jcurr->half, but 2jbest->full, imanotokoro
				copygaplencompactx( gaplen2jcurr[j], gaplen2jbest[j-1], lgth2, 0, *ijppt-1, 0, 1 ); //ryouhou half

			}
			else
#endif
			{
#if DEBUG
				reporterr( "Path: %d-%d->%d-%d, i=%d (%c), j=%d (%c), *ijppt=%d\n", i, j, i-1, j-1, i, seq1[0][i], j, seq2[0][j], *ijppt );
				reporterr( "Copying gaplen1 and gaplen2 (%c%d-%c%d)\n", seq1[0][i], i, seq2[0][j], j );
#endif	
#if USEGAPLENMTX
//				addnewgaplen( gaplen1mtx[i][j], gaplen1mtx[i-1][j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2mtx[i][j], gaplen2mtx[i-1][j-1], gaplen2, lgth2, -1, 0 );
				copygaplenrestricted( gaplen1mtx[i][j], gaplen1mtx[i-1][j-1], lgth1, -1, 0, i, i ); // 0, lgth1
				copygaplenrestricted( gaplen2mtx[i][j], gaplen2mtx[i-1][j-1], lgth2, -1, 0, j, j ); // 0, lgth2
#endif
#if USEGAPLENHALF
				copygaplenrestricted_zurasu( gaplen1half[i][j], gaplen1half[i-1][j-1], lgth1, -1, 0, 0, 0, 1, 1 ); // 0, lgth1
				copygaplenrestricted_zurasu( gaplen2half[i][j], gaplen2half[i-1][j-1], lgth2, -1, 0, 0, 0, 1, 1 ); // 0, lgth2
#endif
//				addnewgaplen( gaplen1jcurr[j], gaplen1jprev[j-1], gaplen1, lgth1, -1, 0 );
//				addnewgaplen( gaplen2jcurr[j], gaplen2jprev[j-1], gaplen2, lgth2, -1, 0 );
				copygaplencompactx( gaplen1jcurr[j], gaplen1jprev[j-1], lgth1, -1, 0, i, i );
//				copygaplenrestricted( gaplen2jcurr[j], gaplen2jprev[j-1], lgth2, -1, 0, j, j );
				copygaplencompactx( gaplen2jcurr[j], gaplen2jprev[j-1], lgth2, -1, 0, 0, 1 ); // half
			}

#if DEBUG
			reporterr( "at the end of j loop, gaplen1jcurr[%d] = \n", j );
			showgaplen( gaplen1jcurr[j], 100 );
			reporterr( "at the end of j loop, gaplen1prev[%d] = \n", j );
			showgaplen( gaplen1jprev[j], 100 );
#endif



#if 1
			freegaplenpartly( gaplen1jcurr[j-1], i-3, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // half!
			freegaplenpartly( gaplen1jbestkamo[j-1], i-3, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbestkamo[j-1], j-3, j-2 ); // -1 dehanaku -2??
			freegaplenpartly( gaplen1jbest[j-1], i-3, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbest[j-1], j-3, j-2 ); // -1 dehanaku -2??
#else
			freegaplenpartly( gaplen1jprev[j-1], 0, i-2 ); // -1 dehanaku -2??
			freegaplenpartly( gaplen1jcurr[j-1], 0, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jcurr[j-1], j-3, j-2 ); // half!
			freegaplenpartly( gaplen1jbestkamo[j-1], 0, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbestkamo[j-1], j-3, j-2 ); // -1 dehanaku -2??
			freegaplenpartly( gaplen1jbest[j-1], 0, i-2 ); // -1 dehanaku -2??
//			freegaplenpartly( gaplen2jbest[j-1], j-3, j-2 ); // -1 dehanaku -2??
#endif


#if USEGAPLENMTX
//			freegaplenpartly( gaplen1mtx[i-1][j-1], 0, i-2 );
//			freegaplenpartly( gaplen2mtx[i-1][j-1], 0, j-2 );
#endif


			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
		lastverticalw[i] = currentw[lgth2-1];

#if 1
//		freegaplenpartly( gaplen1icurr[i-1], i-1, i-1 );
		freegaplenpartly( gaplen1icurr[i-1], 0, lgth1-i );
		freegaplenpartly( gaplen2icurr[i-1], 0, lgth2 );
//		freegaplenpartly( gaplen1ibestkamo[i-1], i-1, i-1 );
		freegaplenpartly( gaplen1ibestkamo[i-1], 0, lgth1-i );
		freegaplenpartly( gaplen2ibestkamo[i-1], 0, lgth2 );
#endif

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
#if 0
		fprintf( stderr, "i=%d, %15.5f \n", i, wm );
#endif
//if( i == 2 ) exit( 1 );
	}

	if( trywarp )
	{
//		fprintf( stderr, "wm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}


#if OUTGAP0TRY
	if( !outgap )
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * ( lgth2 - j ) / 2.0;
		for( i=1; i<lgth1+1; i++ )
			lastverticalw[i] -= offset * ( lgth1 - i  / 2.0);
	}
#endif
		
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/
	if( constraint )
	{
		Atracking_localhom( impmatch, currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, warpis, warpjs, warpbase );
	}
	else
		Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase );

	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );

//	fprintf( stderr, "### impmatch = %f\n", *impmatch );

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );
#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
#endif

//	reporterr( "clearing\n" );
	cleargaplens( gaplens );
	if( masklist1 ) FreeIntMtx( masklist1 ); masklist1 = NULL;
	if( masklist2 ) FreeIntMtx( masklist2 ); masklist2 = NULL;
	if( nmask ) free( nmask ); nmask = NULL;

#if USEGAPLENMTX
/* maikai free */
	reporterr( "Freeing!\n" );
	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen1mtx[i][j] ) FreeGaplenMtx( gaplen1mtx[i][j], 0 );
			gaplen1mtx[i][j] = NULL;
		}
		free( gaplen1mtx[i] );
		gaplen1mtx[i] = NULL;
	}
	free( gaplen1mtx );
	gaplen1mtx = NULL;

	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen2mtx[i][j] ) FreeGaplenMtx( gaplen2mtx[i][j], 0 );
			gaplen2mtx[i][j] = NULL;
		}
		free( gaplen2mtx[i] );
		gaplen2mtx[i] = NULL;
	}
	free( gaplen2mtx );
	gaplen2mtx = NULL;
#endif


#if USEGAPLENHALF
	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen1half[i][j] ) FreeGaplenMtx( gaplen1half[i][j], 0 );
			gaplen1half[i][j] = NULL;
		}
		free( gaplen1half[i] );
		gaplen1half[i] = NULL;
	}
	free( gaplen1half );
	gaplen1half = NULL;

	for( i=0; i<lgth1+1; i++ ) 
	{
		for( j=0; j<lgth2+1; j++ ) 
		{
			if( gaplen2half[i][j] ) FreeGaplenMtx( gaplen2half[i][j], 0 );
			gaplen2half[i][j] = NULL;
		}
		free( gaplen2half[i] );
		gaplen2half[i] = NULL;
	}
	free( gaplen2half );
	gaplen2half = NULL;
#endif
/* maikai free */


#if WMCHECK
	fprintf( stderr, "wm = %f\n", wm - *impmatch);
	fprintf( stderr, "*impmatch = %f\n", *impmatch);

	int kenzan = 0;
	for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
	{
		kenzan += pairgapcount( mseq1[i], mseq2[j] );
	}


	reporterr( "kenzan = %d -> %f\n", kenzan, (double)kenzan /( icyc*jcyc ) );

	double pairscore, nogappairscore, diff;
	char **pseq;
	pseq = AllocateCharMtx( 2, strlen( seq1[0] ) + 1 );
	pairscore = nogappairscore = 0.0;
#if 1
	for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
	{
		strcpy( pseq[0], seq1[i] );
		strcpy( pseq[1], seq2[j] );
		commongappick( 2, pseq );
		c = which[i][j];
		pairscore += eff1[i] * eff2[j] * naivepairscore11_dynmtx( matrices[c], pseq[0], pseq[1], penalty );
		nogappairscore += eff1[i] * eff2[j] * naivepairscore11_dynmtx( matrices[c], pseq[0], pseq[1], 0 );
	}
#else
	for( c=0; c<maxdistclass; c++ )
	{
		for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
		{
			strcpy( pseq[0], seq1[i] );
			strcpy( pseq[1], seq2[j] );
			commongappick( 2, pseq );
			pairscore += eff1s[c][i] * eff2s[c][j] * naivepairscore11_dynmtx( matrices[c], pseq[0], pseq[1], penalty );
			nogappairscore += eff1s[c][i] * eff2s[c][j] * naivepairscore11_dynmtx( matrices[c], pseq[0], pseq[1], 0 );
		}
	}
#endif

	FreeCharMtx( pseq );
	diff = (pairscore - wm + *impmatch ) / (double)strlen( seq1[0] );
	reporterr( "pairscore = %f\n", (double)pairscore );
	reporterr( "pairscore-nogappairscore = %f\n", (double)(pairscore-nogappairscore) );
	reporterr( "pairscore-nogappairscore / penalty = %f\n", (double)(pairscore-nogappairscore)/(double)(fpenalty) );
	reporterr( "diff = %f\n\n", diff );

#if 1
	if( ( !trywarp && fabs( diff ) > 0.01 ) || PFACERROR )
//	if( abs( pairscore - wm +*impmatch ) > 0.01 )
#else
	if( abs( pairscore - wm +*impmatch ) > 0.01 )
#endif
//	if( abs( pairscore - wm +*impmatch ) > 0.01 )
	{
		for( i=0; i<icyc; i++ )
			printf( ">group1\n%s\n", seq1[i] );
		for( j=0; j<jcyc; j++ )
			printf( ">group2\n%s\n", seq2[j] );
		exit( 1 );
	}
#else
	reporterr( "\n" );
#endif

#if 0
//	if( strlen( seq1[0] ) - lgth1 > 100 && icyc > 1 || strlen( seq2[0] ) - lgth2 > 100 & jcyc > 1 )
	if( strstr( seq1[0], "LNDDAT" ) && icyc == 1 || strstr( seq2[0], "LNDDAT" ) && jcyc==1)
	{
		for( i=0; i<icyc; i++ )
			printf( ">group1\n%s\n", seq1[i] );
		for( j=0; j<jcyc; j++ )
			printf( ">group2\n%s\n", seq2[j] );
		exit( 1 );
	}
#endif


	return( wm );
}
