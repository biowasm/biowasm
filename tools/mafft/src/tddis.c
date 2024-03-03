#include "mltaln.h"

#define DEBUG 0
#define USEDISTONTREE 1

#if 0
void mdfymtx( char pair[njob][njob], int s1, double **partialmtx, double **mtx )
#else
void mdfymtx( char **pair, int s1, double **partialmtx, double **mtx )
#endif
{
	int i, j;
	int icount, jcount;
#if DEBUG
	FILE *fp;
	static char name[M][B];

	for( i=0; i<M; i++ ) name[i][0] = 0;
	fprintf( stdout, "s1 = %d\n", s1 );
	for( i=0; i<njob; i++ ) 
	{
		for( j=0; j<njob; j++ ) 
		{
			printf( "%#2d", pair[i][j] );
		}
		printf( "\n" );
	}
#endif

	for( i=0, icount=0; i<njob-1; i++ )
	{
		if( !pair[s1][i] ) continue;
		for( j=i+1, jcount=icount+1; j<njob; j++ ) 
		{
			if( !pair[s1][j] ) continue;
			partialmtx[icount][jcount] = mtx[i][j];
			jcount++;
		}
		icount++;
	}
#if DEBUG
	fp = fopen( "hat2.org", "w" );
	WriteHat2( fp, njob, name, mtx );
	fclose( fp );
	fp = fopen( "hat2.mdf", "w" );
	WriteHat2( fp, icount, name, partialmtx );
	fclose( fp );
#endif
		
}

		
double score_calc( char **seq, int s )  /* method 3  */
{
    int i, j, k, c;
    int len = strlen( seq[0] );
    double score;
    int tmpscore;
    char *mseq1, *mseq2;

    score = 0.0;
    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0;
            c = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                c++;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            /*
            if( mseq1[0] == '-' || mseq2[0] == '-' )
            {
                for( k=0; k<len; k++ )
                {
                    if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                    if( !( mseq1[k] != '-' && mseq2[k] != '-' ) ) 
                    {
                        c--;
                        tmpscore -= penalty;
                        break;
                    }
                    else break;
                }
            }
            if( mseq1[len-1] == '-' || mseq2[len-1] == '-' )
            {
                for( k=0; k<len; k++ )
                {
                    if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                    if( !( mseq1[k] != '-' && mseq2[k] != '-' ) ) 
                    {
                        c--;
                        tmpscore -= penalty;
                        break;
                    }
                    else break;
                }
            }
            */
            score += (double)tmpscore / (double)c;
        }
    }
    score = (double)score / ( ( (double)s * ((double)s-1.0) ) / 2.0 );
	fprintf( stderr, "score in score_calc = %f\n", score );
    return( score );
}

void cpmx_calc( char **seq, double **cpmx, double *eff, int lgth, int clus )
{
	int  i, j, k;
	double totaleff = 0.0;

	for( i=0; i<clus; i++ ) totaleff += eff[i]; 
	for( i=0; i<nalphabets; i++ ) for( j=0; j<lgth; j++ ) cpmx[i][j] = 0.0;
	for( j=0; j<lgth; j++ ) for( k=0; k<clus; k++ )
			cpmx[(int)amino_n[(unsigned char)seq[k][j]]][j] += (double)eff[k] / totaleff;
}


void cpmx_calc_new_bk( char **seq, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0 
{
    int  i, j, k;
    double feff;

    for( i=0; i<nalphabets; i++ ) for( j=0; j<lgth; j++ ) cpmx[i][j] = 0.0;
    for( k=0; k<clus; k++ )
    {
        feff = (double)eff[k];
        for( j=0; j<lgth; j++ ) 
        {
            cpmx[(int)amino_n[(int)seq[k][j]]][j] += feff;
        }
    }
}

void cpmx_calc_add( char **seq, double **cpmx, double *eff, int lgth, int clus ) // lastmem = newmem; summ eff must be 1.0
{
	double neweff, orieff;
	int newmem, i, j;

	newmem = clus-1;
	neweff = eff[clus-1];
	orieff = 1.0-neweff;
#if 1 // TESTING Feb/1/18:00
	for( j=0; j<lgth; j++ )
	{
		for( i=0;i<nalphabets; i++ ) cpmx[i][j] *= orieff;
		cpmx[(unsigned char)amino_n[(unsigned char)seq[newmem][j]]][j] += neweff;
	}
#else // possibly faster?
	for( i=0;i<nalphabets; i++ ) 
	{
		for( j=0; j<lgth; j++ ) cpmx[i][j] *= orieff;
	}
	for( j=0; j<lgth; j++ ) cpmx[(unsigned char)amino_n[(unsigned char)seq[newmem][j]]][j] += neweff;
#endif
}

void cpmx_calc_new( char **seq, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0
{
	int  i, j, k;
	double feff;
	double *cpmxpt, **cpmxptpt;
	char *seqpt;

	j = nalphabets;
	cpmxptpt = cpmx;
	while( j-- )
	{
		cpmxpt = *cpmxptpt++;
		i = lgth;
		while( i-- )
			*cpmxpt++ = 0.0;
	}
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		seqpt = seq[k];
//		fprintf( stderr, "seqpt = %s, lgth=%d\n", seqpt, lgth );
		for( j=0; j<lgth; j++ )
		{
			cpmx[(unsigned char)amino_n[(unsigned char)*seqpt++]][j] += feff;
		}
	}
}
void MScpmx_calc_new( char **seq, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0
{
	int  i, j, k;
	double feff;
	double **cpmxptpt, *cpmxpt;
	char *seqpt;

	j = lgth;
	cpmxptpt = cpmx;
	while( j-- )
	{
		cpmxpt = *cpmxptpt++;
		i = nalphabets;
		while( i-- )
			*cpmxpt++ = 0.0;
	}
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		seqpt = seq[k];
		cpmxptpt = cpmx;
		j = lgth;
		while( j-- )
			(*cpmxptpt++)[(int)amino_n[(unsigned char)*seqpt++]] += feff;
	}
#if 0
	for( j=0; j<lgth; j++ ) for( i=0; i<nalphabets; i++ ) cpmx[j][i] = 0.0;
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		for( j=0; j<lgth; j++ ) 
			cpmx[j][(int)amino_n[(int)seq[k][j]]] += feff;
	}
#endif
}

void cpmx_ribosum( char **seq, char **seqr, char *dir, double **cpmx, double *eff, int lgth, int clus ) // summ eff must be 1.0
{
	int  i, j, k;
	double feff;
	double **cpmxptpt, *cpmxpt;
	char *seqpt, *seqrpt, *dirpt;
	int code, code1, code2;

	j = lgth;
	cpmxptpt = cpmx;
	while( j-- )
	{
		cpmxpt = *cpmxptpt++;
		i = 37;
		while( i-- )
			*cpmxpt++ = 0.0;
	}
	for( k=0; k<clus; k++ )
	{
		feff = (double)eff[k];
		seqpt = seq[k];
		seqrpt = seqr[k];
		dirpt = dir;
		cpmxptpt = cpmx;
		j = lgth;
		while( j-- )
		{
#if 0
			code1 = amino_n[(int)*seqpt];
			if( code1 > 3 ) code = 36;
			else
				code = code1;
#else
			code1 = amino_n[(unsigned char)*seqpt];
			code2 = amino_n[(unsigned char)*seqrpt];
			if( code1 > 3 ) 
			{
				code = 36;
			}
			else if( code2 > 3 )
			{
				code = code1;
			}
			else if( *dirpt == '5' ) 
			{
				code = 4 + code2 * 4  + code1;
			}
			else if( *dirpt == '3' ) 
			{
				code = 20 + code2 * 4  + code1;
			}
			else // if( *dirpt == 'o' ) // nai 
			{
				code = code1;
			}
#endif

//			fprintf( stderr, "%c -> code=%d toa=%d, tog=%d, toc=%d, tot=%d, ton=%d, efee=%f\n", *seqpt, code%4, ribosumdis[code][4+0], ribosumdis[code][4+1], ribosumdis[code][4+2], ribosumdis[code][20+3], ribosumdis[code][36], feff );

			seqpt++;
			seqrpt++;
			dirpt++;
			
			(*cpmxptpt++)[code] += feff;
		}
	}
}

void mseqcat( char **seq1, char **seq2, double **eff, double *effarr1, double *effarr2, char name1[M][B], char name2[M][B], int clus1, int clus2 )
{
	int i, j;
    for( i=0; i<clus2; i++ )
        seq1[clus1+i] = seq2[i];
    for( i=0; i<clus2; i++ ) strcpy( name1[clus1+i], name2[i] );

	for( i=0; i<clus1; i++ ) for( j=0; j<clus1; j++ ) eff[i][j] = effarr1[i]* effarr1[j]; 
	for( i=0; i<clus1; i++ ) for( j=clus1; j<clus1+clus2; j++ ) eff[i][j] = effarr1[i]* effarr2[j-clus1]; 
	for( i=clus1; i<clus1+clus2; i++ ) for( j=0; j<clus1; j++ ) eff[i][j] = effarr2[i-clus1] * effarr1[j]; 
	for( i=clus1; i<clus1+clus2; i++ ) for( j=clus1; j<clus1+clus2; j++ ) eff[i][j] = effarr2[i-clus1] * effarr2[j-clus1]; 
}

	

#if 0
int conjuction( char pair[njob][njob], int s, char **seq, char **aseq, double *peff, double *eff, char name[M][B], char aname[M][B], char *d )
#else
int conjuctionforgaln( int s0, int s1, char **seq, char **aseq, double *peff, double *eff, char **name, char **aname, char *d )
#endif
{
	int m, k;
	char b[B];
	double total;

#if 0
	fprintf( stderr, "s0 = %d, s1 = %d\n", s0, s1 );
#endif

	total = 0.0;
	d[0] = 0;
	for( m=s0, k=0; m<s1; m++ )
	{
		{
			sprintf( b, " %d", m+1 ); 
#if 1
			if( strlen( d ) < 100 ) strcat( d, b );
#else
			strcat( d, b );
#endif
			aseq[k] = seq[m];
			peff[k] = eff[m];
			total += peff[k];
#if 0
			strcpy( aname[k], name[m] );
#endif
			k++;
		}
	}
#if 1
	for( m=0; m<k; m++ )
	{
		peff[m] /= total;
//		fprintf( stderr, "peff[%d] = %f\n", m, peff[m] );
	}
#endif
	return( k );
}

void makegrouprna( RNApair ***group, RNApair ***all, int *memlist )
{
	int k, m;
	for( k=0; (m=*memlist)!=-1; memlist++, k++ )
	{
		group[k] = all[m];
	}
}


int fastconjuction_noweight( int *memlist, char **seq, char **aseq, double *peff, char *d )
{
	int m, k, dln;
	char b[B];
	double total;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 ); 
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		peff[k] = 1.0;
		total += peff[k];
	}
#if 1
	for( m=0; m<k; m++ )
		peff[m] /= total;
#endif
	return( k );
}

int fastconjuction_noname_kozo( int *memlist, char **seq, char **aseq, double *peff, double *eff, double *peff_kozo, double *eff_kozo, char *d )
{
	int m, k, dln;
	char b[B];
	double total;
	double total_kozo;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	total_kozo = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 ); 
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		peff[k] = eff[m];
		peff_kozo[k] = eff_kozo[m];
		total += peff[k];
		total_kozo += peff_kozo[k];
	}
#if 1
	for( m=0; m<k; m++ )
	{
//		fprintf( stderr, "peff[%d] = %f\n", m, peff[m] );
		peff[m] /= total;
	}
	if( total_kozo ) 
	{
		for( m=0; m<k; m++ )
		{
//			fprintf( stderr, "peff[%d] = %f\n", m, peff[m] );
			peff_kozo[m] /= total_kozo;
			if( peff_kozo[m] > 0.0 ) peff_kozo[m] += peff[m];
		}
	}
	else  //iranai
	{
		for( m=0; m<k; m++ )
		{
			peff_kozo[m] = 0.0;
		}
	}
#endif

//	fprintf( stderr, "\n\ntbfast_total_kozo = %f\n\n", total_kozo );
	return( k );
}


#if 0
int fastconjuction_target( int *memlist, char **seq, char **aseq, double *peff, double *eff, char *d, double mineff, int *targetmap )
{
	int m, k, dln;
	char b[B];
	double total;
	int *memlistbk = memlist;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 ); 
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		if( eff[m] < mineff )
			peff[k] = mineff;
		else
			peff[k] = eff[m];

		total += peff[k];
	}
#if 1
	for( m=0; m<k; m++ )
	{
//		fprintf( stderr, "Apr17   peff[%d] = %f\n", m, peff[m] );
		peff[m] /= total;
	}
#endif

	return( k );
}
#endif


int fastconjuction_noname( int *memlist, char **seq, char **aseq, double *peff, double *eff, char *d, double mineff, double *oritotal  )
{
	int m, k, dln;
	char b[B];
	double total;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 ); 
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		if( eff[m] < mineff )
			peff[k] = mineff;
		else
			peff[k] = eff[m];

		total += peff[k];
	}
	if( oritotal ) *oritotal = total;
#if 1
	for( m=0; m<k; m++ )
	{
//		fprintf( stderr, "Apr17   peff[%d] = %20.10f\n", m, peff[m] );
		peff[m] /= total;
	}
#endif
	return( k );
}

int fastconjuction( int *memlist, char **seq, char **aseq, double *peff, double *eff, char name[M][B], char aname[M][B], char *d )
{
	int m, k, dln;
	char b[B];
	double total;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 ); 
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		peff[k] = eff[m];
		total += peff[k];
#if 0
			strcpy( aname[k], name[m] );
#endif
	}
#if 1
	for( m=0; m<k; m++ )
		peff[m] /= total;
#endif
	return( k );
}



int conjuctionfortbfast_old( char **pair, int s, char **seq, char **aseq, double *peff, double *eff, char *d )
{
	int m, k;
	char *b;
	double total;

	b = calloc( B, sizeof( char ) );
#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	for( m=s, k=0; m<njob; m++ )
	{
		if( pair[s][m] != 0 )
		{
			sprintf( b, " %d", m+1 ); 
#if 1
			if( strlen( d ) < 100 ) strcat( d, b );
#else
			strcat( d, b );
#endif
			aseq[k] = seq[m];
			peff[k] = eff[m];
			total += peff[k];
#if 0
			strcpy( aname[k], name[m] );
#endif
			k++;
		}
	}
#if 1
	for( m=0; m<k; m++ )
		peff[m] /= total;
#endif
	free( b );
	return( k );
}

int conjuction( char **pair, int s, char **seq, char **aseq, double *peff, double *eff, char **name, char **aname, char *d )
{
	int m, k;
	char b[B];
	double total;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	for( m=s, k=0; m<njob; m++ )
	{
		if( pair[s][m] != 0 )
		{
			sprintf( b, " %d", m+1 ); 
#if 1
			if( strlen( d ) < 100 ) strcat( d, b );
#else
			strcat( d, b );
#endif
			aseq[k] = seq[m];
			peff[k] = eff[m];
			total += peff[k];
#if 0
			strcpy( aname[k], name[m] );
#endif
			k++;
		}
	}
#if 0
	for( m=0; m<k; m++ )
		peff[m] /= total;
#endif
	return( k );
}
			
void doubledelete( double **cpmx, int d, int len )
{
	int i, j;

	for( i=d; i<len-1; i++ )
	{
		for( j=0; j<nalphabets; j++ )
		{
			cpmx[j][i] = cpmx[j][i+1]; 
		}
	}
}
		
void chardelete( char *seq, int d )
{
	char b[N]; 

		strcpy( b, seq+d+1 );
		strcpy( seq+d, b );
}

int RootBranchNode( int nseq, int ***topol, int step, int branch )
{
	int i, j, s1, s2, value;

	value = 1;
	for( i=step+1; i<nseq-2; i++ ) 
	{
		for( j=0; (s1=topol[i][0][j])>-1; j++ ) 
			if( s1 == topol[step][branch][0] ) value++;
		for( j=0; (s2=topol[i][1][j])>-1; j++ ) 
			if( s2 == topol[step][branch][0] ) value++;
	}
	return( value );
}

void BranchLeafNode( int nseq, int ***topol, int *node, int step, int branch )
{
	int i, j, k, s;

	for( i=0; i<nseq; i++ ) node[i] = 0;
	for( i=0; i<step-1; i++ )
		for( k=0; k<2; k++ ) 
			for( j=0; (s=topol[i][k][j])>-1; j++ ) 
				node[s]++;
	for( k=0; k<branch+1; k++ ) 
		for( j=0; (s=topol[step][k][j])>-1; j++ )
			node[s]++;
}

void RootLeafNode( int nseq, int ***topol, int *node )
{
	int i, j, k, s;
	for( i=0; i<nseq; i++ ) node[i] = 0;
	for( i=0; i<nseq-2; i++ )
		for( k=0; k<2; k++ ) 
			for( j=0; (s=topol[i][k][j])>-1; j++ ) 
				node[s]++;
}
		
void nodeFromABranch( int nseq, int *result, int **pairwisenode, int ***topol, double **len, int step, int num )
{
	int i, s, count;
	int *innergroup;
	int *outergroup1;
#if 0
	int outergroup2[nseq];
	int table[nseq];
#else
	static int *outergroup2 = NULL;
	static int *table = NULL;
	if( outergroup2 == NULL )
	{
		outergroup2 = AllocateIntVec( nseq );
		table = AllocateIntVec( nseq );
	}
#endif
	innergroup = topol[step][num];
	outergroup1 = topol[step][!num];
	for( i=0; i<nseq; i++ ) table[i] = 1;
	for( i=0; (s=innergroup[i])>-1; i++ ) table[s] = 0;
	for( i=0; (s=outergroup1[i])>-1; i++ ) table[s] = 0;
	for( i=0, count=0; i<nseq; i++ ) 
	{
		if( table[i] )
		{
			outergroup2[count] = i;
			count++;
		}
	}
	outergroup2[count] = -1;

	for( i=0; (s=innergroup[i])>-1; i++ )
	{
		result[s] = pairwisenode[s][outergroup1[0]]
				  + pairwisenode[s][outergroup2[0]]
				  - pairwisenode[outergroup1[0]][outergroup2[0]] - 1;
		result[s] /= 2;
	}
	for( i=0; (s=outergroup1[i])>-1; i++ ) 
	{
		result[s] = pairwisenode[s][outergroup2[0]]
				  + pairwisenode[s][innergroup[0]]
				  - pairwisenode[innergroup[0]][outergroup2[0]] + 1;
		result[s] /= 2;
	}
	for( i=0; (s=outergroup2[i])>-1; i++ ) 
	{
		result[s] = pairwisenode[s][outergroup1[0]]
				  + pairwisenode[s][innergroup[0]]
				  - pairwisenode[innergroup[0]][outergroup1[0]] + 1;
		result[s] /= 2;
	}

#if 0
	for( i=0; i<nseq; i++ ) 
		fprintf( stderr, "result[%d] = %d\n", i+1, result[i] );
#endif
}
		




	

		
		
					
					

				
void OneClusterAndTheOther_fast( int locnjob, int *memlist1, int *memlist2, int *s1, int *s2, char *pair, int ***topol, int step, int branch, double **smalldistmtx, double **distmtx, double *distontree )
{
	int i, k, j;
	int r1;
//	char *pair;

//	pair = calloc( locnjob, sizeof( char ) );

	for( i=0; i<locnjob; i++ ) pair[i] = 0;
    for( i=0, k=0; (r1=topol[step][branch][i])>-1; i++ ) 
	{
        pair[r1] = 1;
		memlist1[k++] = r1;
	}
	memlist1[k] = -1;

    for( i=0, k=0; i<locnjob; i++ ) 
    {
        if( !pair[i] ) 
        {
			memlist2[k++] = i;
        }
    }
	memlist2[k] = -1;

	*s1 = memlist1[0];
	*s2 = memlist2[0];

	if( smalldistmtx )
	{
		int im, jm;
#if USEDISTONTREE
		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			smalldistmtx[i][j] = distontree[im] + distontree[jm];
#else
		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			smalldistmtx[i][j] = distmtx[MIN(im,jm)][MAX(im,jm)];
#endif

#if 0
		reporterr( "\n" );
		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			reporterr( "smalldistmtx[%d][%d] = %f\n", i, j, smalldistmtx[i][j] );


		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			smalldistmtx[i][j] = distmtx[MIN(im,jm)][MAX(im,jm)];

		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			reporterr( "old smalldistmtx[%d][%d] = %f\n", i, j, smalldistmtx[i][j] );
if( im > 10 && jm > 10 ) exit( 1 );
#endif
	}
}
		

void makeEffMtx( int nseq, double **mtx, double *vec )
{
	int i, j;
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) 
		mtx[i][j] = vec[i] * vec[j];
}
	
void node_eff( int nseq, double *eff, int *node )
{
	int i;
	extern double ipower( double, int );
	for( i=0; i<nseq; i++ ) 
		eff[i] = ipower( 0.5, node[i] ) + geta2;
	/*
		eff[i] = ipower( 0.5, node[i] ) + 0;
	*/
#if DEBUG
	for( i=0; i<nseq; i++ ) 
#endif
}


int shrinklocalhom( char **pair, int s1, int s2, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int m1, k1, m2, k2;

	for( m1=s1, k1=0; m1<njob; m1++ )
	{
		if( pair[s1][m1] != 0 )
		{
			for( m2=s2, k2=0; m2<njob; m2++ )
			{
				if( pair[s2][m2] != 0 )
				{
					if( localhom[m1][m2].opt == -1 )
						localhomshrink[k1][k2] = NULL;
					else
						localhomshrink[k1][k2] = localhom[m1]+m2;
					k2++;
				}
			}
			k1++;
		}
	}
	return( 0 );
}
int msshrinklocalhom_fast_target( int *memlist1, int *memlist2, LocalHom **localhom, LocalHom ***localhomshrink, char *swaplist, int *targetmap )
{
	int m1, k1, m2, k2, t1, i2;

	for( k1=0; (m1=memlist1[k1])!=-1; k1++ )
	{
		if( targetmap[m1] == -1 )
		{	
			swaplist[k1] = 1;
//			swaplist[k1] = 0; // DAME!!!
			for( k2=0; (m2=memlist2[k2])!=-1; k2++ )
			{
				if( targetmap[m2] == -1 )
				{
					localhomshrink[k1][k2] = NULL;
					continue;
				}
		
				t1 = targetmap[m2]; // start1 <-> start2, end1 <-> end2
				i2 = m1;

				if( localhom[t1][i2].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[t1]+i2;
			}
		}
		else
		{
			swaplist[k1] = 0;
			for( k2=0; (m2=memlist2[k2])!=-1; k2++ )
			{
				t1 = targetmap[m1];
				i2 = m2;
	
				if( localhom[t1][i2].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[t1]+i2;
			}
		}
	}
	return( 0 );
}

int msshrinklocalhom_fast_half( int *memlist1, int *memlist2, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int m1, k1, m2, k2;

	for( k1=0; (m1=memlist1[k1])!=-1; k1++ )
	{
		for( k2=0; (m2=memlist2[k2])!=-1; k2++ )
		{
			if( m1 < m2 )
			{
				if( localhom[m1][m2-m1].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[m1]+m2-m1;
			}
			else
			{
				if( localhom[m2][m1-m2].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[m2]+m1-m2;
			}
		}
	}
	return( 0 );
}

int msshrinklocalhom_fast( int *memlist1, int *memlist2, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int m1, k1, m2, k2;

	for( k1=0; (m1=memlist1[k1])!=-1; k1++ )
	{
		for( k2=0; (m2=memlist2[k2])!=-1; k2++ )
		{
			if( localhom[m1][m2].opt == -1 )
				localhomshrink[k1][k2] = NULL;
			else
				localhomshrink[k1][k2] = localhom[m1]+m2;
		}
	}
	return( 0 );
}
int fastshrinklocalhom_one( int *mem1, int *mem2, int norg, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int k1, k2;
	int *intpt1, *intpt2;

	
	for( intpt1=mem1, k1=0; *intpt1!=-1; intpt1++, k1++ )
	{
		for( intpt2=mem2, k2=0; *intpt2!=-1; intpt2++, k2++ )
		{
			if( *intpt2 != norg ) 
			{
				fprintf( stderr, "ERROR! *intpt2 = %d\n", *intpt2 );
				exit( 1 );
			}
			if( localhom[*intpt1][0].opt == -1 )
				localhomshrink[k1][k2] = NULL;
			else
				localhomshrink[k1][k2] = localhom[*intpt1];
		}
	}
	return( 0 );
}

int fastshrinklocalhom( int *mem1, int *mem2, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int k1, k2;
	int *intpt1, *intpt2;

	
	for( intpt1=mem1, k1=0; *intpt1!=-1; intpt1++, k1++ )
	{
		for( intpt2=mem2, k2=0; *intpt2!=-1; intpt2++, k2++ )
		{
			if( localhom[*intpt1][*intpt2].opt == -1 )
				localhomshrink[k1][k2] = NULL;
			else
				localhomshrink[k1][k2] = localhom[*intpt1]+*intpt2;

//			if( localhomshrink[k1][k2] != NULL )
//				printf( "ori localhomshrink[%d][%d].opt = %f\n", k1, k2, localhomshrink[k1][k2]->opt );
		}
	}
	return( 0 );
}

int fastshrinklocalhom_half_seed( int *mem1, int *mem2, int nseed, int *posinlsh1, int *posinlsh2, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int k1, k2, sk1, sk2;
	int *intpt1, *intpt2;

	for( intpt1=mem1, k1=0, sk1=0; *intpt1!=-1; intpt1++, k1++ )
	{
		if( *intpt1>=nseed ) posinlsh1[k1] = -1;
		else posinlsh1[k1] = sk1++;
	}
	for( intpt2=mem2, k2=0, sk2=0; *intpt2!=-1; intpt2++, k2++ )
	{
		if( *intpt2>=nseed ) posinlsh2[k2] = -1;
		else posinlsh2[k2] = sk2++;
	}
	
	for( intpt1=mem1, sk1=0; *intpt1!=-1; intpt1++ )
	{
		if( *intpt1>=nseed ) continue;
		for( intpt2=mem2, sk2=0; *intpt2!=-1; intpt2++ )
		{
			if( *intpt2>=nseed ) continue;
			if( *intpt1 < *intpt2 )
			{
				if( localhom[*intpt1][*intpt2-*intpt1].opt == -1 )
					localhomshrink[sk1][sk2] = NULL;
				else
					localhomshrink[sk1][sk2] = localhom[*intpt1]+*intpt2-*intpt1;
			}
			else
			{
				if( localhom[*intpt2][*intpt1-*intpt2].opt == -1 )
					localhomshrink[sk1][sk2] = NULL;
				else
					localhomshrink[sk1][sk2] = localhom[*intpt2]+*intpt1-*intpt2;
			}

//			if( localhomshrink[k1][k2] != NULL )
//				printf( "ori localhomshrink[%d][%d].opt = %f, .importance=%f\n", k1, k2, localhomshrink[k1][k2]->opt, localhomshrink[k1][k2]->importance );
			sk2++;
		}
		sk1++;
	}
	return( 0 );
}

int fastshrinklocalhom_half( int *mem1, int *mem2, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int k1, k2;
	int *intpt1, *intpt2;

	
	for( intpt1=mem1, k1=0; *intpt1!=-1; intpt1++, k1++ )
	{
		for( intpt2=mem2, k2=0; *intpt2!=-1; intpt2++, k2++ )
		{
			if( *intpt1 < *intpt2 )
			{
				if( localhom[*intpt1][*intpt2-*intpt1].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[*intpt1]+*intpt2-*intpt1;
			}
			else
			{
				if( localhom[*intpt2][*intpt1-*intpt2].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[*intpt2]+*intpt1-*intpt2;
			}

//			if( localhomshrink[k1][k2] != NULL )
//				printf( "ori localhomshrink[%d][%d].opt = %f, .importance=%f\n", k1, k2, localhomshrink[k1][k2]->opt, localhomshrink[k1][k2]->importance );
		}
	}
	return( 0 );
}


int fastshrinklocalhom_target( int *mem1, int *mem2, LocalHom **localhom, LocalHom ***localhomshrink, char *swaplist, int *targetmap )
{
	int k1, k2;
	int *intpt1, *intpt2;
	int t1, i2;

	
	for( intpt1=mem1, k1=0; *intpt1!=-1; intpt1++, k1++ )
	{
		if( targetmap[*intpt1] == -1 )
		{
			swaplist[k1] = 1;
//			swaplist[k1] = 0; // DAME!!!
			for( intpt2=mem2, k2=0; *intpt2!=-1; intpt2++, k2++ )
			{
				if( targetmap[*intpt2] == -1 )
				{
					localhomshrink[k1][k2] = NULL;
					continue;
				}
	
				t1 = targetmap[*intpt2]; // end1<->end2, start1<->start2
				i2 = *intpt1;
	
				if( localhom[t1][i2].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[t1]+i2;
	
//				if( localhomshrink[k1][k2] != NULL )
//					printf( "localhomshrink[%d][%d].opt = %f\n", k1, k2, localhomshrink[k1][k2]->opt );
//				else
//					printf( "localhomshrink[%d][%d] = NULL\n", k1, k2 );
			}
		}
		else
		{
			swaplist[k1] = 0;
			for( intpt2=mem2, k2=0; *intpt2!=-1; intpt2++, k2++ )
			{
				t1 = targetmap[*intpt1];
				i2 = *intpt2;
	
				if( localhom[t1][i2].opt == -1 )
					localhomshrink[k1][k2] = NULL;
				else
					localhomshrink[k1][k2] = localhom[t1]+i2;
	
//				if( localhomshrink[k1][k2] != NULL )
//					printf( "localhomshrink[%d][%d].opt = %f\n", k1, k2, localhomshrink[k1][k2]->opt );
//				else
//					printf( "localhomshrink[%d][%d] = NULL\n", k1, k2 );
			}
		}
	}
	return( 0 );
}

int msfastshrinklocalhom( int *mem1, int *mem2, LocalHom **localhom, LocalHom ***localhomshrink )
{
	int k1, k2;
	int *intpt1, *intpt2;
	int m1, m2;
	
	for( intpt1=mem1, k1=0; *intpt1!=-1; intpt1++, k1++ )
	{
		for( intpt2=mem2, k2=0; *intpt2!=-1; intpt2++, k2++ )
		{
			m1 = MIN(*intpt1,*intpt2); m2 = MAX(*intpt1,*intpt2);
			if( localhom[m1][m2].opt == -1 )
				localhomshrink[k1][k2] = NULL;
			else
				localhomshrink[k1][k2] = localhom[m1]+m2;
		}
	}
	return( 0 );
}

