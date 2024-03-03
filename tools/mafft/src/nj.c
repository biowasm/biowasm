#include "mltaln.h"
#define DEBUG 0


void topolcpy( int s1[], int s2[], int *mpt1, int *mpt2 )
{
    int i;

    *mpt1 = *mpt2;
    for( i=0; i<*mpt2; i++ )
    {
        s1[i] = s2[i];
    }
}

void topolcat( int s1[], int s2[], int *mpt1, int *mpt2 )
{
    int i;

    for( i=*mpt1; i<*mpt1+*mpt2; i++ )
    {
        s1[i] = s2[i-*mpt1];
    }
    *mpt1 += *mpt2;
}
   
void topolsort( int m, int s[] )
{
    int i, j, im;
    int sm;

    for( j=0; j<m-1; j++ )
    {
        sm = s[j]; im = j;
        for( i=j+1; i<m; i++ )
        {
            if( s[i] < sm )
            {
                sm = s[i];
                im = i;
            }
        }
        s[im] = s[j]; s[j] = sm;
    }
}

void topolswap( int s1[], int s2[], int *mpt1, int *mpt2 )
{
    int i;
    int im;
    int b;
    b = *mpt1; *mpt1 = *mpt2; *mpt2 = b;
    im = MAX(*mpt1,*mpt2);
    for( i=0; i<im; i++ )
    {
        b = s1[i]; s1[i] = s2[i]; s2[i] = b;
    /*
    printf( "s1[%d]=%d\ns2[%d]=%d\n", i, s1[i], i, s2[i] );
    */
    }
}

void reduc( double **mtx, int nseq, int im, int jm )
{
    int i;
    for( i=0; i<nseq; i++ )
    {
        if(    i==im || i==jm
            || mtx[MIN(i,im)][MAX(i,im)] == 9999.9
            || mtx[MIN(i,jm)][MAX(i,jm)] == 9999.9
          ) continue;
        mtx[MIN(i,im)][MAX(i,im)]
        = 0.5 * ( mtx[MIN(i,im)][MAX(i,im)] + mtx[MIN(i,jm)][MAX(i,jm)]
                  - mtx[MIN(im,jm)][MAX(im,jm)] );
        mtx[MIN(i,jm)][MAX(i,jm)] = 9999.9;
    }
    mtx[MIN(im,jm)][MAX(im,jm)] = 9999.9;
}


void  nj( int nseq, double **omtx, int ***topol, double **dis )
{
    int i, j, l, n, m;
    int count;
    double r[M];
    double t;
    double s, sm;
    double totallen = 0.0;
    int im=0, jm=0;
    double len1, len2;
#if 1
	static char **par = NULL;
	static double **mtx = NULL;
	static int **mem = NULL;
	if( par == NULL )
	{
		par = AllocateCharMtx( njob, njob );
		mtx = AllocateDoubleMtx( njob, njob );
		mem = AllocateIntMtx( njob, 2 );
	}
#else
    char par[njob][njob];
	double mtx[njob][njob];
	int mem[njob][2];
#endif
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) mtx[i][j] = omtx[i][j];
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) par[i][j] = 0;
    for( i=0; i<nseq; i++ ) par[i][i] = 1;
//	for( i=0; i<nseq; i++ ) for( j=0; j<2; j++ ) for( l=0; l<nseq+1; l++ ) topol[i][j][l] = -1;
	for( i=0; i<nseq; i++ ) for( j=0; j<2; j++ ) for( l=0; l<nseq; l++ ) topol[i][j][l] = -1;
    for( n=nseq, m=0; n>2; n--, m=nseq-n )
    {
        t = 0.0;
        for( i=0; i<nseq-1; i++ ) for( j=0; j<nseq; j++ ) if( mtx[i][j] < 9999.9 )
            t += mtx[i][j];
        for( i=0; i<nseq; i++ )
        {
            r[i] = 0.0;
            for( l=0; l<nseq; l++ ) 
                if( ( l != i ) && ( mtx[MIN(i,l)][MAX(i,l)] < 9999.9 ) )
                    r[i] += mtx[MIN(i,l)][MAX(i,l)];
        }
        sm = 9999.9;
        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) if( mtx[i][j] < 9999.9)
        {
            s = ( ( 2.0 * t - r[i] - r[j] + (n-2.0)*mtx[i][j] ) ) / ( 2.0*(n-2.0) );
            if ( s < sm )
            {
                sm = s;
                im = i; jm = j;
            }
        }
        len1 = ( (n-2)*mtx[im][jm] + r[im] - r[jm] ) / (2*(n-2));
        len2 = ( (n-2)*mtx[im][jm] - r[im] + r[jm] ) / (2*(n-2));

#if DEBUG
        fprintf( stderr, "STEP-%3d  %3d: L = %5.5f\n", m+1, im+1, len1 );
        fprintf( stderr, "          %3d: L = %5.5f\n",      jm+1, len2 );
#endif

        totallen += len1;
        totallen += len2;

        dis[m][0] = len1;
        dis[m][1] = len2;

        for( l=0, count=0; l<nseq; l++ )
            if( par[im][l] > 0 )
            {
                topol[m][0][count] = l;
                count++;
            }
        mem[m][0] = count;
        for( l=0, count=0; l<nseq; l++ )
            if( par[jm][l] > 0 )
            {
                topol[m][1][count] = l;
                count++;
            }
        mem[m][1] = count;
        for( l=0; l<nseq; l++ )
            par[im][l] += ( par[jm][l] > 0 );
        if( n > 3 ) reduc( mtx, nseq, im, jm );
    }
    for( i=0; i<nseq; i++ )
        if( i!=im && i!=jm && mtx[MIN(i,im)][MAX(i,im)]<9999.9 )
            break;
    len2 = ( mtx[MIN(i,im)][MAX(i,im)] - r[im] + r[i] ) / 2;

/*
    printf("          %3d: L = %5.5f\n", i+1, len2 );
*/
    totallen += len2;

    dis[m][0] = len2;
    dis[m][1] = 0.0;
    for( l=0, count=0; l<nseq; l++ )
        if( par[i][l] > 0 )
        {
            topol[m][0][count] = l;
            count++;
        }
    mem[m][0] = count;
    /*
    printf( " total length == %f\n", totallen );
    */

    topolcpy( topol[nseq-2][1], topol[nseq-3][0], mem[nseq-2]+1, mem[nseq-3] );
    topolcat( topol[nseq-2][1], topol[nseq-3][1], mem[nseq-2]+1, mem[nseq-3]+1 );
    topolsort( mem[nseq-2][1], topol[nseq-2][1] );
	
	if( topol[nseq-2][0][0] > topol[nseq-2][1][0] )
		topolswap( topol[nseq-2][0], topol[nseq-2][1], mem[nseq-2], mem[nseq-2]+1 );

}
