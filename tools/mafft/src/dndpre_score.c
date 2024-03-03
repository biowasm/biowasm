#include "mltaln.h"

#define TEST 0

int main()
{
	int i, j;
	char **seq;
	static char name[M][B];
	static int nlen[M];
	double **mtx;
	FILE *fp;
	int res;

	scoremtx = NOTSPECIFIED;

#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	getnumlen( stdin );
#endif
	rewind( stdin );

	seq = AllocateCharMtx( njob, nlenmax+1 );
	mtx = AllocateDoubleMtx( njob, njob );

#if 0
	FRead( stdin, name, nlen, seq );
#else
	readData( stdin, name, nlen, seq );
#endif

	for( i=0; i<njob-1; i++ ) 
	{
		fprintf( stderr, "%4d/%4d\r", i+1, njob );
		for( j=i+1; j<njob; j++ ) 
			mtx[i][j] = (double)substitution_score( seq[i], seq[j] );
	}
	
#if TEST
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		fprintf( stdout, "i=%d, j=%d, mtx[][] = %f\n", i, j, mtx[i][j] );
#endif

	fp = fopen( "hat2", "w" );
	WriteHat2( fp, njob, name, mtx );
	fclose( fp );
	exit( 0 );
/*
	res = system( ALNDIR "/spgsdl < hat2"  );
	if( res ) exit( 1 );
	else exit( 0 );
*/
}
