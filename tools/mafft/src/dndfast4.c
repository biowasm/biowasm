#include "mltaln.h"
#include <sys/types.h>
#include <unistd.h>
#define DEBUG 0
#define TEST 0


int howmanyx( char *s )
{
	int val = 0;
	if( scoremtx == -1 )
	{
		do
		{
			if( !strchr( "atgcuATGCU", *s ) ) val++;
		} while( *++s );
	}
	else
	{
		do
		{
			if( !strchr( "ARNDCQEGHILKMFPSTWYV", *s ) ) val++;
		} while( *++s );
	}
	return( val );
}

void arguments( int argc, char *argv[] )
{
	int c;

	disopt = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
        while ( c = *++argv[0] )
            switch( c )
            {
				case 'i':
					disopt = 1;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
    if( argc != 0 )
    {
        fprintf( stderr, "options: -i\n" );
        exit( 1 );
    }
}

int main( int argc, char *argv[] )
{
	int ktuple;
	int i, j;
	FILE *hat2p;
	char **seq;
	char **seq1;
	static char name[M][B];
	static char name1[M][B];
	static int nlen1[M];
	double **mtx;
	double **mtx2;
	static int nlen[M];
	char b[B];
	double max;
	char com[B];
	int opt[M];
	int res;
	char *home;
	char queryfile[B];
	char datafile[B];
	char fastafile[B];
	char hat2file[B];
	int pid = (int)getpid();
#if 0
	home = getenv( "HOME" );
#else /* $HOME wo tsukau to fasta ni watasu hikisuu ga afureru */ 
	home = NULL;
#endif

#if DEBUG
	if( home ) fprintf( stderr, "home = %s\n", home );
#endif
	if( !home ) home = "";
	sprintf( queryfile, "%s/tmp/query-%d\0", home, pid );
	sprintf( datafile, "%s/tmp/data-%d\0", home, pid );
	sprintf( fastafile, "%s/tmp/fasta-%d\0", home, pid );
	sprintf( hat2file, "hat2-%d\0", pid );

	arguments( argc, argv );
#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	getnumlen( stdin );
#endif
	rewind( stdin );

	seq = AllocateCharMtx( njob, nlenmax+1 );
	seq1 = AllocateCharMtx( 2, nlenmax+1 );
	mtx = AllocateDoubleMtx( njob, njob );
	mtx2 = AllocateDoubleMtx( njob, njob );

#if 0
	FRead( stdin, name, nlen, seq );
#else
	readData( stdin, name, nlen, seq );
#endif
	if( scoremtx == -1 ) ktuple = 6;
	else                 ktuple = 1;

	for( i=0; i<njob; i++ )
	{
		gappick0( seq1[0], seq[i] ); 
		strcpy( seq[i], seq1[0] );
	}
	for( j=0; j<njob; j++ )
	{
		sprintf( name1[j], "+==========+%d                      \0", j );
		nlen1[j] = nlen[j];
	}
	hat2p = fopen( datafile, "w" );
	if( !hat2p ) ErrorExit( "Cannot open datafile." );
	WriteForFasta( hat2p, njob, name1, nlen1, seq );
	fclose( hat2p );

	for( i=0; i<njob; i++ ) 
	{

		hat2p = fopen( datafile, "w" );
		if( !hat2p ) ErrorExit( "Cannot open datafile." );
		WriteForFasta( hat2p, njob-i, name1+i, nlen1+i, seq+i );
		fclose( hat2p );


		seq1[0] = seq[i];
		nlen1[0] = nlen[i];

		hat2p = fopen( queryfile, "w" );
		if( !hat2p ) ErrorExit( "Cannot open queryfile." );
		WriteForFasta( hat2p, 1, name1+i, nlen1, seq1 ); 
		fclose( hat2p );

		if( scoremtx == -1 )
			sprintf( com, "fasta3 -n -Q -h -b%d -E%d -d%d %s %s %d > %s\0", M, M, 0, queryfile, datafile, ktuple, fastafile );
		else
			sprintf( com, "fasta3 -Q -h -b%d -E%d -d%d %s %s %d > %s\0", M, M, 0, queryfile, datafile, ktuple, fastafile );
		res = system( com );
		if( res ) ErrorExit( "error in fasta" );

		hat2p = fopen( fastafile, "r" );
		if( hat2p == NULL ) 
			ErrorExit( "file 'fasta.$$' does not exist\n" );
		ReadFasta3( hat2p, mtx[i], njob-i, name1 );

		if( i == 0 )
			for( j=0; j<njob; j++ ) opt[j] = (int)mtx[0][j];

		fclose( hat2p );

#if 1
		{
			int ii, jj;
			if( i < njob-1 ) for( jj=i; jj<i+5; jj++ ) 
				fprintf( stdout, "mtx[%d][%d] = %f\n", i+1, jj+1, mtx[i][jj] );
		}
#endif
		fprintf( stderr, "query : %#4d\n", i+1 );
	}

	for( i=0; i<njob; i++ )
	{
		max = mtx[i][i];
		if( max == 0.0 )
		{
			for( j=0; j<njob; j++ )
				mtx2[i][j] = 2.0;
		}
		else
		{
			for( j=0; j<njob; j++ )
			{
				mtx2[i][j] = ( max - mtx[MIN(i,j)][MAX(i,j)] ) / max * 2.0;
//				fprintf( stdout, "max = %f, mtx[%d][%d] = %f -> %f\n", max, i+1, j+1, mtx[i][j], mtx2[i][j] );
			}
		}
	}
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
	{
//		fprintf( stdout, "mtx2[%d][%d] = %f, %f\n", i+1, j+1, mtx2[i][j], mtx2[j][i] );
		mtx2[i][j] = MIN( mtx2[i][j], mtx2[j][i] );
	}
#if 0
	{
		int ii, jj;
		if( i < njob-1 ) for( jj=i+1; jj<njob; jj++ ) 
			fprintf( stderr, "mtx2[][] = %f\n", mtx2[i][jj] );
	}
#endif

	for( i=0; i<njob; i++ ) name[i][0] = '=';

	if( disopt )
	{
		strcpy( b, name[0] );
		sprintf( name[0], "=query====lgth=%#04d-%04d %.*s\0", nlen[0], howmanyx( seq[0] ), B-30, b );
#if 0
		strins(  b, name[0] );
#endif
		for( i=1; i<njob; i++ ) 
		{	
			strcpy( b, name[i] );
			sprintf( name[i], "=opt=%#04d=lgth=%#04d-%04d %.*s\0", opt[i], nlen[i], howmanyx( seq[i] ), B-30, b );
#if 0
			strins( b, name[i] );
#endif
		}
	}

	hat2p = fopen( hat2file, "w" );
	if( !hat2p ) ErrorExit( "Cannot open hat2." );
	WriteHat2( hat2p, njob, name, mtx2 );
	fclose( hat2p );

	sprintf( com, "/bin/rm %s %s %s", queryfile, datafile, fastafile );
	system( com );

#if 0
	sprintf( com, ALNDIR "/supgsdl < %s\0", hat2file );
	res = system( com );
	if( res ) ErrorExit( "error in spgsdl" );
#endif

	sprintf( com, "mv %s hat2", hat2file );
	res = system( com );
	if( res ) ErrorExit( "error in mv" );

	SHOWVERSION;
	exit( 0 );
}
