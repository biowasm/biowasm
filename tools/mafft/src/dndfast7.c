#include "mltaln.h"
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

	inputfile = NULL;
	disopt = 0;
	divpairscore = 0;
	swopt = "";

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'I':
					disopt = 1;
					break;
				case 'A':
					swopt = "-A";
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
	FILE *hat3p;
	FILE *infp;
	char **seq = NULL; // by D.Mathog
	char **seq1;
	char **name;
	char **name1;
	static int nlen1[M];
	double **mtx;
	double **mtx2;
	static int nlen[M];
	static char b[B];
	double max;
	char com[1000];
	int opt[M];
	int res;
	char *home;
	char *fastapath;
	char queryfile[B];
	char datafile[B];
	char fastafile[B];
	char hat2file[B];
	int pid = (int)getpid();
	LocalHom **localhomtable, *tmpptr;
#if 0
	home = getenv( "HOME" );
#else /* $HOME wo tsukau to fasta ni watasu hikisuu ga afureru */ 
	home = NULL;
#endif
	fastapath = getenv( "FASTA_4_MAFFT" );
	if( !fastapath ) 
		fastapath = "fasta34";

#if DEBUG
	if( home ) fprintf( stderr, "home = %s\n", home );
#endif
	if( !home ) home = "";
	sprintf( queryfile, "%s/tmp/query-%d", home, pid );
	sprintf( datafile, "%s/tmp/data-%d", home, pid );
	sprintf( fastafile, "%s/tmp/fasta-%d", home, pid );
	sprintf( hat2file, "hat2-%d", pid );


	arguments( argc, argv );

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



#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	dorp = NOTSPECIFIED;
	getnumlen( infp );
#endif

	if( dorp == 'd' )
	{
		scoremtx = -1;
		pamN = NOTSPECIFIED;
	}
	else
	{
		nblosum = 62;
		scoremtx = 1;
	}
	constants( njob, seq );

	rewind( infp );

	name = AllocateCharMtx( njob, B+1 );
	name1 = AllocateCharMtx( njob, B+1 );
	seq = AllocateCharMtx( njob, nlenmax+1 );
	seq1 = AllocateCharMtx( 2, nlenmax+1 );
	mtx = AllocateDoubleMtx( njob, njob );
	mtx2 = AllocateDoubleMtx( njob, njob );
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
			localhomtable[i][j].opt = -1.0;
			localhomtable[i][j].next = NULL;
		}
    }

#if 0
	FRead( infp, name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
#endif
	fclose( infp );

	if( scoremtx == -1 ) ktuple = 6;
	else                 ktuple = 1;

	for( i=0; i<njob; i++ )
	{
		gappick0( seq1[0], seq[i] ); 
		strcpy( seq[i], seq1[0] );
	}
	for( j=0; j<njob; j++ )
	{
		sprintf( name1[j], "+==========+%d                      ", j );
		nlen1[j] = nlen[j];
	}
	hat2p = fopen( datafile, "w" );
	if( !hat2p ) ErrorExit( "Cannot open datafile." );
	WriteForFasta( hat2p, njob, name1, nlen1, seq );
	fclose( hat2p );

	for( i=0; i<njob; i++ ) 
	{
//		fprintf( stderr, "###  i = %d\n", i );
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
			sprintf( com, "%s %s -z3 -m10  -n -Q  -b%d -E%d -d%d %s %s %d > %s", fastapath, swopt, M, M, M, queryfile, datafile, ktuple, fastafile );
		else
			sprintf( com, "%s %s -z3 -m10  -Q  -b%d -E%d -d%d %s %s %d > %s", fastapath, swopt, M, M, M, queryfile, datafile, ktuple, fastafile );
		res = system( com );
		if( res ) ErrorExit( "error in fasta" );



		hat2p = fopen( fastafile, "r" );
		if( hat2p == NULL ) 
			ErrorExit( "file 'fasta.$$' does not exist\n" );
		if( scoremtx == -1 )
			res = ReadFasta34m10_nuc( hat2p, mtx[i], i, name1, localhomtable[i] );
		else
			res = ReadFasta34m10( hat2p, mtx[i], i, name1, localhomtable[i] );
		fclose( hat2p );

		if( res < njob - i )
		{
			fprintf( stderr, "count (fasta34 -z 3) = %d\n", res );
			exit( 1 );
		}


		if( i == 0 )
			for( j=0; j<njob; j++ ) opt[j] = (int)mtx[0][j];


#if 0
		{
			int ii, jj;
			if( i < njob-1 ) for( jj=i; jj<i+5; jj++ ) 
				fprintf( stdout, "mtx[%d][%d] = %f\n", i+1, jj+1, mtx[i][jj] );
		}
#endif
		fprintf( stderr, "query : %4d / %5d\r", i+1, njob );
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
//				fprintf( stderr, "##### mtx[%d][%d] = %f\n", i, j, mtx[i][j] );
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
		sprintf( name[0], "=query====lgth=%04d-%04d %.*s", nlen[0], howmanyx( seq[0] ), B-30, b );
#if 0
		strins(  b, name[0] );
#endif
		for( i=1; i<njob; i++ ) 
		{	
			strcpy( b, name[i] );
			sprintf( name[i], "=opt=%04d=lgth=%04d-%04d %.*s", opt[i], nlen[i], howmanyx( seq[i] ), B-30, b );
#if 0
			strins( b, name[i] );
#endif
		}
	}

	hat2p = fopen( hat2file, "w" );
	if( !hat2p ) ErrorExit( "Cannot open hat2." );
	WriteHat2_pointer( hat2p, njob, name, mtx2 );
	fclose( hat2p );

#if 1
	fprintf( stderr, "##### writing hat3\n" );
	hat3p = fopen( "hat3", "w" );
	if( !hat3p ) ErrorExit( "Cannot open hat3." );
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
	{
		for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->opt == -1.0 ) continue;
			fprintf( hat3p, "%d %d %d %6.3f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, (void *)tmpptr->next );
		}
	}
	fclose( hat3p );
#endif

	sprintf( com, "/bin/rm %s %s %s", queryfile, datafile, fastafile );
	system( com );

#if 0
	sprintf( com, ALNDIR "/supgsdl < %s", hat2file );
	res = system( com );
	if( res ) ErrorExit( "error in spgsdl" );
#endif

	sprintf( com, "mv %s hat2", hat2file );
	res = system( com );
	if( res ) ErrorExit( "error in mv" );

	SHOWVERSION;
	exit( 0 );
}
