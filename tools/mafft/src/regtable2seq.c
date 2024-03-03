#include "mltaln.h"

#define DEBUG 0

char *regfile;
char *eregfile;

void arguments( int argc, char *argv[] )
{
    int c;

    outnumber = 0;
    inputfile = NULL;
    regfile = NULL;
    eregfile = NULL;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'e':
					eregfile = *++argv;
					fprintf( stderr, "eregfile = %s\n", eregfile );
					--argc;
					goto nextoption;
				case 'r':
					regfile = *++argv;
					fprintf( stderr, "regfile = %s\n", regfile );
					--argc;
					goto nextoption;
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'n' :
					outnumber = 1;
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
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
}

void readereg( FILE *regfp, int **regtable, char **revtable, int *outtable, int *noutpt, int *loutpt )
{
	char gett[1000];
	int j;
	int mem;
	char cmem;
	char reg[5][100];
	char out[100];
	int startpos, endpos;

	*noutpt = 0;
	*loutpt = 0;
	fgets( gett, 999, regfp );
	reg[0][0] = reg[1][0] = reg[2][0] = reg[3][0] = reg[4][0] = 'n';
	reg[0][1] = reg[1][1] = reg[2][1] = reg[3][1] = reg[4][1] = 0;
	sscanf( gett, "%c %s %s %s %s %s", &cmem, reg[0], reg[1], reg[2], reg[3], reg[4] );
	if( cmem != 'e' )
	{
		fprintf( stderr, "Format error\n" );
		exit( 1 );
	}
	for( j=0; j<5; j++ )
	{
//		reporterr( "reg[j]=%s\n", reg[j] );
		sscanf( reg[j], "%d-%d-%c", regtable[0]+(j*2), regtable[0]+(j*2)+1, revtable[0]+j );
		fprintf( stderr, "%d %d-%d\n", 0, regtable[0][j*2], regtable[0][j*2+1] );
		startpos = regtable[0][j*2];
		endpos   = regtable[0][j*2+1];
//		reporterr( "startpod=%d, endpos=%d, *loutpt=%d\n", startpos, endpos, *loutpt );
		if( startpos > endpos )
		{
			endpos   = regtable[0][j*2];
			startpos = regtable[0][j*2+1];
		}
		if( startpos != -1 && endpos != -1 )
			*loutpt += endpos - startpos + 1;
	}

	while( 1 )
	{
		fgets( gett, 999, regfp );
		if( feof( regfp ) ) break;
		sscanf( gett, "%d o=%s", &mem, out );
		if( mem >= njob )
		{
			fprintf( stderr, "Out of range\n" );
			exit( 1 );
		}
		outtable[mem] = atoi( out );
		if( outtable[mem] ) *noutpt += 1;
	}
}

void readreg( FILE *regfp, int **regtable, char **revtable, int *outtable )
{
	char gett[1000];
	int j;
	int mem;
	char reg[5][100];
	char out[100];

	while( 1 )
	{
		fgets( gett, 999, regfp );
		if( feof( regfp ) ) break;
		sscanf( gett, "%d %s %s %s %s %s o=%s", &mem, reg[0], reg[1], reg[2], reg[3], reg[4], out );
		if( mem >= njob )
		{
			fprintf( stderr, "Out of range\n" );
			exit( 1 );
		}
		for( j=0; j<5; j++ )
		{
			sscanf( reg[j], "%d-%d-%c", regtable[mem]+(j*2), regtable[mem]+(j*2)+1, revtable[mem]+j );
			fprintf( stderr, "%d %d-%d\n", mem, regtable[mem][j*2], regtable[mem][j*2+1] );
		}
		outtable[mem] = atoi( out );
	}
}

int main( int argc, char *argv[] )
{
	FILE *infp;
	FILE *regfp;
	int nlenmin;
	int **regtable;
	char **revtable;
	int *outtable;
	int i, nout, lout;
	char **outseq;
	char **name;

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

	dorp = NOTSPECIFIED;
	getnumlen_nogap( infp, &nlenmin );

	if( regfile )
	{
		regfp = fopen( regfile, "r" );
		if( !regfp )
		{
			fprintf( stderr, "Cannot open %s\n", regfile );
			exit( 1 );
		}
		regtable = AllocateIntMtx( njob, 5*2 );
		revtable = AllocateCharMtx( njob, 5 );
		outtable = AllocateIntVec( njob );
		readreg( regfp, regtable, revtable, outtable );
		cutData( infp, regtable, revtable, outtable );
	}
	else if( eregfile )
	{
		regfp = fopen( eregfile, "r" );
		if( !regfp )
		{
			fprintf( stderr, "Cannot open %s\n", eregfile );
			exit( 1 );
		}
		regtable = AllocateIntMtx( 1, 5*2 );
		revtable = AllocateCharMtx( 1, 5 );
		outtable = AllocateIntVec( njob );
		readereg( regfp, regtable, revtable, outtable, &nout, &lout );
		fprintf( stderr, "nout = %d, lout = %d\n", nout, lout );

		outseq = AllocateCharMtx( nout, lout+1 );
		name = AllocateCharMtx( nout, B );

		cutAlignment( infp, regtable, revtable, outtable, name, outseq );
		fprintf( stderr, "gappick! nout = %d\n", nout );
		commongappick( nout, outseq );
		for( i=0; i<nout; i++ )
		{
			fprintf( stdout, "%s\n", name[i] );
			fprintf( stdout, "%s\n", outseq[i] );
		}
	}
	else
	{
		catData( infp );
	}

	fprintf( stderr, "Strategy:\n" );
	fprintf( stderr, " Not-Aligned\n" );

//	fprintf( stdout, "%d x %d - %d %c\n", njob, nlenmax, nlenmin, dorp );
	return( 0 );
}
