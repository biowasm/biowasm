#include "mltaln.h"

#define DEBUG 0

char *directionfile;
static int show_R_ = 1;

static int subalignment;
static int subalignmentoffset;

void arguments( int argc, char *argv[] )
{
    int c;

	inputfile = NULL;
	directionfile = NULL;
	subalignment = 0;
	subalignmentoffset = 0;
	show_R_ = 1;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'd':
					directionfile = *++argv;
					fprintf( stderr, "directionfile = %s\n", directionfile );
					--argc;
					goto nextoption;
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'r':
					show_R_ = 0;
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



int main( int argc, char *argv[] )
{
	FILE *infp;
	FILE *difp;
	int nlenmin;
	char **name;
	char **seq;
	char *tmpseq;
	char line[100];
	int *nlen;
	int i, j;
	int nsubalignments, maxmem;
	int **subtable = NULL;
	int *preservegaps = NULL;
	char firstdir;
	char *directions;

	arguments( argc, argv );

	reporterr( "subalignment = %d\n", subalignment );
	reporterr( "subalignmentoffset = %d\n", subalignmentoffset );


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

	if( directionfile )
	{
		difp = fopen( directionfile, "r" );
		if( !difp )
		{
			fprintf( stderr, "Cannot open %s\n", directionfile );
			exit( 1 );
		}
	}
	else
	{
		fprintf( stderr, "Give directionfile!\n" );
	}


	dorp = NOTSPECIFIED;
	getnumlen_casepreserve( infp, &nlenmin );

	fprintf( stderr, "%d x %d - %d %c\n", njob, nlenmax, nlenmin, dorp );

	seq = AllocateCharMtx( njob, nlenmax+1 );
	tmpseq = AllocateCharVec( MAX( B, nlenmax )+1 );
	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	directions = calloc( njob, sizeof( int ) );

	readData_pointer_casepreserve( infp, name, nlen, seq );





	for( i=0; i<njob; i++ )
	{
		fgets( line, 99, difp );
		if( line[0] != '_' )
		{
			fprintf( stderr, "Format error!\n" );
			exit( 1 );
		}
		if( line[1] == 'R' )
		{

			directions[i] = 'R';
			sreverse( tmpseq, seq[i] );
			strcpy( seq[i], tmpseq );

			strncpy( tmpseq, name[i]+1, B-3 );
			tmpseq[B-3] = 0;
			if( show_R_ )
			{
				strcpy( name[i]+1, "_R_" );
				strcpy( name[i]+4, tmpseq );
			}
			else
			{
				strcpy( name[i]+1, tmpseq );
			}
		}
		else if( line[1] == 'F' )
		{
			directions[i] = 'F';
		}
		else
		{
			fprintf( stderr, "Format error!\n" );
			exit( 1 );
		}
	}

	if( subalignment )
	{
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem );
		reporterr(       "nsubalignments = %d\n", nsubalignments );
		reporterr(       "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		preservegaps = AllocateIntVec( njob );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );

		for( j=0; j<nsubalignments; j++ ) 
		{
			reporterr( "Checking directions of sequences in subalignment%d\n", j );
			firstdir = directions[subtable[j][0]];
			reporterr( "firstdir = %c\n", firstdir );
			for( i=0; subtable[j][i]>-1; i++ )
			{
				if( directions[subtable[j][i]] != firstdir )
				{
					reporterr( "\n\n#############################################################################\n" );
					reporterr( "\nDirection of nucleotide sequences seems to be inconsistent.\n" );
					reporterr( "Please check the following two sequences:\n" );
					reporterr( "	Sequece no.%d (%s)\n", subtable[j][0]+1, name[subtable[j][0]] );
					reporterr( "	Sequece no.%d (%s)\n", subtable[j][i]+1, name[subtable[j][i]] );
					reporterr( "\nThese sequences are in sub alignment no.%d in your setting of --merge,\nbut their directions seem to be different.\n\n", j+1 );
					reporterr( "#############################################################################\n\n\n\n" );
					exit( 1 );
				}
			}
			reporterr( "OK!\n" );
		}
	}


	for( i=0; i<njob; i++ )
	{
		fprintf( stdout, ">%s\n", name[i]+1 );
		fprintf( stdout, "%s\n", seq[i] );
	}

	free( nlen );
	FreeCharMtx( seq );
	FreeCharMtx( name );
	free( tmpseq );

	return( 0 );
}
