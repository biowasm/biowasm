#include "mltaln.h"

#define DEBUG 0

double maxunusual;

static double count_unusual( char *seq, char *usual )
{
	int i;
	char *pt;
	int count, len;
	count = 0;
	pt = seq;
	while( *pt )
	{
		if( !strchr( usual, *pt ) ) 
			count++;
		pt++;
	}
//	reporterr( "%d/%d=%f\n", count, pt-seq, ((double)count/(pt-seq)) );
	return( (double)count / (pt-seq) );
}

static void shortenN( char *seq, char unknown )
{
	int i;
	int status;
	char *out = seq;
	int unknownU = toupper(unknown);
	status = 0;
	while( *seq )
	{
		if( unknownU != toupper(*seq) )  // hikouritsu?
		{
			*out++ = *seq++;
			status = 0;
		}
		else if( status == 0 )
		{
			*out++ = unknown;
			seq++;
			status = 1;
		}
		else
			seq++;
	}
	*out = 0;
}


void arguments( int argc, char *argv[] )
{
    int c;

	maxunusual = 0.05;
	inputfile = NULL;
	dorp = NOTSPECIFIED;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'm':
					maxunusual = myatof( *++argv );
					fprintf( stderr, "maxunusual = %f\n", maxunusual );
					--argc;
					goto nextoption;
				case 'i':
					inputfile = *++argv;
//					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
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
	int nlenmin;
	char **name;
	char **seq;
	int *nlen;
	int i;
	char *usual;
	char unknown;
	int nout;
	char *tmpseq;

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


//	dorp = NOTSPECIFIED;
	getnumlen_casepreserve( infp, &nlenmin );

//	fprintf( stderr, "%d x %d - %d %c\n", njob, nlenmax, nlenmin, dorp );

	seq = AllocateCharMtx( njob, nlenmax+1 );
	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	tmpseq = AllocateCharVec( nlenmax+1 );

//	readData_pointer( infp, name, nlen, seq );
	readData_pointer_casepreserve( infp, name, nlen, seq );
	fclose( infp );

//	for( i=0; i<njob; i++ ) gappick_samestring( seq[i] );

#if 0
	FILE *origfp;
	origfp = fopen( "_original", "w" );
	if( !origfp )
	{
		fprintf( stderr, "Cannot open _original\n" );
		exit( 1 );
	}
	for( i=0; i<njob; i++ )
	{
		nlen[i] = strlen( seq[i] );
		fprintf( origfp, ">%s\n", name[i]+1 );
		if( seq[i][nlen[i]-1] == '\n' ) seq[i][nlen[i]-1] = 0;
		fprintf( origfp, "%s\n", seq[i] );
	}
	fclose( origfp );
#endif

	if( dorp == 'p' )
	{
		usual = "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-";
		unknown = 'X';
	}
	else
	{
		usual = "ATGCUatgcu-";
		unknown = 'n';
	}
	nout = 0;
	for( i=0; i<njob; i++ )
	{
		gappick0( tmpseq, seq[i] );
		if( count_unusual( tmpseq, usual ) <= maxunusual ) 
		{
			shortenN( tmpseq, unknown ); // 2022/Apr

			fprintf( stdout, ">%s\n", name[i]+1 );
//			fprintf( stdout, "%s\n", seq[i] );
			fprintf( stdout, "%s\n", tmpseq ); // 2022/Apr
			nout++;
		}
	}

	if( nout < njob )
	{
		if( dorp == 'p' )
			fprintf( stderr, "\n\nRemoved %d sequence(s) where the frequency of ambiguous amino acids > %5.3f\n\n\n", njob-nout, maxunusual );
		else
			fprintf( stderr, "\n\nRemoved %d sequence(s) where the frequency of ambiguous bases > %5.3f\n\n\n", njob-nout, maxunusual );
	}
	
	free( nlen );
	free( tmpseq );
	FreeCharMtx( seq );
	FreeCharMtx( name );

	return( 0 );
}
