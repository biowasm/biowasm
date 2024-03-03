#include "mltaln.h"

#define DEBUG 0

void arguments( int argc, char *argv[] )
{
    int c;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
//					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
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
	double nfreq;

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
	getnumlen_nogap_countn( infp, &nlenmin, &nfreq );

	fprintf( stdout, "%d x %d - %d %c nfreq=%f\n", njob, nlenmax, nlenmin, dorp, nfreq );

	fclose( infp );
	return( 0 );
}
