#include "mltaln.h"

#define DEBUG 0

char *weboutfile = NULL;


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
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'w':
					weboutfile = *++argv;
					fprintf( stderr, "weboutfile = %s\n", weboutfile );
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
	FILE *weboutfp;
	int nlenmin;
	int isaligned = 0;

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

	if( weboutfile )
	{
		weboutfp = fopen( weboutfile, "w" );
		if( !weboutfp )
		{
			fprintf( stderr, "Cannot open %s\n", weboutfile );
			exit( 1 );
		}
	}

	dorp = NOTSPECIFIED;
	if( weboutfile )
	{
		getnumlen_nogap_outallreg_web( infp, weboutfp, &nlenmin, &isaligned );
		if( isaligned ) fprintf( stdout, "Aligned\n" );
		else fprintf( stdout, "Not aligned\n" );
	}
	else
		getnumlen_nogap_outallreg( infp, &nlenmin );

	return( 0 );

}
