#include "mltaln.h"

#define DEBUG 0

void arguments( int argc, char *argv[] )
{
    int c;

	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	scoremtx = NOTSPECIFIED;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					fprintf( stderr, "kimuraR = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'D':
					scoremtx = -1;
					break;
				case 'P':
					scoremtx = 0;
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
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
}


int main( int ac, char **av )
{
	int *nlen;
	static char **name, **seq;
	double score;
	extern double score_calc_for_score( int, char ** );

	arguments( ac, av );

	getnumlen( stdin );
	rewind( stdin );

	nlen = AllocateIntVec( njob );
	name = AllocateCharMtx( njob, B+1 );
	seq = AllocateCharMtx( njob, nlenmax+2 );

	readData_pointer( stdin, name, nlen, seq );

	if( !isaligned( njob, seq ) ) ErrorExit( "Not aligned." );

	constants( njob, seq );

	score = score_calc_for_score( njob, seq );
	if( scoremtx == 0 ) score += offset;

	fprintf( stdout, "score = %f\n", score );
	if     ( scoremtx ==  0 ) fprintf( stdout, "JTT %dPAM\n", pamN );
	else if( scoremtx ==  1 ) fprintf( stdout, "Dayhoff( machigai ga aru )\n" );
	else if( scoremtx ==  2 ) fprintf( stdout, "M-Y\n" );
	else if( scoremtx == -1 ) fprintf( stdout, "DNA 1:%d\n", kimuraR );

	fprintf( stdout, "gap penalty = %+6.2f, %+6.2f, %+6.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
	exit( 0 );
}
