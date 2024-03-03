#include "mltaln.h"
#include <errno.h>

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SHISHAGONYU 0 // for debug


// from tbfast
static int treein;
static int treeout;


// from pairlocalalign
static int stdout_dist;


static void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	nadd = 0;
	inputfile = NULL;
	fftkeika = 0;
	pslocal = -1000.0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'x';
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	divpairscore = 0;
	stdout_dist = 0;
//	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_OP = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	ppenalty_EX = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	specificityconsideration = 0.0;
	usenaivescoreinsteadofalignmentscore = 0;
	specifictarget = 0;
	nwildcard = 0;
	compacttree = 2; // tsuneni!
	treein = 0;
	treeout = 0;
	fastathreshold = 2.7;
	constraint = 2;

//	localhomfile = 0; // tbfast.c no wo tsukaunode comment out

//	reporterr( "argc=%d\n", argc );
//	reporterr( "*argv=%s\n", *argv );
//	reporterr( "(*argv)[0]=%c\n", (*argv)[0] );
    while( --argc > 0 && (*++argv)[0] == '-' )
	{
//		reporterr( "(*argv)[0] in while loop = %s\n", (*argv) );
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
//					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'O':
					ppenalty_OP = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'E':
					ppenalty_EX = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
//					fprintf( stderr, "kimuraR = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					fprintf( stderr, "blosum %d\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
//					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
//					fprintf( stderr, "TM %d\n", pamN );
					--argc;
					goto nextoption;
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					--argc;
					goto nextoption;
#if 0
				case 'l':
					ppslocal = (int)( atof( *++argv ) * 1000 + 0.5 );
					pslocal = (int)( 600.0 / 1000.0 * ppslocal + 0.5);
//					fprintf( stderr, "ppslocal = %d\n", ppslocal );
//					fprintf( stderr, "pslocal = %d\n", pslocal );
					--argc;
					goto nextoption;
#else
#endif
				case 'C':
					nthread = myatoi( *++argv );
					if( nthread == 0 ) nthread = 1;
//					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
#ifndef enablemultithread
					nthread = 1;
#endif
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
//					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'u':
					specificityconsideration = (double)myatof( *++argv );
//					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
				case 'c':
					stdout_dist = 1;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
#if 0
				case 'r':
					fmodel = -1;
					break;
#endif
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
#if 0
				case 'Q':
					calledByXced = 1;
					break;
				case 'x':
					disp = 1;
					break;
				case 'a':
					alg = 'a';
					break;
				case 'S':
					alg = 'S';
					break;
#endif
				case 'N':
					alg = 'N';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'L':
					alg = 'L';
					break;
				case 'Z':
					usenaivescoreinsteadofalignmentscore = 1;
					break;
				case 'B': // hitsuyou! memopt -M -B no tame
					break;
#if 0
				case 'Y':
					alg = 'Y'; // nadd>0 no toki nomi. moto no hairetsu to atarashii hairetsuno alignmnt -> L;
					break;
				case 's':
					alg = 's';
					break;
				case 'G':
					alg = 'G';
					break;
				case 'B': // hitsuyou! memopt -M -B no tame
					break;
				case 'T':
					alg = 'T';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'R':
					alg = 'R';
					break;
				case 'r':
					alg = 'r'; // nadd>0 no toki nomi. moto no hairetsu to atarashii hairetsuno alignmnt -> R, last
					break;
				case 'V':
					alg = 'V';
					break;
#endif
				case 'T': // tbfast.c no noalign ni taiou
					break;

				case 'F':
					use_fft = 1;
					break;
				case 'U':
					treein = 1;
					break;
				case 't':
					treeout = 1;
					break;
				case 'y':
					divpairscore = 1;
					break;
				case '=':
					specifictarget = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
				case 'q':
					lhlimit = myatoi( *++argv );
					--argc; 
					goto nextoption;
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
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
    if( argc != 0 ) 
    {
        fprintf( stderr, "pairlocalalign options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}

int main( int argc, char *argv[] )
{
	static int  *nlen = NULL;	
	static int *selfscore = NULL;
	static char **name = NULL, **seq = NULL;
	static double *eff = NULL;
	int i;
	static int ***topol = NULL;
	static Treedep *dep = NULL;
	static double **len = NULL;
	FILE *infp = NULL;

	char c;

	arguments( argc, argv );

	if( alg != 'A' && alg != 'L' && alg != 'N' )
	{
		reporterr( "alg %c is not yet supported\n", alg );
		exit( 1 );
	}
	if( alg != 'N' && usenaivescoreinsteadofalignmentscore == 1 )
	{
		reporterr( "The combination of usenaivescoreinsteadofalignmentscore and alg %c is not yet supported\n", alg );
		exit( 1 );
	}

	if( fastathreshold < 0.0001 )
	{
		constraint = 0; 
		lhlimit = 0;
	}

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

	getnumlen( infp );
	rewind( infp );


	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}


#if !defined(mingw) && !defined(_MSC_VER)
	setstacksize( 200 * njob ); // topolorder() de ookime no stack wo shiyou.
#endif


	seq = AllocateCharMtx( njob, nlenmax+1 );

	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	selfscore = AllocateIntVec( njob );

	topol = AllocateIntCub( njob, 2, 0 );
	len = AllocateFloatMtx( njob, 2 );
	eff = AllocateDoubleVec( njob );


	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );


#if 0
	readData( infp, name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
	fclose( infp );
#endif

	constants( njob, seq );

#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();
	initFiles();

//	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

//	writePre( njob, name, nlen, seq, 0 );

	if( treein ) loadtree( njob, topol, len, name, nlen, dep, treeout );
	pairalign_node( njob, nlenmax, name, seq, topol, len, dep, treein, treeout );


	FreeCharMtx( seq ); seq = NULL;

	FreeCharMtx( name ); name = NULL;
	free( nlen ); nlen = NULL;
	free( selfscore ); selfscore = NULL;

	for( i=0; i<njob; i++ ) 
	{
		free( topol[i][0] ); topol[i][0] = NULL;
		free( topol[i][1] ); topol[i][1] = NULL;
		free( topol[i] ); topol[i] = NULL;
		free( len[i] ); len[i] = NULL;
	}
	free( topol ); topol = NULL;
	free( len ); len = NULL;
	free( eff ); eff = NULL;
	free( dep ); dep = NULL;


	freeconstants();
	closeFiles();
	FreeCommonIP();

	SHOWVERSION; // koko?

	return( 0 );

#if 0 // gui mitaiou
chudan:
	if( seq ) FreeCharMtx( seq ); seq = NULL;

	if( name ) FreeCharMtx( name ); name = NULL;
	if( nlen ) free( nlen ); nlen = NULL;
	if( selfscore ) free( selfscore ); selfscore = NULL;



	if( topol ) FreeIntCub( topol ); topol = NULL;
	if( dep ) free( dep ); dep = NULL;

	freeconstants();
	closeFiles();
	FreeCommonIP();

	if( topol )
	{
		for( i=0; i<njob; i++ ) 
		{
			if( topol[i] )
			{
				if( topol[i][0] ) free( topol[i][0] ); topol[i][0] = NULL;
				if( topol[i][1] ) free( topol[i][1] ); topol[i][1] = NULL;
			}
			free( topol[i] ); topol[i] = NULL;
		}
		free( topol ); topol = NULL;
	}
	if( len ) 
	{
		for( i=0; i<njob; i++ ) 
		{
			if( len[i] ) free( len[i] ); len[i] = 0;
		}
		free( len ); len = NULL;
	}
	if( eff ) free( eff ); eff = NULL;
	if( dep ) free( dep ); dep = NULL;

	return( 0 );
#endif

}
