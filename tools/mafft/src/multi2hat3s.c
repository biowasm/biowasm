#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 1
#define TSUYOSAFACTOR 100


static int nhomologs;
static int seedoffset;

void strip( char *s )
{
	char *pt = s;
	while( *++pt )
		if( *pt == '\n' ) *pt = 0;
}


void arguments( int argc, char *argv[] )
{
    int c;

	seedoffset = 0;
	nhomologs = 1;
	inputfile = NULL;
	fftkeika = 0;
	pslocal = -1000.0;
	constraint = 0;
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
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_OP = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	ppenalty_EX = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "seed = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 't':
					nhomologs = myatoi( *++argv );
					fprintf( stderr, "nhomologs = %d\n", nhomologs );
					--argc;
					goto nextoption;
				case 'o':
					seedoffset = myatoi( *++argv );
					fprintf( stderr, "seedoffset = %d\n", seedoffset );
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
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
	if( alg == 'C' && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : C, o\n" );
		exit( 1 );
	}
}

int countamino( char *s, int end )
{
	int val = 0;
	while( end-- )
		if( *s++ != '-' ) val++;
	return( val );
}

static void pairalign( char **name, int nlen[M], char **seq, double *effarr, int alloclen )
{
	int i, j;
	FILE *hat3p;
	double pscore = 0.0; // by D.Mathog
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	char *aseq;
	static char **pseq;
	LocalHom **localhomtable, *tmpptr;
	double tsuyosa;

	if( nhomologs < 1 ) nhomologs = 1; // tsuyosa=0.0 wo sakeru
	tsuyosa = (double)nhomologs * nhomologs * TSUYOSAFACTOR;
	fprintf( stderr, "tsuyosa = %f\n", tsuyosa );
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

	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		pseq = AllocateCharMtx( 2, 0 );
		aseq = AllocateCharVec( nlenmax*9+1 );
#if 0
#else
#endif
	}

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif


//	writePre( njob, name, nlen, aseq, 0 );

	hat3p = fopen( "hat3", "w" );
	if( !hat3p ) ErrorExit( "Cannot open hat3." );
	fprintf( stderr, "\n" );
	for( i=0; i<njob-1; i++ )
	{
		for( j=i+1; j<njob; j++ )
		{
			pseq[0] = seq[i];
			pseq[1] = seq[j];

			if( strlen( pseq[0] ) != strlen( pseq[1] ) )
			{
				fprintf( stderr, "## ERROR  ###\n" );
				fprintf( stderr, "Not aligned,  %s - %s\n", name[i], name[j] );
				fprintf( stderr, "## ERROR  ###\n" );
				exit( 1 );
			}


			fprintf( stderr, "adding %d-%d\r", i, j );
			putlocalhom2( pseq[0], pseq[1], localhomtable[i]+j, 0, 0, (int)pscore, strlen( pseq[0] ), 'k' );
			for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				if( tmpptr->start1 == -1 ) continue;
				fprintf( hat3p, "%d %d %d %6.3f %d %d %d %d k\n", i+seedoffset, j+seedoffset, tmpptr->overlapaa, tmpptr->opt * tsuyosa, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2 ); 
			}
		}
	}
	fprintf( stderr, "\n" );
	fclose( hat3p );

#if DEBUG
	fprintf( stderr, "calling FreeLocalHomTable\n" );
#endif
	FreeLocalHomTable( localhomtable, njob );
#if DEBUG
	fprintf( stderr, "done. FreeLocalHomTable\n" );
#endif
}

static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    fprintf( stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else if( alg == 'S' ) 
		fprintf( fp, "Apgorithm S\n" );
	else if( alg == 'C' ) 
		fprintf( fp, "Apgorithm A+/C\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'x' )
		fprintf( fp, "Tree = UPGMA (3).\n" );
	else if( treemethod == 's' )
		fprintf( fp, "Tree = UPGMA (2).\n" );
	else if( treemethod == 'p' )
		fprintf( fp, "Tree = UPGMA (1).\n" );
	else
		fprintf( fp, "Unknown tree.\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
        else
        {
            if( fftscore )
                fprintf( fp, "Basis : Polarity and Volume\n" );
            else
                fprintf( fp, "Basis : 20 amino acids\n" );
        }
        fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
        fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
    }
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}
	 

int main( int argc, char *argv[] )
{
	static int  nlen[M];	
	static char **name, **seq;
	static char **bseq;
	static double *eff;
	int i;
	char c;
	int alloclen;
	FILE *infp;

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

	getnumlen( infp );
	rewind( infp );

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	name = AllocateCharMtx( njob, B+1 );
	seq = AllocateCharMtx( njob, nlenmax*9+1 );
	bseq = AllocateCharMtx( njob, nlenmax*9+1 );
	alloclen = nlenmax*9;

	eff = AllocateDoubleVec( njob );

#if 0
	Read( name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
#endif
	fclose( infp );

	constants( njob, seq );

#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illeagal character %c\n", c );
		exit( 1 );
	}

//	writePre( njob, name, nlen, seq, 0 );

	for( i=0; i<njob; i++ ) eff[i] = 1.0;


	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );


//	for( i=0; i<njob; i++ ) fprintf( stdout, ">_seed_%s\n%s\n", name[i]+1, bseq[i] ); // CHUUI!!
	for( i=0; i<njob; i++ ) fprintf( stdout, ">_seed_%s\n%s\n", name[i]+1, seq[i] );

	pairalign( name, nlen, seq, eff, alloclen );

	fprintf( trap_g, "done.\n" );
#if DEBUG
	fprintf( stderr, "closing trap_g\n" );
#endif
	fclose( trap_g );

#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif
	SHOWVERSION;
	return( 0 );
}
