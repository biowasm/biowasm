#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 1

double corethr;
int coreext;

void arguments( int argc, char *argv[] )
{
    int c;

	fftkeika = 1;
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
	scoremtx = 0;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	corethr = .5;
	coreext = 0;

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
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
					fprintf( stderr, "blosum %d\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					fprintf( stderr, "weighti = %f\n", fastathreshold );
					--argc;
					goto nextoption;
				case 'i':
					corethr = atof( *++argv );
					fprintf( stderr, "corethr = %f\n", corethr );
					--argc;
					goto nextoption;
				case 'm':
					fmodel = 1;
					break;
				case 'c':
					coreext = 1;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
				case 'Q':
					calledByXced = 1;
					break;
				case 's':
					treemethod = 's';
					break;
				case 'x':
					treemethod = 'x';
					break;
				case 'p':
					treemethod = 'p';
					break;
				case 'a':
					alg = 'a';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'S':
					alg = 'S';
					break;
				case 'C':
					alg = 'C';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					disp = 1;
					break;
				case 'o':
					outgap = 0;
					break;
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
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
	static char **oseq;
	static double **pscore;
	static double *eff;
	static double **node0, **node1;
	static double *gapc;
	static double *avgap;
	double tmpavgap;
	int i, j, m, goffset;
	static int ***topol;
	static double **len;
	FILE *prep;
	char c;
	int corestart, coreend;
	int alloclen;
	int winsize;
	char *pt, *ot;
	double gapmin;

	arguments( argc, argv );

	getnumlen( stdin );
	rewind( stdin );

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	seq = AllocateCharMtx( njob, nlenmax*9+1 );
	name = AllocateCharMtx( njob, B+1 );
	oseq = AllocateCharMtx( njob, nlenmax*9+1 );
	alloclen = nlenmax*9;

	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	pscore = AllocateDoubleMtx( njob, njob );
	eff = AllocateDoubleVec( njob );
	node0 = AllocateDoubleMtx( njob, njob );
	node1 = AllocateDoubleMtx( njob, njob );
	gapc = AllocateDoubleVec( alloclen );
	avgap = AllocateDoubleVec( alloclen );

#if 0
	Read( name, nlen, seq );
#else
	readData_pointer( stdin, name, nlen, seq );
#endif

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

	writePre( njob, name, nlen, seq, 0 );

	if( tbutree == 0 )
	{
		for( i=1; i<njob; i++ ) 
		{
			if( nlen[i] != nlen[0] ) 
			{
				fprintf( stderr, "Input pre-aligned seqences or make hat2.\n" );
				exit( 1 );
			}
		}
		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		{
		/*
			pscore[i][j] = (double)score_calc1( seq[i], seq[j] );
		*/
			pscore[i][j] = (double)substitution_hosei( seq[i], seq[j] );
		}
	}
	else
	{
		fprintf( stderr, "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2_pointer( prep, njob, name, pscore );
		fclose( prep );
		fprintf( stderr, "done.\n" );

#if 0
		prep = fopen( "hat2_check", "w" );
		WriteHat2( prep, njob, name, pscore );
		fclose( prep );
#endif

	}

	fprintf( stderr, "Constructing dendrogram ... " );
	if( treemethod == 'x' )
		supg( njob, pscore, topol, len );
	else if( treemethod == 's' )
		spg( njob, pscore, topol, len );
	else if( treemethod == 'p' )
		upg2( njob, pscore, topol, len );
	else 
		ErrorExit( "Incorrect tree\n" );
	fprintf( stderr, "done.\n" );

	countnode( njob, topol, node0 );
	if( tbrweight )
	{
		weight = 3; 
#if 0
		utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
		counteff_simple( njob, topol, len, eff );
#endif
	}
	else
	{
		for( i=0; i<njob; i++ ) eff[i] = 1.0;
	}


	for( i=0; i<nlenmax; i++ )
	{
		gapc[i] = 0.0;
		for( j=0; j<njob; j++ )
		{
			if( seq[j][i] == '-' ) gapc[i] += eff[j];
		}
	}

	gapmin = 1.0;
	winsize = fftWinSize;
	goffset = winsize/2;
	tmpavgap = 0.0;
	corestart = coreend = -1;
	for( i=0; i<winsize; i++ )
	{
		tmpavgap += gapc[i];
	}
	for( i=winsize; i<nlenmax; i++ )
	{
		m = i - goffset;
		avgap[m] = tmpavgap / winsize;
//		fprintf( stdout, "%d %f %f\n", m, avgap[m], gapc[i] );
		if( avgap[m] < corethr )
		{
			if( corestart == -1 )
				corestart = i - winsize;
//			fprintf( stdout, "ok, gapmin = %f, corestart = %d, coreend = %d\n", gapmin, corestart, coreend );
			if( avgap[m] < gapmin )
			{ 
				gapmin = avgap[m];
			}
			coreend = i;
		}
		tmpavgap -= gapc[i-winsize];
		tmpavgap += gapc[i];
	}
	if( corestart == -1 || coreend == -1 )
	{
		corestart = 0;
		coreend = nlenmax-1;
	}

	for( i=0; i<njob; i++ )
	{
		pt = oseq[i];
		m = winsize;
		while( m-- ) *pt++ = '-';
		for( j=corestart; j<=coreend; j++ )
			*pt++ = seq[i][j];
		m = winsize;
		while( m-- ) *pt++ = '-';
		*pt = 0;

		ot = oseq[i]+winsize-1;
		pt = seq[i]+corestart-1;
		if( coreext ) m = winsize;
		else m = 0;
		while( m && --pt > seq[i] )
			if( *pt != '-' )
			{
				*ot-- = *pt;
				m--;
			}

		ot = oseq[i]+winsize+coreend-corestart+1;
		pt = seq[i]+coreend;
		if( coreext ) m = winsize;
		else m = 0;
		while( m && *(++pt) )
		{
			if( *pt != '-' ) 
			{
				*ot++ = *pt;
				m--;
			}
		}
		fprintf( stdout, ">%s\n", name[i] );
		fprintf( stdout, "%s\n", oseq[i] );
	}

	exit( 1 );

	SHOWVERSION;
	return( 0 );
}
