#include "mltaln.h"

#define DEBUG 0 

#if DEBUG
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
double getrusage_sec()
{
    struct rusage t;
    struct timeval tv;
    getrusage(RUSAGE_SELF, &t);
    tv = t.ru_utime;
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}
#endif


int intcmp( int *str1, int *str2 )
{
	while( *str1 != -1 && *str2 != -1 )
		if( *str1++ != *str2++ ) return( 1 );
	if( *str1 != *str2 ) return( 1 );
	return( 0 );
}

char **arguments( int argc, char *argv[] )
{
    int c = 0;
	
	fmodel = 0;
	nblosum = 62;
	calledByXced = 0;
	devide = 0;
	fftscore = 1;
	use_fft = 1;
	nevermemsave = 0;
	alg = 'A';
    weight = 0;
    utree = 1;
	tbutree = 0;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 0;
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
    scoremtx = 1;
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = 0; // chokusetsu yobareru kara
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
	treemethod = 'x';


    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'P':
					dorp = 'p';
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'N':
					use_fft = 0;
					break;
				case 'n':
					nevermemsave = 1;
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'Q':
					alg = 'Q';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'd':
					disp = 1;
					break;
				case 'O':
					outgap = 0;
					break;
				case 'a':
					fmodel = 1;
					break;
				case 'u':
					tbrweight = 0;
					break;
				case 'U':
					tbrweight = -1;
					break;
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
                    fprintf( stderr, "kappa = %d\n", kimuraR );
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
					TMorJTT = JTT;
					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "tm %d\n", pamN );
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
    if( argc != 2 ) 
    {
        fprintf( stderr, "options: Check source file ! %c ?\n", c );
        exit( 1 );
    }
	fprintf( stderr, "tbitr = %d, tbrweight = %d, tbweight = %d\n", tbitr, tbrweight, tbweight );
//	readOtherOptions( &ppid, &fftThreshold, &fftWinSize );
	return( argv ); 

}

void GroupAlign( int nseq1, int nseq2, char **name, int *nlen, char **seq, char **aseq, char **mseq1, char **mseq2, int ***topol, double **len, double *eff, int alloclen )
{
	int i;
	int clus1, clus2;
	int s1, s2;
	double pscore;
	static char **name1, **name2;
	double *effarr = eff;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	static char *indication1, *indication2;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int intdum;
#if DEBUG
	double time1, time2;
#endif


//	fprintf( stderr, "in GroupAlign fftWinSize   = %d\n", fftWinSize );
//	fprintf( stderr, "in GroupAlign fftThreshold = %d\n", fftThreshold );

	if( effarr1 == NULL ) 
	{
		name1 = AllocateCharMtx( nseq1, B );
		name2 = AllocateCharMtx( nseq2, B );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
#if 0
#else
#endif
	}

	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


		
	s1 = 0;
	s2 = nseq1;
//	fprintf( stdout, "nseq1 = %d\n", nseq1 );



	clus1 = conjuctionforgaln( 0, nseq1, aseq, mseq1, effarr1, effarr, name, name1, indication1 );
	clus2 = conjuctionforgaln( nseq1, njob, aseq, mseq2, effarr2, effarr, name, name2, indication2 );
/*
	fprintf( stderr, "before align all\n" );
	display( aseq, njob );
	fprintf( stderr, "\n" );
	fprintf( stderr, "before align 1 %s \n", indication1 );
	display( mseq1, clus1 );
	fprintf( stderr, "\n" );
	fprintf( stderr, "before align 2 %s \n", indication2 );
	display( mseq2, clus2 );
	fprintf( stderr, "\n" );
*/

	commongappick( nseq1, mseq1 );
	commongappick( nseq2, mseq2 );


#if DEBUG
	time1 = getrusage_sec();
	fprintf( stderr, "entering Falign\n" );
#endif
	if( use_fft )
	{
		if( alg == 'M' )
			pscore = Falign_udpari_long( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, alloclen, &intdum );
		else
			pscore = Falign( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, alloclen, &intdum, NULL, 0, NULL );
	}
	else
	{
		if( alg == 'M' )
			pscore = MSalignmm( n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, NULL, NULL, NULL, 0.0, 0.0 );
		else
			pscore = A__align( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	}
#if DEBUG
		time2 = getrusage_sec();
		fprintf( stdout, "### %d - %d,  %f\n", clus1, clus2, time2-time1 );
		fflush( stdout );
#endif

	
/*
	fprintf( stderr, "after align 1 %s \n", indication1 );
	display( mseq1, clus1 );
	fprintf( stderr, "\n" );
	fprintf( stderr, "after align 2 %s \n", indication2 );
	display( mseq2, clus2 );
	fprintf( stderr, "\n" );
*/

	fprintf( stderr, "group-to-group %s /%s     %f\n", indication1, indication2, pscore );
	if( disp ) display( aseq, njob );
	fprintf( stderr, "\n" );

/*
	trap = fopen( "pre", "r+" );
	if( !trap ) ErrorExit( 1 );
	WriteGapFill( trap, njob, name, nlen, aseq );
	fclose( trap );
	fprintf( stdout, "nseq1 = %d\n", nseq1 );
*/
}

static void WriteOptions( FILE *fp )
{
	fprintf( fp, "tree-base method\n" );
	if( tbweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbweight == 3 ) fprintf( fp, "reversely weighted\n" );
	if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
	else if( scoremtx ==  1 ) fprintf( fp, "Dayhoff( machigai ga aru )\n" );
	else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	else if( scoremtx == -1 ) fprintf( fp, "DNA\n" );

	if( scoremtx )
		fprintf( fp, "Gap Penalty = %d, %d\n", penalty, offset );
	else
		fprintf( fp, "Gap Penalty = %d\n", penalty );
}
	 

int main( int argc, char *argv[] )
{
	char **argv2;
	static int  *nlen;	
	static char **name, **seq;
	static char **seq1, **seq2;
	static char **mseq1, **mseq2;
	static char **aseq;
	static char **bseq;
	static double **pscore;
	static double *eff;
	int i, j, len1, len2;
	static int ***topol;
	static double **len;
	FILE *gp1, *gp2;
	char c;
	int nlenmax1, nlenmax2, nseq1, nseq2;
	int alloclen;

	argv2 = arguments( argc, argv );

	fprintf( stderr, "####### in galn\n" );

	initFiles();

	fprintf( stderr, "file1 = %s\n", argv2[0] );
	fprintf( stderr, "file2 = %s\n", argv2[1] );

	gp1 = fopen( argv2[0], "r" ); if( !gp1 ) ErrorExit( "cannot open file1" );
	gp2 = fopen( argv2[1], "r" ); if( !gp2 ) ErrorExit( "cannot open file2" );

#if 0
	PreRead( gp1, &nseq1, &nlenmax1 );
	PreRead( gp2, &nseq2, &nlenmax2 );
#else
    getnumlen( gp1 );
	nseq1 = njob; nlenmax1 = nlenmax;
    getnumlen( gp2 );
	nseq2 = njob; nlenmax2 = nlenmax;
#endif

	njob = nseq1 + nseq2;
	nlenmax = MAX( nlenmax1, nlenmax2 );

	rewind( gp1 );
	rewind( gp2 );


	name = AllocateCharMtx( njob, B );
	nlen = AllocateIntVec( njob );
	seq1 = AllocateCharMtx( nseq1, nlenmax*3 );
	seq2 = AllocateCharMtx( nseq2, nlenmax*3 );
	seq  = AllocateCharMtx( njob, 1 );
	aseq = AllocateCharMtx( njob, nlenmax*3 );
	bseq = AllocateCharMtx( njob, nlenmax*3 );
	mseq1 = AllocateCharMtx( njob, 1 );
	mseq2 = AllocateCharMtx( njob, 1 );
	alloclen = nlenmax * 3;

	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	if( tbrweight > 0 ) pscore = AllocateDoubleMtx( njob, njob );
	eff = AllocateDoubleVec( njob );

#if 0
    njob=nseq2; FRead( gp2, name+nseq1, nlen+nseq1, seq2 );
	njob=nseq1; FRead( gp1, name, nlen, seq1 );
#else
    njob=nseq2; readDataforgaln( gp2, name+nseq1, nlen+nseq1, seq2 );
	njob=nseq1; readDataforgaln( gp1, name, nlen, seq1 );
#endif
	njob = nseq1 + nseq2;


#if 0  // CHUUI
	commongappick( nseq1, seq1 );
	commongappick( nseq2, seq2 );
#endif

	for( i=0; i<nseq1; i++ ) seq[i] = seq1[i];
	for( i=nseq1; i<njob; i++ ) seq[i] = seq2[i-nseq1];
/*
	Write( stdout, njob, name, nlen, seq );
*/

    constants( njob, seq );

    WriteOptions( trap_g );

    c = seqcheck( seq );
    if( c )
    {
        fprintf( stderr, "Illeagal character %c\n", c );
        exit( 1 );
    }
    for( i=1; i<nseq1; i++ ) 
    {
        if( nlen[i] != nlen[0] ) 
            ErrorExit( "group1 is not aligned." );
    }
    for( i=nseq1+1;  i<njob; i++ ) 
    {
        if( nlen[i] != nlen[nseq1] ) 
            ErrorExit( "group2 is not aligned." );
    }


    if( tbutree == 0 )
	{
		if( tbrweight > 0 )
		{
			for( i=0; i<nseq1; i++ ) 
			{
				for( j=i+1; j<nseq1; j++ )
				{
					pscore[i][j] = (double)substitution_hosei( seq[i], seq[j] );
//					fprintf( stderr, "%d-%d, %5.1f \n", i, j, pscore[i][j] );
				}
				for( j=nseq1; j<njob; j++ )
				{
					pscore[i][j] = 3.0;
//					fprintf( stderr, "%d-%d, %5.1f \n", i, j, pscore[i][j] );
				}
			}
			for( i=nseq1; i<njob-1; i++ ) 
			{
				for( j=i+1; j<njob; j++ )
				{
					pscore[i][j] = (double)substitution_hosei( seq[i], seq[j] );
//					fprintf( stderr, "%d-%d, %5.1f \n", i, j, pscore[i][j] );
				}
			}
//			fprintf( stderr, "\n" );

		}

    }
   	else
	{
		fprintf( stderr, "Not supported\n" );
		exit( 1 );
#if 0
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2( prep, njob, name, pscore );
		fclose( prep );
#endif
	}

	if( tbrweight > 0 )
	{
		fprintf( stderr, "Constructing dendrogram ... " );
		if( treemethod == 'x' )
			veryfastsupg( njob, pscore, topol, len );
		else
			ErrorExit( "Incorrect tree\n" );
		fprintf( stderr, "done.\n" );

		weight = 3;
		counteff_simple( njob, topol, len, eff );
//		for( i=0; i<njob; i++ ) fprintf( stderr, "eff[%d] = %f\n", i, eff[i] );
	}
	else if ( tbrweight == 0 )
	{
		for( i=0; i<njob; i++ ) eff[i] = 1.0;
	}
	else
	{
		getweightfromname( njob, eff, name );
	}

	//for( i=0; i<njob; i++ ) reporterr( "weight[%d] = %f\n", i, eff[i] );
	//exit( 1 );

	len1 = strlen( seq[0] );
	len2 = strlen( seq[nseq1] );
	if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
//	if( len1 > 30000 || len2 > 30000 )
	{       
		fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2 );
		alg = 'M';
	}       
        


	reporterr( "GroupAglin..\n" );

	GroupAlign( nseq1, nseq2, name, nlen, seq, aseq, mseq1, mseq2, topol, len, eff, alloclen );

#if 0
	writePre( njob, name, nlen, aseq, 1 );
#else
	writeDataforgaln( stdout, njob, name, nlen, aseq );
#endif

	SHOWVERSION;
	return( 0 );
}
