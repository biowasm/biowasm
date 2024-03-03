#include "mltaln.h"


#define TREE 1
#define PICKSIZE 50 // must be >= 3
#define WEIGHT 0
#define TOKYORIPARA 0.70 // 0.70
#define TOKYORIPARA_A 0.70  // changed
#define LENFAC 1
#define HUKINTOTREE 1
#define DIANA 0
#define MAX6DIST 10.0

// kouzoutai ni sasareru pointer ha static

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

#define END_OF_VEC -1

static char *fastapath;
static int doalign;
static int fromaln;
static int uselongest;
static int treeout;
static int classsize;
static int picksize;
static int reorder;
static int pid;
static int maxdepth = 0;
static double tokyoripara;

#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define DLENFACA 0.01
#define DLENFACB 2500
#define DLENFACC 2500
#define DLENFACD 0.1

static char datafile[1000];
static char queryfile[1000];
static char resultfile[1000];

typedef struct _scores
{
	double score;
	int selfscore;
	int orilen;
	int *pointt;
	int numinseq;
	char *name;
//	char *seq; // reallo
//	char **seqpt;
	int shimon;
} Scores;

int intcompare( const int *a, const int *b )
{
	return( *a - *b );
}

int lcompare( const Scores *a, const Scores *b )
{
	if( a->orilen < b->orilen ) return 1;
	else if( a->orilen > b->orilen ) return -1;
	else return 0;
}

int dcompare( const Scores *a, const Scores *b )
{
	if( a->score > b->score ) return 1;
	else if( a->score < b->score ) return -1;
	else
	{
		if( a->selfscore < b->selfscore ) return 1;
		else if( a->selfscore > b->selfscore ) return -1;
		else 
		{
			if( a->orilen < b->orilen ) return 1;
			else if( a->orilen > b->orilen ) return -1;
			else return 0;
		}
	}
}


static void getfastascoremtx( int **tmpaminodis )
{
	FILE *qfp;
	FILE *dfp;
	FILE *rfp;
	int i, j;
	char aa;
	int slen;
	int res;
	char com[10000];
	static char *tmpseq;
	static char *tmpname;
	double *resvec;

	if( scoremtx == -1 )
	{
		tmpaminodis['a']['a'] = 5;
		tmpaminodis['g']['g'] = 5;
		tmpaminodis['c']['c'] = 5;
		tmpaminodis['t']['t'] = 5;
		tmpaminodis['n']['n'] = -1;

		return;
	}


	tmpseq = calloc( 2000, sizeof( char ) );
	tmpname = calloc( B, sizeof( char ) );
	resvec = calloc( 1, sizeof( double ) );

//	fprintf( stderr, "xformatting .. " );
	dfp = fopen( datafile, "w" );
	if( !dfp ) ErrorExit( "Cannot open datafile." );
	sprintf( tmpname, ">+===========+%d                      ", 0 );
	strcpy( tmpseq, "AAAAAAXXXXXX" );
	strcat( tmpseq, "CCCCCCXXXXXX" );
	strcat( tmpseq, "DDDDDDXXXXXX" );
	strcat( tmpseq, "EEEEEEXXXXXX" );
	strcat( tmpseq, "FFFFFFXXXXXX" );
	strcat( tmpseq, "GGGGGGXXXXXX" );
	strcat( tmpseq, "HHHHHHXXXXXX" );
	strcat( tmpseq, "IIIIIIXXXXXX" );
	strcat( tmpseq, "KKKKKKXXXXXX" );
	strcat( tmpseq, "LLLLLLXXXXXX" );
	strcat( tmpseq, "MMMMMMXXXXXX" );
	strcat( tmpseq, "NNNNNNXXXXXX" );
	strcat( tmpseq, "PPPPPPXXXXXX" );
	strcat( tmpseq, "QQQQQQXXXXXX" );
	strcat( tmpseq, "RRRRRRXXXXXX" );
	strcat( tmpseq, "SSSSSSXXXXXX" );
	strcat( tmpseq, "TTTTTTXXXXXX" );
	strcat( tmpseq, "VVVVVVXXXXXX" );
	strcat( tmpseq, "WWWWWWXXXXXX" );
	strcat( tmpseq, "YYYYYYXXXXXX" );
	slen = strlen( tmpseq );
	writeData_pointer( dfp, 1, &tmpname, &slen, &tmpseq );
	fclose( dfp );
	fprintf( stderr, "done.\n" );

	for( i=0; i<20; i++ )
	{
		aa = amino[i];
//		fprintf( stderr, "checking %c\n", aa );
		*tmpseq = 0;
		sprintf( tmpname, ">+===========+%d                      ", 0 );
		for( j=0; j<6; j++ )
			sprintf( tmpseq+strlen( tmpseq ), "%c", aa );
		qfp = fopen( queryfile, "w" );
		if( !qfp ) ErrorExit( "Cannot open queryfile." );
		writeData_pointer( qfp, 1, &tmpname, &slen, &tmpseq );
		fclose( qfp );

		if( scoremtx == -1 ) 
			sprintf( com, "%s -z3 -m10  -n -Q -H -b%d -E%d -d%d %s %s %d > %s", fastapath,  M, M, 0, queryfile, datafile, 6, resultfile );
		else
			sprintf( com, "%s -z3 -m10  -p -Q -H -b%d -E%d -d%d %s %s %d > %s", fastapath,  M, M, 0, queryfile, datafile, 2, resultfile );
		res = system( com );
		if( res )
		{
			fprintf( stderr, "error in %s", fastapath );
			exit( 1 );
		}

		rfp = fopen( resultfile, "r" );
		if( rfp == NULL )  
			ErrorExit( "file 'fasta.$$' does not exist\n" );
		res = ReadFasta34m10_scoreonly( rfp, resvec, 1 );
		fprintf( stderr, "%c: %f\n", 'A'+i, *resvec/6 );
		fclose( rfp );
		if( ( (int)*resvec % 6 ) > 0.0 )
		{
			fprintf( stderr, "Error in blast, *resvec=%f\n", *resvec );
			fprintf( stderr, "Error in blast, *resvec/6=%f\n", *resvec/6 );
			exit( 1 );
		}
		tmpaminodis[(int)aa][(int)aa] = (int)( *resvec / 6 );
//		fprintf( stderr, "*resvec=%f, tmpaminodis[aa][aa] = %d\n", *resvec, tmpaminodis[aa][aa] );
	}
	tmpaminodis['X']['X'] = -1;
	free( tmpname );
	free( tmpseq );
	free( resvec );
}

#if 0
static void getblastscoremtx( int **tmpaminodis )
{
	FILE *qfp;
	FILE *dfp;
	FILE *rfp;
	int i, j;
	char aa;
	int slen;
	int res;
	char com[10000];
	static char *tmpseq;
	static char *tmpname;
	double *resvec;

	if( scoremtx == -1 )
	{
		tmpaminodis['a']['a'] = 1;
		tmpaminodis['g']['g'] = 1;
		tmpaminodis['c']['c'] = 1;
		tmpaminodis['t']['t'] = 1;

		return;
	}


	tmpseq = calloc( 2000, sizeof( char ) );
	tmpname = calloc( B, sizeof( char ) );
	resvec = calloc( 1, sizeof( double ) );

//	fprintf( stderr, "xformatting .. " );
	dfp = fopen( datafile, "w" );
	if( !dfp ) ErrorExit( "Cannot open datafile." );
	sprintf( tmpname, "\0", i ); // BUG!!
	strcpy( tmpseq, "AAAAAAXXXXXX" );
	strcat( tmpseq, "CCCCCCXXXXXX" );
	strcat( tmpseq, "DDDDDDXXXXXX" );
	strcat( tmpseq, "EEEEEEXXXXXX" );
	strcat( tmpseq, "FFFFFFXXXXXX" );
	strcat( tmpseq, "GGGGGGXXXXXX" );
	strcat( tmpseq, "HHHHHHXXXXXX" );
	strcat( tmpseq, "IIIIIIXXXXXX" );
	strcat( tmpseq, "KKKKKKXXXXXX" );
	strcat( tmpseq, "LLLLLLXXXXXX" );
	strcat( tmpseq, "MMMMMMXXXXXX" );
	strcat( tmpseq, "NNNNNNXXXXXX" );
	strcat( tmpseq, "PPPPPPXXXXXX" );
	strcat( tmpseq, "QQQQQQXXXXXX" );
	strcat( tmpseq, "RRRRRRXXXXXX" );
	strcat( tmpseq, "SSSSSSXXXXXX" );
	strcat( tmpseq, "TTTTTTXXXXXX" );
	strcat( tmpseq, "VVVVVVXXXXXX" );
	strcat( tmpseq, "WWWWWWXXXXXX" );
	strcat( tmpseq, "YYYYYYXXXXXX" );
	slen = strlen( tmpseq );
	writeData_pointer( dfp, 1, &tmpname, &slen, &tmpseq );
	fclose( dfp );
	if( scoremtx == -1 )
		sprintf( com, "formatdb  -p f -i %s -o F", datafile );
	else
		sprintf( com, "formatdb  -i %s -o F", datafile );
	system( com );
	fprintf( stderr, "done.\n" );

	for( i=0; i<20; i++ )
	{
		aa = amino[i];
		fprintf( stderr, "checking %c\n", aa );
		*tmpseq = 0;
		for( j=0; j<6; j++ )
			sprintf( tmpseq+strlen( tmpseq ), "%c", aa );
		qfp = fopen( queryfile, "w" );
		if( !qfp ) ErrorExit( "Cannot open queryfile." );
		writeData_pointer( qfp, 1, &tmpname, &slen, &tmpseq );
		fclose( qfp );

		sprintf( com, "blastall -b %d -G 10 -E 1 -e 1e10 -p blastp -m 7  -i %s -d %s >  %s\0", 1, queryfile, datafile, resultfile );
		res = system( com );
		if( res )
		{
			fprintf( stderr, "error in %s", "blastall" );
			exit( 1 );
		}

		rfp = fopen( resultfile, "r" );
		if( rfp == NULL )  
			ErrorExit( "file 'fasta.$$' does not exist\n" );
		res = ReadBlastm7_scoreonly( rfp, resvec, 1 );
		fprintf( stdout, "%c: %f\n", 'A'+i, *resvec/6 );
		fclose( rfp );
		if( ( (int)*resvec % 6 ) > 0.0 )
		{
			fprintf( stderr, "Error in blast, *resvec=%f\n", *resvec );
			fprintf( stderr, "Error in blast, *resvec/6=%f\n", *resvec/6 );
			exit( 1 );
		}
		tmpaminodis[aa][aa] = (int)( *resvec / 6 );
	}
	tmpaminodis['X']['X'] = 0;
	free( tmpname );
	free( tmpseq );
	free( resvec );
	
}
#endif

static double *callfasta( char **seq, Scores *scores, int nin, int *picks, int query, int rewritedata )
{
	double *val;
	FILE *qfp;
	FILE *dfp;
	FILE *rfp;
	int i;
	char com[10000];
	static char datafile[1000];
	static char queryfile[1000];
	static char resultfile[1000];
	static int pid;
	static char *tmpseq;
	static char *tmpname;
	int slen;
	int res;
	static Scores *scoresbk = NULL;
	static int ninbk = 0;

	if( pid == 0 )
	{
		pid = (int)getpid();
		sprintf( datafile, "/tmp/data-%d", pid );
		sprintf( queryfile, "/tmp/query-%d", pid );
		sprintf( resultfile, "/tmp/fasta-%d", pid );

		tmpseq = calloc( nlenmax+1, sizeof( char ) );
		tmpname = calloc( B+1, sizeof( char ) );
	}

	val = calloc( nin, sizeof( double ) );
//	fprintf( stderr, "nin=%d, q=%d\n", nin, query );

	if( rewritedata )
	{
		scoresbk = scores;
		ninbk = nin;
//		fprintf( stderr, "\nformatting .. " );
		dfp = fopen( datafile, "w" );
		if( !dfp ) ErrorExit( "Cannot open datafile." );
		if( picks == NULL ) for( i=0; i<nin; i++ )
		{
//			fprintf( stderr, "i=%d / %d / %d\n", i,  nin, njob );
//			fprintf( stderr, "nlenmax = %d\n", nlenmax );
//			fprintf( stderr, "scores[i].orilen = %d\n", scores[i].orilen );
//			fprintf( stderr, "strlen( seq[scores[i].numinseq] = %d\n", strlen( seq[scores[i].numinseq] ) );
			gappick0( tmpseq, seq[scores[i].numinseq] );
			sprintf( tmpname, ">+===========+%d                      ", i );
			slen = scores[i].orilen;
			writeData_pointer( dfp, 1, &tmpname, &slen, &tmpseq );
		}
		else for( i=0; i<nin; i++ )
		{
			gappick0( tmpseq, seq[scores[picks[i]].numinseq] );
			sprintf( tmpname, ">+===========+%d                      ", i );
			slen = scores[picks[i]].orilen;
			writeData_pointer( dfp, 1, &tmpname, &slen, &tmpseq );
		}
		fclose( dfp );
	}


	gappick0( tmpseq, seq[scores[query].numinseq] );
	sprintf( tmpname, ">+==========+%d                      ", 0 );
	slen = scores[query].orilen;
	qfp = fopen( queryfile, "w" );
	if( !qfp ) ErrorExit( "Cannot open queryfile." );
	writeData_pointer( qfp, 1, &tmpname, &slen, &tmpseq );
	fclose( qfp );

//	fprintf( stderr, "calling fasta, nin=%d\n", nin );

	if( scoremtx == -1 ) 
		sprintf( com, "%s  -z3 -m10  -n -Q -H -b%d -E%d -d%d %s %s %d > %s",  fastapath, nin, nin, 0, queryfile, datafile, 6, resultfile );
	else
		sprintf( com, "%s  -z3 -m10  -p -Q -H -b%d -E%d -d%d %s %s %d > %s",  fastapath, nin, nin, 0, queryfile, datafile, 2, resultfile );
	res = system( com );
	if( res )
	{
		fprintf( stderr, "error in %s", fastapath );
		exit( 1 );
	}
//	fprintf( stderr, "fasta done\n" );

//exit( 1 );

	rfp = fopen( resultfile, "r" );
	if( rfp == NULL )  
		ErrorExit( "file 'fasta.$$' does not exist\n" );

//	fprintf( stderr, "reading fasta\n" );
	if( scoremtx == -1 ) 
		res = ReadFasta34m10_scoreonly_nuc( rfp, val, nin );
	else
		res = ReadFasta34m10_scoreonly( rfp, val, nin );
//	fprintf( stderr, "done. val[0] = %f\n", val[0] );


	fclose( rfp );

#if 0
	for( i=0; i<nin; i++ )
		fprintf( stderr, "r[%d-%d] = %f\n", 0, i, val[i] );
	exit( 1 );
#endif

	return( val );
}
#if 0
static double *callblast( char **seq, Scores *scores, int nin, int query, int rewritedata )
{
	double *val;
	FILE *qfp;
	FILE *dfp;
	FILE *rfp;
	int i, j;
	char com[10000];
	static char datafile[1000];
	static char queryfile[1000];
	static char resultfile[1000];
	static int pid;
	static char *tmpseq;
	static char *tmpname;
	char *seqptr;
	int slen;
	int res;
	static Scores *scoresbk = NULL;
	static int ninbk = 0;

	if( pid == 0 )
	{
		pid = (int)getpid();
		sprintf( datafile, "/tmp/data-%d\0", pid );
		sprintf( queryfile, "/tmp/query-%d\0", pid );
		sprintf( resultfile, "/tmp/fasta-%d\0", pid );

		tmpseq = calloc( nlenmax+1, sizeof( char ) );
		tmpname = calloc( B+1, sizeof( char ) );
	}

	val = calloc( nin, sizeof( double ) );
//	fprintf( stderr, "nin=%d, q=%d\n", nin, query );

	if( rewritedata )
	{
		scoresbk = scores;
		ninbk = nin;
		fprintf( stderr, "\nformatting .. " );
		dfp = fopen( datafile, "w" );
		if( !dfp ) ErrorExit( "Cannot open datafile." );
		for( i=0; i<nin; i++ )
		{
//			fprintf( stderr, "i=%d / %d / %d\n", i,  nin, njob );
//			fprintf( stderr, "nlenmax = %d\n", nlenmax );
//			fprintf( stderr, "scores[i].orilen = %d\n", scores[i].orilen );
//			fprintf( stderr, "strlen( seq[scores[i].numinseq] = %d\n", strlen( seq[scores[i].numinseq] ) );
			gappick0( tmpseq, seq[scores[i].numinseq] );
			sprintf( tmpname, "+===========+%d                      \0", i );
			slen = scores[i].orilen;
			writeData_pointer( dfp, 1, &tmpname, &slen, &tmpseq );
		}
		fclose( dfp );
			
		if( scoremtx == -1 )
			sprintf( com, "formatdb  -p f -i %s -o F", datafile );
		else
			sprintf( com, "formatdb  -i %s -o F", datafile );
		system( com );
//		fprintf( stderr, "done.\n" );
	}


	gappick0( tmpseq, seq[scores[query].numinseq] );
	sprintf( tmpname, "+==========+%d                      \0", 0 );
	slen = scores[query].orilen;
	qfp = fopen( queryfile, "w" );
	if( !qfp ) ErrorExit( "Cannot open queryfile." );
	writeData_pointer( qfp, 1, &tmpname, &slen, &tmpseq );
	fclose( qfp );
//	fprintf( stderr, "q=%s\n", tmpseq );

	fprintf( stderr, "\ncalling blast .. \n" );
	if( scoremtx == -1 ) 
		sprintf( com, "blastall -b %d -e 1e10 -p blastn -m 7  -i %s -d %s >  %s\0", nin, queryfile, datafile, resultfile );
	else
		sprintf( com, "blastall -b %d -G 10 -E 1 -e 1e10 -p blastp -m 7  -i %s -d %s >  %s\0", nin, queryfile, datafile, resultfile );
	res = system( com );
	if( res ) ErrorExit( "error in blast" );

	rfp = fopen( resultfile, "r" );
	if( rfp == NULL )  
		ErrorExit( "file 'fasta.$$' does not exist\n" );
	res = ReadBlastm7_scoreonly( rfp, val, nin );
	fclose( rfp );

#if 0
	for( i=0; i<nin; i++ )
		fprintf( stderr, "r[%d-%d] = %f\n", 0, i, val[i] );
#endif

	return( val );
}
#endif

#if 0
static void selhead( int *ar, int n )
{
	int min = *ar;
	int *minptr = ar;
	int *ptr = ar;
	int tmp;
	n--;
	ar++;
	while( n-- )
	{
		if( ( tmp = *ptr++ ) < min )
		{
			min = tmp;
			minptr = ptr;
		}
	}
	if( minptr != ar )
	{
		tmp = *ar;
		*ar = min;
		*minptr = tmp;
	}
	return;
}
#endif

void arguments( int argc, char *argv[] )
{
    int c;

	doalign = 0;
	fromaln = 0;
	treeout = 0;
	uselongest = 1;
	reorder = 1;
	nevermemsave = 0;
	inputfile = NULL;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
	force_fft = 0;
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
	treemethod = 'X';
	sueff_global = 0.1;
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = -1530;
	ppenalty_ex = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
	classsize = NOTSPECIFIED;
	picksize = NOTSPECIFIED;
	tokyoripara = NOTSPECIFIED;
	legacygapcost = 0;
	nwildcard = 0;
	outnumber = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'p':
					picksize = myatoi( *++argv );
					fprintf( stderr, "picksize = %d\n", picksize );
					--argc;
					goto nextoption;
				case 's':
					classsize = myatoi( *++argv );
					fprintf( stderr, "groupsize = %d\n", classsize );
					--argc;
					goto nextoption;
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "poffset = %d\n", poffset );
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
//					fprintf( stderr, "blosum %d\n", nblosum );
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
				case 'T':
					tokyoripara = (double)atof( *++argv );
					--argc;
					goto nextoption;
				case 'l':
					uselongest = 0;
					break;
				case 'n' :
					outnumber = 1;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'S':
					doalign = 'f';
					break;
				case 'Z':
					fromaln = 1;
					break;
				case 'U':
					doalign = 1;
					break;
				case 'x':
					reorder = 0;
					break;
				case 't':
					treeout = 1;
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
				case 'L':
					legacygapcost = 1;
					break;
#if 0
				case 'R':
					fftRepeatStop = 1;
					break;
				case 'Q':
					calledByXced = 1;
					break;
				case 'a':
					alg = 'a';
					break;
				case 'R':
					alg = 'R';
					break;
				case 'Q':
					alg = 'Q';
					break;
#endif
				case 'A':
					alg = 'A';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'M':
					alg = 'M';
					break;
				case 'C':
					alg = 'C';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					use_fft = 1;
					force_fft = 1;
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
				case 'J':
					tbutree = 0;
					break;
				case 'X':
					treemethod = 'X'; // tsukawareteiru ????
					sueff_global = atof( *++argv );
					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
				case 'E':
					treemethod = 'E'; // upg (average)
					break;
				case 'q':
					treemethod = 'q'; // minimum
					break;
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case ':':
					nwildcard = 1;
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

static int nunknown = 0;

int seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	return( grp-grpbk );
}

int seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	return( grp-grpbk );
}

void makecompositiontable_p( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
	{
#if 0
		table[point]++;
#else
		if( (unsigned int)table[point]++ >= INT_MAX )
		{
			reporterr( "Overflow. table[point]=%d>INT_MAX(%d).\n", table[point], INT_MAX );
			exit( 1 );
		}
#endif
	}
}

static int localcommonsextet_p( int *table, int *pointt )
{
	int value = 0;
	unsigned int tmp;
	int point;
	static int *memo = NULL;
	static int *ct = NULL;
	static int *cp;

	if( !memo )
	{
		memo = (int *)calloc( tsize, sizeof( int ) );
		if( !memo ) ErrorExit( "Cannot allocate memo\n" );
		ct = (int *)calloc( MIN( maxl, tsize )+1, sizeof( int ) );
		if( !ct ) ErrorExit( "Cannot allocate memo\n" );
	}

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
#if 0
		if( tmp >= INT_MAX )
		{
			reporterr( "Overflow.\n" );
			reporterr( "cp-ct=%d, point=%d, tmp=%d, memo[point]=%d>INT_MAX(%d)\n", cp-ct, point, tmp, memo[point], INT_MAX );
			exit( 1 );
		}
#endif
	}
	*cp = END_OF_VEC;
	
	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	return( value );
}

void makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable( int *pointt, int *n )
{
	int point;
	register int *p;

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

#if 1
static void pairalign( int nseq, int *nlen, char **seq, int *mem1, int *mem2, double *weight, int *alloclen )
{
	int l, len1, len2;
	int clus1, clus2;
	double pscore, tscore;
	static int *fftlog;
	static char *indication1, *indication2;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static char **mseq1, **mseq2;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
#if 0
	int i, j;
#endif


	if( effarr1 == NULL ) 
	{
		fftlog = AllocateIntVec( nseq );
		effarr1 = AllocateDoubleVec( nseq );
		effarr2 = AllocateDoubleVec( nseq );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		mseq1 = AllocateCharMtx( nseq, 0 );
		mseq2 = AllocateCharMtx( nseq, 0 );
		for( l=0; l<nseq; l++ ) fftlog[l] = 1;
	}

	tscore = 0.0;
	m1 = mem1[0];
	m2 = mem2[0];
	len1 = strlen( seq[m1] );
	len2 = strlen( seq[m2] );
	if( *alloclen < len1 + len2 )
	{
		fprintf( stderr, "\nReallocating.." );
		*alloclen = ( len1 + len2 ) + 1000;
		ReallocateCharMtx( seq, nseq, *alloclen + 10 ); 
		fprintf( stderr, "done. *alloclen = %d\n", *alloclen );
	}

#if WEIGHT
	clus1 = fastconjuction_noname( mem1, seq, mseq1, effarr1, weight, indication1, 0.0 );
	clus2 = fastconjuction_noname( mem2, seq, mseq2, effarr2, weight, indication2, 0.0 );
#else
	clus1 = fastconjuction_noweight( mem1, seq, mseq1, effarr1, indication1 );
	clus2 = fastconjuction_noweight( mem2, seq, mseq2, effarr2, indication2 );
#endif

#if 0
	for( i=0; i<clus1; i++ )
		fprintf( stderr, "in p seq[%d] = %s\n", mem1[i], seq[mem1[i]] );
	for( i=0; i<clus2; i++ )
		fprintf( stderr, "in p seq[%d] = %s\n", mem2[i], seq[mem2[i]] );
#endif

#if 0
	fprintf( stderr, "group1 = %.66s", indication1 );
	if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
	fprintf( stderr, "\n" );
	fprintf( stderr, "group2 = %.66s", indication2 );
	if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
	fprintf( stderr, "\n" );
#endif

//	fprintf( stdout, "mseq1 = %s\n", mseq1[0] );
//	fprintf( stdout, "mseq2 = %s\n", mseq2[0] );

	if( !nevermemsave && ( alg != 'M' && ( len1 > 10000 || len2 > 10000  ) ) )
	{
		fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
		alg = 'M';
		if( commonIP ) FreeIntMtx( commonIP );
		commonIP = 0;
		commonAlloc1 = 0;
		commonAlloc2 = 0;
	}

	if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
	else	    				   ffttry = 0;

	if( force_fft || ( use_fft && ffttry ) )
	{
		fprintf( stderr, "\bf" );
		if( alg == 'M' )
		{
			fprintf( stderr, "\bm" );
//			fprintf( stderr, "%d-%d", clus1, clus2 );
			pscore = Falign_udpari_long( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
		}
		else
		{
//			fprintf( stderr, "%d-%d", clus1, clus2 );
			pscore = Falign( NULL, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
		}
	}
	else
	{
		fprintf( stderr, "\bd" );
		fftlog[m1] = 0;
		switch( alg )
		{
			case( 'a' ):
				pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
				break;
			case( 'M' ):
				fprintf( stderr, "\bm" );
//				fprintf( stderr, "%d-%d", clus1, clus2 );
				pscore = MSalignmm( n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, NULL, NULL, NULL, 0.0, 0.0 );
				break;
			case( 'A' ):
				if( clus1 == 1 && clus2 == 1 )
				{
//					fprintf( stderr, "%d-%d", clus1, clus2 );
					pscore = G__align11( n_dis_consweight_multi, mseq1, mseq2, *alloclen, outgap, outgap );
				}
				else
				{
//					fprintf( stderr, "%d-%d", clus1, clus2 );
					pscore = A__align( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
				}
				break;
			default:
				ErrorExit( "ERROR IN SOURCE FILE" );
		}
	}
#if SCOREOUT
	fprintf( stderr, "score = %10.2f\n", pscore );
#endif
	nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );
	return;
}
#endif

#if 0
static void treebase( int nseq, int *nlen, char **aseq, double *eff, int nalign, int ***topol, int *alloclen ) // topol
{
	int i, l;
	int nlim;
	int clus1, clus2;

	nlim = nalign-1;
	for( l=0; l<nlim; l++ )
	{
		fprintf( stderr, "in treebase, l = %d\n", l );
		fprintf( stderr, "aseq[0] = %s\n", aseq[0] );
		fprintf( stderr, "aseq[topol[l][0][0]] = %s\n", aseq[topol[l][0][0]] );
		pairalign( nseq, nlen, aseq, topol[l][0], topol[l][1], eff, alloclen );
		free( topol[l][0] );
		free( topol[l][1] );
	}
}
#endif

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

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = .UPGMA (average)\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree = Minimum linkage.\n" );
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
	 
#if 1
static int splitseq_mq( Scores *scores, int nin, int *nlen, char **seq, char **orialn, char **name, char *inputfile, int uniform, char **tree, int *alloclen, int *order, int *whichgroup, double *weight, int *depthpt, int qinoya )
{
	int val;
	int ii, jj;
	int *mptr;
	int treelen = 0; // by Mathog, a guess
	static int groupid = 0;
	static int branchid = 0;
	int i, j;
	int selfscore0;
	double **dfromc;
	double **dfromcp;
	double **pickmtx;
	double **yukomtx;
	static int *table1;
	Scores **outs, *ptr;
	int *numin;
	int *tsukau;
	int belongto;
	char **children = NULL; // by Mathog, a guess
	char *tmptree;
	static int *orderpos = NULL;
	int rn;
	int npick;
	int nyuko;
	int *picks;
	int *yukos;
	int *s_p_map;
	int *p_o_map;
	int *s_y_map;
	int *y_o_map;
	int *closeh;
	int nkouho;
	int *pickkouho;
	int *iptr;
	int *jptr;
	int aligned;
	int ***topol;
	int *treeorder;
	int picktmp;
	double **len;
	double minscore;
//	double *minscoreinpick;
	double *hanni;
	double lenfac;
	double longer;
	double shorter;
	static char **mseq1 = NULL;
	static char **mseq2 = NULL;
	double *blastresults = NULL; // by Mathog, a guess
	static int palloclen = 0;
	double maxdist;

	if( orderpos == NULL )
		orderpos = order;
	if( palloclen == 0 )
		palloclen = *alloclen * 2;
	if( mseq1 == NULL && doalign == 1 )
	{
		mseq1 = AllocateCharMtx( 1, palloclen );
		mseq2 = AllocateCharMtx( 1, palloclen );
	}

	if( nin == 0 ) 
	{
#if TREE
		if( treeout )
		{
			*tree = (char *)calloc( 1, sizeof( char ) );
			**tree = 0;
		}
#endif
		return 1;
	}

	if( nin < 2 || uniform == -1 ) // kokodato muda deha nai ga
	{
		fprintf( stderr, "\nLeaf  %d / %d                ", ++branchid, njob );
#if 0
		outputfile = AllocateCharVec( strlen( inputfile ) + 100 );
		sprintf( outputfile, "%s-%d", inputfile, branchid );
		if( uniform > 0 )
//			sprintf( outputfile, "%su%d", outputfile, uniform );
			sprintf( outputfile + strlen(outputfile), "u%d", uniform );
		fprintf( stderr, "GROUP %d: %d member(s) (%d) %s\n", branchid, nin, scores[0].numinseq, outputfile );
		outfp = fopen( outputfile, "w" );
		free( outputfile );
		if( outfp == NULL )
		{
			fprintf( stderr, "Cannot open %s\n", outputfile );
			exit( 1 );
		}
		for( j=0; j<nin; j++ )
			fprintf( outfp, ">G%d %s\n%s\n", branchid, scores[j].name+1, seq[scores[j].numinseq] );
		fclose( outfp );
#endif


#if TREE
		if( treeout )
		{
			treelen = 0;
			tmptree = calloc( 100, sizeof( char ) );
			for( j=0; j<nin; j++ )
			{
				treelen += sprintf( tmptree, "%d", scores[j].numinseq+1 );
			}
			free( tmptree );
	
			*tree = (char *)calloc( treelen + nin + 15, sizeof( char ) );
			**tree = '\n';
			if( nin > 1 ) 
			{
				*(*tree+1) = '(';
				*(*tree+2) = '\0';
			}
			else
			{
				*(*tree+1) = '\0';
			}
			for( j=0; j<nin-1; j++ )
			{
				sprintf( *tree+strlen( *tree ), "%d,", scores[j].numinseq+1 );
			}
			sprintf( *tree+strlen( *tree ), "%d", scores[j].numinseq+1 );
			if( nin > 1 ) strcat( *tree, ")\n" );
			else strcat( *tree, "\n" );
//			fprintf( stdout, "*tree = %s\n", *tree );
		}

#endif
		for( j=0; j<nin; j++ )
		{
			*orderpos++ = scores[j].numinseq;
//			fprintf( stderr, "*order = %d\n", scores[j].numinseq );
		}

		return 1;
	}



	if( uselongest )
	{
		i = nin;
		ptr = scores;
		selfscore0 = scores->selfscore;
		belongto = 0;
		while( i-- )
		{
//			fprintf( stderr, "ptr-scores=%d, numinseq = %d, score = %f\n", ptr-scores, ptr->numinseq+1, ptr->score );
			if( ptr->selfscore > selfscore0 )
			{
				selfscore0 = ptr->selfscore;
				belongto = ptr-scores;
			}
			ptr++;
		} 
#if 1 
		if( belongto != 0 )
		{
//			fprintf( stderr, "swap %d %s\n<->\n%d %s\n", 0, scores->name, belongto, (scores+belongto)->name );
			ptr = calloc( 1, sizeof( Scores ) );
			*ptr = scores[belongto];
			scores[belongto] = *scores;
			*scores = *ptr;
			free( ptr );
		}
#endif
	}
	else
	{
		qsort( scores, nin, sizeof( Scores ), (int (*)())lcompare );
		belongto = (int)( 0.5 * nin );
//		fprintf( stderr, "lengths = %d, %d, %d\n", scores->orilen, scores[belongto].orilen, scores[nin-1].orilen );
		if( belongto != 0 )
		{
//			fprintf( stderr, "swap %d %s\n<->\n%d %s\n", 0, scores->name, belongto, (scores+belongto)->name );
			ptr = calloc( 1, sizeof( Scores ) );
			*ptr = scores[belongto];
			scores[belongto] = *scores;
			*scores = *ptr;
			free( ptr );
		}
	}

	if( qinoya != scores->numinseq )
//	if( 1 || qinoya != scores->numinseq )
	{
//		fprintf( stdout, "### scores->numinseq = %d, qinoya=%d, depth=%d\n", scores->numinseq, qinoya, *depthpt );


		if( doalign )
		{
			if( doalign == 'f' )
			{
				blastresults = callfasta( seq, scores, nin, NULL, 0, 1 );
				if( scores->selfscore != (int)blastresults[0] )
				{
					fprintf( stderr, "\n\nWARNING1: selfscore\n" );
					fprintf( stderr, "scores->numinseq = %d\n", scores->numinseq+1 );
					fprintf( stderr, "scores->orilen = %d\n", scores->orilen );
					fprintf( stderr, "scores->selfscore = %d, but blastresults[0] = %f\n", scores->selfscore, blastresults[0] );
//					if( abs( scores->selfscore - (int)blastresults[0] ) > 2 )
//						exit( 1 );
//					scores->selfscore = (int)blastresults[0]; //iinoka?
	
//					exit( 1 );
				}
			}
			else
				gappick0( mseq1[0], seq[scores->numinseq] );
		}
		else
		{
			table1 = (int *)calloc( tsize, sizeof( int ) );
			if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
			makecompositiontable_p( table1, scores[0].pointt );
		}
	
		selfscore0 = scores[0].selfscore;
		for( i=0; i<nin; i++ ) 
		{
			if( scores->orilen > scores[i].orilen )
			{
				longer = (double)scores->orilen;
				shorter = (double)scores[i].orilen;
			}
			else
			{
				longer = (double)scores[i].orilen; // nai
				shorter = (double)scores->orilen; //nai
			}

#if LENFAC
			lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//			lenfac = 1.0 / ( (double)LENFACA + (double)LENFACB / ( (double)longer + (double)LENFACC ) + (double)shorter / (double)longer * LENFACD );
//			fprintf( stderr, "lenfac = %f l=%d,%d\n", lenfac,scores->orilen, scores[i].orilen );
#else
			lenfac = 1.0;
#endif

			if( doalign )
			{
				if( doalign == 'f' )
				{
					scores[i].score = ( 1.0 - blastresults[i] / MIN( scores->selfscore, scores[i].selfscore ) ) * 1;
					if( scores[i].score < 0.0 ) scores[i].score = 0.0;
				}
				else
				{
					if( fromaln )
					{
//						scores[i].score = ( 1.0 - (double)G__align11_noalign( n_disLN, mseq1, mseq2, palloclen ) / MIN( selfscore0, scores[i].selfscore ) ) * 1;
						scores[i].score = ( 1.0 - (double)naivepairscore11( orialn[scores[i].numinseq], orialn[scores->numinseq], penalty ) / MIN( selfscore0, scores[i].selfscore ) ) * 1;
					}
					else
					{
						if( *depthpt == 0 ) fprintf( stderr, "\r%d / %d   ", i, nin );
						gappick0( mseq2[0], seq[scores[i].numinseq] );
//						fprintf( stdout, "### before calc scores[%d] = %f (%c)\n", i, scores[i].score, qinoya == scores->numinseq?'o':'x' );
						scores[i].score = ( 1.0 - (double)G__align11_noalign( n_disLN, -1200, -60, mseq1, mseq2, palloclen ) / MIN( selfscore0, scores[i].selfscore ) ) * 1;
//						fprintf( stderr, "scores[i] = %f\n", scores[i].score );
//						fprintf( stderr, "m1=%s\n", seq[scores[0].numinseq] );
//						fprintf( stderr, "m2=%s\n", seq[scores[i].numinseq] );
//						fprintf( stdout, "### before calc scores[%d] = %f (%c)\n", i, scores[i].score, qinoya == scores->numinseq?'o':'x' );
					}
				}
			}
			else
			{
				scores[i].score = ( 1.0 - (double)localcommonsextet_p( table1, scores[i].pointt ) / MIN( selfscore0, scores[i].selfscore ) ) * lenfac;
				if( scores[i].score > MAX6DIST ) scores[i].score = MAX6DIST;
			}
//			if( i ) fprintf( stderr, "%d-%d d %4.2f len %d %d\n", 1, i+1, scores[i].score, scores->orilen, scores[i].orilen );
		}
		if( doalign == 'f' ) free( blastresults );
		if( doalign == 0 ) free( table1 );
//exit( 1 );
	}

//	fprintf( stderr, "sorting .. " );
	qsort( scores, nin, sizeof( Scores ), (int (*)())dcompare );
//	fprintf( stderr, "done.\n" );


	maxdist = scores[nin-1].score;
	if( fromaln ) // kanzen itch ga misalign sareteiru kamoshirenai.
	{
		if( scores[0].shimon == scores[nin-1].shimon && !strcmp( seq[scores[0].numinseq], seq[scores[nin-1].numinseq] ) ) 
		{
			maxdist = 0.0;
		}
//		fprintf( stderr, "maxdist?? = %f, nin=%d, %d inori\n", scores[nin-1].score, nin, scores[nin-1].numinseq+1 );
	}

//	fprintf( stderr, "maxdist? = %f, nin=%d\n", scores[nin-1].score, nin );

	if( nin == 1 ) fprintf( stderr, "nin=1, scores[0].score = %f\n", scores[0].score );

// kokoni if( nin < 2 || ... )

	picks = AllocateIntVec( nin+1 );
	s_p_map = AllocateIntVec( nin+1 );
	s_y_map = AllocateIntVec( nin+1 );
	pickkouho = AllocateIntVec( nin+1 );
	closeh = AllocateIntVec( nin+1 );

//	nkouho = getkouho( pickkouho, (picksize+100)/nin, nin, scores, seq );
//	nkouho = getkouho( pickkouho, 1.0, nin, scores, seq ); // zenbu
//	fprintf( stderr, "selecting kouhos phase 2\n"  );
//	if( nkouho == 0 )
//	{
//		fprintf( stderr, "selecting kouhos, phase 2\n"  );
//		nkouho = getkouho( pickkouho, 1.0, nin, scores, seq );
//	}
//	fprintf( stderr, "\ndone\n\n"  );
	for( i=0; i<nin; i++ ) pickkouho[i] = i+1; nkouho = nin-1; // zenbu



	iptr = picks;
	*iptr++ = 0;
	npick = 1;
	if( nkouho > 0 )
	{
//		fprintf( stderr, "pickkouho[0] = %d\n", pickkouho[0] );
//		fprintf( stderr, "pickkouho[nin-1] = %d\n", pickkouho[nin-1] );
		picktmp = pickkouho[nkouho-1];
//		fprintf( stderr, "\nMOST DISTANT kouho=%d, nin=%d, nkouho=%d\n", picktmp, nin, nkouho );
		nkouho--;
		if( ( scores[picktmp].shimon == scores[0].shimon ) && ( !strcmp( seq[scores[0].numinseq], seq[scores[picktmp].numinseq] ) ) )
		{
//			fprintf( stderr, "known, j=%d (%d inori)\n", 0, scores[picks[0]].numinseq );
//			fprintf( stderr, "%s\n%s\n", seq[scores[picktmp].numinseq], seq[scores[picks[0]].numinseq] );
			;
		}
		else
		{
			*iptr++ = picktmp;
			npick++;
//			fprintf( stderr, "ok, %dth pick = %d (%d inori)\n", npick, picktmp, scores[picktmp].numinseq );
		}
	}
	i = 1;
	while( npick<picksize && nkouho>0 )
	{
		if( i )
		{
			i = 0;
			rn = nkouho * 0.5;
//			fprintf( stderr, "rn = %d\n", rn );
		}
		else
		{
			rn = rnd() * (nkouho);
		}
		picktmp = pickkouho[rn];
//		fprintf( stderr, "rn=%d/%d (%d inori), kouho=%d, nin=%d, nkouho=%d\n", rn, nkouho, scores[pickkouho[rn]].numinseq, pickkouho[rn], nin, nkouho );

//		fprintf( stderr, "#kouho before swap\n" );
//		for( i=0; i<nkouho; i++ ) fprintf( stderr, "%d ",  pickkouho[i] ); fprintf( stderr, "\n" );

		nkouho--;
		pickkouho[rn] = pickkouho[nkouho];
#if 1
//		fprintf( stderr, "#kouho after swap\n" ); 
//		for( i=0; i<nkouho; i++ ) fprintf( stderr, "%d ",  pickkouho[i] ); fprintf( stderr, "\n" );
		for( j=0; j<npick; j++ )
		{
			if( scores[picktmp].shimon == scores[picks[j]].shimon && !strcmp( seq[scores[picks[j]].numinseq], seq[scores[picktmp].numinseq] ) ) 
				break;
		}
		if( j == npick )
#endif
		{
//			fprintf( stderr, "ok, %dth pick = %d (%d inori)\n", npick, picktmp, scores[picktmp].numinseq );
			npick++;
			*iptr++ = picktmp;
		}
		else
		{
//			fprintf( stderr, "known, j=%d (%d inori)\n", j, scores[picks[j]].numinseq );
		}
	}
#if 0
	for( i=0; i<nin; i++ )
	{
		fprintf( stderr, "i=%d/%d, scores[%d].score = %f, inori=%d\n", i, nin, i, scores[i].score, scores[i].numinseq );
	}
	fprintf( stderr, "range:nin=%d scores[%d].score <= %f\n", nin, npick, scores[nin-1].score);
	for( i=0; i<npick; i++ )
	{
		fprintf( stderr, "i=%d/%d, scores[%d].score = %f, inori=%d\n", i, npick, picks[i], scores[picks[i]].score, scores[picks[i]].numinseq );
	}
exit( 1 );
#endif

//	fprintf( stderr, "\nnkouho=%d, defaultq2 = %d (%d inori)\n", nkouho, picks[npick-1], scores[picks[npick-1]].numinseq );

	qsort( picks, npick, sizeof( int ), (int (*)())intcompare );

//	fprintf( stderr, "allocating..\n" );

//	fprintf( stderr, "allocating outs, npick = %d\n", npick );
	numin = calloc( npick, sizeof( int ) );
	tsukau = calloc( npick, sizeof( int ) );
	outs = calloc( npick, sizeof( Scores * ) );
	for( i=0; i<npick; i++ ) outs[i] = NULL;
	topol = AllocateIntCub( npick, 2, 0 );
	treeorder = AllocateIntVec( npick + 1 );
	len = AllocateFloatMtx( npick, 2 );
	pickmtx = AllocateFloatHalfMtx( npick );
//	yukomtx = AllocateFloatHalfMtx( npick );
//	minscoreinpick = AllocateDoubleVec( npick );
	yukos = AllocateIntVec( npick );
	p_o_map = AllocateIntVec( npick+1 );
	y_o_map = AllocateIntVec( npick+1 );
	hanni = AllocateFloatVec( npick );

	for( i=0; i<nin; i++ ) s_p_map[i] = -1;
//	fprintf( stderr, "npick = %d\n", npick );
//	fprintf( stderr, "picks =" );
	for( i=0; i<npick; i++ )
	{
		s_p_map[picks[i]] = i;
		p_o_map[i] = scores[picks[i]].numinseq;
//		fprintf( stderr, " %d (%dinori)\n", picks[i], scores[picks[i]].numinseq+1 );
	}
//	fprintf( stderr, "\n" );

#if 0
	fprintf( stderr, "p_o_map =" );
	for( i=0; i<npick; i++ )
	{
		fprintf( stderr, " %d", p_o_map[i]+1 );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "picks =" );
	for( i=0; i<npick; i++ )
	{
		fprintf( stderr, " %d", picks[i] );
	}
	fprintf( stderr, "\n" );
#endif

	for( j=0; j<nin; j++ )
	{
		if( s_p_map[j] != -1 )
		{
			pickmtx[0][s_p_map[j]] = (double)scores[j].score;
//			fprintf( stderr, "pickmtx[0][%d] = %f\n", s_p_map[j], pickmtx[0][s_p_map[j]] );
		}
	}

	for( j=1; j<npick; j++ )
	{
		if( doalign )
		{
			if( doalign == 'f' )
			{
//				blastresults = callfasta( seq, scores, npick-j+1, picks+j-1, picks[j], 1 );
				blastresults = callfasta( seq, scores, npick, picks, picks[j], (j==1) );
				if( scores[picks[j]].selfscore != (int)blastresults[j] )
				{
					fprintf( stderr, "\n\nWARNING2: selfscore j=%d/%d\n", j, npick );
					fprintf( stderr, "scores[picks[j]].numinseq = %d\n", scores[picks[j]].numinseq+1 );
					fprintf( stderr, "scores[picks[j]].orilen = %d\n", scores[picks[j]].orilen );
					fprintf( stderr, "scores[picks[j]].selfscore = %d, but blastresults[j] = %f\n", scores[picks[j]].selfscore, blastresults[j] );
//					if( abs( scores[picks[j]].selfscore - (int)blastresults[j] ) > 2 )
//						exit( 1 );
//					scores->selfscore = (int)blastresults[0]; //iinoka?
				}
			}
			else
				gappick0( mseq1[0], seq[scores[picks[j]].numinseq] );
		}
		else
		{
			table1 = (int *)calloc( tsize, sizeof( int ) );
			if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
			makecompositiontable_p( table1, scores[picks[j]].pointt );
		}
	
		selfscore0 = scores[picks[j]].selfscore;
		pickmtx[j][0] = 0.0;
	    for( i=j+1; i<npick; i++ ) 
		{
			if( scores[picks[j]].orilen > scores[picks[i]].orilen )
			{
				longer = (double)scores[picks[j]].orilen;
				shorter = (double)scores[picks[i]].orilen;
			}
			else
			{
				longer = (double)scores[picks[i]].orilen;
				shorter = (double)scores[picks[j]].orilen;
			}
	
	#if LENFAC
			lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
	//		lenfac = 1.0 / ( (double)LENFACA + (double)LENFACB / ( (double)longer + (double)LENFACC ) + (double)shorter / (double)longer * LENFACD );
	//		fprintf( stderr, "lenfac = %f l=%d,%d\n", lenfac,scores->orilen, scores[i].orilen );
	#else
			lenfac = 1.0;
	#endif
	
			if( doalign )
			{
				if( doalign == 'f' )
				{
					pickmtx[j][i-j] = ( 1.0 - blastresults[i] / MIN( selfscore0, scores[picks[i]].selfscore ) ) * 1;
					if( pickmtx[j][i-j] < 0.0 ) pickmtx[j][i-j] = 0.0;
				}
				else
				{
					if( fromaln )
					{
						fprintf( stderr, "%d-%d/%d\r", j, i, npick );
						pickmtx[j][i-j] = ( 1.0 - (double)naivepairscore11(  orialn[scores[picks[i]].numinseq], orialn[scores[picks[j]].numinseq], penalty ) / MIN( selfscore0, scores[picks[i]].selfscore ) ) * 1;
					}
					else
					{
//						fprintf( stderr, "\r%d / %d   ", i, nin );
						gappick0( mseq2[0], seq[scores[picks[i]].numinseq] );
						pickmtx[j][i-j] = ( 1.0 - (double)G__align11_noalign( n_disLN, -1200, -60, mseq1, mseq2, palloclen ) / MIN( selfscore0, scores[picks[i]].selfscore ) ) * 1;
	//					fprintf( stderr, "scores[picks[i]] = %f\n", scores[picks[i]].score );
					}
				}
			}
			else
			{
				pickmtx[j][i-j] = ( 1.0 - (double)localcommonsextet_p( table1, scores[picks[i]].pointt ) / MIN( selfscore0, scores[picks[i]].selfscore ) ) * lenfac;
				if( pickmtx[j][i-j] > MAX6DIST ) pickmtx[j][i-j] = MAX6DIST;
			}

		}
		if( doalign == 'f' ) free( blastresults );
		if( doalign == 0 ) free( table1 );
	}

	dfromcp = AllocateDoubleMtx( npick, nin );
	dfromc = AllocateDoubleMtx( npick, 0 );

	for( i=0; i<npick; i++ ) for( j=0; j<nin; j++ )
		dfromcp[i][j] = -0.5;
	for( j=0; j<nin; j++ )
	{
		dfromcp[0][j] = ( scores[j].score );
//		fprintf( stderr, "j=%d, s_p_map[j]=%d\n", j, s_p_map[j] );
	}

	for( i=0; i<npick; i++ ) for( j=i; j<npick; j++ )
	{
		dfromcp[i][picks[j]] = dfromcp[j][picks[i]] = pickmtx[i][j-i];
	}

#if 0
	fprintf( stderr, "pickmtx = \n" );
	for( i=0; i<npick; i++ )
	{
		for( j=i; j<npick; j++ )
		{
			fprintf( stderr, "pickmtx[%d][%d] = %f\n", p_o_map[i]+1, p_o_map[j]+1, pickmtx[i][j-i] );
		}
	}
	exit( 1 );
#endif



//	for( i=0; i<npick-1; i++ ) for( j=i; j<npick; j++ )
//		fprintf( stderr, "dist[%d][%d] = %f\n", p_o_map[i]+1, p_o_map[j]+1, pickmtx[i][j-i] );

	for( i=0; i<npick; i++ ) tsukau[i] = 1;
	for( i=0; i<nin; i++ ) closeh[i] = -1;
	for( i=0; i<npick; i++ ) 
	{
		closeh[picks[i]] = picks[i];
//		fprintf( stderr, "i=%d/%d, picks[i]=%d, %d inori, closeh[%d] = %d \n", i, npick, picks[i], p_o_map[i]+1, picks[i], closeh[picks[i]] );
	}
#if 0
	fprintf( stderr, "closeh = \n" );
	for( i=0; i<nin; i++ )
	{
		fprintf( stderr, "%d ", closeh[i] );
	}
	fprintf( stderr, "\n" );
#endif
#if DIANA
	for( i=0; i<npick-1; i++ ) for( j=i; j<npick; j++ )
		fprintf( stderr, "dist[%d][%d] = %f\n", p_o_map[i]+1, p_o_map[j]+1, pickmtx[i][j-i] );
	fprintf( stderr, "DIANA!!\n" );
	if( npick > 2 )
	{
		double avdist;
		double avdist1;
		double avdist2;
		double maxavdist;
		int splinter;
		int count;
		int dochokoho;
		splinter = 0;
		int *docholist;
		int *docholistbk;
		maxavdist = 0.0;
		for( i=0; i<npick; i++ )
		{
			avdist = 0.0;
			for( j=i+1; j<npick; j++ )
			{
				avdist += pickmtx[i][j-i];
			}
			for( j=0; j<i; j++ )
			{
				avdist += pickmtx[j][i-j];
			}
			avdist /= (npick-1);
			fprintf( stderr, "avdist[%d] = %f\n", p_o_map[i] + 1, avdist );
			if( maxavdist < avdist ) 
			{
				maxavdist = avdist;
				splinter = i;
			}
		}
		fprintf( stderr, "splinter = %d (%d inori), maxavdist = %f\n", splinter, p_o_map[splinter]+1, maxavdist );

		docholist = AllocateIntVec( npick );
		docholistbk = AllocateIntVec( npick );
		for( i=0; i<npick; i++ ) docholist[i] = 0;
		docholist[splinter] = 1;
		while( 1 )
		{
			for( i=0; i<npick; i++ ) docholistbk[i] = docholist[i]; 
			for( dochokoho = 0; dochokoho<npick; dochokoho++ )
			{
				fprintf( stderr, "dochokoho=%d\n", dochokoho );
				if( docholist[dochokoho] ) continue;
				count = 0;
				avdist1 = 0.0;
				i=dochokoho;
				{
					for( j=i+1; j<npick; j++ )
					{
						if( docholist[j] || j == dochokoho ) continue;
						avdist1 += pickmtx[i][j-i];
						count++;
					}
					for( j=0; j<i; j++ )
					{
						if( docholist[j] || j == dochokoho ) continue;
						avdist1 += pickmtx[j][i-j];
						count++;
					}
				}
				if( count < 1 ) avdist1 = 0.0;
				else avdist1 /= (double)count;
				fprintf( stderr, "docho %d (%dinori), avdist1 = %f\n", dochokoho, p_o_map[dochokoho] + 1, avdist1 );

				count = 0;
				avdist2 = 0.0;
				i=dochokoho;
				{
					for( j=i+1; j<npick; j++ )
					{
						if( !docholist[j] || j == dochokoho ) continue;
						avdist2 += pickmtx[i][j-i];
						count++;
					}
					for( j=0; j<i; j++ )
					{
						if( !docholist[j] || j == dochokoho ) continue;
						avdist2 += pickmtx[j][i-j];
						count++;
					}
				}
				if( count < 1 ) avdist2 = 0.0;
				else avdist2 /= (double)count;
				fprintf( stderr, "docho %d (%dinori), avdist2 = %f\n", dochokoho, p_o_map[dochokoho] + 1, avdist2 );

				if( avdist2 < avdist1 ) 
				{
					docholist[dochokoho] = 1;
					hanni[dochokoho] = avdist2;
				}
				else
				{
					docholist[dochokoho] = 0;
					hanni[dochokoho] = avdist1;
				}
				fprintf( stderr, "avdist1=%f, avdist2=%f\n", avdist1, avdist2 );

			}
			for( i=0; i<npick; i++ ) if( docholist[i] != docholistbk[i] ) break;
			if( i == npick ) break;

			fprintf( stderr, "docholist = \n" );
			for( i=0; i<npick; i++ ) fprintf( stderr, "%d ", docholist[i] );
			fprintf( stderr, "\n" );
		}
		fprintf( stderr, "docholist = \n" );
		for( i=0; i<npick; i++ ) fprintf( stderr, "%d ", docholist[i] );
		fprintf( stderr, "\n" );

		for( i=0; i<npick; i++ ) if( docholist[i] == 0 ) break;
		yukos[0] = picks[i];
		for( i=0; i<npick; i++ ) if( docholist[i] == 1 ) break;
		yukos[1] = picks[splinter];

		for( i=0; i<npick; i++ ) 
		{
			if( docholist[i] == 0 ) closeh[picks[i]] = yukos[0];
			if( docholist[i] == 1 ) closeh[picks[i]] = yukos[1];
		}
//		for( i=0; i<npick; i++ ) closeh[picks[i]] = -1; // CHUUI !! iminai
		nyuko = 2;
		free( docholist );
		free( docholistbk );
	}
	else if( npick > 1 )
	{
		nyuko = 2;
		yukos[0] = picks[0]; yukos[1] = picks[1];
		closeh[picks[0]] = yukos[0];
		closeh[picks[1]] = yukos[1];
	}
	else
	{
		nyuko = 1;
		yukos[0] = picks[0];
		closeh[picks[0]] = yukos[0];
	}
#elif HUKINTOTREE
	if( npick > 2 )
	{
#if 0
		double avdist;
		double maxavdist;
		int count;
		int splinter;
		maxavdist = 0.0;
		splinter=0;
		for( i=0; i<npick; i++ )
		{
			avdist = 0.0;
			for( j=i+1; j<npick; j++ )
			{
				avdist += pickmtx[i][j-i];
			}
			for( j=0; j<i; j++ )
			{
				avdist += pickmtx[j][i-j];
			}
			avdist /= (npick-1);
			fprintf( stderr, "avdist[%d] = %f\n", p_o_map[i] + 1, avdist );
			if( maxavdist < avdist ) 
			{
				maxavdist = avdist;
				splinter = i;
			}
		}
		fprintf( stderr, "splinter = %d (%d inori), maxavdist = %f\n", splinter, p_o_map[splinter]+1, maxavdist );
#endif


//		fprintf( stderr, "check kaishi =>, npick=%d members = \n", npick );
//		for( i=0; i<npick; i++ ) fprintf( stderr, "%d (%d)", p_o_map[i]+1, picks[i] );
//		fprintf( stderr, "\n" );
		for( i=0; i<npick-1; i++ ) 
		{
			if( tsukau[i] == 0 ) continue;
			for( j=i+1; j<npick; j++ )
			{
//				double kijun = maxdist *  1/(npick-2);
//				double kijun = maxavdist * tokyoripara;
				double kijun;
				kijun = maxdist * tokyoripara;  // atode kakunin
//				fprintf( stderr, "%d-%d\n", i, j );
//				fprintf( stderr, "maxdist = %f\n", maxdist );
//				if( i==0 && j == 1 ) continue; // machigai!! CHUUI!!
//				if( maxdist == pickmtx[i][j-i] ) continue;
				if( tsukau[j] == 0 ) continue;
//				fprintf( stderr, "checking %d-%d (%d-%d) %f, kijun=%f\n", p_o_map[i]+1, p_o_map[j]+1, i, j, pickmtx[i][j-i], kijun );
				if( pickmtx[i][j-i] < kijun )
				{
//					fprintf( stderr, "dame!! %d => %d, because %f < %f\n", p_o_map[j]+1, p_o_map[i]+1, pickmtx[i][j-i], kijun );
#if 0
					if( scores[picks[i]].orilen > scores[picks[j]].orilen )
					{
						fprintf( stderr, "%d => %d\n", p_o_map[j]+1, p_o_map[i]+1 );
						tsukau[j] = 0;
					}
					else
					{
						fprintf( stderr, "%d => %d\n", p_o_map[i]+1, p_o_map[j]+1 );
						tsukau[i] = 0;
					}
					if( 0 && j == npick-1 ) tsukau[i] = 0;
					else  			   tsukau[j] = 0;
					fprintf( stderr, "tsukau[%d] = %d (%d inori)\n", j, tsukau[j], p_o_map[j]+1 );
#else
					tsukau[j] = 0;
					closeh[picks[j]] = closeh[picks[i]];
//					fprintf( stderr, "%d => tsukawanai\n", j );
#endif
				}
			}
		}
	}
	for( ii=0,i=0; i<npick; i++ )
	{
		if( tsukau[i] )
		{
			dfromc[ii] = dfromcp[i];
			ii++;
		}
		else
		{
			free( dfromcp[i] );
			dfromcp[i] = NULL;
		}
	}
	dfromc[ii] = NULL;

	for( ii=0,i=0; i<npick; i++ )
	{
		if( tsukau[i] )
		{
			for( jj=ii,j=i; j<npick; j++ )
			{
				if( tsukau[j] )
				{
					pickmtx[ii][jj-ii] = pickmtx[i][j-i];
					jj++;
				}
			}
			ii++;
		}
	}
	for( ; ii<npick; ii++ )
	{
		free( pickmtx[ii] ); 
		pickmtx[ii] = NULL;
	}

	for( ii=0,i=0; i<npick; i++ )
	{
		if( tsukau[i] )
		{
			yukos[ii++] = picks[i];
		}
	}


	nyuko = ii;
	yukomtx = pickmtx;
	pickmtx = NULL;

#endif
#if 0
	for( i=0; i<npick; i++ ) for( j=i; j<npick; j++ )
	{
		if( tsukau[i] == 1 && tsukau[j] == 1 )
			fprintf( stderr, "dist[%d][%d] = %f (ok)\n", p_o_map[i]+1, p_o_map[j]+1, pickmtx[i][j-i]  );
		else if( tsukau[i] == 0 && tsukau[j] == 0 )
			fprintf( stderr, "dist[%d][%d] = %f (xx)\n", p_o_map[i]+1, p_o_map[j]+1, pickmtx[i][j-i]  );
		else	
			fprintf( stderr, "%d-%d, okashii\n", p_o_map[i]+1, p_o_map[j]+1 );
	}
#endif
//	FreeFloatHalfMtx( pickmtx, npick );


	for( i=0; i<nin; i++ ) s_y_map[i] = -1;
//	fprintf( stderr, "npick = %d\n", npick );
//	fprintf( stderr, "yukos =" );
	for( i=0; i<nyuko; i++ )
	{
		s_y_map[yukos[i]] = i;
		y_o_map[i] = scores[yukos[i]].numinseq;
//		fprintf( stderr, " %d\n", yukos[i] );
	}
//	fprintf( stderr, "\n" );
#if 0
	for( i=0; i<nyuko; i++ )
	{
		fprintf( stderr, "y_o_map[%d] = %d\n", i, y_o_map[i]+1 );
	}
	for( i=0; i<nyuko; i++ ) for( j=i; j<nyuko; j++ )
	{
		fprintf( stderr, "yukodist[%d][%d] = %f (ok)\n", y_o_map[i]+1, y_o_map[j]+1, yukomtx[i][j-i]  );
	}
#endif

	for( i=0; i<nin; i++ )
	{
		if( closeh[i] != -1 )
		{
//			fprintf( stderr, "closeh[%d,%dinori] = %d,%dinori\n", i, scores[i].numinseq+1, closeh[i], scores[closeh[i]].numinseq+1 );
		}
	}

#if 0
		for( i=0; i<nyuko; i++ )
		{
			minscoreinpick[i] = 99.9;
			for( j=i+1; j<nyuko; j++ )
			{
				if( minscoreinpick[i] > yukomtx[i][j-i] )
					minscoreinpick[i] = yukomtx[i][j-i];
			}
			for( j=0; j<i; j++ )
			{
				if( minscoreinpick[i] > yukomtx[j][i-j] )
					minscoreinpick[i] = yukomtx[j][i-j];
			}
			fprintf( stderr, "minscoreinpick[%d(%dinori)] = %f\n", i, y_o_map[i]+1, minscoreinpick[i] );
		}
#endif


#if TREE
	if( treeout )
	{
		children = calloc( nyuko+1, sizeof( char * ) );
		for( i=0; i<nyuko+1; i++ ) children[i] = NULL;
	}
#endif
//	fprintf( stderr, "done..\n" );
	
//	fprintf( stderr, "classifying, nyuko=%d \n", nyuko );
	if( nyuko == 1 )
	{
		if( npick != 1 )
		{
			fprintf( stderr, "okashii, nyuko = 1, shikashi npick = %d\n", npick );
			exit( 1 );
		}
//		fprintf( stderr, "### itchi suru hazu, nazenara scores[nin-1].score=%f, selfscores=%d,%d\n", scores[nin-1].score, scores[nin-1].selfscore, scores->selfscore );
//		fprintf( stderr, "seq[%d] = scores->seq = \n%s\n", scores->numinseq, seq[scores->numinseq] );

		uniform = -1;
		for( j=0; j<nin; j++ ) 
		{
			belongto = 0;
			outs[belongto] = realloc( outs[belongto], sizeof( Scores ) * ( numin[belongto] + 1 ) );
			outs[belongto][numin[belongto]] = scores[j];
			numin[belongto]++;
		}
	}
	else
	{

#if 0
		fprintf( stderr, "yukos = \n" );
		for( i=0; i<nyuko; i++ ) fprintf( stderr, "%d ", y_o_map[i] + 1 );
		fprintf( stderr, "\n" );
#endif
		fprintf( stderr, "\n\n%dx%d distance matrix\n", nyuko, nin );

		for( i=1; i<nyuko; i++ )
		{
			fprintf( stderr, "%d / %d \r", i, nyuko );

			if( doalign )
			{
				if( doalign == 'f' )
				{
					blastresults = callfasta( seq, scores, nin, NULL, yukos[i], (i==1) );
				}
				else
					gappick0( mseq1[0], seq[scores[yukos[i]].numinseq] );
			}
			else
			{
				table1 = (int *)calloc( tsize, sizeof( int ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				makecompositiontable_p( table1, scores[yukos[i]].pointt );
			}
		
			selfscore0 = scores[yukos[i]].selfscore;
			for( j=0; j<nin; j++ ) 
			{
				if( scores[yukos[i]].orilen > scores[j].orilen )
				{
					longer = scores[yukos[i]].orilen;
					shorter = scores[j].orilen;
				}
				else
				{
					shorter = scores[yukos[i]].orilen;
					longer = scores[j].orilen;
				}

#if LENFAC
//				lenfac = 1.0 / ( (double)LENFACA + (double)LENFACB / ( (double)longer + (double)LENFACC ) + (double)shorter / (double)longer * LENFACD );
				lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//				lenfac = 1.0 / ( shorter / longer * LENFACD + LENFACB / ( longer + LENFACC ) + LENFACA );
//				fprintf( stderr, "lenfac = %f, l=%d, %d\n", lenfac, scores[yukos[i]].orilen, scores[j].orilen );
#else
				lenfac = 1.0;
#endif
#if 0 // iihazu -> dame
				ii = s_y_map[j]; jj=s_y_map[yukos[i]];
				if( ii != -1 && jj != -1 )
				{
					if( dfromc[ii][yukos[jj]] != -0.5 )
					{
						dfromc[i][j] = dfromc[ii][yukos[jj]];
					}
					else
					{
						if( ii > jj )
						{
							kk = jj;
							jj = ii;
							ii = kk;
						}
						dfromc[ii][yukos[jj]] = 
						dfromc[i][j] = yukomtx[ii][jj-ii];
					}
				}
				else
#else
				if( dfromc[i][j] == -0.5 )
#endif
				{
					if( doalign )
					{
						if( doalign == 'f' )
						{
							dfromc[i][j] = 
							( 1.0 - blastresults[j] / MIN( selfscore0, scores[j].selfscore ) ) * 1;
							if( dfromc[i][j] < 0.0 ) dfromc[i][j] = 0.0;
						}
						else
						{
							if( fromaln )
							{
								dfromc[i][j] = ( 1.0 - (double)naivepairscore11( orialn[scores[j].numinseq], orialn[scores[yukos[i]].numinseq], penalty ) / MIN( selfscore0, scores[j].selfscore ) ) * 1;
							}
							else
							{
								gappick0( mseq2[0], seq[scores[j].numinseq] );
								dfromc[i][j] = ( 1.0 - (double)G__align11_noalign( n_disLN, -1200, -60, mseq1, mseq2, palloclen ) / MIN( selfscore0, scores[j].selfscore ) ) * 1;
							}
						}
					}
					else
					{
						dfromc[i][j] = ( 1.0 - (double)localcommonsextet_p( table1, scores[j].pointt ) / MIN( selfscore0, scores[j].selfscore ) ) * lenfac;
						if( dfromc[i][j] > MAX6DIST ) dfromc[i][j] = MAX6DIST;
					}
				}
//				fprintf( stderr, "i,j=%d,%d (%d,%d)/ %d,%d, dfromc[][]=%f \n", i, j, scores[yukos[i]].numinseq+1, scores[j].numinseq+1, nyuko, nin, dfromc[i][j] );

//				if( i == 1 )
//					fprintf( stdout, "&&& dfromc[%d][%d] (%d,%d) = %f\n", i, j, p_o_map[i], scores[j].numinseq, dfromc[i][j] );
			}
//			fprintf( stderr, "i=%d, freeing\n", i );
			if( !doalign ) free( table1 );
			if( doalign && doalign == 'f' ) free( blastresults );
		}
		fprintf( stderr, "                \r" );




		for( i=0; i<nyuko; i++ ) numin[i] = 0;
//		fprintf( stderr, "### itchi shinai hazu, nazenara scores[nin-1].score=%f, selfscores=%d,%d, len=%d,%d, nin=%d\n", scores[nin-1].score, scores[nin-1].selfscore, scores->selfscore, scores->orilen, scores[nin-1].orilen, nin );
		for( j=0; j<nin; j++ ) 
		{
#if 0
			belongto = s_y_map[j];
			{
				fprintf( stderr, "belongto = %d (%dinori)\n", belongto, y_o_map[belongto]+1 ); 
			}
			if( belongto == -1 && closeh[j] != -1 )
#endif
#if 0
			if( closeh[j] != -1 )
			{
				belongto = s_y_map[closeh[j]];
//				if( belongto != -1 )
//					fprintf( stderr, "known, %d(%dinori)->%d(%dinori)\n", j, scores[j].numinseq+1, belongto, y_o_map[belongto]+1 );
			}
			else
//			if( belongto == -1 )
#else
			belongto = s_y_map[j];
			if( belongto == -1 )
#endif
			{
				belongto = 0;  // default ha horyu
				minscore = dfromc[0][j];
				for( i=0; i<nyuko; i++ )
				{
//					fprintf( stderr, "checking %d/%d,%d/%d (%d-%d inori) minscore=%f, dfromc[0][j]=%f, dfromc[i][j]=%f\n", i, nyuko, j, nin, y_o_map[i], scores[j].numinseq, minscore, dfromc[0][j], dfromc[i][j] );
					if( scores[j].shimon == scores[yukos[i]].shimon && !strcmp( seq[scores[j].numinseq], seq[y_o_map[i]] ) ) 
					{
//						fprintf( stderr, "yuko-%d (%d in ori) to score-%d (%d inori) ha kanzen itch\n", i, y_o_map[i], j, scores[j].numinseq );
						belongto = i;
						break;
					}
					if( dfromc[i][j] < minscore )
//					if( dfromc[i][j] < minscore && minscore-dfromc[i][j] > ( minscoreinpick[yukos[i]] + minscoreinpick[j] ) * 1.0 )
//					if( rnd() < 0.5 ) // CHUUI !!!!!
					{
//						fprintf( stderr, "yuko-%d (%d in ori) to score-%d (%d inori) ha tikai, %f>%f\n", i, y_o_map[i]+1, j, scores[j].numinseq+1, minscore, dfromc[i][j] );
						minscore = dfromc[i][j];
						belongto = i;
					}
				}
			}
#if 0
			if( dfromc[belongto][j] > minscoreinpick[belongto] )
			{
				fprintf( stderr, "dame, %f > %f\n", dfromc[belongto][j], minscoreinpick[belongto] );
				belongto = npick;
			}
			else
				fprintf( stderr, "ok, %f < %f\n", dfromc[belongto][j], minscoreinpick[belongto] );
#endif
//			fprintf( stderr, "j=%d (%d inori) -> %d (%d inori) d=%f\n", j, scores[j].numinseq+1, belongto, y_o_map[belongto]+1, dfromc[belongto][j] );
//			fprintf( stderr, "numin = %d\n", numin[belongto] );
			outs[belongto] = realloc( outs[belongto], sizeof( Scores ) * ( numin[belongto] + 1 ) );
			outs[belongto][numin[belongto]] = scores[j];
			numin[belongto]++;

		}
		free( dfromcp );
		FreeDoubleMtx( dfromc );

//		fprintf( stderr, "##### npick = %d\n", npick );
//		fprintf( stderr, "##### nyuko = %d\n", nyuko );


		if( nyuko > 2 )
		{
			fprintf( stderr, "upgma       " );
//			veryfastsupg_double_realloc_nobk_halfmtx( nyuko, yukomtx, topol, len );
			fixed_musclesupg_double_realloc_nobk_halfmtx( nyuko, yukomtx, topol, len, NULL, 1, 1 );
			fprintf( stderr, "\r                      \r" );
		}
		else
		{
			topol[0][0] = (int *)realloc( topol[0][0], 2 * sizeof( int ) );
			topol[0][1] = (int *)realloc( topol[0][1], 2 * sizeof( int ) );
			topol[0][0][0] = 0;
			topol[0][0][1] = -1;
			topol[0][1][0] = 1;
			topol[0][1][1] = -1;
		}
		FreeFloatHalfMtx( yukomtx, npick );

#if 0
		ii = nyuko-1;
		fprintf( stderr, "nyuko = %d, topol[][] = \n", nyuko );
		for( j=0; j<nyuko-1; j++ )
		{
			fprintf( stderr, "STEP%d \n", j );
			for( i=0; ; i++ )
			{
				fprintf( stderr, "%d ", ( topol[j][0][i] )+0 );
				if( topol[j][0][i] == -1 ) break;
			}
			fprintf( stderr, "\n" );
			for( i=0; ; i++ )
			{
				fprintf( stderr, "%d ", ( topol[j][1][i] )+0 );
				if( topol[j][1][i] == -1 ) break;
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "\n" );
		}
		exit( 1 );
#endif
		jptr = treeorder; 
		iptr = topol[nyuko-2][0]; while( *iptr != -1 ) *jptr++ = *iptr++;
		iptr = topol[nyuko-2][1]; while( *iptr != -1 ) *jptr++ = *iptr++;
		*jptr++ = -1;
		for( j=0; j<nyuko; j++ )
		{
//			fprintf( stderr, "treeorder[%d] = %d\n", j, treeorder[j] );
			if( treeorder[j] == -1 ) break;
		}
	}

	aligned = 1;
	for( i=0; i<nyuko; i++ )
	{
		ii = treeorder[i];
#if 0
		if( numin[ii] > 1 )
		{
			fprintf( stderr, "\ncalling a child, pick%d (%d inori): # of mem=%d\n", i, p_o_map[ii]+1, numin[ii] );
			for( j=0; j<numin[ii]; j++ )
			{
				fprintf( stderr, "%d ", outs[ii][j].numinseq+1 );
			}
			fprintf( stderr, "\n" );
		}
#endif
		aligned *= splitseq_mq( outs[ii], numin[ii], nlen, seq, orialn, name, inputfile, uniform, children+ii, alloclen, order, whichgroup, weight, depthpt, scores->numinseq );
	}


	for( i=0; i<nyuko; i++ )
	{
		if( !numin[i] )
		{
			fprintf( stderr, "i=%d/%d, ERROR!\n", i, nyuko );
			for( j=0; j<nyuko; j++ )
				fprintf( stderr, "numin[%d] = %d (rep=%d inori)\n", j, numin[j], y_o_map[j] );
			exit( 1 );
		}
	}

#if TREE
	if( treeout )
	{
		treelen = 0;
		for( i=0; i<nyuko; i++ )
			treelen += strlen( children[i] );
		*tree = calloc( treelen + nin * 3, sizeof ( char ) );
	}
#endif


	if( nin >= classsize || !aligned )
		val = 0;
	else
		val = 1;

	if( nyuko > 1 )
	{
		int *mem1p, *mem2p;
		int mem1size, mem2size;
		int v1 = 0, v2 = 0, v3 = 0;
		int nlim;
		int l;
		static int *mem1 = NULL;
		static int *mem2 = NULL;
		char **parttree = NULL; // by Mathog

#if TREE
		if( treeout )
		{
			parttree = (char **)calloc( nyuko, sizeof( char * ) );
			for( i=0; i<nyuko; i++ )
			{
//				fprintf( stderr, "allocating parttree, size = %d\n", treelen + nin * 5 );
				parttree[i] = calloc( treelen + nin * 5, sizeof ( char ) );
				strcpy( parttree[i], children[i] );
				free( children[i] );
			}
			free( children );
		}
#endif
		if( mem1 == NULL )
		{
			mem1 = AllocateIntVec( njob+1 );
			mem2 = AllocateIntVec( njob+1 );
		}

//		veryfastsupg_double_realloc_nobk_halfmtx( nyuko, yukomtx, topol, len );
	
//		counteff_simple_double( nyuko, topol, len, eff );


		nlim = nyuko-1;
		for( l=0; l<nlim; l++ )
		{
			mem1p = topol[l][0];
			mptr = mem1;
			mem1size = 0;
			while( *mem1p != -1 )
			{
//				fprintf( stderr, "*mem1p = %d (%d inori), numin[]=%d\n", *mem1p, p_o_map[*mem1p], numin[*mem1p] );
				i = numin[*mem1p]; ptr = outs[*(mem1p++)];
				mem1size += i;
				while( i-- )
				{
					*mptr++ = (ptr++)->numinseq;
				}
			}
			*mptr = -1;

			mem2p = topol[l][1];
			mptr = mem2;
			mem2size = 0;
			while( *mem2p != -1 )
			{
//				fprintf( stderr, "*mem2p = %d (%d inori), numin[]=%d\n", *mem2p, p_o_map[*mem2p], numin[*mem2p] );
				i = numin[*mem2p]; ptr = outs[*(mem2p++)];
				mem2size += i;
				while( i-- )
				{
					*mptr++ = (ptr++)->numinseq;
				}
			}
			*mptr = -1;

			qsort( mem1, mem1size, sizeof( int ), (int (*)())intcompare );
			qsort( mem2, mem2size, sizeof( int ), (int (*)())intcompare );
//			selhead( mem1, numin[0] );
//			selhead( mem2, numin[1] );


#if 0
			fprintf( stderr, "\n" );
			fprintf( stderr, "mem1 (nin=%d) = \n", nin );
			for( i=0; ; i++ )
			{
				fprintf( stderr, "%d ", mem1[i]+1 );
				if( mem1[i] == -1 ) break;
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "mem2 (nin=%d) = \n", nin );
			for( i=0; ; i++ )
			{
				fprintf( stderr, "%d ", mem2[i]+1 );
				if( mem2[i] == -1 ) break;
			}
			fprintf( stderr, "\n" );
#endif

#if 0
			fprintf( stderr, "before pairalign, l = %d, nyuko=%d, mem1size=%d, mem2size=%d\n", l, nyuko, mem1size, mem2size );
			fprintf( stderr, "before alignment\n" );
			for( j=0; j<mem1size; j++ )
				fprintf( stderr, "%s\n", seq[mem1[j]] );
			fprintf( stderr, "----\n" );
			for( j=0; j<mem2size; j++ )
				fprintf( stderr, "%s\n", seq[mem2[j]] );
			fprintf( stderr, "----\n\n" );
#endif

			if( val )
			{
				fprintf( stderr, "\r  Alignment %d-%d                                 \r", mem1size, mem2size );
				if( *mem1 < *mem2 )
					pairalign( njob, nlen, seq, mem1, mem2, weight, alloclen );
				else
					pairalign( njob, nlen, seq, mem2, mem1, weight, alloclen );
			}

#if TREE
			if( treeout )
			{
				v1 = topol[l][0][0];
				v2 = topol[l][1][0];
	
//				fprintf( stderr, "nyuko=%d, v1=%d, v2=%d\n", nyuko, v1, v2 );
				if( v1 > v2 )
				{
					v3 = v1;
					v1 = v2;
					v2 = v3;
				}
//				fprintf( stderr, "nyuko=%d, v1=%d, v2=%d after sort\n", nyuko, v1, v2 );
//				fprintf( stderr, "nyuko=%d, v1=%d, v2=%d\n", nyuko, v1, v2 );
//				fprintf( stderr, "v1=%d, v2=%d, parttree[v1]=%s, parttree[v2]=%s\n", v1, v2, parttree[v1], parttree[v2] );
				sprintf( *tree, "(%s,%s)", parttree[v1], parttree[v2] );
				strcpy( parttree[v1], *tree );
//				fprintf( stderr, "parttree[%d] = %s\n", v1, parttree[v1] );
//				fprintf( stderr, "*tree = %s\n", *tree );
				free( parttree[v2] ); parttree[v2] = NULL;
			}
#endif

#if 0
			fprintf( stderr, "after alignment\n" );
			for( j=0; j<mem1size; j++ )
				fprintf( stderr, "%s\n", seq[mem1[j]] );
			fprintf( stderr, "----\n" );
			for( j=0; j<mem2size; j++ )
				fprintf( stderr, "%s\n", seq[mem2[j]] );
			fprintf( stderr, "----\n\n" );
#endif
		}
#if TREE
		if( treeout )
		{
			free( parttree[v1] ); parttree[v1] = NULL;
//			fprintf( stderr, "*tree = %s\n", *tree );
//			FreeCharMtx( parttree );
			free( parttree ); parttree = NULL;
		}
#endif

#if 0
		fprintf( stderr, "after alignment\n" );
		for( i=0; i<nyuko; i++ )
		{
			for( j=0; j<numin[i]; j++ )
				fprintf( stderr, "%s\n", seq[outs[i][j].numinseq] );
		}
#endif

		if( val )
		{
			groupid++;
			mptr = mem1; while( *mptr != -1 ) 
			{
#if 0
				fprintf( stdout, "==g1-%d \n", *mptr+1 );
				fprintf( stdout, "%s \n", seq[*mptr] );
#endif
				whichgroup[*mptr] = groupid;
				weight[*mptr++] *= 0.5;
			}
	
			mptr = mem2; while( *mptr != -1 ) 
			{
#if 0
				fprintf( stdout, "=g2-%d ", *mptr+1 );
				fprintf( stdout, "%s \n", seq[*mptr] );
#endif
				whichgroup[*mptr] = groupid;
				weight[*mptr++] *= 0.5;
			}
	
			if( numin[1] == 0 )
			{
				mptr = mem1; while( *mptr != -1 ) 
				{
					whichgroup[*mptr] = groupid;
					weight[*mptr++] *= (double)2.0/numin[0];
				}
			}
		}
		{
			if( *depthpt > maxdepth ) maxdepth = *depthpt;
			(*depthpt)++;
		}
	}
	else
	{
#if TREE
		if( treeout )
		{
			sprintf( *tree, "%s", children[0] );
			free( children[0] );
			free( children );
		}
#endif
	}
	for( i=0; i<npick; i++ ) free( (void *)outs[i] );
//	FreeFloatHalfMtx( pickmtx, npick );
//	FreeFloatHalfMtx( yukomtx, npick );
	FreeFloatMtx( len );
	FreeIntCub( topol );
	FreeIntVec( treeorder );
	free( outs );
	free( numin );
	free( picks );
	free( yukos );
	free( s_p_map );
	free( s_y_map );
	free( p_o_map );
	free( y_o_map );
	free( hanni );
	free( closeh );
	free( pickkouho );
	free( tsukau );
//	free( minscoreinpick );
	return val;
}
#endif

static void alignparaphiles( int nseq, int *nlen, double *weight, char **seq, int nmem, int *members, int *alloclen )
{
	int i, ilim;
	int *mem1 = AllocateIntVec( nmem );
	int *mem2 = AllocateIntVec( 2 );

	mem2[1] = -1;
	ilim = nmem-1;
	for( i=0; i<ilim; i++ )
	{
		mem1[i] = members[i];
		mem1[i+1] = -1;
		mem2[0] = members[i+1];
		pairalign( nseq, nlen, seq, mem1, mem2, weight, alloclen );
	}
	free( mem1 );
	free( mem2 );
}











int main( int argc, char *argv[] )
{
	static char **name, **seq, **orialn;
	static int *grpseq;
	static char *tmpseq;
	static int  **pointt;
	static int *nlen;
	int i, st, en;
	FILE *infp;
	FILE *treefp;
	char *treefile = NULL; //by Mathog
	char c;
	int alloclen;
	static int *order;
	static int *whichgroup;
	static double *weight;
	static char tmpname[B+100];
	int groupnum;
	int groupid;
	int pos;
	int *paramem;
	int npara;
	int completed;
	int orilen;
	int pscore;
	char *pt;
	int **tmpaminodis;
	static char com[1000];
	int depth;
	int aan;

	static Scores *scores;
	static int *table1;
	static char **tree;



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


	if( picksize == NOTSPECIFIED || picksize < 2 )
		picksize = PICKSIZE;

	if( classsize == NOTSPECIFIED || classsize < 0 )
	{
		classsize = njob + 1;
	}
	else
	{
//		picksize = MIN( picksize, (int)sqrt( classsize ) + 1);
	}

	if( tokyoripara == NOTSPECIFIED )
	{
		if( doalign ) tokyoripara = TOKYORIPARA_A;
		else tokyoripara = TOKYORIPARA;
	}
	if( picksize > njob )
		tokyoripara = 0.0;


	alloclen = nlenmax * 2;
	name = AllocateCharMtx( njob, B+1 );

	if( classsize == 1 )
		seq = AllocateCharMtx( njob, 0 );
	else
		seq = AllocateCharMtx( njob, alloclen+1 );


	nlen = AllocateIntVec( njob ); 
	tmpseq = calloc( nlenmax+1, sizeof( char )  );
	pointt = AllocateIntMtx( njob, 0 );
	grpseq = AllocateIntVec( nlenmax + 1 );
	order = (int *)calloc( njob + 1, sizeof( int ) );
	whichgroup = (int *)calloc( njob, sizeof( int ) );
	weight = (double *)calloc( njob, sizeof( double ) );

	fprintf( stderr, "alloclen = %d in main\n", alloclen );

	for( i=0; i<njob; i++ ) whichgroup[i] = 0;
	for( i=0; i<njob; i++ ) weight[i] = 1.0;
	for( i=0; i<njob; i++ ) order[i] = -1;

	if( classsize == 1 )
		readData_varlen( infp, name, nlen, seq );
	else
		readData_pointer( infp, name, nlen, seq );

	fclose( infp );

	if( fromaln ) doalign = 1;

	if( fromaln )
	{
		orialn = AllocateCharMtx( njob, alloclen+1 );
		for( i=0; i<njob; i++ )
		{
			if( strlen( seq[i] ) != nlenmax )
			{
				fprintf( stderr, "Input sequences must be aligned\n" );
				exit( 1 );
			}
			strcpy( orialn[i], seq[i] );
		}
	}

	constants( njob, seq );

	if( dorp == 'd' ) tsize = (int)pow( 4, 6 );
	else              tsize = (int)pow( 6, 6 );

	if( dorp == 'd' )
	{
		lenfaca = DLENFACA;
		lenfacb = DLENFACB;
		lenfacc = DLENFACC;
		lenfacd = DLENFACD;
	}
	else
	{
		lenfaca = PLENFACA;
		lenfacb = PLENFACB;
		lenfacc = PLENFACC;
		lenfacd = PLENFACD;
	}

	maxl = 0;
	for( i=0; i<njob; i++ ) 
	{
		gappick0( tmpseq, seq[i] );
		nlen[i] = strlen( tmpseq );
		strcpy( seq[i], tmpseq );
		pointt[i] = AllocateIntVec( nlen[i]+1 ); // ??

		if( nlen[i] < 6 )
		{
			fprintf( stderr, "Seq %d, too short, %d characters\n", i+1, nlen[i] );
			fprintf( stderr, "name = %s\n", name[i] );
			fprintf( stderr, "seq = %s\n", seq[i] );
			exit( 1 );
//			continue;
		}
		if( nlen[i] > maxl ) maxl = nlen[i];
		if( dorp == 'd' ) /* nuc */
		{
			if( seq_grp_nuc( grpseq, tmpseq ) < 6 )
			{
				fprintf( stderr, "Seq %d, too short.\n", i+1 );
				fprintf( stderr, "name = %s\n", name[i] );
				fprintf( stderr, "seq = %s\n", seq[i] );
				exit( 1 );
//				continue;
			}
			makepointtable_nuc( pointt[i], grpseq );
		}
		else                 /* amino */
		{
			if( seq_grp( grpseq, tmpseq ) < 6 )
			{
				fprintf( stderr, "Seq %d, too short.\n", i+1 );
				fprintf( stderr, "name = %s\n", name[i] );
				fprintf( stderr, "seq = %s\n", seq[i] );
				exit( 1 );
//				continue;
			}
			makepointtable( pointt[i], grpseq );
		}
//		fprintf( stdout, ">%s\n", name[i] );
//		fprintf( stdout, "%s\n", seq[i] );
	}
	if( nunknown ) fprintf( stderr, "\nThere are %d ambiguous characters\n", nunknown );
//	exit( 1 );

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

	pid = (int)getpid();
	sprintf( datafile, "/tmp/data-%d", pid );
	sprintf( queryfile, "/tmp/query-%d", pid );
	sprintf( resultfile, "/tmp/fasta-%d", pid );

	scores = (Scores *)calloc( njob, sizeof( Scores ) );

//	fprintf( stderr, "\nCalculating i-i scores ... \n" );
	for( i=0; i<njob; i++ ) 
	{
		orilen = strlen( seq[i] );
		scores[i].numinseq = i; // irimasu
		scores[i].orilen = orilen;
	}

	if( doalign == 'f' )
	{
		fastapath = getenv( "FASTA_4_MAFFT" );
		if( !fastapath ) fastapath = "fasta34";
		fprintf( stderr, "fastapath=%s\n", fastapath );
		tmpaminodis = AllocateIntMtx( 0x80, 0x80 );
		getfastascoremtx( tmpaminodis );
	}
	else
		tmpaminodis = NULL;
	
	for( i=0; i<njob; i++ )
	{
		scores[i].pointt = pointt[i];
		scores[i].shimon = (int)pointt[i][0] + (int)pointt[i][1] + (int)pointt[i][scores[i].orilen-6];
		scores[i].name = name[i];
		if( doalign )
		{
			fprintf( stderr, "\r %05d/%05d   ", i, njob );
			free( scores[i].pointt );
			if( doalign == 'f' )
			{
#if 0
#define KIZAMI 100
				int ipos = (int)( i / KIZAMI ) * KIZAMI;
				int iposamari = i % KIZAMI;

				fprintf( stderr, "%d / %d\r", i, njob );
//				fprintf( stderr, "ipos = %d\n", ipos );
//				fprintf( stderr, "iposamari = %d\n", iposamari );

//				fprintf( stderr, " calling blast, i=%d\n", i );
//				blastresults = callfasta( seq, scores+i, 1, 0, 1 );
				blastresults = callfasta( seq, scores+ipos, MIN(KIZAMI,njob-ipos), NULL, iposamari, (iposamari==0)  );
//				fprintf( stderr, "done., i=%d\n\n", i );
				scores[i].selfscore = (int)blastresults[iposamari]; 
#if 0
				for( j=0; j<100; j++ )
				{
					fprintf( stderr, "res[%d] = %f\n", j, blastresults[j] );
				}
#endif
//				fprintf( stderr, "%d->selfscore = %d\n", i, scores[i].selfscore );
				free( blastresults );
#else
				pscore = 0;
				if( scoremtx == -1 )
				{
					st = 1;
					en = 0;
					for( pt=seq[i]; *pt; pt++ )
					{
						if( *pt == 'u' ) *pt = 't';
						aan = amino_n[(int)*pt];
						if( aan<0 || aan >= 4 ) *pt = 'n';

						if( *pt == 'n' ) 
						{
							en++;
							if( st ) continue;
							else pscore += tmpaminodis[(int)*pt][(int)*pt];
						}
						else
						{
							st = 0;
							en = 0;
							pscore += tmpaminodis[(int)*pt][(int)*pt];
						}
					}
					scores[i].selfscore = pscore - en * tmpaminodis['n']['n']; 
				}
				else
				{
					st = 1;
					en = 0;
					for( pt=seq[i]; *pt; pt++ )
					{
						aan = amino_n[(int)*pt];
						if( aan<0 || aan >= 20 ) *pt = 'X';
						if( *pt == 'X' ) 
						{
							en++;
							if( st ) continue;
							else pscore += tmpaminodis[(int)*pt][(int)*pt];
						}
						else
						{
							st = 0;
							en = 0;
							pscore += tmpaminodis[(int)*pt][(int)*pt];
						}
					}
					scores[i].selfscore = pscore - en * tmpaminodis['X']['X']; 
				}
#endif
			}
			else
			{
				pscore = 0;
				for( pt=seq[i]; *pt; pt++ )
				{
//					pscore += amino_dis[(int)*pt][(int)*pt];
					pscore += amino_dis[(int)*pt][(int)*pt];
				}
				scores[i].selfscore = pscore; 
			}
//			fprintf( stderr, "selfscore[%d] = %d\n", i+1, scores[i].selfscore );
		}
		else
		{
			table1 = (int *)calloc( tsize, sizeof( int ) );
			if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
			makecompositiontable_p( table1, pointt[i] );
			scores[i].selfscore = localcommonsextet_p( table1, pointt[i] );
			free( table1 );
		}
	}
	if( tmpaminodis ) FreeIntMtx( tmpaminodis );

	depth = 0;
#if TREE
	if( treeout )
	{
		tree = (char **)calloc( 1, sizeof( char *) );
		*tree = NULL;
//		splitseq_bin( scores, njob, nlen, seq, name, inputfile, 0, tree, &alloclen, order, whichgroup, weight );
		completed = splitseq_mq( scores, njob, nlen, seq, orialn, name, inputfile, 0, tree, &alloclen, order, whichgroup, weight, &depth, -1 );
		treefile = (char *)calloc( strlen( inputfile ) + 10, sizeof( char ) );
		if( inputfile )
			sprintf( treefile, "%s.tree", inputfile );
		else
			sprintf( treefile, "splittbfast.tree" );
		treefp = fopen( treefile, "w" );
		fprintf( treefp, "%s\n", *tree );
		fclose( treefp );
	}
	else
		completed = splitseq_mq( scores, njob, nlen, seq, orialn, name, inputfile, 0, tree, &alloclen, order, whichgroup, weight, &depth, -1 );
#else
	completed = splitseq_mq( scores, njob, nlen, seq, orialn, name, inputfile, 0, tree, &alloclen, order, whichgroup, weight, &depth, -1 );
#endif

	fprintf( stderr, "\nDone.\n\n" );

#if 1
	groupnum = 0;
	groupid = -1;
	paramem = NULL;
	npara = 0;
	for( i=0; i<njob; i++ )
	{
		pos = order[i];
		if( whichgroup[pos] != groupid )
		{
			groupnum++;
			groupid = whichgroup[pos];
		}
		if( whichgroup[pos] )
		{
			if( paramem )
			{
				paramem[npara] = -1;
				if( npara > 1 && classsize > 2 ) 
				{
					qsort( paramem, npara, sizeof( int ), (int (*)(const void *, const void*))intcompare );
//					selhead( paramem, npara );
					alignparaphiles( njob, nlen, weight, seq, npara, paramem, &alloclen );
				}
				free( paramem ); paramem = NULL; npara = 0;
			}
			sprintf( tmpname, "Group-%d %s", groupnum, name[pos]+1 );
		}
		else
		{
			paramem = realloc( paramem, sizeof( int) * ( npara + 2 ) );
			paramem[npara++] = pos;
			sprintf( tmpname, "Group-para %s", name[pos]+1 );
		}
		tmpname[B-1] = 0;
		if( classsize > 1 && classsize <= njob )
			strcpy( name[pos]+1, tmpname );
	}
	if( paramem )
	{
		paramem[npara] = -1;
		if( npara > 1 && classsize > 2 ) 
		{
			qsort( paramem, npara, sizeof( int ), (int (*)(const void *, const void*))intcompare );
//			selhead( paramem, npara );
			alignparaphiles( njob, nlen, weight, seq, npara, paramem, &alloclen );
		}
		free( paramem ); paramem = NULL; npara = 0;
	}
#else
	for( i=0; i<njob; i++ )
	{
		sprintf( tmpname, "Group-%d %s", whichgroup[i], name[i]+1 );
		strcpy( name[i]+1, tmpname );
	}
#endif


//	maketanni( name, seq,  njob, nlenmax, nlen );

	fclose( trap_g );

#if DEBUG
	fprintf( stderr, "writing alignment to stdout\n" );
#endif
	if( reorder )
		writeData_reorder_pointer( stdout, njob, name, nlen, seq, order );
	else
		writeData_pointer( stdout, njob, name, nlen, seq );
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif
	if( classsize == 1 )
	{
		fprintf( stderr, "\n\n" );
		fprintf( stderr, "----------------------------------------------------------------------------\n" );
		fprintf( stderr, "\n" );
		fprintf( stderr, "nseq = %d\n", njob );
		fprintf( stderr, "groupsize = %d, picksize=%d\n", classsize, picksize );
		fprintf( stderr, "The input sequences have been sorted so that similar sequences are close.\n" );
		if( reorder )
			fprintf( stderr, "The order of sequences has been changed according to estimated similarity.\n" );
#if TREE
		if( treeout )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "A guide tree is in the '%s' file.\n", treefile );
		}
//		else
//		{
//			fprintf( stderr, "To output guide tree,\n" );
//			fprintf( stderr, "%% %s -t  -i %s\n", progName( argv[0] ), "inputfile" );
//		}
#endif
		if( !doalign )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "mafft --dpparttree might give a better result, although slow.\n" );
			fprintf( stderr, "mafft --fastaparttree is also available if you have fasta34.\n" );
		}
		fprintf( stderr, "\n" );
		fprintf( stderr, "----------------------------------------------------------------------------\n" );
	}
	else if( groupnum > 1 )
	{
		fprintf( stderr, "\n\n" );
		fprintf( stderr, "----------------------------------------------------------------------------\n" );
		fprintf( stderr, "\n" );
		fprintf( stderr, "groupsize = %d, picksize=%d\n", classsize, picksize );
		fprintf( stderr, "The input sequences have been classified into %d groups + some paraphyletic groups\n", groupnum );
		fprintf( stderr, "Note that the alignment is not completed.\n" );
		if( reorder )
			fprintf( stderr, "The order of sequences has been changed according to estimated similarity.\n" );
#if TREE
		if( treeout )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "A guide tree is in the '%s' file.\n", treefile );
		}
//		else
//		{
//			fprintf( stderr, "To output guide tree,\n" );
//			fprintf( stderr, "%% %s -t  -i %s\n", progName( argv[0] ), "inputfile" );
//		}
#endif
		if( !doalign )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "mafft --dpparttree might give a better result, although slow.\n" );
			fprintf( stderr, "mafft --fastaparttree is also available if you have fasta34.\n" );
		}
		fprintf( stderr, "\n" );
		fprintf( stderr, "----------------------------------------------------------------------------\n" );
	}			
	else
	{
		fprintf( stderr, "\n\n" );
		fprintf( stderr, "----------------------------------------------------------------------------\n" );
		fprintf( stderr, "\n" );
		fprintf( stderr, "nseq = %d\n", njob );
		fprintf( stderr, "groupsize = %d, partsize=%d\n", classsize, picksize );
//		fprintf( stderr, "A single alignment containing all the input sequences has been computed.\n" );
//		fprintf( stderr, "If the sequences are highly diverged and you feel there are too many gaps,\n" );
//		fprintf( stderr, "please try \n" );
//		fprintf( stderr, "%% mafft --parttree --groupsize 100 inputfile\n" );
//		fprintf( stderr, "which classifies the sequences into several groups with <~ 100 sequences\n" );
//		fprintf( stderr, "and performs only intra-group alignments.\n" );
		if( reorder )
			fprintf( stderr, "The order of sequences has been changed according to estimated similarity.\n" );
#if TREE
		if( treeout )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "A guide tree is in the '%s' file.\n", treefile );
		}
//		else
//		{
//			fprintf( stderr, "To output guide tree,\n" );
//			fprintf( stderr, "%% %s -t  -i %s\n", progName( argv[0] ), "inputfile" );
//		}
#endif
		if( !doalign || fromaln )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "mafft --dpparttree might give a better result, although slow.\n" );
			fprintf( stderr, "mafft --fastaparttree is also available if you have fasta34.\n" );
		}
		fprintf( stderr, "\n" );
		fprintf( stderr, "----------------------------------------------------------------------------\n" );
	}
#if TREE
	if( treeout ) free( treefile );
#endif

#if 0
	fprintf( stdout, "weight =\n" );
	for( i=0; i<njob; i++ )
		fprintf( stdout, "%d: %f\n", i+1, weight[i] );
#endif

	if( doalign == 'f' )
	{
		strcpy( com, "rm -f" );
		strcat( com, " " );
		strcat( com, datafile );
		strcat( com, "*  " );
		strcat( com, queryfile );
		strcat( com, " " );
		strcat( com, resultfile );
		fprintf( stderr, "%s\n", com );
		system( com );
	}

	SHOWVERSION;

	return( 0 );
}
