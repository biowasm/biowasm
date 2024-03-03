#include "mltaln.h"


#define REPORTCOSTS 0

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SKIP 1

#define ITERATIVECYCLE 2

#define END_OF_VEC -1

static int treein;
static int topin;
static int treeout;
static int noalign;
static int distout;
static int tuplesize;
static int subalignment;
static int subalignmentoffset;
static int nguidetree;
static int sparsepickup;
static int keeplength;
static int ndeleted;
static int mapout;
static int smoothing;
static double maxdistmtxsize;
static int nthreadtb;
static int useexternalanchors;
static int oneiteration;
static double maxanchorseparation;

#if 0
#define PLENFACA 0.0123
#define PLENFACB 10252
#define PLENFACC 10822
#define PLENFACD 0.5
#define DLENFACA 0.01
#define DLENFACB 2445
#define DLENFACC 2412
#define DLENFACD 0.1
#else
#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define D6LENFACA 0.01
#define D6LENFACB 2500
#define D6LENFACC 2500
#define D6LENFACD 0.1
#define D10LENFACA 0.01
#define D10LENFACB 1000000
#define D10LENFACC 1000000
#define D10LENFACD 0.0
#endif

typedef struct _jobtable
{
    int i;  
    int j;  
} Jobtable;

typedef struct _msacompactdistmtxthread_arg
{
	int njob;
	int thread_no;
	int *selfscore;
	double **partmtx;
	char **seq;
	int **skiptable;
	double *mindist;
	int *mindistfrom;
 	int *jobpospt;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msacompactdistmtxthread_arg_t;

typedef struct _compactdistmtxthread_arg
{
	int njob;
	int thread_no;
	int *nogaplen;
	int **pointt;
	int *selfscore;
	double **partmtx;
	int *jobpospt;
	double *mindist;
	int *mindistfrom;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} compactdistmtxthread_arg_t;

typedef struct _msadistmtxthread_arg
{
	int njob;
	int thread_no;
	int *selfscore;
	double **iscore;
	double **partmtx;
	char **seq;
	int **skiptable;
	Jobtable *jobpospt;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msadistmtxthread_arg_t;

#ifdef enablemultithread
// ue futatsu ha singlethread demo tsukau
typedef struct _treebasethread_arg
{
	int thread_no;
	int njob;
	int *nrunpt;
	int *nlen;
	int *jobpospt;
	int ***topol;
	Treedep *dep;
	double ***cpmxhist;
	int **memhist;
	char **aseq;
	double *effarr;
	int *alloclenpt;
	int *fftlog;
	char *mergeoralign;
	double **newdistmtx;
	int *selfscore;
	ExtAnch *extanch;
	int **anchindex;
	pthread_mutex_t *mutex;
	pthread_cond_t *treecond;
} treebasethread_arg_t;

typedef struct _distancematrixthread_arg
{
	int thread_no;
	int njob;
	int *jobpospt;
	int **pointt;
	double **mtx;
	pthread_mutex_t *mutex;
} distancematrixthread_arg_t;
#endif


void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	nthreadpair = 1;
	nthreadtb = 1;
	outnumber = 0;
	topin = 0;
	treein = 0;
	treeout = 0;
	distout = 0;
	noalign = 0;
	nevermemsave = 0;
	inputfile = NULL;
	nadd = 0;
	addprofile = 1;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
	useexternalanchors = 0;
	oneiteration = 0;
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
	ppenalty_dist = NOTSPECIFIED;
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
	scoreout = 0;
	spscoreout = 0;
	tuplesize = 6;
	subalignment = 0;
	subalignmentoffset = 0;
	legacygapcost = 0;
	specificityconsideration = 0.0;
	nguidetree = 1;
	sparsepickup = 0;
	keeplength = 0;
	mapout = 0;
	smoothing = 0;
	nwildcard = 0;
	maxanchorseparation = 1000.0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					reporterr(       "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
					reporterr(       "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'V':
					ppenalty_dist = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					reporterr(       "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					reporterr(       "blosum %d / kimura 200 \n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					reporterr(       "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					reporterr(       "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'C':
					nthreadpair = nthread = myatoi( *++argv );
					reporterr(       "nthread = %d\n", nthread );
					reporterr(       "nthreadpair = %d\n", nthread );
					if( strchr( *argv, '-' ) )
						nthreadtb = myatoi( strchr( *argv, '-' )+1 );
					else
						nthreadtb = nthread;
					reporterr(       "nthreadtb = %d\n", nthreadtb );
					--argc; 
					goto nextoption;
				case 's':
					specificityconsideration = (double)myatof( *++argv );
//					reporterr(       "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
				case 'y':
					distout = 1;
					break;
				case 't':
					treeout = 1;
					break;
				case '^':
					treeout = 2;
					break;
				case 'T':
					noalign = 1;
					break;
				case 'r':
					oneiteration = 1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'L':
					legacygapcost = 1;
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'x':
					maxanchorseparation = myatof( *++argv );
					--argc; 
					goto nextoption;
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'R':
					fftRepeatStop = 1;
					break;
#endif
				case 'n' :
					outnumber = 1;
					break;
#if 0
				case 's':
					treemethod = 's';
					break;
				case 'q':
					treemethod = 'q'; // minimum
					break;
#endif
				case 'q':
					sparsepickup = myatoi( *++argv );
//					reporterr(       "sparsepickup = %d\n", sparsepickup );
					--argc; 
					goto nextoption;
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
//					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
				case 'E':
					nguidetree = myatoi( *++argv );
//					reporterr(       "nguidetree = %d\n", nguidetree );
					--argc; 
					goto nextoption;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'R':
					alg = 'R';
					break;
#endif
				case 'A':
					alg = 'A';
					break;
				case '&':
					alg = 'a';
					break;
				case '@':
					alg = 'd';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'M':
					alg = 'M';
					break;
#if 0
				case 'S' :
					scoreout = 1; // for checking parallel calculation
					break;
#else
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
#endif
				case 'B': // hitsuyou! memopt -M -B no tame
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'l':
					useexternalanchors = 1;
				case 'G':
					use_fft = 1;
					force_fft = 1;
					break;
#if 0
				case 'V':
					topin = 1;
					break;
#endif
				case 'U':
					treein = 1;
					break;
				case 'u':
					weight = 0;
					tbrweight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
#if 1
				case 'd':
					disp = 1;
					break;
#endif
#if 1
				case 'O':
					outgap = 0;
					break;
#else
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
				case 'J':
					tbutree = 0;
					break;
#if 0
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
#endif
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'W':
					tuplesize = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'Z':
					checkC = 1;
					break;
#endif
				case 'Y':
					keeplength = 1;
					break;
				case 'z':
					mapout = 2;
					break;
				case 'Z':
					mapout = 1;
					break;
				case 'p':
					smoothing = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
                default:
                    reporterr(       "illegal option %c\n", c );
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
        reporterr(       "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		reporterr(       "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}

static int varpairscore( int nseq, int npick, int nlenmax, char **seq, int seed )
{
	int i, j, npair;
	int *slist;
	char **pickseq;
	double score;
	double scoreav;
	double scoreav2;
	double scorestd;
	double scorevar;
	slist = calloc( nseq, sizeof( int ) );
	pickseq = AllocateCharMtx( npick, nlenmax );
	reporterr( "nseq = %d, nlenmax=%d, seed=%d\n", nseq, nlenmax, seed );

	srand( seed );

	for( i=0; i<nseq; i++ ) slist[i] = i;
//	for( i=0; i<nseq; i++ ) reporterr( "slist[%d] = %d\n", i, slist[i] );

	stringshuffle( slist, nseq );
	for( i=0; i<npick; i++ ) gappick0( pickseq[i], seq[slist[i]] );

	scoreav = 0.0;
	scoreav2 = 0.0;
	npair = npick * (npick-1) / 2;
	for( i=1; i<npick; i++ ) 
	{
		reporterr( "%d / %d\r", i, npick );
		for( j=0; j<i; j++ ) 
		{
			score = G__align11_noalign( n_dis_consweight_multi, -1200, -60, pickseq+i, pickseq+j, nlenmax );
			scoreav += score;
			scoreav2 += score * score;
			printf( "score = %d\n", (int)score );
		}
	}

	scoreav /= (double)npair;
	scoreav2 /= (double)npair;
	scorevar = ( scoreav2 - scoreav * scoreav )*npair/(npair-1);
	scorestd = sqrt( scorevar );
	printf( "av = %f\n", scoreav );
	printf( "stddev = %f\n", scorestd );
	printf( "cv = %f\n", scorestd/scoreav );

	FreeCharMtx( pickseq );

	if( scorestd/scoreav < 0.2 ) return( 's' );
	else return( 't' );
}

static void pickup( int n, int *seqlen, int ***topol, char **name, char **seq ) // memsave ni mitaiou
{
	int i, j, k, m;
	int **longestseq;
	int **longestlen;
	int *select;
	char **nameout, **seqout;
	int *nlenout;
	char **namenotused, **seqnotused;
	int *nlennotused;
	FILE *notusedfp;

	longestseq = AllocateIntMtx( n-1, 2 );
	longestlen = AllocateIntMtx( n-1, 2 );
	select = AllocateIntVec( n );
	for( i=0; i<n; i++ ) select[i] = 0;
	nameout = AllocateCharMtx( n, 0 );
	seqout = AllocateCharMtx( n, 0 );
	nlenout = AllocateIntVec( n );
	namenotused = AllocateCharMtx( n, 0 );
	seqnotused = AllocateCharMtx( n, 0 );
	nlennotused = AllocateIntVec( n );

	for( i=0; i<n-1; i++ )
	{
//		reporterr( "STEP %d\n", i );
		longestlen[i][0] = -1;
		longestseq[i][0] = -1;
		for( j=0; (m=topol[i][0][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][0] )
			{
				longestlen[i][0] = seqlen[m];
				longestseq[i][0] = m;
			}
//			reporterr( "%d ", topol[i][0][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][0], longestseq[i][0] );


		longestlen[i][1] = -1;
		longestseq[i][1] = -1;
		for( j=0; (m=topol[i][1][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][1] )
			{
				longestlen[i][1] = seqlen[m];
				longestseq[i][1] = m;
			}
//			reporterr( "%d ", topol[i][1][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][1], longestseq[i][1] );
	}

	m = 1;
	for( i=n-2; i>-1; i-- )
	{
//		reporterr( "longest[%d][0] = %d (%d)\n", i, longestlen[i][0], longestseq[i][0] );
//		reporterr( "longest[%d][1] = %d (%d)\n", i, longestlen[i][1], longestseq[i][1] );
		select[longestseq[i][0]] = 1;
		select[longestseq[i][1]] = 1;
		m += 1;
		if( m >= sparsepickup ) break;
	}
	for( i=0, k=0, j=0; i<n; i++ ) 
	{
		if( select[i] )
		{
			nameout[k] = name[i];
			seqout[k] = seq[i];
			nlenout[k] = strlen( seqout[k] );
			k++;
		}
		else
		{
			namenotused[j] = name[i];
			seqnotused[j] = seq[i];
			nlennotused[j] = strlen( seqnotused[j] );
			j++;
		}
	}
	writeData_pointer( stdout, m, nameout, nlenout, seqout );

	notusedfp = fopen( "notused", "w" );
	writeData_pointer( notusedfp, n-m, namenotused, nlennotused, seqnotused );
	fclose( notusedfp );


	free( nameout );
	free( nlenout );
	free( seqout );
	free( namenotused );
	free( nlennotused );
	free( seqnotused );
	FreeIntMtx( longestseq );
	FreeIntMtx( longestlen );
	free( select );
}


static int nunknown = 0;

void seq_grp_nuc( int *grp, char *seq )
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
	if( grp - grpbk < tuplesize )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(unsigned char)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void makecompositiontable_p( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
	{
#if 1
		table[point]++;
#else // kakunin shinai
		if( (unsigned int)table[point]++ >= INT_MAX )
		{
			reporterr( "Overflow. table[point]=%d>INT_MAX(%d).\n", table[point], INT_MAX );
			exit( 1 );
		}
#endif
	}
}


void makepointtable_nuc_dectet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *262144;
	point += *n++ * 65536;
	point += *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ *262144;
		point *= 4;
		point += *n++;
		*pointt++ = point;

	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc_octet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 16384;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

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

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

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

static double preferenceval( int ori, int pos, int max ) // for debug
{
	pos -= ori;
	if( pos < 0 ) pos += max;
	return( 0.00000000000001 * pos );
}

static void *compactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistmtxthread_arg_t *targ = (compactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	int *nogaplen = targ->nogaplen;
	int **pointt = targ->pointt;
 	int *jobpospt = targ->jobpospt;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
	int i, j;
	double tmpdist, preference, tmpdistx; //, tmpdisty;
	int *table1;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == -1 )
			{
				pthread_mutex_unlock( targ->mutex );
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i-1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif 
		{
			i = *jobpospt;
			if( i == -1 )
			{
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i-1;
		}

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			if( nthreadpair )
				reporterr(       "\r% 5d / %d (thread %4d)", njob-i, njob, thread_no );
			else
				reporterr(       "\r% 5d / %d", njob-i, njob );
		}
		makecompositiontable_p( table1, pointt[i] );

//		for( j=i+1; j<njob; j++ ) 
		for( j=i-1; j>-1; j-- ) 
		{

			tmpdist = distcompact( nogaplen[i], nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j] );
			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

//			preference = preferenceval( j, i, njob );
//			tmpdisty = tmpdist + preference;
//			if( tmpdisty < mindist[j] )
//			{
//				mindist[j] = tmpdisty;
//				mindistfrom[j] = i;
//			}

			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		} 
		free( table1 );
	}
}

static void *ylcompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistmtxthread_arg_t *targ = (compactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	int *nogaplen = targ->nogaplen;
	int **pointt = targ->pointt;
 	int *jobpospt = targ->jobpospt;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
	int i, j;
	double tmpdist, preference, tmpdistx, tmpdisty;
	int *table1;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i+1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif 
		{
			i = *jobpospt;
			if( i == njob-1 )
			{
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			*jobpospt = i+1;
		}

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			if( nthreadpair )
				reporterr(       "\r% 5d / %d (thread %4d)", i+1, njob, thread_no );
			else
				reporterr(       "\r% 5d / %d", i+1, njob );
		}
		makecompositiontable_p( table1, pointt[i] );

		for( j=i+1; j<njob; j++ ) 
//		for( j=i-1; j>-1; j-- ) 
		{

			tmpdist = distcompact( nogaplen[i], nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j] );
			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

			preference = preferenceval( j, i, njob );
			tmpdisty = tmpdist + preference;
			if( tmpdisty < mindist[j] )
			{
				mindist[j] = tmpdisty;
				mindistfrom[j] = i;
			}

			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		} 
		free( table1 );
	}
}


static void *msacompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msacompactdistmtxthread_arg_t *targ = (msacompactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
 	int *jobpospt = targ->jobpospt;
	double tmpdist, preference, tmpdistx; //, tmpdisty;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == -1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			*jobpospt = i-1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			i = *jobpospt;
			if( i == -1 )
			{
				return( NULL );
			}
			*jobpospt = i-1;
		}

		if( i % 100 == 0 ) 
		{
			if( nthreadpair )
				fprintf( stderr, "\r% 5d / %d (thread %4d)", njob-i, njob, thread_no );
			else
				fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		for( j=i-1; j>-1; j-- ) 
//		for( j=i+1; j<njob; j++ ) 
		{
			tmpdist = distcompact_msa( seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j] ); // osoikedo,

			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

//			preference = preferenceval( j, i, njob );
//			tmpdisty = tmpdist + preference;
//			if( tmpdisty < mindist[j] )
//			{
//				mindist[j] = tmpdisty;
//				mindistfrom[j] = i;
//			}
			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		}
	}
}

static void *ylmsacompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msacompactdistmtxthread_arg_t *targ = (msacompactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
 	int *jobpospt = targ->jobpospt;
	double tmpdist, preference, tmpdistx, tmpdisty;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			*jobpospt = i+1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			i = *jobpospt;
			if( i == njob-1 )
			{
				return( NULL );
			}
			*jobpospt = i+1;
		}

		if( i % 100 == 0 ) 
		{
			if( nthreadpair )
				fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
			else
				fprintf( stderr, "\r% 5d / %d", i, njob );
		}

//		for( j=i-1; j>-1; j-- ) 
		for( j=i+1; j<njob; j++ ) 
		{
			tmpdist = distcompact_msa( seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j] ); // osoikedo,

			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

			preference = preferenceval( j, i, njob );
			tmpdisty = tmpdist + preference;
			if( tmpdisty < mindist[j] )
			{
				mindist[j] = tmpdisty;
				mindistfrom[j] = i;
			}
			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		}
	}
}

#if 1
static void *msadistmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msadistmtxthread_arg_t *targ = (msadistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;


	double ssi, ssj, bunbo, iscoretmp;
	int i, j;
	int nlim = njob-1;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = jobpospt->i; // (jobpospt-i)++ dato, shuuryou hantei no mae ni ++ surunode, tomaranakunaru.

			if( i == nlim )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			jobpospt->i += 1;
			pthread_mutex_unlock( targ->mutex );
			if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		}
		else
#endif
		{
			i = (jobpospt->i)++;
			if( i == nlim ) return( NULL );
			if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		ssi = selfscore[i];
		for( j=i+1; j<njob; j++ )
		{
			ssj = selfscore[j];
			bunbo = MIN( ssi, ssj );
//fprintf( stderr, "bunbo = %f\n", bunbo );
//fprintf( stderr, "naivepairscorefast() = %f\n", naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) );
			if( bunbo == 0.0 )
				iscoretmp = 2.0; // 2013/Oct/17
			else
			{
				iscoretmp = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast 
				if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17
	
			}
			if( iscoretmp < 0.0 ) 
			{
				reporterr( "WARNING: negative distance, iscoretmp = %f\n", iscoretmp );
				iscoretmp = 0.0;
			}
			iscore[i][j-i] = iscoretmp;
//			printf( "i,j=%d,%d, iscoretmp=%f\n", i, j, iscoretmp );

		}
	}
}
#else
static void *msadistmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msadistmtxthread_arg_t *targ = (msadistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;


	double ssi, ssj, bunbo, iscoretmp;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthreadpair ) pthread_mutex_lock( targ->mutex );
#endif
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
#ifdef enablemultithread
				if( nthreadpair ) pthread_mutex_unlock( targ->mutex );
#endif
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
#ifdef enablemultithread
		if( nthreadpair ) pthread_mutex_unlock( targ->mutex );
#endif


		if( nthreadpair )
		{
			if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		}
		else
		{
			if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d", i, njob );
		}
		ssi = selfscore[i];
		ssj = selfscore[j];
		bunbo = MIN( ssi, ssj );
//fprintf( stderr, "bunbo = %f\n", bunbo );
//fprintf( stderr, "naivepairscorefast() = %f\n", naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) );
		if( bunbo == 0.0 )
			iscoretmp = 2.0; // 2013/Oct/17
		else
		{
			iscoretmp = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast 
			if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17

		}
		iscore[i][j-i] = iscoretmp;


	}
}
#endif

#ifdef enablemultithread
static void *distancematrixthread( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int *jobpospt = targ->jobpospt;
	int **pointt = targ->pointt;
	double **mtx = targ->mtx;

	int *table1;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = *jobpospt;
		if( i == njob )
		{
			pthread_mutex_unlock( targ->mutex );
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		*jobpospt = i+1;
		pthread_mutex_unlock( targ->mutex );

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 100 == 0 )
		{
			reporterr(       "\r% 5d / %d (thread %4d)", i+1, njob, thread_no );
		}
		makecompositiontable_p( table1, pointt[i] );

		for( j=i; j<njob; j++ ) 
		{
			mtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] );
		} 
		free( table1 );
	}
}

static void recountpositions( ExtAnch *pairanch, int n1, int n2, char **seq1, char **seq2 ) // loop no junban kentou
{
	int i, j, k, len, pos;
	int *map;

	len = strlen( seq1[0] )+1;
	map = calloc( sizeof( int ), len );

	for( k=0; k<n1; k++ )
	{
		pos = 0;
		for( i=0; i<len; i++ )
		{
			if( seq1[k][i] != '-' )
			{
				map[pos] = i;
				pos++;
			}
		}

		for( j=0; pairanch[j].i>-1; j++ )
		{
			if( pairanch[j].i == k ) 
			{
//				reporterr( "pairanch[%d].endi: %d->%d\n", j, pairanch[j].endi, map[pairanch[j].endi] );
				pairanch[j].starti = map[pairanch[j].starti];
				pairanch[j].endi = map[pairanch[j].endi];
			}
		}
	}
	free( map );

	len = strlen( seq2[0] )+1;
	map = calloc( sizeof( int ), len );
	for( k=0; k<n2; k++ )
	{
		pos = 0;
		for( i=0; i<len; i++ )
		{
			if( seq2[k][i] != '-' )
			{
				map[pos] = i;
				pos++;
			}
		}
		for( j=0; pairanch[j].i>-1; j++ )
		{
			if( pairanch[j].j == k )
			{
//				reporterr( "pairanch[%d].endj: %d->%d\n", j, pairanch[j].endj, map[pairanch[j].endj] );
				pairanch[j].startj = map[pairanch[j].startj];
				pairanch[j].endj = map[pairanch[j].endj];
			}
		}
	}
	free( map );
}

static int anchidcomp( const void *p, const void *q )
{
	if ( ((ExtAnch *)q)->i != ((ExtAnch *)p)->i )
		return ((ExtAnch *)p)->i - ((ExtAnch *)q)->i;
	return ((ExtAnch *)p)->j - ((ExtAnch *)q)->j;
}

static int anchcomp( const void *p, const void *q )
{
	if ( ((ExtAnch *)q)->starti != ((ExtAnch *)p)->starti )
		return ((ExtAnch *)p)->starti - ((ExtAnch *)q)->starti;
	return (int)((void *)p - (void *)q);
}

static int anchscorecomp( const void *p, const void *q )
{
	if ( ((ExtAnch *)q)->score != ((ExtAnch *)p)->score )
		return ((ExtAnch *)q)->score - ((ExtAnch *)p)->score;
	return (int)((void *)q - (void *)p);
}


static void indexanchors( ExtAnch *a, int **idx )
{
	int n;
	for( n=0; a[n].i>-1; n++ )
		;

	qsort( a, n, sizeof( ExtAnch ), anchidcomp );


	for( n=0; a[n].i>-1; n++ )
	{
//		reporterr( "%d, %dx%d, %d-%d x %d-%d\n", n, a[n].i, a[n].j, a[n].starti, a[n].endi, a[n].startj, a[n].endj );
		if( idx[a[n].i][a[n].j] == -1 ) idx[a[n].i][a[n].j] = n;
	}
#if 0
	int m;
	for( n=0; n<njob; n++ ) for( m=n+1; m<njob; m++ )
		reporterr( "%dx%d -> %d\n", n, m, idx[n][m] );
	exit( 1 );
#endif
}


static void checkanchors_strongestfirst_considerseparation( ExtAnch *a, int s, double gapratio1, double gapratio2 )
{
	int p, q;
	qsort( a, s, sizeof( ExtAnch ), anchscorecomp );

	for( p=0; a[p].i>-1; p++ )
	{
		if( a[p].starti == -1 ) continue;

		int nearest, mindist;
		double zurei, zurej;
		if( p )
		{
			mindist = 999999999;
			for( q=0; q<p; q++ )
			{
				if( a[q].starti == -1 ) continue;
				if( abs( a[p].starti - a[q].starti ) < mindist )
				{
					nearest = q;
					mindist = abs( a[p].starti - a[q].starti );
				}
			}
			//reporterr( "nearest=%d\n", nearest );
			if( a[nearest].starti < a[p].starti )
			{
				zurei = (double)( a[p].starti - a[nearest].endi )/(1.0+gapratio1);
				zurej = (double)( a[p].startj - a[nearest].endj )/(1.0+gapratio2);
			}
			else
			{
				zurei = (double)( a[nearest].starti - a[p].endi )/(1.0+gapratio1);
				zurej = (double)( a[nearest].startj - a[p].endj )/(1.0+gapratio2);
			}
		}
		else
			zurei = zurej = 0.0;
		if( fabs( zurei - zurej ) > maxanchorseparation )
//		if( fabs( zurei - zurej ) > maxanchorseparation || zurei > maxanchorseparation || zurej > maxanchorseparation ) // test
		{
//			reporterr( "warning: long internal gaps in %d-%d, |%5.2f-%5.2f - %5.2f| = %5.2f > %5.2f\n", a[p].i, a[p].j, nogaplenestimation1, nogaplenestimation2, zureij, fabs( zureij - ( nogaplenestimation1, nogaplenestimation2 ) ), maxanchorseparation );
			a[p].starti = a[p].startj = a[p].startj = a[p].endj = -1;
			continue;
		}

//		reporterr( "P score=%d, %d-%d, %d-%d\n", a[p].score, a[p].starti, a[p].endi, a[p].startj, a[p].endj );
		for( q=p+1; a[q].i>-1; q++ )
		{
			if( a[q].starti == -1 ) continue;
//			reporterr( "Q score=%d, %d-%d, %d-%d\n", a[q].score, a[q].starti, a[q].endi, a[q].startj, a[q].endj );


			if( a[p].endi < a[q].starti && a[p].endj < a[q].startj ) 
			{
//				reporterr( "consistent\n" );
				;
			}
			else if( a[p].endi == a[q].starti && a[p].endj < a[q].startj && a[q].starti<a[q].endi ) 
			{
				a[q].starti += 1; // 1 zai overlap
			}
			else if( a[p].endi < a[q].starti && a[p].endj == a[q].startj && a[q].startj<a[q].endj ) 
			{
				a[q].startj += 1; // 1 zai overlap
			}
			else if( a[q].endi < a[p].starti && a[q].endj < a[p].startj )
			{
//				reporterr( "consistent\n" );
				;
			}
			else if( a[q].endi == a[p].starti && a[q].endj < a[p].startj && a[q].starti<a[q].endi ) // bug in v7.442
			{
				a[q].endi -= 1; // 1 zai overlap
			}
			else if( a[q].endi < a[p].starti && a[q].endj == a[p].startj && a[q].startj<a[q].endj )
			{
				a[q].endj -= 1; // 1 zai overlap
			}
			else 
			{
//				reporterr( "INconsistent\n" );
				a[q].starti = a[q].startj = a[q].startj = a[q].endj = -1;
			}
		}
		if( p % 1000 == 0 ) reporterr( "%d/%d\r", p, s );
	}

	qsort( a, s, sizeof( ExtAnch ), anchcomp );
#if 0
	reporterr( "after filtering and sorting\n" );
	for( p=0; a[p].i>-1; p++ )
	{
		reporterr( "a[%d].starti,j=%d,%d, score=%d\n", p, a[p].starti, a[p].startj, a[p].score );
	}
#endif
}

static void checkanchors_strongestfirst( ExtAnch *a, int s )
{
	int p, q;
	qsort( a, s, sizeof( ExtAnch ), anchscorecomp );

	for( p=0; a[p].i>-1; p++ )
	{
		if( a[p].starti == -1 ) continue;

//		reporterr( "P score=%d, %d-%d, %d-%d\n", a[p].score, a[p].starti, a[p].endi, a[p].startj, a[p].endj );
		for( q=p+1; a[q].i>-1; q++ )
		{
			if( a[q].starti == -1 ) continue;
//			reporterr( "Q score=%d, %d-%d, %d-%d\n", a[q].score, a[q].starti, a[q].endi, a[q].startj, a[q].endj );


			if( a[p].endi < a[q].starti && a[p].endj < a[q].startj ) 
			{
//				reporterr( "consistent\n" );
				;
			}
			else if( a[p].endi == a[q].starti && a[p].endj < a[q].startj && a[q].starti<a[q].endi ) 
			{
				a[q].starti += 1; // 1 zai overlap
			}
			else if( a[p].endi < a[q].starti && a[p].endj == a[q].startj && a[q].startj<a[q].endj ) 
			{
				a[q].startj += 1; // 1 zai overlap
			}
			else if( a[q].endi < a[p].starti && a[q].endj < a[p].startj )
			{
//				reporterr( "consistent\n" );
				;
			}
			else if( a[q].endi == a[p].starti && a[q].endj < a[p].startj && a[q].starti<a[q].endi ) // bug in v7.442
			{
				a[q].endi -= 1; // 1 zai overlap
			}
			else if( a[q].endi < a[p].starti && a[q].endj == a[p].startj && a[q].startj<a[q].endj )
			{
				a[q].endj -= 1; // 1 zai overlap
			}
			else 
			{
//				reporterr( "INconsistent\n" );
				a[q].starti = a[q].startj = a[q].startj = a[q].endj = -1;
			}
		}
		if( p % 1000 == 0 ) reporterr( "%d/%d\r", p, s );
	}

	qsort( a, s, sizeof( ExtAnch ), anchcomp );
#if 0
	reporterr( "after filtering and sorting\n" );
	for( p=0; a[p].i>-1; p++ )
	{
		reporterr( "a[%d].starti,j=%d,%d, score=%d\n", p, a[p].starti, a[p].startj, a[p].score );
	}
#endif
}


static double gapnongapratio( int n, char **s )
{
	int i, j, len;
	char *seq, *pt1, *pt2;
	double fv, ng;

	len = strlen( s[0] );
	seq = calloc( len+1, sizeof( char ) );

	fv = 0.0;
	ng = 0.0;
	for( i=0; i<n; i++  )
	{
		pt1 = s[i];
		while( *pt1 == '-' ) pt1++;
		pt2 = seq;
		while( *pt1 != 0 ) *pt2++ = *pt1++;
		*pt2 = *pt1; // 0
		pt1 = pt2-1;
		while( *pt1 == '-' ) pt1--;
		*(pt1+1) = 0;	
//		reporterr( "seq[i]=%s\n", s[i] );
//		reporterr( "seq=%s\n", seq );
		len = pt1-seq+1;
		for( j=0; j<len; j++ )
			if( seq[j] == '-' ) 
				fv+=1.0;
			else
				ng+=1.0;
	}
	free( seq );
	return( fv/ng );
}

static void	pickpairanch( ExtAnch **pairanch, ExtAnch *extanch, int **anchindex, int n1, int n2, int *m1, int *m2, char **seq1, char **seq2 ) // loop no junban wo kaeta hou ga iikamo
{
	int i, j, k, s;
	s = 0;
#if 0
	reporterr( "m1,m2=\n" );
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		reporterr( "%d,%d\n", m1[i], m2[j] );
	}
#endif
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
#if 1
		if( m1[i] < m2[j] ) 
		{
//			reporterr( "%dx%d, %dx%d -> jump to %d\n", i, j, m1[i], m2[j], anchindex[m1[i]][m2[j]] );
			k = anchindex[m1[i]][m2[j]];
			while( ( k!=-1 ) && ( extanch[k].i == m1[i] && extanch[k].j == m2[j] ) )
			{
				s++;
				k++;
			}
		}
		else
		{
//			reporterr( "%dx%d, %dx%d -> jump to %d\n", j, i, m1[i], m2[j], anchindex[m2[j]][m1[i]] );
			k = anchindex[m2[j]][m1[i]];
			while( ( k!=-1 ) && ( extanch[k].i == m2[j] && extanch[k].j == m1[i] ) )
			{
				s++;
				k++;
			}
		}
#else
		k = 0;
		while( extanch[k].i > -1 ) // kanari muda
		{
			//reporterr( "m1[i],m2[j]=%d,%d ? extanch[k].i,j=%d,%d k=%d\n", m1[i], m2[j], extanch[k].i, extanch[k].j, k );
			if( ( extanch[k].i == m1[i] && extanch[k].j == m2[j] ) || ( extanch[k].i == m2[j] && extanch[k].j == m1[i] ) ) 
			{
				//reporterr( "hit, extanch[k].startj=%d\n", extanch[k].startj );
				s++;
			}
			k++;
		}
#endif
	}
	*pairanch = calloc( sizeof( ExtAnch ), s+1 );

	s = 0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
#if 1
		if( m1[i] < m2[j] ) 
		{
			k = anchindex[m1[i]][m2[j]];
			while( ( k!=-1 ) && ( extanch[k].i == m1[i] && extanch[k].j == m2[j] ) )
			{
//				if( extanch[k].starti + 1 < extanch[k].endi )
				{
					(*pairanch)[s].i = i;
					(*pairanch)[s].j = j;
					(*pairanch)[s].starti = extanch[k].starti; // map mae
					(*pairanch)[s].endi = extanch[k].endi; // map mae
					(*pairanch)[s].startj = extanch[k].startj; // map mae 
					(*pairanch)[s].endj = extanch[k].endj; // map mae
					(*pairanch)[s].score = extanch[k].score;
					s++;
				}
				k++;
			}
		}
		else
		{
			k = anchindex[m2[j]][m1[i]];
			while( ( k!=-1 ) && ( extanch[k].i == m2[j] && extanch[k].j == m1[i] ) )
			{
//				if( extanch[k].starti + 1 < extanch[k].endi )
				{
					(*pairanch)[s].i = i;
					(*pairanch)[s].j = j;
					(*pairanch)[s].starti = extanch[k].startj; // map mae
					(*pairanch)[s].endi = extanch[k].endj; // map mae
					(*pairanch)[s].startj = extanch[k].starti; // map mae 
					(*pairanch)[s].endj = extanch[k].endi; // map mae
					(*pairanch)[s].score = extanch[k].score;
					s++;
				}
				k++;
			}
		}
#else
		k = 0;
		while( extanch[k].i > -1 ) // kanari muda
		{
			if( extanch[k].i == m1[i] && extanch[k].j == m2[j] )
			{
				(*pairanch)[s].i = i;
				(*pairanch)[s].j = j;
				(*pairanch)[s].starti = extanch[k].starti; // map mae
				(*pairanch)[s].endi = extanch[k].endi; // map mae
				(*pairanch)[s].startj = extanch[k].startj; // map mae 
				(*pairanch)[s].endj = extanch[k].endj; // map mae
				(*pairanch)[s].score = extanch[k].score;
				s++;
			}
			if( extanch[k].j == m1[i] && extanch[k].i == m2[j] )
			{
				(*pairanch)[s].i = i;
				(*pairanch)[s].j = j;
				(*pairanch)[s].starti = extanch[k].startj; // map mae
				(*pairanch)[s].endi = extanch[k].endj; // map mae
				(*pairanch)[s].startj = extanch[k].starti; // map mae 
				(*pairanch)[s].endj = extanch[k].endi; // map mae
				(*pairanch)[s].score = extanch[k].score;
				s++;
			}
			k++;
		}
#endif
	}
	(*pairanch)[s].i = (*pairanch)[s].j = -1;

	recountpositions( *pairanch, n1, n2, seq1, seq2 );
//	truncateseq_group( *pairanch, seq1, seq2, n1, n2 );
//	copybackanchors( *pairanch, ddn1, n2, seq1, seq2 ); // tabun dame

#if 0
	reporterr( "Before check\n" );
	for( k=0; (*pairanch)[k].i>-1; k++ )
	{
		if( (*pairanch)[k].starti!=-1)
			reporterr( "seq1-%d,seq2-%d %d-%d,%d-%d\n", (*pairanch)[k].i, (*pairanch)[k].j, (*pairanch)[k].starti, (*pairanch)[k].endi, (*pairanch)[k].startj, (*pairanch)[k].endj );
	}
#endif

#if 0
	reporterr( "\ngroup1=\n" );
	for( i=0; m1[i]>-1; i++ )
		reporterr( "%d ", m1[i] );
	reporterr( "\n" );
	reporterr( "\ngroup2=\n" );
	for( i=0; m2[i]>-1; i++ )
		reporterr( "%d ", m2[i] );
	reporterr( "\n" );
#endif

	reporterr( "Checking external anchors\n" );
	if( maxanchorseparation != -1.0 )
		checkanchors_strongestfirst_considerseparation( *pairanch, s, gapnongapratio( n1, seq1 ), gapnongapratio( n2, seq2 ) );
	else
		checkanchors_strongestfirst( *pairanch, s );


#if 0
	reporterr( "After check\n" );
	for( k=0; (*pairanch)[k].i>-1; k++ )
	{
		if( (*pairanch)[k].starti!=-1)
			reporterr( "seq1-%d,seq2-%d %d-%d,%d-%d\n", (*pairanch)[k].i, (*pairanch)[k].j, (*pairanch)[k].starti, (*pairanch)[k].endi, (*pairanch)[k].startj, (*pairanch)[k].endj );
	}
#endif
}

static void *treebasethread( void *arg )
{
	treebasethread_arg_t *targ = (treebasethread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int *nrunpt = targ->nrunpt;
	int njob = targ->njob;
	int *nlen = targ->nlen;
	int *jobpospt = targ->jobpospt;
	int ***topol = targ->topol;
	Treedep *dep = targ->dep;
	double ***cpmxhist = targ->cpmxhist;
	int **memhist = targ->memhist;
	char **aseq = targ->aseq;
	double *effarr = targ->effarr;
	int *alloclen = targ->alloclenpt;
	int *fftlog = targ->fftlog;
	char *mergeoralign = targ->mergeoralign;
	double **newdistmtx = targ->newdistmtx;
	int *selfscore = targ->selfscore;
	ExtAnch *extanch = targ->extanch;
	int **anchindex = targ->anchindex;

	char **mseq1, **mseq2;
	char **localcopy;
	int i, m, j, l;
	int immin, immax;
	int len1, len2;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	double **dynamicmtx;
	int ssi, ssm, bunbo;
	int tm, ti;
	int **localmem = NULL;
	double ***cpmxchild0, ***cpmxchild1;
	double orieff1, orieff2;
	ExtAnch *pairanch = NULL;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif



	tscore = 0;
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	localcopy = calloc( njob, sizeof( char * ) );
	for( i=0; i<njob; i++ ) localcopy[i] = NULL;
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	if( specificityconsideration )
		dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
	localmem = (int **)calloc( sizeof( int * ), 2 );



#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		l = *jobpospt;
		if( l == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
			Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
			D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
			free( mseq1 );
			free( mseq2 );
			free( localcopy );
			free( effarr1 );
			free( effarr2 );
			free( indication1 );
			free( indication2 );
			if( specificityconsideration )
				FreeDoubleMtx( dynamicmtx );
			free( localmem );
			return( NULL );
		}
		*jobpospt = l+1;

//		reporterr( "l=%d, child0=%d, child1=%d\n", l, dep[l].child0, dep[l].child1 );

		if( dep[l].child0 != -1 )
		{
			while( dep[dep[l].child0].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
		if( dep[l].child1 != -1 ) 
		{
			while( dep[dep[l].child1].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}




//		while( *nrunpt >= nthread ) // bug
		while( *nrunpt >= nthreadtb ) // tabun iranai
			pthread_cond_wait( targ->treecond, targ->mutex ); // tabun iranai
		(*nrunpt)++;

		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
#if 0
		localmem[0][0] = -1;
		posinmem=topolorderz( localmem[0], topol, dep, l, 0 ) - localmem[0];
		localmem[1][0] = -1;
		posinmem=topolorderz( localmem[1], topol, dep, l, 1 ) - localmem[1];
#else
		if( dep[l].child0 == -1 ) 
		{
			localmem[0] = calloc( sizeof( int ), 2 );
			localmem[0][0] = m1;
			localmem[0][1] = -1;
			clus1 = 1;
		}
		else
		{
			localmem[0] = memhist[dep[l].child0];
			clus1 = intlen( localmem[0] );
		}
		if( dep[l].child1 == -1 ) 
		{
			localmem[1] = calloc( sizeof( int ), 2 );
			localmem[1][0] = m2;
			localmem[1][1] = -1;
			clus2 = 1;
		}
		else
		{
			localmem[1] = memhist[dep[l].child1];
			clus2 = intlen( localmem[1] );
		}

		if( l != njob-2 )
		{
			memhist[l] = calloc( sizeof( int ), clus1+clus2+1 );
			intcpy( memhist[l], localmem[0] );
			intcpy( memhist[l]+clus1, localmem[1] );
			memhist[l][clus1+clus2] = -1;
		}
#endif

// moved, 2018/Mar/10. Must be after changing memhist[l]
		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
			dep[l].done = 1;
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
//			free( topol[l][0] ); topol[l][0] = NULL;
//			free( topol[l][1] ); topol[l][1] = NULL;
//			free( topol[l] ); topol[l] = NULL;
			pthread_mutex_unlock( targ->mutex );
			free( localmem[0] );
			free( localmem[1] );
			continue;
		}


//		reporterr( "l=%d, dep[l].child0=%d, dep[l].child1=%d\n", l, dep[l].child0, dep[l].child1 );
		if( dep[l].child0 == -1 ) cpmxchild0 = NULL; else cpmxchild0 = cpmxhist+dep[l].child0;
		if( dep[l].child1 == -1 ) cpmxchild1 = NULL; else cpmxchild1 = cpmxhist+dep[l].child1;
//		reporterr( "cpmxchild0=%p, cpmxchild1=%p\n", cpmxchild0, cpmxchild1 );


//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		if( specificityconsideration )
			makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
		else
			dynamicmtx = n_dis_consweight_multi;
//		reporterr( "dynamicmtx[0][1] = %f\n", dynamicmtx[0][1] );

		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen <= len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}


		for( i=0; (j=localmem[0][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
//			localcopy[j] = aseq[j];	
		}
		for( i=0; (j=localmem[1][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
//			localcopy[j] = aseq[j];	
		}

		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
		}

		if( alg == 'M' ) // hoka no thread ga M ni shitakamo shirenainode
		{
//			reporterr(       "Freeing commonIP (thread %d)\n", thread_no );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

		pthread_mutex_unlock( targ->mutex );

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( localmem[0], localcopy, mseq1, effarr1, effarr, indication1, 0.0, &orieff1 );
		clus2 = fastconjuction_noname( localmem[1], localcopy, mseq2, effarr2, effarr, indication2, 0.0, &orieff2 );
#else
		clus1 = fastconjuction_noweight( topol[l][0], localcopy, mseq1, effarr1,  indication1 );
		clus2 = fastconjuction_noweight( topol[l][1], localcopy, mseq2, effarr2,  indication2 );
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		if( l < 500 || l % 100 == 0 ) reporterr(       "\rSTEP % 5d / %d (thread %4d)", l+1, njob-1, thread_no );
#if 0
		reporterr( "\nclus1=%d, clus2=%d\n", clus1, clus2 );
#endif

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/



//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		else						   ffttry = 0;
		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
//		reporterr( "fftlog=%d,%d, ffttry=%d\n", fftlog[m1], fftlog[m2], ffttry );

		if( useexternalanchors )
		{
//			reporterr( "%%%% %d vs %d\n", m1, m2 );
			pickpairanch( &pairanch, extanch, anchindex, clus1, clus2, localmem[0], localmem[1], mseq1, mseq2 );
//			reporterr( "pairanch: %d:%d\n", pairanch[0].starti, pairanch[0].startj );
			pscore = Falign_givenanchors( pairanch, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			free( pairanch );
			pairanch = NULL;
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " f\b\b" );
			if( alg == 'M' )
			{
				if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
			}
		}
		else
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " d\b\b" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
					if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
						pscore = A__align( dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

		if( disp ) display( localcopy, njob );

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
				ti = localmem[0][i];
				ssi = selfscore[localmem[0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[localmem[1][m]];
					tm = localmem[1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}






		pthread_mutex_lock( targ->mutex );
		dep[l].done = 1;
		(*nrunpt)--;
		pthread_cond_broadcast( targ->treecond );

		for( i=0; (j=localmem[0][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		for( i=0; (j=localmem[1][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );

//		reporterr( "at step %d\n", l );
//		use_getrusage();

		pthread_mutex_unlock( targ->mutex );



		for( i=0; (j=localmem[0][i])!=-1; i++ )
		{
			if(localcopy[j] ) free( localcopy[j] );
			localcopy[j] = NULL;
		}
		for( i=0; (j=localmem[1][i])!=-1; i++ )
		{
			if( localcopy[j] ) free( localcopy[j] );
			localcopy[j] = NULL;
		}


//		if( topol[l][0] ) free( topol[l][0] );
//		topol[l][0] = NULL;
//		if( topol[l][1] ) free( topol[l][1] );
//		topol[l][1] = NULL;
//		if( topol[l] ) free( topol[l] );
//		topol[l] = NULL;


//		reporterr(       "\n" );

		free( localmem[0] );
		free( localmem[1] );
	}
#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
}
#endif

static int dooneiteration( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, int **memhist, double ***cpmxhist, double *effarr, double **newdistmtx, int *selfscore, ExtAnch *extanch, int **anchindex, int *alloclen, int (*callback)(int, int, char*) )
{
	int l, ll, len1, len2, i, j;
	int clus1, clus2;
	double pscore;
	char *indication1 = NULL, *indication2 = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	int *fftlog = NULL; // fixed at 2006/07/26
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	int *alreadyaligned = NULL;
	double **dynamicmtx = NULL;
	int **localmem = NULL;
	double ***cpmxchild0, ***cpmxchild1;
	double orieff1, orieff2;
	double oscore, nscore;
	ExtAnch *pairanch;
	char **oseq1, **oseq2;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif


	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
		alreadyaligned = AllocateIntVec( njob );
		if( specificityconsideration )
			dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
		localmem = calloc( sizeof( int * ), 2 );
	}
	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	if( callback && callback( 0, 50, "Progressive alignment" ) ) goto chudan_tbfast;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0 // chain you
	localmem[0][0] = -1;
	localmem[1][0] = -1;
	clus1 = 1;// chain ni hitsuyou
#endif

#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	for( ll=0; ll<njob*ITERATIVECYCLE; ll++ )
	{
		l = ll % njob;
		cpmxchild0 = NULL;
		cpmxchild1 = NULL;

		localmem[0] = calloc( sizeof( int ), 2 );
		localmem[0][0] = l;
		localmem[0][1] = -1;
		clus1 = 1;
		m1 = localmem[0][0];

		localmem[1] = calloc( sizeof( int ), njob );
		for( i=0,j=0; i<njob; i++ )
			if( i != l ) localmem[1][j++] = i;
		localmem[1][j] = -1;
		clus2 = njob-1;
		m2 = localmem[1][0];

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		if( specificityconsideration )
			makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
		else
			dynamicmtx = n_dis_consweight_multi;
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );

		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
		}

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( localmem[0], aseq, mseq1, effarr1, effarr, indication1, 0.0, &orieff1 );
		clus2 = fastconjuction_noname( localmem[1], aseq, mseq2, effarr2, effarr, indication2, 0.0, &orieff2 );
#else
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
//		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
//		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif

		intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, strlen( mseq1[0] ), &oscore );

		oseq1 = AllocateCharMtx( 1, len1+1 );
		oseq2 = AllocateCharMtx( njob-1, len2+1 );
		for( i=0; i<clus1; i++ ) strcpy( oseq1[i], mseq1[i] );
		for( i=0; i<clus2; i++ ) strcpy( oseq2[i], mseq2[i] );

		newgapstr = "-";
		commongappick( clus2, mseq2 );
		commongappick( clus1, mseq1 );

		if( l < 500 || l % 100 == 0 ) reporterr(       "\rIteration % 5d / %d ", ll+1, njob*ITERATIVECYCLE );
		if( callback && callback( 0, 50+50*l/(njob-1), "Progressive alignment" ) ) goto chudan_tbfast;
#if 0
		reporterr( "\nclus1=%d, clus2=%d\n", clus1, clus2 );
		for( i=0; i<clus1; i++ ) reporterr( "effarr1[%d]=%f\n", i, effarr1[i] );
		for( i=0; i<clus2; i++ ) reporterr( "effarr2[%d]=%f\n", i, effarr2[i] );
#endif

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/

		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );

		if( useexternalanchors )
		{
			pickpairanch( &pairanch, extanch, anchindex, clus1, clus2, localmem[0], localmem[1], mseq1, mseq2 );
//			reporterr( "pairanch: %d:%d\n", pairanch[0].starti, pairanch[0].startj );
			pscore = Falign_givenanchors( pairanch, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			free( pairanch );
			pairanch = NULL;
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " f\b\b" );
			if( alg == 'M' )
			{
				if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
//				reporterr(       "######### mseq1[0] = %s\n", mseq1[0] );
			}
		}
		else
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " d\b\b" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
					if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
						if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//						reporterr(       "\n\n %d - %d (%d x %d) : \n", topol[l][0][0], topol[l][1][0], clus1, clus2 );
						pscore = A__align( dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}

					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
		intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, strlen( mseq1[0] ), &nscore );
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		if( nscore < oscore )
		{
			for( i=0; i<clus1; i++ ) strcpy( mseq1[i], oseq1[i] );
			for( i=0; i<clus2; i++ ) strcpy( mseq2[i], oseq2[i] );
		}
		FreeCharMtx( oseq1 );
		FreeCharMtx( oseq2 );

		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
//		reporterr(       "\n" );



#if 0
		if( localmem[1][0] == 13 ) 
		{
			reporterr( "OUTPUT!\n" );
			for( i=0; i<clus1; i++ ) reporterr( ">g1\n%s\n", mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr( ">g2\n%s\n", mseq2[i] );
			exit( 1 );
		}
#endif

//		free( topol[l][0] ); topol[l][0] = NULL;
//		free( topol[l][1] ); topol[l][1] = NULL;
//		free( topol[l] ); topol[l] = NULL;


//		reporterr(       ">514\n%s\n", aseq[514] );
		free( localmem[0] );
		free( localmem[1] );

	}

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( fftlog );
	if( specificityconsideration )
		FreeDoubleMtx( dynamicmtx );
	free( alreadyaligned );
	free( localmem );
	effarr1 = NULL;
	return( 0 );

	chudan_tbfast:

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	if( effarr1 ) free( effarr1 ); effarr1 = NULL;
	if( effarr2 ) free( effarr2 ); effarr2 = NULL;
	if( indication1 ) free( indication1 ); indication1 = NULL;
	if( indication2 ) free( indication2 ); indication2 = NULL;
	if( fftlog ) free( fftlog ); fftlog = NULL;
	if( alreadyaligned ) free( alreadyaligned ); alreadyaligned = NULL;
	if( specificityconsideration )
	{
		if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); dynamicmtx = NULL;
	}
	if( localmem ) free( localmem ); localmem = NULL;
#if SKIP
	if( skiptable1 ) FreeIntMtx( skiptable1 ); skiptable1 = NULL;
	if( skiptable2 ) FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif

	return( 1 );
}
static int treebase( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, int **memhist, double ***cpmxhist, double *effarr, double **newdistmtx, int *selfscore, ExtAnch *extanch, int **anchindex, int *alloclen, int (*callback)(int, int, char*) )
{
	int l, len1, len2, i, m, immin, immax;
	int len1nocommongap, len2nocommongap;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1 = NULL, *indication2 = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	int *fftlog = NULL; // fixed at 2006/07/26
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	int m1, m2;
	int *gaplen = NULL;
	int *gapmap = NULL;
	int *alreadyaligned = NULL;
	double **dynamicmtx = NULL;
	double ssi, ssm, bunbo;
	int tm, ti;
	int gapmaplen;
	int **localmem = NULL;
	double ***cpmxchild0, ***cpmxchild1;
	double orieff1, orieff2;
	ExtAnch *pairanch;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif


	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
		gaplen = AllocateIntVec( *alloclen+10 );
		gapmap = AllocateIntVec( *alloclen+10 );
		alreadyaligned = AllocateIntVec( njob );
		if( specificityconsideration )
			dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
		localmem = calloc( sizeof( int * ), 2 );
	}
	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	if( callback && callback( 0, 50, "Progressive alignment" ) ) goto chudan_tbfast;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0 // chain you
	localmem[0][0] = -1;
	localmem[1][0] = -1;
	clus1 = 1;// chain ni hitsuyou
#endif

#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{
		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
//		reporterr( " at the beginning of the loop, clus1,clus2=%d,%d\n", clus1, clus2 );

//		reporterr( "l=%d, dep[l].child0=%d, dep[l].child1=%d\n", l, dep[l].child0, dep[l].child1 );
		if( dep[l].child0 == -1 ) cpmxchild0 = NULL; else cpmxchild0 = cpmxhist+dep[l].child0;
		if( dep[l].child1 == -1 ) cpmxchild1 = NULL; else cpmxchild1 = cpmxhist+dep[l].child1;
//		reporterr( "cpmxchild0=%p, cpmxchild1=%p\n", cpmxchild0, cpmxchild1 );

#if 0
		if(  l > 0 && dep[l].child0 == l-1 && dep[l].child1 == -1 && dep[dep[l].child0].child1 == -1 )
		{
			localmem[0][clus1] = topol[l-1][1][0];
			localmem[0][clus1+1] = -1;

			localmem[1][0] = topol[l][1][0];
			localmem[1][1] = -1;
		}
		else
		{
			localmem[0][0] = -1;
			posinmem = topolorderz( localmem[0], topol, dep, l, 0 ) - localmem[0];
			localmem[1][0] = -1;
			posinmem = topolorderz( localmem[1], topol, dep, l, 1 ) - localmem[1];
		}
#else
		if( dep[l].child0 == -1 ) 
		{
			localmem[0] = calloc( sizeof( int ), 2 );
			localmem[0][0] = m1;
			localmem[0][1] = -1;
			clus1 = 1;
		}
		else
		{
			localmem[0] = memhist[dep[l].child0];
			clus1 = intlen( localmem[0] );
		}
		if( dep[l].child1 == -1 ) 
		{
			localmem[1] = calloc( sizeof( int ), 2 );
			localmem[1][0] = m2;
			localmem[1][1] = -1;
			clus2 = 1;
		}
		else
		{
			localmem[1] = memhist[dep[l].child1];
			clus2 = intlen( localmem[1] );
		}

		if( l != njob-2 )
		{
			memhist[l] = calloc( sizeof( int ), clus1+clus2+1 );
			intcpy( memhist[l], localmem[0] );
			intcpy( memhist[l]+clus1, localmem[1] );
			memhist[l][clus1+clus2] = -1;
		}
#endif

		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
//			free( topol[l][0] ); topol[l][0] = NULL;
//			free( topol[l][1] ); topol[l][1] = NULL;
//			free( topol[l] ); topol[l] = NULL;
			free( localmem[0] );
			free( localmem[1] );
			continue;
		}

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		if( specificityconsideration )
			makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
		else
			dynamicmtx = n_dis_consweight_multi;
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );


		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				reporterr(       "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				reporterr(       "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( localmem[0], aseq, mseq1, effarr1, effarr, indication1, 0.0, &orieff1 );
		clus2 = fastconjuction_noname( localmem[1], aseq, mseq2, effarr2, effarr, indication2, 0.0, &orieff2 );
#else
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
//		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
//		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif











		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' )
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";

		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap );
			commongappick( clus2, mseq2 );
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap );
			commongappick( clus1, mseq1 );
			len1nocommongap = strlen( mseq1[0] );
		}

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		if( l < 500 || l % 100 == 0 ) reporterr(       "\rSTEP % 5d / %d ", l+1, njob-1 );
		if( callback && callback( 0, 50+50*l/(njob-1), "Progressive alignment" ) ) goto chudan_tbfast;
#if 0
		reporterr( "\nclus1=%d, clus2=%d\n", clus1, clus2 );
#endif

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/


		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		else						   ffttry = 0;
		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
//		reporterr( "fftlog=%d,%d, ffttry=%d\n", fftlog[m1], fftlog[m2], ffttry );

		if( useexternalanchors )
		{
			pickpairanch( &pairanch, extanch, anchindex, clus1, clus2, localmem[0], localmem[1], mseq1, mseq2 );
//			reporterr( "pairanch: %d:%d\n", pairanch[0].starti, pairanch[0].startj );
			pscore = Falign_givenanchors( pairanch, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			free( pairanch );
			pairanch = NULL;
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " f\b\b" );
			if( alg == 'M' )
			{
				if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
//				reporterr(       "######### mseq1[0] = %s\n", mseq1[0] );
			}
		}
		else
		{
			if( l < 500 || l % 100 == 0 ) reporterr(       " d\b\b" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					if( l < 500 || l % 100 == 0 ) reporterr(       "m" );
					if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					break;
				case( 'd' ):
					if( 1 && clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = D__align_ls( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
						if( l < 500 || l % 100 == 0 ) if( ( cpmxchild1 && *cpmxchild1 ) || ( cpmxchild0 && *cpmxchild0 ) ) reporterr(       " h" );
//						reporterr(       "\n\n %d - %d (%d x %d) : \n", topol[l][0][0], topol[l][1][0], clus1, clus2 );
						pscore = A__align( dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist+l, orieff1, orieff2 );
					}

					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
//		reporterr(       "\n" );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
			reporterr( "Check source!!!\n" );
			exit( 1 );
		}
		if( mergeoralign[l] == '2' )
		{
//			if( localkeeplength ) ndeleted += deletenewinsertions( clus1, clus2, mseq1, mseq2, NULL );
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP0 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP0 mseq2[%d] = \n%s\n", i, mseq2[i] );
			gapmaplen = strlen( mseq1[0] )-len1nocommongap+len1;
			adjustgapmap( gapmaplen, gapmap, mseq1[0] );
#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP1 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP1 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif
//			if( clus1 + clus2 < njob ) restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			if( smoothing )
			{
				restorecommongapssmoothly( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps_bothorders( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, gapmaplen, *alloclen, alg, '-' );
			}
			else
			{
				restorecommongaps( njob, njob-(clus1+clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps( njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, *alloclen, alg, '-' );
			}

#if 0
			reporterr( "\n" );
			for( i=0; i<clus1; i++ ) reporterr(       ">STEP3 mseq1[%d] = \n%s\n", i, mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr(       ">STEP3 mseq2[%d] = \n%s\n", i, mseq2[i] );
#endif

#if 0
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; i<clus1; i++ ) 
			{
				reporterr( "mseq1[%d] bef change = %s\n", i, mseq1[i] );
				eq2dash( mseq1[i] );
				reporterr( "mseq1[%d] aft change = %s\n", i, mseq1[i] );
			}
			for( i=0; i<clus2; i++ ) 
			{
				reporterr( "mseq2[%d] bef change = %s\n", i, mseq2[i] );
				eq2dash( mseq2[i] );
				reporterr( "mseq2[%d] aft change = %s\n", i, mseq2[i] );
			}
			for( i=0; i<clus1; i++ ) eq2dash( mseq1[i] );
			for( i=0; i<clus2; i++ ) eq2dash( mseq2[i] );
#endif


			eq2dashmatometehayaku( mseq1, clus1 );
			eq2dashmatometehayaku( mseq2, clus2 );

			for( i=0; (m=localmem[1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
#if SKIP
//				makeskiptable( 1, skiptable1, mseq1+i ); // allocate suru.
#endif
				ti = localmem[0][i];
				ssi = selfscore[localmem[0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[localmem[1][m]];
					tm = localmem[1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0.0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}

#if 0
		if( localmem[1][0] == 13 ) 
		{
			reporterr( "OUTPUT!\n" );
			for( i=0; i<clus1; i++ ) reporterr( ">g1\n%s\n", mseq1[i] );
			for( i=0; i<clus2; i++ ) reporterr( ">g2\n%s\n", mseq2[i] );
			exit( 1 );
		}
#endif

//		free( topol[l][0] ); topol[l][0] = NULL;
//		free( topol[l][1] ); topol[l][1] = NULL;
//		free( topol[l] ); topol[l] = NULL;


//		reporterr(       ">514\n%s\n", aseq[514] );
		free( localmem[0] );
		free( localmem[1] );
	}

#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( fftlog );
	free( gaplen );
	free( gapmap );
	if( specificityconsideration )
		FreeDoubleMtx( dynamicmtx );
	free( alreadyaligned );
	free( localmem );
	effarr1 = NULL;
	return( 0 );

	chudan_tbfast:

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	Falign_givenanchors( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0 );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	if( effarr1 ) free( effarr1 ); effarr1 = NULL;
	if( effarr2 ) free( effarr2 ); effarr2 = NULL;
	if( indication1 ) free( indication1 ); indication1 = NULL;
	if( indication2 ) free( indication2 ); indication2 = NULL;
	if( fftlog ) free( fftlog ); fftlog = NULL;
	if( gaplen ) free( gaplen ); gaplen = NULL;
	if( gapmap ) free( gapmap ); gapmap = NULL;
	if( alreadyaligned ) free( alreadyaligned ); alreadyaligned = NULL;
	if( specificityconsideration )
	{
		if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); dynamicmtx = NULL;
	}
	if( localmem ) free( localmem ); localmem = NULL;
#if SKIP
	if( skiptable1 ) FreeIntMtx( skiptable1 ); skiptable1 = NULL;
	if( skiptable2 ) FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif

	return( 1 );
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
    reporterr(       "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
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
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = UPGMA (average).\n" );
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

static double **preparepartmtx( int nseq )
{
	int i;
	double **val;
	double size;

	val = (double **)calloc( nseq, sizeof( double *) );;
	size = 0;

	if( compacttree == 1 )
	{
		for( i=0; i<nseq; i++ )
		{
			size += (double)sizeof( double ) * nseq;
			if( size > maxdistmtxsize ) 
			{
				reporterr( "\n\nThe size of full distance matrix is estimated to exceed %.2fGB.\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "Will try the calculation using a %d x %d matrix.\n", nseq, i );
				reporterr( "This calculation will be slow due to the limited RAM space.\n", i, nseq );
				reporterr( "To avoid the slowdown, please try '--initialramusage xGB' (x>>%.2f),\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "if larger RAM space is available.\n" );
				reporterr( "Note that xGB is NOT the upper limit of RAM usage.\n" );
				reporterr( "Two to three times larger space may be used for building a guide tree.\n" );
				reporterr( "Memory usage of the MSA stage depends on similarity of input sequences.\n\n" );
//				reporterr( "If the RAM is small, try '--initialramusage xGB' with a smaller x value.\n" );
				reporterr( "The '--memsavetree' option uses smaller RAM space.\n" );
				reporterr( "If tree-like relationship can be ignored, try '--pileup' or '--randomchain'.\n\n" );
				reporterr( "The result of --initialramusage xGB is almost identical to the default, except for rounding differences.\n" );

				reporterr( "In the cases of --memsavetree, --pileup and --randomchain, the result differs from the default.\n\n" );
				break;
			}
			val[i] = (double *)calloc( nseq, sizeof( double ) );
		}
		if( i == nseq ) reporterr( "The full matrix will be used.\n" );

		for( ;i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	else
	{
		for( i=0; i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	return( val );
}

int disttbfast( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*))
{
	int  *nlen = NULL;	
	int  *nogaplen = NULL;	
	char **name = NULL, **seq = NULL;
	char **mseq1 = NULL, **mseq2 = NULL;
	char **bseq = NULL;
	double *eff = NULL;
	int i, j;
	int ***topol = NULL;
	int *addmem = NULL;
	Treedep *dep = NULL;
	int **memhist = NULL;
	double ***cpmxhist = NULL;
	double **len = NULL;
	FILE *infp = NULL;
//	FILE *adfp;
	char c;
	int alloclen;
	double longer, shorter;
	double lenfac;
	double bunbo;

	FILE *orderfp = NULL, *hat2p = NULL;
	int *grpseq = NULL;
	char *tmpseq = NULL;
	int  **pointt = NULL;
	double **mtx = NULL; // by D. Mathog
	int *table1 = NULL;
	char b[B];
	int ien, nlim;
	int includememberres0, includememberres1;
	double unweightedspscore;
	int alignmentlength;
	char *mergeoralign = NULL;
	int foundthebranch;
	int nsubalignments = 0, maxmem;
	int **subtable = NULL;
	int *insubtable = NULL;
	int *preservegaps = NULL;
	char ***subalnpt = NULL;
	int val;
	char **tmpargv = NULL;
	int iguidetree;
	int *selfscore = NULL;
	int calcpairdists;
	int	**skiptable = NULL;
	char algbackup;
	char *originalgaps = NULL;
	char **addbk = NULL;
	GapPos **deletelist = NULL;
	FILE *dlf = NULL;
	int randomseed;
	int **localmem = NULL;
	int posinmem;
// for compacttree
	int *mindistfrom = NULL;
	double *mindist = NULL;
	double **partmtx = NULL;
// for compacttree
	ExtAnch *extanch = NULL;
	int **anchindex = NULL;


	if( ngui )
	{
		initglobalvariables();
		njob = ngui;
		nlenmax = 0;
		for( i=0; i<njob; i++ )
		{
			ien = strlen( seqgui[i] );
			if( ien > nlenmax ) nlenmax = ien;
		}
		infp = NULL;
//		stderr = fopen( "/dev/null", "a" ); // Windows????
		tmpargv = AllocateCharMtx( argc, 0 );
		for( i=0; i<argc; i++ ) tmpargv[i] = argv[i];
		gmsg = 1;
	}
	else
		gmsg = 0; // iranai

	arguments( argc, argv );
	algbackup = alg; // tbfast wo disttbfast ni ketsugou shitatame.
#ifndef enablemultithread
	nthreadpair = nthread = 0;
#endif


	if( ngui )
	{
		for( i=0; i<argc; i++ ) 
		{
//			free( tmpargv[i] );
			argv[i] = tmpargv[i];
		}
		free( tmpargv );
	}
	else
	{
		if( inputfile )
		{
			infp = fopen( inputfile, "rb" );
			if( !infp )
			{
				reporterr(       "Cannot open %s\n", inputfile );
				exit( 1 );
			}
		}
		else
			infp = stdin;
	
		getnumlen( infp );
		rewind( infp );
	}
	
	if( njob > 10000000 )
	{
		reporterr(       "The number of sequences must be < %d\n", 10000000 );
//		reporterr(       "Please try the --parttree option for such large data.\n" );
		exit( 1 );
	}

	if( njob < 2 )
	{
		seq = AllocateCharMtx( 2, nlenmax*1+1 );
    	name = AllocateCharMtx( 2, B+1 );
	    nlen = AllocateIntVec( 2 ); 
		readData_pointer( infp, name, nlen, seq );
		fclose( infp );
		gappick0( seq[1], seq[0] );
		writeData_pointer( stdout, njob, name, nlen, seq+1 );
		reporterr( "Warning: Only %d sequence found.\n", njob ); 
		FreeCharMtx( seq );
		FreeCharMtx( name );
		free( nlen );
		exit( 0 );
	}

	if( specificityconsideration != 0.0 && nlenmax)
	{
		if( nlenmax > 100000 )
		{
			reporterr( "\n" );
			reporterr( "Too long to apply --allowshift or --unalignlevel>0\n" );
			reporterr( "Please use the normal mode.\n" );
			reporterr( "Please also note that MAFFT does not assume genomic rearrangements.\n" );
			reporterr( "\n" );
			exit( 1 );
		}
	}


#if !defined(mingw) && !defined(_MSC_VER)
	setstacksize( 200 * njob ); // topolorder() de ookime no stack wo shiyou.
#endif

	if( subalignment )
	{
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem );
		reporterr(       "nsubalignments = %d\n", nsubalignments );
		reporterr(       "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		insubtable = AllocateIntVec( njob );
		preservegaps = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		for( i=0; i<njob; i++ ) preservegaps[i] = 0;
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );
	}



	seq = AllocateCharMtx( njob, nlenmax*1+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	eff = AllocateDoubleVec( njob );
	mergeoralign = AllocateCharVec( njob );

	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );

	if( nadd ) addmem = AllocateIntVec( nadd+1 );

	localmem = AllocateIntMtx( 2, njob+1 );

#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob ); 
    nogaplen = AllocateIntVec( njob ); 
	if( ngui )
	{
		if( copydatafromgui( namegui, seqgui, name, nlen, seq ) )
			exit( 1 );
	}
	else
	{
		readData_pointer( infp, name, nlen, seq );
		fclose( infp );
	}
#endif

	if( useexternalanchors ) // nogaplen ha ato de uwagaki sareru kamo
	{
		char *tmpseq = calloc( nlenmax+5, sizeof( char ) );
		for( i=0; i<njob; i++ )
		{
//			reporterr( "i=%d, nlenmax=%d, len=%d\n", i, nlenmax, strlen( seq[i] ) );
			gappick0( tmpseq, seq[i] );
			nogaplen[i] = strlen( tmpseq );
		}
		extanch = calloc( sizeof( ExtAnch ), 1 );
		extanch[0].i=-1;
		extanch[0].j=-1;
		reporterr( "reading anchors\n" );
		readexternalanchors( &extanch, njob, nogaplen ); // allocate sareru
		anchindex = AllocateIntMtx( njob, njob ); // sukoshi muda
		for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) anchindex[i][j] = -1;
		reporterr( "sorting anchors\n" );
		indexanchors( extanch, anchindex );
		//checkanchors_internal( extanch ); // comment out -> equivalent to v7.448
		free( tmpseq );
	}

	constants( njob, seq );


#if 0
	reporterr(       "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		reporterr(       "Illegal character %c\n", c );
		exit( 1 );
	}

	reporterr(       "\n" );

//	reporterr(       "tuplesize = %d, dorp = %c\n", tuplesize, dorp );
	if( dorp == 'p' && tuplesize != 6 )
	{
		reporterr(       "tuplesize must be 6 for aa sequence\n" );
		exit( 1 );
	}
	if( dorp == 'd' && tuplesize != 6 && tuplesize != 10 )
	{
		reporterr(       "tuplesize must be 6 or 10 for dna sequence\n" );
		exit( 1 );
	}

	if( treein )
	{
		int npickx;
		treein = check_guidetreefile( &randomseed, &npickx, &maxdistmtxsize );
		if( treein == 't' )
		{
			varpairscore( njob, npickx, nlenmax, seq, randomseed );
			exit( 1 );
		}
		else if( treein == 'c' )
		{
			compacttree = 1;
			treein = 0;
//			use_fft = 0; // kankeinai?
//			maxdistmtxsize = 5 * 1000 * 1000; // 5GB. ato de kahen ni suru.
//			maxdistmtxsize =  1.0 * 1000 * 1000 * 1000; // 5GB. ato de kahen ni suru.
		}
		else if( treein == 'Y' )
		{
			compacttree = 4; // youngest linkage, 3 ha tbfast de tsukaunode ichiou sakeru
			treein = 0;
//			use_fft = 0; // kankeinai?
		}
		else if( treein == 'S' || treein == 'C' )
		{
			compacttree = 2; // 3 ha tbfast de tsukaunode ichiou sakeru
			treein = 0;
//			use_fft = 0; // kankeinai?
		}
		else if( treein == 'a' )
		{
//			reporterr( "Compute pairwise scores\n" );
			if( njob > 200000 )
			{
				reporterr( "Chain?\n" );
				treein = 's';
				nguidetree = 1;
			}
			else if( njob < 100 || 't' == varpairscore( njob, npickx, nlenmax, seq, randomseed ) )
			{
				if( treein == 'c' ) exit( 1 );
				reporterr( "Tree!\n" );
				treein = 0;
				nguidetree = 2;
			}
			else
			{
				reporterr( "Chain!\n" );
				treein = 's';
				nguidetree = 1;
			}
		}
		else if ( treein != 0 ) // auto no toki arieru
			nguidetree = 1;
	}

# if 0 // tameshini
	if( sueff_global < 0.0001 || compacttree == 2 )
	{
		nthread = 0;
		nthreadtb = 0;
	}
#endif
//	if( njob > 10000 ) nthreadtb = 0; 
	if( njob > 20000 ) nthreadtb = 0; 
// 2018/Jan.  Hairetsu ga ooi toki
// 1. topolorder_lessargs no stack ga tarinakunaru
// 2. localcopy no tame kouritsu warui

	if( compacttree == 1 )
	{
		if( maxdistmtxsize > (double)njob * (njob-1) * sizeof( double ) / 2 ) 
		{
			reporterr( "Use conventional tree.\n" );
			compacttree = 0;
		}
	}

	if( !treein )
	{
		reporterr(       "\n\nMaking a distance matrix ..\n" );
		if( callback && callback( 0, 0, "Distance matrix" ) ) goto chudan;

		tmpseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		pointt = AllocateIntMtx( njob, nlenmax+1 );
		if( !compacttree ) mtx = AllocateFloatHalfMtx( njob ); 
		if( dorp == 'd' ) tsize = (int)pow( 4, tuplesize );
		else              tsize = (int)pow( 6, 6 );

		if( dorp == 'd' && tuplesize == 6 )
		{
			lenfaca = D6LENFACA;
			lenfacb = D6LENFACB;
			lenfacc = D6LENFACC;
			lenfacd = D6LENFACD;
		}
		else if( dorp == 'd' && tuplesize == 10 )
		{
			lenfaca = D10LENFACA;
			lenfacb = D10LENFACB;
			lenfacc = D10LENFACC;
			lenfacd = D10LENFACD;
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
			nogaplen[i] = strlen( tmpseq );
			if( nogaplen[i] < 6 )
			{
//				reporterr(       "Seq %d, too short, %d characters\n", i+1, nogaplen[i] );
//				reporterr(       "Please use mafft-ginsi, mafft-linsi or mafft-ginsi\n\n\n" );
//				exit( 1 );
			}
			if( nogaplen[i] > maxl ) maxl = nogaplen[i];
			if( dorp == 'd' ) /* nuc */
			{
				seq_grp_nuc( grpseq, tmpseq );
//				makepointtable_nuc( pointt[i], grpseq );
//				makepointtable_nuc_octet( pointt[i], grpseq );
				if( tuplesize == 10 )
					makepointtable_nuc_dectet( pointt[i], grpseq );
				else if( tuplesize == 6 )
					makepointtable_nuc( pointt[i], grpseq );
				else
				{
					reporterr(       "tuplesize=%d: not supported\n", tuplesize );
					exit( 1 );
				}
			}
			else                 /* amino */
			{
				seq_grp( grpseq, tmpseq );
				makepointtable( pointt[i], grpseq );
			}
		}
		if( nunknown ) reporterr(       "\nThere are %d ambiguous characters.\n", nunknown );


		if( compacttree )
		{

			reporterr( "Compact tree, step 1\n" );
			mindistfrom = (int *)calloc( njob, sizeof( int ) );
			mindist = (double *)calloc( njob, sizeof( double ) );
			selfscore = (int *)calloc( njob, sizeof( int ) );
			partmtx = preparepartmtx( njob );


			for( i=0; i<njob; i++ )
			{
				table1 = (int *)calloc( tsize, sizeof( int ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				makecompositiontable_p( table1, pointt[i] );
				selfscore[i] = commonsextet_p( table1, pointt[i] );
				free( table1 );
				table1 = NULL;
			}
			commonsextet_p( NULL, NULL );

#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				compactdistmtxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				double **mindistthread;
				int **mindistfromthread;
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;
				targ = calloc( nthreadpair, sizeof( compactdistmtxthread_arg_t ) );
				handle = calloc( nthreadpair, sizeof( pthread_t ) );
				mindistthread = AllocateDoubleMtx( nthreadpair, njob );
				mindistfromthread = AllocateIntMtx( nthreadpair, njob );
				pthread_mutex_init( &mutex, NULL );

		
				for( j=0; j<nthreadpair; j++ )
				{
					for( i=0; i<njob; i++ )
					{
						mindistthread[j][i] = 999.9;
						mindistfromthread[j][i] = -1;
					}
					targ[j].thread_no = j;
					targ[j].nogaplen = nogaplen;
					targ[j].pointt = pointt;
					targ[j].selfscore = selfscore;
					targ[j].partmtx = partmtx;
					targ[j].njob = njob;
					targ[j].mindist = mindistthread[j];
					targ[j].mindistfrom = mindistfromthread[j];
					targ[j].jobpospt = &jobpos;
					targ[j].mutex = &mutex;
	
					if( compacttree == 4 )
						pthread_create( handle+j, NULL, ylcompactdisthalfmtxthread, (void *)(targ+j) );
					else
						pthread_create( handle+j, NULL, compactdisthalfmtxthread, (void *)(targ+j) );
				}
		
				for( j=0; j<nthreadpair; j++ ) pthread_join( handle[j], NULL );

				for( i=0; i<njob; i++ )
				{
					mindist[i] = 999.9;
					mindistfrom[i] = -1;
					for( j=0; j<nthreadpair; j++ )
					{
						if( mindistthread[j][i] < mindist[i] )
						{
							mindist[i] = mindistthread[j][i];
							mindistfrom[i] = mindistfromthread[j][i];
						}
					}
				}
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug


				pthread_mutex_destroy( &mutex );
				FreeDoubleMtx( mindistthread );
				FreeIntMtx( mindistfromthread );
				free( handle );
				free( targ );
	
			}
			else
#endif
			{
				compactdistmtxthread_arg_t *targ;
				int jobpos;
		
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;
				targ = calloc( 1, sizeof( compactdistmtxthread_arg_t ) );
		
				{
					for( i=0; i<njob; i++ )
					{
						mindist[i] = 999.9;
						mindistfrom[i] = -1;
					}
					targ[0].thread_no = 0;
					targ[0].nogaplen = nogaplen;
					targ[0].pointt = pointt;
					targ[0].selfscore = selfscore;
					targ[0].partmtx = partmtx;
					targ[0].njob = njob;
					targ[0].mindist = mindist;
					targ[0].mindistfrom = mindistfrom;
					targ[0].jobpospt = &jobpos;
	
					if( compacttree == 4 )
						ylcompactdisthalfmtxthread( targ );
					else
						compactdisthalfmtxthread( targ );
				}

				free( targ );

				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug
			}

//			reporterr( "\n" );
//			for( i=0; i<njob; i++ ) reporterr( "mindist[%d] = %f, mindistfrom[%d] = %d\n", i, mindist[i], i, mindistfrom[i] );
			reporterr( "\ndone.\n" );

#if 0
			reporterr( "\npartmtx = .\n" );
			for( i=0; i<njob; i++ )
			{
				reporterr( "i=%d\n", i );
				if( partmtx[i] ) for( j=0; j<njob; j++ ) reporterr( "%f ", partmtx[i][j]);
				else reporterr( "nil" );
				reporterr( "\n", i );
			}
#endif
		}
		else
		{
#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				distancematrixthread_arg_t *targ; 
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
	
				jobpos = 0; 
				targ = calloc( nthreadpair, sizeof( distancematrixthread_arg_t ) ); 
				handle = calloc( nthreadpair, sizeof( pthread_t ) ); 
				pthread_mutex_init( &mutex, NULL );
	
				for( i=0; i<nthreadpair; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].jobpospt = &jobpos;
					targ[i].pointt = pointt;
					targ[i].mtx = mtx;
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, distancematrixthread, (void *)(targ+i) );
				}
			
				for( i=0; i<nthreadpair; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
				for( i=0; i<njob; i++ )
				{
					table1 = (int *)calloc( tsize, sizeof( int ) );
					if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
					if( i % 100 == 0 )
					{
						reporterr(       "\r% 5d / %d", i+1, njob );
						if( callback && callback( 0, i*25/njob, "Distance matrix" ) ) goto chudan;
					}
					makecompositiontable_p( table1, pointt[i] );
			
					for( j=i; j<njob; j++ ) 
					{
						mtx[i][j-i] = (double)commonsextet_p( table1, pointt[j] );
					} 
					free( table1 ); table1 = NULL;
				}
			}
			reporterr(       "\ndone.\n\n" );
			ien = njob-1;
	
			for( i=0; i<ien; i++ )
			{
				for( j=i+1; j<njob; j++ ) 
				{
					if( nogaplen[i] > nogaplen[j] )
					{
						longer=(double)nogaplen[i];
						shorter=(double)nogaplen[j];
					}
					else
					{
						longer=(double)nogaplen[j];
						shorter=(double)nogaplen[i];
					}
//					if( tuplesize == 6 )
					lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//					else
//						lenfac = 1.0;
//					reporterr(       "lenfac = %f (%.0f,%.0f)\n", lenfac, longer, shorter );
					bunbo = MIN( mtx[i][0], mtx[j][0] );
					if( bunbo == 0.0 )
						mtx[i][j-i] = 2.0; // 2013/Oct/17 -> 2bai
					else
						mtx[i][j-i] = ( 1.0 - mtx[i][j-i] / bunbo ) * lenfac * 2.0; // 2013/Oct/17 -> 2bai
//					reporterr(       "##### mtx = %f, mtx[i][0]=%f, mtx[j][0]=%f, bunbo=%f\n", mtx[i][j-i], mtx[i][0], mtx[j][0], bunbo );
				}
			}
			if( disopt )
			{
				for( i=0; i<njob; i++ ) 
				{
					sprintf( b, "=lgth = %04d", nogaplen[i] );
					strins( b, name[i] );
				}
			}
			FreeIntMtx( pointt ); pointt = NULL;
			commonsextet_p( NULL, NULL );
		}
		free( grpseq ); grpseq = NULL;
		free( tmpseq ); tmpseq = NULL;

#if 0 // writehat2 wo kakinaosu -> iguidetree loop nai ni idou
		if( distout )
		{
			hat2p = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
			fclose( hat2p );
		}
#endif

	}
#if 0 
	else 
	{
		reporterr(       "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2_double( prep, njob, name, mtx ); // name chuui
		fclose( prep );
		reporterr(       "done.\n" );
	}
#endif

//	reporterr( "after computing distance matrix," );
//	use_getrusage();

	if( nadd && keeplength )
	{
		originalgaps = (char *)calloc( nlenmax+1, sizeof( char) );
		recordoriginalgaps( originalgaps, njob-nadd, seq );

		if( mapout )
		{
			addbk = (char **)calloc( nadd+1, sizeof( char * ) );
			for( i=0; i<nadd; i++ )
			{
				ien = strlen( seq[njob-nadd+i] );
				addbk[i] = (char *)calloc( ien + 1, sizeof( char ) );
				gappick0( addbk[i], seq[njob-nadd+i] );
			}
			addbk[nadd] = NULL;
		}
		else
			addbk = NULL;
	}
	else
	{
		originalgaps = NULL;
		addbk = NULL;
	}
	


	for( iguidetree=0; iguidetree<nguidetree; iguidetree++ )
//	for( iguidetree=0; ; iguidetree++ )
	{

		alg = algbackup; // tbfast wo disttbfast ni ketsugou shitatame.


		topol = AllocateIntCub( njob, 2, 0 );
		len = AllocateFloatMtx( njob, 2 );

		if( iguidetree == nguidetree - 1 ) calcpairdists = 0;
		else                               calcpairdists = 1;
	
		if( treein )
		{
			nguidetree = 1; //  iranai
			calcpairdists = 0; // iranai
			if( treein == (int)'l' )
			{
				loadtree( njob, topol, len, name, nogaplen, dep, treeout );
			}
			else if( treein == (int)'s' )
			{
				createchain( njob, topol, len, name, nogaplen, dep, treeout, 1, randomseed );
				nthreadpair = nthread = 0;
				weight = 0; // mafft.tmpl kara idou
				tbrweight = 0; // mafft.tmpl kara idou
				nthreadtb = 0; // 2017/Nov
			}
			else if( treein == (int)'p' )
			{
				createchain( njob, topol, len, name, nogaplen, dep, treeout, 0, randomseed );
				nthreadpair = nthread = 0;
				weight = 0; // mafft.tmpl kara idou
				tbrweight = 0; // mafft.tmpl kara idou
				nthreadtb = 0; // 2017/Nov
			}
			else
			{
				reporterr( "Error. treein = %d or %c\n", treein, treein );
				exit( 1 );
			}
		}
		else if( topin )
		{
			reporterr(       "Loading a topology ... " );
			reporterr(       "--topin has been disabled\n" );
			exit( 1 );
//			loadtop( njob, mtx, topol, len );
//			FreeFloatHalfMtx( mtx, njob );
		}
		else
		{
			if( distout )
			{
				hat2p = fopen( "hat2", "w" );
				WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
				// writehat2 wo kakinaosu
				fclose( hat2p );

				if( !treeout && noalign )  // 2016Jul31
				{
					writeData_pointer( stdout, njob, name, nlen, seq );
					reporterr(       "\n" );
					SHOWVERSION;
					goto chudan;
//					return( 0 );
				}
			}

			if( subalignment ) // merge ha localmem ni mitaiou
			{
				reporterr(       "Constructing a UPGMA tree ... " );
				fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained( njob, mtx, topol, len, name, nlen, dep, nsubalignments, subtable, !calcpairdists );
				if( !calcpairdists ) 
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else if( compacttree )
			{
				reporterr(       "Constructing a tree ... " );
// 2017Oct26
				if( compacttree == 4 )
					compacttree_memsaveselectable( njob, partmtx, mindistfrom, mindist, pointt, selfscore, bseq, skiptable, topol, len, name, nogaplen, dep, treeout, 2, 1 ); // howcompact == 2
				else
					compacttreegivendist( njob, mindist, mindistfrom, topol, len, name, dep, treeout );
				if( mindistfrom ) free( mindistfrom ); mindistfrom = NULL;
				if( mindist ) free( mindist );; mindist = NULL;
				if( selfscore ) free( selfscore ); selfscore = NULL;
				if( bseq ) FreeCharMtx( bseq ); bseq = NULL; // nikaime dake
				if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL; // nikaime dake
				if( pointt ) FreeIntMtx( pointt ); pointt = NULL; // ikkaime dake.
				free( partmtx );
			}
			else if( treeout )
			{
				reporterr(       "Constructing a UPGMA tree (treeout, efffree=%d) ... ", !calcpairdists );
				fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave( njob, mtx, topol, len, name, nogaplen, dep, !calcpairdists, treeout );
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else
			{
				reporterr(       "Constructing a UPGMA tree (efffree=%d) ... ", !calcpairdists );
				fixed_musclesupg_double_realloc_nobk_halfmtx_memsave( njob, mtx, topol, len, dep, 1, !calcpairdists );
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
		}
//		else 
//			ErrorExit( "Unknown tree method\n" );




		if( calcpairdists ) selfscore = AllocateIntVec( njob );


		if( callback && callback( 0, 25, "Guide tree" ) ) goto chudan;
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 50, "Guide tree" ) ) goto chudan;

		if( sparsepickup && iguidetree == nguidetree-1 )
		{
			reporterr(       "Sparsepickup! \n" );
			pickup( njob, nogaplen, topol, name, seq );
			reporterr(       "done. \n" );
			SHOWVERSION;
			goto chudan;
		}
//		reporterr( "after tree building" );
//		use_getrusage();


		if( treein == 's' || treein == 'p' )
		{
			localmem[0][0] = topol[0][0][0];
			for( i=1; i<njob; i++ )
				localmem[0][i] = topol[i-1][1][0];
		}
		else
		{
			localmem[0][0] = -1;
			posinmem = topolorderz( localmem[0], topol, dep, njob-2, 2 ) - localmem[0];
		}
	
		orderfp = fopen( "order", "w" );
		if( !orderfp )
		{
			reporterr(       "Cannot open 'order'\n" );
			exit( 1 );
		}
#if 0
		for( i=0; (j=topol[njob-2][0][i])!=-1; i++ )
		{
			fprintf( orderfp, "%d\n", j );
		}
		for( i=0; (j=topol[njob-2][1][i])!=-1; i++ )
		{
			fprintf( orderfp, "%d\n", j );
		}
#else
		for( i=0; i<njob; i++ )
			fprintf( orderfp, "%d\n", localmem[0][i] );
#endif
		fclose( orderfp );

	
		if( ( treeout || distout )  && noalign ) 
		{
			writeData_pointer( stdout, njob, name, nlen, seq );
			reporterr(       "\n" );
			SHOWVERSION;
			goto chudan;
//			return( 0 );
		}
		
		if( tbrweight )
		{
			weight = 3; 
#if 0
			utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
			counteff_simple_double_nostatic_memsave( njob, topol, len, dep, eff );
//			counteff_simple_double_nostatic( njob, topol, len, eff );
#endif
		}
		else
		{
			for( i=0; i<njob; i++ ) eff[i] = 1.0;
		}
	
#if 0
		for( i=0; i<njob; i++ )
			reporterr(       "eff[%d] = %20.16f\n", i, eff[i] );
		exit( 1 );
#endif
	
	
		FreeFloatMtx( len ); len = NULL;
	
		bseq = AllocateCharMtx( njob, nlenmax*2+1 );
		alloclen = nlenmax*2+1;


	
		if( nadd )
		{
			alignmentlength = strlen( seq[0] );
			for( i=0; i<njob-nadd; i++ )
			{
				if( alignmentlength != strlen( seq[i] ) )
				{
					reporterr(       "#################################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# The original %d sequences must be aligned\n", njob-nadd );
					reporterr(       "# alignmentlength = %d, but strlen(seq[%d])=%d\n", alignmentlength, i, (int)strlen( seq[i] ) );
					reporterr(       "#################################################################################\n" );
					goto chudan; // TEST!!
					//exit( 1 );
				}
			}
			if( addprofile )
			{
				alignmentlength = strlen( seq[njob-nadd] );
				for( i=njob-nadd; i<njob; i++ )
				{
					if( alignmentlength != strlen( seq[i] ) )
					{
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# The %d additional sequences must be aligned\n", nadd );
						reporterr(       "# Otherwise, try the '--add' option, instead of '--addprofile' option.\n" );
						reporterr(       "###############################################################################\n" );
						exit( 1 );
					}
				}
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				foundthebranch = 0;
				for( i=0; i<njob-1; i++ )
				{
					localmem[0][0] = -1;
					posinmem = topolorderz( localmem[0], topol, dep, i, 0 ) - localmem[0];
					localmem[1][0] = -1;
					posinmem = topolorderz( localmem[1], topol, dep, i, 1 ) - localmem[1];

					if( samemember( localmem[0], addmem ) ) // jissainiha nai
					{
						mergeoralign[i] = '1';
						foundthebranch = 1;
					}
					else if( samemember( localmem[1], addmem ) ) // samemembern ni henkou kanou
					{
						mergeoralign[i] = '2';
						foundthebranch = 1;
					}
					else
					{
						mergeoralign[i] = 'n';
					}
				}
				if( !foundthebranch )
				{
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# There is no appropriate position to add the %d sequences in the guide tree.\n", nadd );
					reporterr(       "# Check whether the %d sequences form a monophyletic cluster.\n", nadd );
					reporterr(       "# If not, try the '--add' option, instead of the '--addprofile' option.\n" );
					reporterr(       "############################################################################### \n" );
					exit( 1 );
				}
				commongappick( nadd, seq+njob-nadd );
				for( i=njob-nadd; i<njob; i++ ) strcpy( bseq[i], seq[i] );
			}
			else
			{
				for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'n';
#if 0
				for( j=njob-nadd; j<njob; j++ )
				{
					addmem[0] = j;
					addmem[1] = -1;
					for( i=0; i<njob-1; i++ )
					{
						reporterr( "Looking for samemember, %d-%d/%d\n", j, i, njob );
						localmem[0][0] = -1;
						posinmem = 0;
						topolorder( localmem[0], &posinmem, topol, dep, i, 0 );
						localmem[1][0] = -1;
						posinmem = 0;
						topolorder( localmem[1], &posinmem, topol, dep, i, 1 );

						if( samemembern( localmem[0], addmem, 1 ) ) // arieru
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '1';
						}
						else if( samemembern( localmem[1], addmem, 1 ) )
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '2';
						}
					}
				}
#else
				for( i=0; i<njob-1; i++ )
				{
//					reporterr( "Looking for samemember, %d-%d/%d\n", j, i, njob );
					localmem[0][0] = -1;
					posinmem = topolorderz( localmem[0], topol, dep, i, 0 ) - localmem[0];
					localmem[1][0] = -1;
					posinmem = topolorderz( localmem[1], topol, dep, i, 1 ) - localmem[1];

					for( j=njob-nadd; j<njob; j++ )
					{
						addmem[0] = j;
						addmem[1] = -1;

						if( samemembern( localmem[0], addmem, 1 ) ) // arieru
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '1';
						}
						else if( samemembern( localmem[1], addmem, 1 ) )
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '2';
						}
					}
				}
#endif
		
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				nlim = njob-1;
//				for( i=0; i<njob-1; i++ )
				for( i=0; i<nlim; i++ )
				{
					localmem[0][0] = -1;
					posinmem = topolorderz( localmem[0], topol, dep, i, 0 ) - localmem[0];
					localmem[1][0] = -1;
					posinmem = topolorderz( localmem[1], topol, dep, i, 1 ) - localmem[1];

					includememberres0 = includemember( localmem[0], addmem );
					includememberres1 = includemember( localmem[1], addmem );
//					if( includemember( topol[i][0], addmem ) && includemember( topol[i][1], addmem ) )
					if( includememberres0 && includememberres1 )
					{
						mergeoralign[i] = 'w';
					}
					else if( includememberres0 )
					{
						mergeoralign[i] = '1';
					}
					else if( includememberres1 )
					{
						mergeoralign[i] = '2';
					}
				}
#if 0
				for( i=0; i<njob-1; i++ )
				{
					reporterr(       "mem0 = " );
					for( j=0; topol[i][0][j]>-1; j++ )	reporterr(       "%d ", topol[i][0][j] );
					reporterr(       "\n" );
					reporterr(       "mem1 = " );
					for( j=0; topol[i][1][j]>-1; j++ )	reporterr(       "%d ", topol[i][1][j] );
					reporterr(       "\n" );
					reporterr(       "i=%d, mergeoralign[] = %c\n", i, mergeoralign[i] );
				}
#endif
				for( i=njob-nadd; i<njob; i++ ) gappick0( bseq[i], seq[i] );
			}
	
//			if( !keeplength ) commongappick( njob-nadd, seq );
			commongappick( njob-nadd, seq );

			for( i=0; i<njob-nadd; i++ ) strcpy( bseq[i], seq[i] );

		}
//--------------- kokokara ----
		else if( subalignment )
		{
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
			for( i=0; i<nsubalignments; i++ )
			{
				reporterr(       "Checking subalignment %d:\n", i+1 );
				alignmentlength = strlen( seq[subtable[i][0]] );
//				for( j=0; subtable[i][j]!=-1; j++ )
//					reporterr(       " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
				for( j=0; subtable[i][j]!=-1; j++ )
				{
					if( subtable[i][j] >= njob ) // check sumi
					{
						reporterr(       "No such sequence, %d.\n", subtable[i][j]+1 );
						exit( 1 );
					}
					if( alignmentlength != strlen( seq[subtable[i][j]] ) )
					{
						reporterr(       "\n" );
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# Subalignment %d must be aligned.\n", i+1 );
						reporterr(       "# Please check the alignment lengths of following sequences.\n" );
						reporterr(       "#\n" );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][0]+1, name[subtable[i][0]]+1, alignmentlength );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j]+1, name[subtable[i][j]]+1, (int)strlen( seq[subtable[i][j]] ) );
						reporterr(       "#\n" );
						reporterr(       "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" );
						if( subalignmentoffset )
						{
							reporterr(       "#\n" );
							reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
							reporterr(       "# In this case, the rule of numbering is:\n" );
							reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
							reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
						}
						reporterr(       "###############################################################################\n" );
						reporterr(       "\n" );
						goto chudan; // TEST!!
						//exit( 1 );
					}
					insubtable[subtable[i][j]] = 1;
				}
				for( j=0; j<njob-1; j++ )
				{
					if( includemember( topol[j][0], subtable[i] ) && includemember( topol[j][1], subtable[i] ) )
					{
						mergeoralign[j] = 'n';
					}
				}
				foundthebranch = 0;
				for( j=0; j<njob-1; j++ )
				{
					if( samemember( topol[j][0], subtable[i] ) || samemember( topol[j][1], subtable[i] ) )
					{
						foundthebranch = 1;
						reporterr(       " -> OK\n" );
						break;
					}
				}
				if( !foundthebranch )
				{
					system( "cp infile.tree GuideTree" ); // tekitou
					reporterr(       "\n" );
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# Subalignment %d does not seem to form a monophyletic cluster\n", i+1 );
					reporterr(       "# in the guide tree ('GuideTree' in this directory) internally computed.\n" );
					reporterr(       "# If you really want to use this subalignment, pelase give a tree with --treein \n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/treein.html\n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/merge.html\n" );
					if( subalignmentoffset )
					{
						reporterr(       "#\n" );
						reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
						reporterr(       "# In this case, the rule of numbering is:\n" );
						reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
						reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
					}
					reporterr(       "############################################################################### \n" );
					reporterr(       "\n" );
					goto chudan; // TEST!!
					//exit( 1 );
				}
//				commongappick( seq[subtable[i]], subalignment[i] ); // irukamo
			}
#if 0
			for( i=0; i<njob-1; i++ )
			{
				reporterr(       "STEP %d\n", i+1 );
				reporterr(       "group1 = " );
				for( j=0; topol[i][0][j] != -1; j++ )
					reporterr(       "%d ", topol[i][0][j]+1 );
				reporterr(       "\n" );
				reporterr(       "group2 = " );
				for( j=0; topol[i][1][j] != -1; j++ )
					reporterr(       "%d ", topol[i][1][j]+1 );
				reporterr(       "\n" );
				reporterr(       "%d -> %c\n\n", i, mergeoralign[i] );
			}
#endif
	
			for( i=0; i<njob; i++ ) 
			{
				if( insubtable[i] ) strcpy( bseq[i], seq[i] );
				else gappick0( bseq[i], seq[i] );
			}
	
			for( i=0; i<nsubalignments; i++ ) 
			{
				for( j=0; subtable[i][j]!=-1; j++ ) subalnpt[i][j] = bseq[subtable[i][j]];
				if( !preservegaps[i] ) commongappick( j, subalnpt[i] );
			}
	
#if 0 // --> iguidetree loop no soto he
			FreeIntMtx( subtable );
			free( insubtable );
			for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
			free( subalnpt );
			free( preservegaps );
#endif
		}
//--------------- kokomade ----
		else
		{
			for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
		}

		if( calcpairdists ) for( i=0; i<njob; i++ ) selfscore[i] = naivepairscore11( seq[i], seq[i], penalty_dist ); // (int)?
	
		reporterr(       "Progressive alignment %d/%d... \n", iguidetree+1, nguidetree );

//		reporterr( "\nbefore treebase" );
//		use_getrusage();

		cpmxhist = (double ***)calloc( njob-1, sizeof( double ** ) );
		for( i=0; i<njob-1; i++ ) cpmxhist[i] = NULL; 

		memhist = (int **)calloc( njob-1, sizeof( int * ) );
		for( i=0; i<njob-1; i++ ) memhist[i] = NULL; 

#if REPORTCOSTS
		time_t starttime, startclock;
		starttime = time(NULL);
		startclock = clock();
#endif
	
#ifdef enablemultithread
		if( nthreadtb > 0 && nadd == 0 ) // nthreadpair ha minai
		{
			treebasethread_arg_t *targ; 
			int jobpos;
			pthread_t *handle;
			pthread_mutex_t mutex;
			pthread_cond_t treecond;
			int *fftlog;
			int nrun;
			int nthread_yoyu;
	
			nthread_yoyu = nthreadtb * 1;
			nrun = 0;
			jobpos = 0; 
			targ = calloc( nthread_yoyu, sizeof( treebasethread_arg_t ) ); 
			fftlog = AllocateIntVec( njob );
			handle = calloc( nthread_yoyu, sizeof( pthread_t ) ); 
			pthread_mutex_init( &mutex, NULL );
			pthread_cond_init( &treecond, NULL );
	
			for( i=0; i<njob; i++ ) dep[i].done = 0; 
			for( i=0; i<njob; i++ ) fftlog[i] = 1; 
	
			for( i=0; i<nthread_yoyu; i++ )
			{
				targ[i].thread_no = i;
				targ[i].njob = njob;
				targ[i].nrunpt = &nrun;
				targ[i].nlen = nlen;
				targ[i].jobpospt = &jobpos;
				targ[i].topol = topol;
				targ[i].dep = dep;
				targ[i].cpmxhist = cpmxhist;
				targ[i].memhist = memhist;
				targ[i].aseq = bseq;
				targ[i].effarr = eff;
				targ[i].alloclenpt = &alloclen;
				targ[i].fftlog = fftlog;
				targ[i].mergeoralign = mergeoralign;
				targ[i].extanch = extanch;
				targ[i].anchindex = anchindex;
#if 1 // tsuneni SEPARATELYCALCPAIRDISTS
				targ[i].newdistmtx = NULL;
				targ[i].selfscore = NULL;
#else
				if( calcpairdists ) // except for last cycle
				{
					targ[i].newdistmtx = mtx;
					targ[i].selfscore = selfscore;
				}
				else
				{
					targ[i].newdistmtx = NULL;
					targ[i].selfscore = NULL;
				}
#endif
				targ[i].mutex = &mutex;
				targ[i].treecond = &treecond;
	
				pthread_create( handle+i, NULL, treebasethread, (void *)(targ+i) );
			}

			for( i=0; i<nthread_yoyu; i++ )
			{
				pthread_join( handle[i], NULL );
			}
			pthread_mutex_destroy( &mutex );
			pthread_cond_destroy( &treecond );
			free( handle );
			free( targ );
			free( fftlog );

		}
		else
#endif
		{
#if 0
			if( calcpairdists ) // except for last
			{
				if( treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, mtx, selfscore, &alloclen, callback ) ) goto chudan;
			}
			else
#endif
			{
//				if( treebase( keeplength && (iguidetree==nguidetree-1), nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, NULL, NULL, deletemap, deletelag, &alloclen, callback ) ) goto chudan;
				if( treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, memhist, cpmxhist, eff, NULL, NULL, extanch, anchindex, &alloclen, callback ) ) goto chudan;
			}
		}


#if REPORTCOSTS
//		use_getrusage();
		reporterr( "\ntb %d, real = %f min\n", iguidetree, (float)(time(NULL) - starttime)/60.0 );
		reporterr( "tb %d, user = %f min\n", iguidetree, (float)(clock()-startclock)/CLOCKS_PER_SEC/60);
#endif
//		reporterr( "after treebase, " );
//		use_getrusage();
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 100, "Progressive alignment" ) ) goto chudan;
//		free( topol[njob-1][0] ); topol[njob-1][0]=NULL;
//		free( topol[njob-1][1] ); topol[njob-1][1]=NULL;
//		free( topol[njob-1] ); topol[njob-1]=NULL;
//		free( topol ); topol=NULL;
		FreeIntCub( topol ); topol = NULL;
#if 1 // 2021/Jun/24
		if( cpmxhist ) // nakutemo yoi
		{
			for( i=0; i<njob-1; i++ )
			{
				if( cpmxhist[i] ) 
				{
//					reporterr( "freeing cpmxhist[%d]\n", i );
					FreeDoubleMtx( cpmxhist[i] ); cpmxhist[i] = NULL;
				}
			}
			free( cpmxhist ); cpmxhist = NULL;
		}
#else
		if( cpmxhist[njob-2] ) 
		{
//			reporterr( "freeing cpmxhist[njob-2]\n" );
			FreeDoubleMtx( cpmxhist[njob-2] ); cpmxhist[njob-2] = NULL;
		}
		free( cpmxhist ); cpmxhist = NULL;
#endif

		free( memhist ); memhist = NULL;
//		reporterr( "after freeing topol, " );
//		use_getrusage();


//		reporterr( "compacttree = %d, calcpairdist = %d\n", compacttree, calcpairdists );


//		reporterr( "\nbseq[njob-3] = %s\n", bseq[njob-3] );
//		reporterr( "bseq[njob-2] = %s\n", bseq[njob-2] );
//		reporterr( "bseq[njob-1] = %s\n", bseq[njob-1] );



// Distance matrix from MSA SEPARATELYCALCPAIRDISTS
//		if( iguidetree < nguidetree-1 )
#ifdef enablemultithread
//		if( nthread>0 && nadd==0 ) if( calcpairdists )
		if( calcpairdists && !compacttree )
#else
//		if( 0 && nadd==0 ) if( calcpairdists ) // zettai nai
		if( calcpairdists && !compacttree )
#endif
		{
			reporterr( "Making a distance matrix from msa.. \n" );
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, bseq ); // allocate suru.
#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
	
				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( nthreadpair, sizeof( msadistmtxthread_arg_t ) );
				handle = calloc( nthreadpair, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
	
				for( i=0; i<nthreadpair; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].iscore = mtx;
					targ[i].seq = bseq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, msadistmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthreadpair; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
//				reporterr( "Check source!\n" );
//				exit( 1 );

#if 1
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;

				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( 1, sizeof( msadistmtxthread_arg_t ) );
	
				{
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].iscore = mtx;
					targ[0].seq = bseq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
	
					msadistmtxthread( targ );
				}
	
				free( targ );
#endif
			}
			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;
			reporterr(       "\ndone.\n\n" );
			free( selfscore ); selfscore = NULL;
			FreeCharMtx( bseq ); bseq = NULL;
		}
		else if( calcpairdists && compacttree )
		{
			reporterr( "Making a compact tree from msa, step 1.. \n" );
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, bseq ); // allocate suru.
			mindistfrom = (int *)calloc( njob, sizeof( int ) );
			mindist = (double *)calloc( njob, sizeof( double ) );
			partmtx = preparepartmtx( njob );
#ifdef enablemultithread
			if( nthreadpair > 0 )
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				double **mindistthread;
				int **mindistfromthread;

				mindistthread = AllocateDoubleMtx( nthreadpair, njob );
				mindistfromthread = AllocateIntMtx( nthreadpair, njob );
				targ = calloc( nthreadpair, sizeof( msacompactdistmtxthread_arg_t ) );
				handle = calloc( nthreadpair, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;

				for( i=0; i<nthreadpair; i++ )
				{
					for( j=0; j<njob; j++ )
					{
						mindistthread[i][j] = 999.9;
						mindistfromthread[i][j] = -1;
					}
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].partmtx = partmtx;
					targ[i].seq = bseq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mindistfrom = mindistfromthread[i];
					targ[i].mindist = mindistthread[i];
					targ[i].mutex = &mutex;
	
					if( compacttree == 4 )
						pthread_create( handle+i, NULL, ylmsacompactdisthalfmtxthread, (void *)(targ+i) );
					else
						pthread_create( handle+i, NULL, msacompactdisthalfmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthreadpair; i++ ) pthread_join( handle[i], NULL );
				pthread_mutex_destroy( &mutex );

				for( i=0; i<njob; i++ )
				{
					mindist[i] = 999.9;
					mindistfrom[i] = -1;
					for( j=0; j<nthreadpair; j++ )
					{
						if( mindistthread[j][i] < mindist[i] )
						{
							mindist[i] = mindistthread[j][i];
							mindistfrom[i] = mindistfromthread[j][i];
						}
					}
				}
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug

				free( handle );
				free( targ );
				FreeDoubleMtx( mindistthread );
				FreeIntMtx( mindistfromthread );
			}
			else
#endif
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				if( compacttree == 4 )
					jobpos = 0;
				else
					jobpos = njob-1;
				targ = calloc( 1, sizeof( msacompactdistmtxthread_arg_t ) );

				{
					for( j=0; j<njob; j++ )
					{
						mindist[j] = 999.9;
						mindistfrom[j] = -1;
					}
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].partmtx = partmtx;
					targ[0].seq = bseq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
					targ[0].mindistfrom = mindistfrom;
					targ[0].mindist = mindist;
	
					if( compacttree == 4 )
						ylmsacompactdisthalfmtxthread( targ );
					else
						msacompactdisthalfmtxthread( targ );
//					msacompactdistmtxthread( targ );
				}
				free( targ );
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug
			}
			reporterr( "\n" );
//			free( selfscore ); selfscore = NULL; // mada tsukau
//			FreeCharMtx( bseq ); bseq = NULL; // mada tsukau
//			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;

//			reporterr( "\n" );
//			for( i=0; i<njob; i++ ) reporterr( "mindist[%d] = %f\n", i, mindist[i] );
//			exit( 1 );
		}
// Distance matrix from MSA end
//		reporterr( "at the end of guidetree loop, " );
//		use_getrusage();

	}


	if( oneiteration )
	{
		reporterr( "Iterative refinement (one vs others)\n" );
		dooneiteration( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, memhist, cpmxhist, eff, NULL, NULL, extanch, anchindex, &alloclen, callback );
	}

#if DEBUG
	reporterr(       "closing trap_g\n" );
#endif
//	fclose( trap_g );
//	reporterr( "after guidetree loop, " );
//	use_getrusage();

	if( keeplength )
	{

		dlf = fopen( "_deletelist", "w" );
		deletelist = (GapPos **)calloc( nadd+1, sizeof( GapPos * ) );
		for( i=0; i<nadd; i++ )
		{
			deletelist[i] = calloc( 1, sizeof( GapPos ) );
			deletelist[i][0].pos = -1;
			deletelist[i][0].len = 0;
		}
		deletelist[nadd] = NULL;
		ndeleted = deletenewinsertions_whole( njob-nadd, nadd, bseq, bseq+njob-nadd, deletelist );

		for( i=0; i<nadd; i++ )
		{
			if( deletelist[i] )
				for( j=0; deletelist[i][j].pos!=-1; j++ )
//					fprintf( dlf, "%d %d\n", njob-nadd+i, deletelist[i][j] ); // 0origin
					fprintf( dlf, "%d %d %d\n", njob-nadd+i, deletelist[i][j].pos, deletelist[i][j].len ); // 0origin
		}
		fclose( dlf );

		restoreoriginalgaps( njob, bseq, originalgaps );
		free( originalgaps ); originalgaps = NULL; // 2017/Nov/15

		if( mapout )
		{
			dlf = fopen( "_deletemap", "w" );
			if( mapout == 1 )
				reconstructdeletemap( nadd, addbk, deletelist, bseq+njob-nadd, dlf, name+njob-nadd );
			else
				reconstructdeletemap_compact( nadd, addbk, deletelist, seq+njob-nadd, dlf, name+njob-nadd );
			FreeCharMtx( addbk );
			addbk = NULL;
			fclose( dlf );
		}

//		FreeIntMtx( deletelist );
//		deletelist = NULL;
		for( i=0; deletelist[i] != NULL; i++ ) free( deletelist[i] );
		free( deletelist );
		deletelist = NULL;
	}

	if( scoreout )
	{
		unweightedspscore = plainscore( njob, bseq );
		reporterr(       "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore );
		reporterr(       "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
		reporterr(       "\n\n" );
	}

#if DEBUG
	reporterr(       "writing alignment to stdout\n" );
#endif


	val = 0;
	if( ngui ) 
	{
		ien = strlen( bseq[0] );
		if( ien > lgui )
		{
			reporterr( "alignmentlength = %d, gui allocated %d", ien, lgui );
			val = GUI_LENGTHOVER;
		}
		else
		{
			for( i=0; i<njob; i++ ) 
			{
#if 1
				strcpy( seqgui[i], bseq[i] );
#else
				free( seqgui[i] );
				seqgui[i] =  bseq[i];
#endif
			}
		}
	}
	else
	{
		writeData_pointer( stdout, njob, name, nlen, bseq );
	} 

	if( spscoreout ) reporterr( "Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore( njob, bseq ) );
	nthread = nthreadtb; // toriaezu
	SHOWVERSION;
	if( ndeleted > 0 )
	{
		reporterr( "\nTo keep the alignment length, %d letters were DELETED.\n", ndeleted );
		if( mapout )
			reporterr( "The deleted letters are shown in the (filename).map file.\n" );
		else
			reporterr( "To know the positions of deleted letters, rerun the same command with the --mapout option.\n" );
	}



	if( subalignment )
	{
		FreeIntMtx( subtable );
		free( insubtable );
		for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
		free( subalnpt );
		free( preservegaps );
	}
	if( useexternalanchors )
	{
		free( extanch );
		FreeIntMtx( anchindex );
	}


#if 1 // seqgui[i] =  bseq[i] no toki bseq ha free shinai
	FreeCharMtx( bseq );
#endif
	FreeCharMtx( name );
    free( nlen );

	free( mergeoralign );
	FreeCharMtx( seq );
    free( nogaplen );

	free( mseq1 );
	free( mseq2 );
//	FreeIntCub( topol ); // 
//	FreeFloatMtx( len ); //
//	free( mergeoralign ); //
	free( dep );

	if( nadd ) free( addmem );
	FreeIntMtx( localmem );
	free( eff );
	freeconstants();
	closeFiles();
	FreeCommonIP();
	if( originalgaps ) free( originalgaps ); originalgaps = NULL;
	if( deletelist ) 
	{
		for( i=0; deletelist[i] != NULL; i++ ) free( deletelist[i] );
		free( deletelist );
		deletelist = NULL;
	}

//	use_getrusage();

	return( val );

chudan:

	if( nlen ) free( nlen ); nlen = NULL;
	if( seq ) FreeCharMtx( seq ); seq = NULL;
	if( mseq1 ) free( mseq1 ); mseq1 = NULL;
	if( mseq2 ) free( mseq2 ); mseq2 = NULL;
//	if( topol ) 
//	{
//		for( i=0; i<njob; i++ )
//		{
//			if( topol[i] && topol[i][0] ) 
//			{
//				free( topol[i][0] ); topol[i][0] = NULL;
//			}
//			if( topol[i] && topol[i][1] ) 
//			{
//				free( topol[i][1] ); topol[i][1] = NULL;
//			}
//			if( topol[i] ) free( topol[i] ); topol[i] = NULL;
//		}
//		free( topol ); topol = NULL;
//	}
	if( topol ) FreeIntCub( topol ); topol = NULL;
	if( len ) FreeFloatMtx( len ); len = NULL;
	if( eff ) free( eff ); eff = NULL;
	if( mergeoralign ) free( mergeoralign ); mergeoralign = NULL;
	if( dep ) free( dep ); dep = NULL;
	if( addmem ) free( addmem ); addmem = NULL;
	if( localmem ) FreeIntMtx( localmem ); localmem = NULL;
	if( name ) FreeCharMtx( name ); name = NULL;
	if( nogaplen ) free( nogaplen ); nogaplen = NULL;

	if( tmpseq ) free( tmpseq ); tmpseq = NULL;
	if( grpseq ) free( grpseq ); grpseq = NULL;
	if( pointt ) FreeIntMtx( pointt ); pointt = NULL;
	if( mtx ) FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
	if( table1 ) free( table1 ); table1 = NULL;

	if( bseq ) FreeCharMtx( bseq ); bseq = NULL;
	if( selfscore ) free( selfscore ); selfscore = NULL;
	if( skiptable ) FreeIntMtx( skiptable ); skiptable = NULL;
	if( originalgaps ) free( originalgaps ); originalgaps = NULL;
	if( deletelist ) 
	{
		for( i=0; deletelist[i] != NULL; i++ ) free( deletelist[i] );
		free( deletelist );
		deletelist = NULL;
	}


	if( subtable ) FreeIntMtx( subtable ); subtable = NULL;
	if( insubtable ) free( insubtable ); insubtable = NULL;
	for( i=0; i<nsubalignments; i++ ) 
	{
		if( subalnpt[i] ) free( subalnpt[i] ); subalnpt[i] = NULL;
	}
	if( subalnpt ) free( subalnpt ); subalnpt = NULL;
	if( preservegaps ) free( preservegaps ); preservegaps = NULL;


	if( mindistfrom ) free( mindistfrom ); mindistfrom = NULL;
	if( mindist ) free( mindist ); mindist = NULL;

	if( cpmxhist )
	{
		for( i=0; i<njob-1; i++ )
		{
			if( cpmxhist[i] ) FreeDoubleMtx( cpmxhist[i] ); cpmxhist[i] = NULL;
		}
		free( cpmxhist ); cpmxhist = NULL;
	}

	if( memhist )
	{
		for( i=0; i<njob-1; i++ )
		{
			if( memhist[i] ) free( memhist[i] ); memhist[i] = NULL;
		}
		free( memhist ); memhist = NULL;
	}

	freeconstants();
	closeFiles();
	FreeCommonIP();

	return( GUI_CANCEL );
}

int main( int argc, char **argv )
{
	int res = disttbfast( 0, 0, NULL, NULL, argc, argv, NULL );
	if( res == GUI_CANCEL ) res = 0; // treeout de goto chudan wo riyousuru
	return res;
}
