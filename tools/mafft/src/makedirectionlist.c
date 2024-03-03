#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

#define GLOBAL 0

#define END_OF_VEC -1

//int nadd;
double thresholdtorev;
int dodp;
int addfragment;
int mode = '2';
int reflim = 1000;
int contrastsort = 1;

typedef struct _thread_arg
{
	int iend; 
	char **seq;
	int *map;
	char *tmpseq;
	int *res;
	int **spointt;
	int *table1;
	int iq;
#ifdef enablemultithread
	int *jshare;
	int thread_no;
	pthread_mutex_t *mutex_counter;
#endif
} thread_arg_t;

typedef struct _selfdpthread_arg
{
	int iend;
	char **seq;
	double *res;
#ifdef enablemultithread
	int *jshare;
	int thread_no;
	pthread_mutex_t *mutex_counter;
#endif
} selfdpthread_arg_t;

typedef struct _contrast
{
	int pos; 
	double dif;
} contrastarr;

static void	*selfdpthread( void *arg )
{
	selfdpthread_arg_t *targ = (selfdpthread_arg_t *)arg;
	int iend = targ->iend;
	char **seq = targ->seq;
	double *res = targ->res;
#ifdef enablemultithread
	int thread_no = targ->thread_no;
	int *jshare = targ->jshare; 
#endif
	int j;
	char **revseq;

	revseq = AllocateCharMtx( 1, nlenmax+1 );

	j = -1;
	while( 1 )
	{
#ifdef enablemultithread
		if( nthread )
		{
			pthread_mutex_lock( targ->mutex_counter );
			j = *jshare;
			if( j%100 == 0 ) reporterr( "%d / %d (thread %d)   \r", j, iend, thread_no );
			if( j == iend )
			{
				pthread_mutex_unlock( targ->mutex_counter );
				break;
			}
			++(*jshare);
			pthread_mutex_unlock( targ->mutex_counter );
		}
		else
#endif
		{
			j++;
			if( j%100 == 0 ) reporterr( "%d / %d      \r", j, iend );
			if( j == iend ) 
			{
				break;
			}
		}

		sreverse( revseq[0], seq[j] );
#if GLOBAL
		res[j] =  G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, seq+j, seq+j, 0 );
		res[j] -= G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, seq+j, revseq, 0 );
#else
		res[j] =  L__align11_noalign( n_dis_consweight_multi, seq+j, seq+j );
		res[j] -= L__align11_noalign( n_dis_consweight_multi, seq+j, revseq );
#endif
	}

	creverse( 0 );
	FreeCharMtx( revseq );
#if GLOBAL
	G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
#else
	L__align11_noalign( NULL, NULL, NULL );
#endif
	return( NULL );
}

#if 0
static void partshuffle( int size, int outsize, int *ary )
{
	int i;

//	reporterr( "ary before shuffle = \n" );
 //   for(i=0;i<size;i++) reporterr( "%d ", ary[i] );
//	reporterr( "\n" );

    for(i=0;i<outsize;i++)
    {
        int j = rand()%size;
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }

//	reporterr( "ary after shuffle = \n" );
 //   for(i=0;i<outsize;i++) reporterr( "%d ", ary[i] );
//	reporterr( "|" );
 //   for(i=outsize;i<size;i++) reporterr( "%d ", ary[i] );
//	reporterr( "\n" );
}
#endif

void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	inputfile = NULL;
	nadd = 0;
	dodp = 0;
	alg = 'a';
	alg = 'm';
	dorp = NOTSPECIFIED;
	fmodel = 0;
//	ppenalty = (int)( -2.0 * 1000 - 0.5 );
//	ppenalty_ex = (int)( -0.1 * 1000 - 0.5 );
//	poffset = (int)( 0.1 * 1000 - 0.5 ); 
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = 2;
	pamN = 200;
	thresholdtorev = 0.0;
	addfragment = 0;


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
				case 'I':
					nadd = myatoi( *++argv );
					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
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
					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 't':
					thresholdtorev = atof( *++argv );
					fprintf( stderr, "thresholdtorev = %f\n", thresholdtorev );
					--argc; 
					goto nextoption;
				case 'o':
					mode = *(*++argv);
					fprintf( stderr, "mode = %c\n", mode );
					--argc; 
					goto nextoption;
				case 'r':
					reflim = myatoi(*++argv);
					fprintf( stderr, "reflim = %d\n", reflim );
					--argc; 
					goto nextoption;
				case 'c':
					contrastsort = 0;
					break;
				case 'd':
					dodp = 1;
					break;
				case 'F':
					addfragment = 1;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'S':
					alg = 'S';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'm':
					alg = 'm';
					break;
				case 'G':
					alg = 'G';
					break;
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
}





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
//			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
			;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
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
		tmp = amino_grp[(int)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
//			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
			;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
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
#else
		if( (unsigned int)table[point]++ >= INT_MAX )
		{
			reporterr( "Overflow. table[point]=%d>INT_MAX(%d).\n", table[point], INT_MAX );
			exit( 1 );
		}
#endif
	}
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

static int localcommonsextet_p2( int *table, int *pointt )
{
	int value = 0;
	unsigned int tmp;
	int point;
	int *memo;
	int *ct;
	int *cp;

	if( *pointt == -1 )
		return( 0 );

	memo = (int *)calloc( tsize, sizeof( int ) );
	if( !memo ) ErrorExit( "Cannot allocate memo\n" );
	ct = (int *)calloc( MIN( maxl, tsize )+1, sizeof( int ) );
	if( !ct ) ErrorExit( "Cannot allocate memo\n" );

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
#if 0 // kakunin shinai
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

	free( memo );
	free( ct );
	return( value );
}

static int compfunc( const void *a, const void *b )
{
	return ((contrastarr *)b)->dif - ((contrastarr *)a)->dif; // correct
//	return ((contrastarr *)a)->dif - ((contrastarr *)b)->dif; // incorrect!
} 

static void makecontrastorder6mer( int *order, int **pointt, int **pointt_rev, char **seq, int iend, int shift )
{
	int i;
	double *res;
	contrastarr *arr;
	int *table1, *table1_rev;


	arr = calloc( iend, sizeof( contrastarr ) );
	res = calloc( iend, sizeof( double ) );

	for( i=0; i<iend; i++ )
	{
		if( i % 100 == 1 ) reporterr( "%d   \r", i );
		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		makecompositiontable_p( table1, pointt[i] );
		res[i] = localcommonsextet_p2( table1, pointt[i] );
		free( table1 );

		table1_rev = (int *)calloc( tsize, sizeof( int ) );
		if( !table1_rev ) ErrorExit( "Cannot allocate table1\n" );
		makecompositiontable_p( table1_rev, pointt_rev[i] );
		res[i] -= localcommonsextet_p2( table1_rev, pointt[i] );
		free( table1_rev );

	}

	for( i=0; i<iend; i++ )
	{
		arr[i].pos = i;
		arr[i].dif = res[i];
	}

	qsort( arr, iend, sizeof( contrastarr ), compfunc );

	for( i=0; i<iend; i++ )
		order[i] = arr[i].pos + shift;

//	for( i=0; i<iend; i++ ) reporterr( "%f\n", arr[i].dif );
//	reporterr( "highest contrast, %s\n", seq[order[0]] );
//	reporterr( "lowest contrast, %s\n", seq[order[iend-1]] );

	free( arr );
	free( res );

}
static void makecontrastorder( int *order, char **seq, int iend, int shift )
{
	int i;
	double *res;
	contrastarr *arr;

	arr = calloc( iend, sizeof( contrastarr ) );
	res = calloc( iend, sizeof( double ) );

#ifdef enablemultithread
	if( nthread )
	{
		int j;
		pthread_t *handle;
		pthread_mutex_t mutex_counter;
		selfdpthread_arg_t *targ;
		int *jsharept;
		
		targ = calloc( nthread, sizeof( selfdpthread_arg_t ) );
		handle = calloc( nthread, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex_counter, NULL );
		jsharept = calloc( 1, sizeof(int) );
		*jsharept = 0;
		
		for( j=0; j<nthread; j++ )
		{
			targ[j].iend = iend;
			targ[j].seq = seq;
			targ[j].res = res; 
			targ[j].jshare = jsharept;
			targ[j].mutex_counter = &mutex_counter;
			targ[j].thread_no = j;
			pthread_create( handle+j, NULL, selfdpthread, (void *)(targ+j) );
		}
		for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
		pthread_mutex_destroy( &mutex_counter );
		free( handle );
		free( targ );
		free( jsharept );
	}
	else
#endif
	{
		selfdpthread_arg_t *targ;
		targ = calloc( 1, sizeof( selfdpthread_arg_t ) );
		targ[0].iend = iend;
		targ[0].seq = seq;
		targ[0].res = res; 
		selfdpthread( targ );
		free( targ );
	}

	for( i=0; i<iend; i++ )
	{
		arr[i].pos = i;
		arr[i].dif = res[i];
	}

	qsort( arr, iend, sizeof( contrastarr ), compfunc );

	for( i=0; i<iend; i++ )
		order[i] = arr[i].pos + shift;

//	for( i=0; i<iend; i++ ) reporterr( "%f\n", arr[i].dif );
//	reporterr( "highest contrast, %s\n", seq[order[0]] );
//	reporterr( "lowest contrast, %s\n", seq[order[iend-1]] );

	free( arr );
	free( res );

}


static void	*directionthread( void *arg )
{
	thread_arg_t *targ = (thread_arg_t *)arg;
	int iend = targ->iend;
	char **seq = targ->seq;
	int *map = targ->map;
	char *tmpseq = targ->tmpseq;
	int *res = targ->res;
	int **spointt = targ->spointt;
	int *table1 = targ->table1;
//	int iq = targ->iq;
#ifdef enablemultithread
//	int thread_no = targ->thread_no;
	int *jshare = targ->jshare; 
#endif
	int j;
	char **mseq1, **mseq2;


	if( dodp ) // nakuserukamo
	{
		mseq1 = AllocateCharMtx( 1, 0 );
		mseq2 = AllocateCharMtx( 1, 0 );
	}

	j = -1;
	while( 1 )
	{
#ifdef enablemultithread
		if( nthread )
		{
			pthread_mutex_lock( targ->mutex_counter );
			j = *jshare;
			if( j == iend )
			{
				pthread_mutex_unlock( targ->mutex_counter );
				break;
			}
			++(*jshare);
			pthread_mutex_unlock( targ->mutex_counter );
		}
		else
#endif
		{
			j++;
			if( j == iend ) 
			{
//				if( iq%100==1 ) fprintf( stderr, "\r %d / %d  \r", iq, njob );
				break;
			}
		}


		if( dodp )
		{
//			strcpy( mseq1[0], tmpseq );
//			strcpy( mseq2[0], seq[j] );
			mseq1[0] = tmpseq;
			mseq2[0] = seq[map[j]];
#if GLOBAL
			res[j] = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, 0 );
#else
			res[j] = L__align11_noalign( n_dis_consweight_multi, mseq1, mseq2 );
#endif
		}
		else
		{
//			reporterr( "\n\nj=%d, map[j]=%d\n\n", j, map[j] );
			res[j] = localcommonsextet_p2( table1, spointt[map[j]] );
		}
	}
	if( dodp ) // nakuserukamo
	{
		free( mseq1 );
		free( mseq2 );
#if GLOBAL
		G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
#else
		L__align11_noalign( NULL, NULL, NULL );
#endif
	}
//	else
//		if( nthread )  // inthread == 0 no toki free suru to, error. nazeda
//			localcommonsextet_p( NULL, NULL );
	return( NULL );
}

int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static int  *nogaplen;	
	static char **name, **seq;
	int i, j, istart, iend, ic;
	FILE *infp;
//	FILE *adfp;
	char c;

	int *grpseq;
	char *tmpseq, *revseq;
	int  **pointt, **pointt_rev, **spointt;
	double res_forward, res_reverse, res_max;
	int ires, mres, mres2;
	int *res, *resr, *resf;
	int *map;
	static int *table1, *table1_rev;
	static char **mseq1f, **mseq1r, **mseq2;
	int *contrastorder;

	arguments( argc, argv );
#ifndef enablemultithread
	nthread = 0;
#endif

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

	if( alg == 'a' )
	{
		if( nlenmax < 10000 )
			alg = 'G';
		else
			alg = 'S';
	}

	seq = AllocateCharMtx( njob, nlenmax*1+1 );

#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob ); 
    nogaplen = AllocateIntVec( njob ); 
	readData_pointer( infp, name, nlen, seq );
	fclose( infp );

	if( dorp != 'd' )
	{
		fprintf( stderr, "Not necessary!\n" );
		for( i=0; i<njob; i++ ) 
			fprintf( stdout, "_F_%-10.10s\n", name[i]+1 );
		exit( 1 );
	}
#endif

	constants( njob, seq );


#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

	fprintf( stderr, "\n" );
	if( alg == 'G' ) // dp to the first sequence
	{
		mseq1f = AllocateCharMtx( 1, nlenmax+nlenmax );
		mseq1r = AllocateCharMtx( 1, nlenmax+nlenmax );
		mseq2 = AllocateCharMtx( 1, nlenmax+nlenmax );
	    tmpseq = AllocateCharVec( MAX( nlenmax, B ) +1 );

		gappick0( mseq1f[0], seq[0] );
		sreverse( mseq1r[0], mseq1f[0] );
		strcpy( seq[0], mseq1f[0] );

		if( nadd )
			istart = njob - nadd;
		else
			istart = 1;

		fprintf( stderr, "\n" );

		for( i=0; i<istart; i++ )
		{
			gappick0( tmpseq, seq[i] );
			strcpy( seq[i], tmpseq );
			strcpy( tmpseq, name[i] );
			strcpy( name[i], "_F_" );
			strncpy( name[i]+3, tmpseq+1, 10 );
			name[i][13] = 0;
		}
		for( i=istart; i<njob; i++ ) 
		{
			gappick0( mseq2[0], seq[i] );

#if GLOBAL
			res_forward = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1f, mseq2, 0 );
			res_reverse = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1r, mseq2, 0 );
#else
			res_forward = L__align11_noalign( n_dis_consweight_multi, mseq1f, mseq2 );
			res_reverse = L__align11_noalign( n_dis_consweight_multi, mseq1r, mseq2 );
#endif
#if 0

			strcpy( mseq2[0], seq[i] );
			strcpy( mseq1f[0], seq[0] );
			res_forward = G__align11( n_dis_consweight_multi, mseq1f, mseq2, nlenmax*2, 0, 0 );
			fprintf( stdout, "%s\n", mseq1f[0] );
			fprintf( stdout, "%s\n", mseq2[0] );

			strcpy( mseq2[0], seq[i] );
			sreverse( mseq1r[0], seq[0] );
			res_reverse = G__align11( n_dis_consweight_multi, mseq1r, mseq2, nlenmax*2, 0, 0 );
			fprintf( stdout, "%s\n", mseq1r[0] );
			fprintf( stdout, "%s\n", mseq2[0] );
#endif

//			fprintf( stdout, "\nscore_for(%d,%d) = %f\n", 0, i, res_forward );
//			fprintf( stdout, "score_rev(%d,%d) = %f\n", 0, i, res_reverse );
			res_max = MAX(res_reverse,res_forward);
			if( (res_reverse-res_forward)/res_max > thresholdtorev ) // tekitou
			{
//				fprintf( stderr, "REVERSE!!!\n" );
				sreverse( seq[i], mseq2[0] );

				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_R_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
			}
			else
			{
				strcpy( seq[i], mseq2[0] );

				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_F_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
			}
		}
		FreeCharMtx( mseq1f );
		FreeCharMtx( mseq1r );
		FreeCharMtx( mseq2 );
		free( tmpseq );
	}
	else if( alg == 'm' )
	{

		if( dodp ) // nakuserukamo
		{
			mseq1f = AllocateCharMtx( 1, nlenmax+1);
			mseq1r = AllocateCharMtx( 1, nlenmax+1 );
			mseq2 = AllocateCharMtx( 1, nlenmax+1 );
		}
		else
		{
//			nthread = 0; // heiretsu keisan no kouritsu ha warui node
			spointt = AllocateIntMtx( njob, 0 ); 
			pointt = AllocateIntMtx( njob, nlenmax+1 );
			pointt_rev = AllocateIntMtx( njob, nlenmax+1 );
		}
	    tmpseq = AllocateCharVec( MAX( nlenmax, B ) +1 );
	    revseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		res = AllocateIntVec( njob );
		resr = AllocateIntVec( njob );
		resf = AllocateIntVec( njob );
		map = AllocateIntVec( njob );
		contrastorder = AllocateIntVec( njob );
		if( dorp == 'd' ) tsize = (int)pow( 4, 6 );
		else              tsize = (int)pow( 6, 6 ); // iranai

		maxl = 0;
		for( i=0; i<njob; i++ )
		{
			gappick0( tmpseq, seq[i] );
			nogaplen[i] = strlen( tmpseq );
			if( nogaplen[i] > maxl ) maxl = nogaplen[i];
		}

		reporterr( "Step 1/2\n" );


		if( !dodp )
		{
			if( nadd )
				iend = njob - nadd;
			else
				iend = 0; // keisan shinai
	
			for( i=0; i<iend; i++ )
			{
				gappick0( tmpseq, seq[i] );
				strcpy( seq[i], tmpseq );
				seq_grp_nuc( grpseq, tmpseq );
				makepointtable_nuc( pointt[i], grpseq );
				spointt[i] = pointt[i];
			}
	
			if( nadd )
				istart = njob - nadd;
			else
				istart = 0;
			for( i=istart; i<njob; i++ ) 
			{
				gappick0( tmpseq, seq[i] );
				strcpy( seq[i], tmpseq );
				sreverse( revseq, tmpseq );

				seq_grp_nuc( grpseq, tmpseq );
				makepointtable_nuc( pointt[i], grpseq );
//				makecompositiontable_p( table1, pointt[i] ); -> moto no basho ni modosu
				seq_grp_nuc( grpseq, revseq );
				makepointtable_nuc( pointt_rev[i], grpseq );
//				makecompositiontable_p( table1_rev, pointt_rev[i] ); -> moto no basho ni modosu
				spointt[i] = pointt[i];


//				reporterr( "pointt[i] = %p\n", pointt[i] );
//				reporterr( "pointt[i][0] = %p\n", pointt[i][0] );

			}	
		}


		if( contrastsort ) // sukoshi chuui
		{


			if( nadd )
			{
				iend = njob-nadd;
				for( i=0; i<iend; i++ ) contrastorder[i] = i;
				istart = njob-nadd;
				iend = nadd;
			}
			else
			{
				istart = 0;
				iend = njob;
			}

			if( dodp )
				makecontrastorder( contrastorder+istart, seq+istart, iend, istart );
			else
				makecontrastorder6mer( contrastorder+istart, pointt+istart, pointt_rev+istart, seq+istart, iend, istart );
		}
		else
		{
			for( i=0; i<njob; i++ ) contrastorder[i] = i;
		}


//		reporterr( "contrastorder = \n" );
//		for( i=0; i<njob; i++ )
//			reporterr( "%d ", contrastorder[i] );
//		reporterr( "\n" );



		if( nadd )
			iend = njob - nadd;
		else
			iend = 1;
		for( i=0; i<iend; i++ )
		{
			ic = contrastorder[i];
//			fprintf( stdout, "%d, SKIP\n", i );
			gappick0( tmpseq, seq[ic] );
			strcpy( seq[ic], tmpseq );
//			if( !nadd ) strcpy( seq[i], tmpseq ); // seq ha tsukawanaikara ii.

#if 0 // -> makecontrastorder() no mae ni idou
			if( !dodp )
			{
				seq_grp_nuc( grpseq, tmpseq );
				makepointtable_nuc( pointt[ic], grpseq );
				spointt[ic] = pointt[ic];
			}
#endif

			strcpy( tmpseq, name[ic] );
			strcpy( name[ic], "_F_" );
			strncpy( name[ic]+3, tmpseq+1, 10 );
			name[ic][13] = 0;
		}

		reporterr( "\n\nStep 2/2\n" );

		if( nadd )
			istart = njob - nadd;
		else
			istart = 1;
		for( i=istart; i<njob; i++ ) 
		{
//			fprintf( stderr, "\r %d / %d ", i, njob );
			ic = contrastorder[i];
			gappick0( tmpseq, seq[ic] );
			strcpy( seq[ic], tmpseq );
			sreverse( revseq, tmpseq );

#if 0 // -> makecontrastorder() no mae ni idou
			if( !dodp )
			{
				table1 = (short *)calloc( tsize, sizeof( short ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				table1_rev = (short *)calloc( tsize, sizeof( short ) );
				if( !table1_rev ) ErrorExit( "Cannot allocate table1_rev\n" );
				seq_grp_nuc( grpseq, tmpseq );
				makepointtable_nuc( pointt[ic], grpseq );
				makecompositiontable_p( table1, pointt[ic] );
				seq_grp_nuc( grpseq, revseq );
				makepointtable_nuc( pointt_rev[ic], grpseq );
				makecompositiontable_p( table1_rev, pointt_rev[ic] );
			}
#else
			if( !dodp )
			{
				table1 = (int *)calloc( tsize, sizeof( int ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				table1_rev = (int *)calloc( tsize, sizeof( int ) );
				if( !table1_rev ) ErrorExit( "Cannot allocate table1_rev\n" );
				makecompositiontable_p( table1, pointt[ic] );
				makecompositiontable_p( table1_rev, pointt_rev[ic] );
			}
#endif

			if( nadd && addfragment )
				iend = njob-nadd;
			else
				iend = i;


			if( iend > reflim ) 
			{
//				reporterr( "iend = %d -> %d\n", iend, reflim );
#if 0
				for( j=0; j<iend; j++ ) map[j] = j;
				partshuffle( iend, reflim, map );
#else
				for( j=0; j<iend; j++ ) map[j] = contrastorder[j];
#endif
				iend = reflim; // approximation
			}
			else
			{
#if 0
				for( j=0; j<iend; j++ ) map[j] = j;
#else
				for( j=0; j<iend; j++ ) map[j] = contrastorder[j];
#endif
			}

//			reporterr( "reflim = %d, seq[%d] = %s\n", reflim, contrastorder[0], seq[contrastorder[0]] );

#ifdef enablemultithread
			if( nthread )
			{
				pthread_t *handle;
				pthread_mutex_t mutex_counter;
				thread_arg_t *targ;
				int *jsharept;
		
				targ = calloc( nthread, sizeof( thread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex_counter, NULL );
				jsharept = calloc( 1, sizeof(int) );
				*jsharept = 0;
		
				if( i%100==1 ) fprintf( stderr, " %d / %d (%d threads)   \r", i, njob, nthread );
				for( j=0; j<nthread; j++ )
				{
					targ[j].iend = iend;
					targ[j].map = map;
					targ[j].seq = seq;
					targ[j].tmpseq = tmpseq; 
					targ[j].res = resf; 
					targ[j].spointt = spointt; 
					targ[j].table1 = table1; 
					targ[j].jshare = jsharept;
					targ[j].iq = i; // iranai
					targ[j].mutex_counter = &mutex_counter;
					targ[j].thread_no = j;
					pthread_create( handle+j, NULL, directionthread, (void *)(targ+j) );
				}
				for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
				pthread_mutex_destroy( &mutex_counter );
				free( handle );
				free( targ );
				free( jsharept );
			}
			else
#endif
			{
				thread_arg_t *targ;

				if( i%100==1 ) fprintf( stderr, " %d / %d   \r", i, njob );
				targ = calloc( 1, sizeof( thread_arg_t ) );
				targ[0].iend = iend;
				targ[0].map = map;
				targ[0].seq = seq;
				targ[0].tmpseq = tmpseq; 
				targ[0].res = resf; 
				targ[0].spointt = spointt; 
				targ[0].table1 = table1; 
				targ[0].iq = i;  // iranai
				directionthread( targ );
				free( targ );
			}



#ifdef enablemultithread
			if( nthread )
			{
				pthread_t *handle;
				pthread_mutex_t mutex_counter;
				thread_arg_t *targ;
				int *jsharept;
		
				targ = calloc( nthread, sizeof( thread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex_counter, NULL );
				jsharept = calloc( 1, sizeof(int) );
				*jsharept = 0;
		
				for( j=0; j<nthread; j++ )
				{
					targ[j].iend = iend;
					targ[j].seq = seq;
					targ[j].map = map;
					targ[j].tmpseq = revseq; 
					targ[j].res = resr; 
					targ[j].spointt = spointt; 
					targ[j].table1 = table1_rev; 
					targ[j].jshare = jsharept;
					targ[j].iq = i; // iranai
					targ[j].mutex_counter = &mutex_counter;
					targ[j].thread_no = j;
					pthread_create( handle+j, NULL, directionthread, (void *)(targ+j) );
				}
				for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
				pthread_mutex_destroy( &mutex_counter );
				free( handle );
				free( targ );
				free( jsharept );
			}
			else
#endif
			{
				thread_arg_t *targ;
				targ = calloc( 1, sizeof( thread_arg_t ) );
				targ[0].iend = iend;
				targ[0].seq = seq;
				targ[0].map = map;
				targ[0].tmpseq = revseq; 
				targ[0].res = resr; 
				targ[0].spointt = spointt;
				targ[0].table1 = table1_rev; 
				targ[0].iq = i;  // iranai
				directionthread( targ );
				free( targ );
			}

			if( mode == '2' )
			{
				mres = mres2 = 0;
				for( j=0; j<iend; j++ )
				{
					ires = resf[j];
//					fprintf( stdout, "ires (%d,%d) = %d\n", i, j, ires );
//					fflush( stdout );
					if( ires>mres2 ) 
					{
						if( ires>mres ) 
						{
							mres2 = mres;
							mres = ires;
						}
						else
							mres2 = ires;
					}
				}
				res_forward = (double)( mres + mres2 ) / 2;
				mres = mres2 = 0;
				for( j=0; j<iend; j++ )
				{
					ires = resr[j];
					if( ires>mres2 )
					{
						if( ires>mres ) 
						{
							mres2 = mres;
							mres = ires;
						}
						else
							mres2 = ires;
					}
				}
				res_reverse = (double)( mres + mres2 ) / 2;
				res_max = MAX(res_reverse,res_forward);
			}
//			reporterr( "i=%d, res_reverse = %f\n", i, res_reverse );
			else if( mode == '1' )
			{
				res_reverse = 0.0;
				for( j=0; j<iend; j++ ) if( res_reverse < (double)resr[j] ) res_reverse = (double)resr[j];
				res_forward = 0.0;
				for( j=0; j<iend; j++ ) if( res_forward < (double)resf[j] ) res_forward = (double)resf[j];
				res_max = 1.0;
			}

			else if( mode == 'd' )
			{
				res_reverse = 0.0;
				for( j=0; j<iend; j++ ) if( res_reverse < (double)(resr[j]-resf[j]) ) res_reverse = (double)(resr[j]-resf[j]);
				res_forward = 0.0;
				for( j=0; j<iend; j++ ) if( res_forward < (double)(resf[j]-resr[j]) ) res_forward = (double)(resf[j]-resr[j]);
				res_max = 1.0;
			}

			else if( mode == 'a' )
			{
				res_reverse = 0.0;
				for( j=0; j<iend; j++ ) res_reverse += (double)resr[j];
				res_reverse /= (double)iend;
				res_forward = 0.0;
				for( j=0; j<iend; j++ ) res_forward += (double)resf[j];
				res_forward /= (double)iend;
				res_max = 1.0;
			}
			else
			{
				reporterr( "Unknown mode!\n" );
				exit( 1 );
			}


			if( (res_reverse>res_forward) ) // tekitou
//			if( (res_reverse-res_forward)/res_max > thresholdtorev ) // tekitou
			{
				strcpy( seq[ic], revseq );

				strcpy( tmpseq, name[ic] );
				strcpy( name[ic], "_R_" );
				strncpy( name[ic]+3, tmpseq+1, 10 );
				name[ic][13] = 0;
				if( !dodp ) spointt[ic] = pointt_rev[ic];
			}
			else
			{
				strcpy( tmpseq, name[ic] );
				strcpy( name[ic], "_F_" );
				strncpy( name[ic]+3, tmpseq+1, 10 );
				name[ic][13] = 0;
				if( !dodp ) spointt[ic] = pointt[ic];
			}

			if( !dodp )
			{
				free( table1 );
				free( table1_rev );
			}
		}

		if( name[0][1] == 'R' )
		{
			for( j=0; j<njob; j++ ) 
			{
				if( name[j][1] == 'R' ) 
					name[j][1] = 'F';
				else
					name[j][1] = 'R';
			}
		}

		creverse( 0 );
		free( tmpseq );
		free( revseq );
		free( grpseq );
		free( res );
		free( resr );
		free( resf );
		free( map );
		free( nlen );
		free( nogaplen );
		free( contrastorder );
		if( dodp )
		{
			FreeCharMtx( mseq1f );
			FreeCharMtx( mseq1r );
			FreeCharMtx( mseq2 );
		}
		else
		{
			FreeIntMtx( pointt );
			FreeIntMtx( pointt_rev );
			free( spointt );
		}
	}
	else
	{
		fprintf( stderr, "Unknown alg %c\n", alg );
		exit( 1 );
	}
//	writeData_pointer( stdout, njob, name, nlen, seq );
	for( i=0; i<njob; i++ ) 
	{
//		fprintf( stdout, ">%s\n", name[i] );
//		fprintf( stdout, "%s\n", seq[i] );
		fprintf( stdout, "%s\n", name[i] );
	}

	FreeCharMtx( seq );
	FreeCharMtx( name );
	freeconstants();
	closeFiles();

	fprintf( stderr, "\n" );
	SHOWVERSION;
	return( 0 );
}

