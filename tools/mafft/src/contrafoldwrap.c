#include "mltaln.h"

#define DEBUG 0

static char *whereiscontrafold;

void unknown_n( char *out, char *in )
{
	while( *in )
	{
		if( *in == 'a' || *in == 'A' )
			*out = 'A';
		else if( *in == 't' || *in == 'T' || *in == 'u' || *in == 'U' )
			*out = 'U';
		else if( *in == 'g' || *in == 'G' )
			*out = 'G';
		else if( *in == 'c' || *in == 'C' )
			*out = 'C';
		else if( *in == '-' )
			*out = '-';
		else
			*out = 'N';

		out++;
		in++;
	}
	*out = 0;
}

void outcontrafold( FILE *fp, RNApair **pairprob, int length )
{
	int i;
	RNApair *pt;
	for( i=0; i<length; i++ ) for( pt=pairprob[i]; pt->bestpos!=-1; pt++ )
	{
		if( pt->bestpos > i ) 
			fprintf( fp, "%d %d %f\n", i, pt->bestpos, pt->bestscore );
	}
}

#if 1
static void readcontrafold( FILE *fp, RNApair **pairprob, int length )
{
	char gett[10000];
	int *pairnum;
	char *pt;
	int i;
	int left, right;
	double prob;

	pairnum = (int *)calloc( length, sizeof( int ) );
	for( i=0; i<length; i++ ) pairnum[i] = 0;

	while( 1 )
	{
		if( feof( fp ) ) break;
		fgets( gett, 9999, fp );

//		fprintf( stderr, "gett=%s\n", gett );

		pt = gett;

		sscanf( gett, "%d ", &left );
		left--;

//		fprintf( stderr, "left=%d\n", left );
		pt = strchr( pt, ' ' ) + 1;
//		fprintf( stderr, "pt=%s\n", pt );

		while( (pt = strchr( pt, ' ' ) ) )
		{
			pt++;
//			fprintf( stderr, "pt=%s\n", pt );
			sscanf( pt, "%d:%lf", &right, &prob );
			right--;

//			fprintf( stderr, "%d-%d, %f\n", left, right, prob );

			pairprob[left] = (RNApair *)realloc( pairprob[left], (pairnum[left]+2) * sizeof( RNApair ) );
			pairprob[left][pairnum[left]].bestscore = prob;
			pairprob[left][pairnum[left]].bestpos = right;
			pairnum[left]++;
			pairprob[left][pairnum[left]].bestscore = -1.0;
			pairprob[left][pairnum[left]].bestpos = -1;
//			fprintf( stderr, "%d-%d, %f\n", left, right, prob );

			pairprob[right] = (RNApair *)realloc( pairprob[right], (pairnum[right]+2) * sizeof( RNApair ) );
			pairprob[right][pairnum[right]].bestscore = prob;
			pairprob[right][pairnum[right]].bestpos = left;
			pairnum[right]++;
			pairprob[right][pairnum[right]].bestscore = -1.0;
			pairprob[right][pairnum[right]].bestpos = -1;
//			fprintf( stderr, "%d-%d, %f\n", right, left, prob );
		}
	}
	free( pairnum );
}
#endif

void arguments( int argc, char *argv[] )
{
    int c;
	inputfile = NULL;
	dorp = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	whereiscontrafold = NULL;

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
				case 'd':
					whereiscontrafold = *++argv;
					fprintf( stderr, "whereiscontrafold = %s\n", whereiscontrafold );
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
	static char com[10000];
	static int  *nlen;	
	int left, right;
	int res;
	static char **name, **seq, **nogap;
	static int **gapmap;
	static int *order;
	int i, j;
	FILE *infp;
	RNApair ***pairprob;
	RNApair **alnpairprob;
	RNApair *pairprobpt;
	RNApair *pt;
	int *alnpairnum;
	double prob;
	int adpos;

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

	if( !whereiscontrafold )
		whereiscontrafold = "";

	getnumlen( infp );
	rewind( infp );

	if( dorp != 'd' )
	{
		fprintf( stderr, "nuc only\n" );
		exit( 1 );
	}

	seq = AllocateCharMtx( njob, nlenmax*2+1 );
	nogap = AllocateCharMtx( njob, nlenmax*2+1 );
	gapmap = AllocateIntMtx( njob, nlenmax*2+1 );
	order = AllocateIntVec( njob );
	name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob );
	pairprob = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	alnpairprob = (RNApair **)calloc( nlenmax, sizeof( RNApair * ) );
	alnpairnum = AllocateIntVec( nlenmax );

	for( i=0; i<nlenmax; i++ ) alnpairnum[i] = 0;

	readData_pointer( infp, name, nlen, seq );
	fclose( infp );

	for( i=0; i<njob; i++ )
	{
		pairprob[i] = (RNApair **)calloc( nlenmax, sizeof( RNApair * ) );
		for( j=0; j<nlenmax; j++ )
		{
			pairprob[i][j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
			pairprob[i][j][0].bestpos = -1;
			pairprob[i][j][0].bestscore = -1.0;
		}
		unknown_n( nogap[i], seq[i] );
		order[i] = i;
	}
	for( j=0; j<nlenmax; j++ )
	{
		alnpairprob[j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
		alnpairprob[j][0].bestpos = -1;
		alnpairprob[j][0].bestscore = -1.0;
	}


	constants( njob, seq );

	fprintf( stderr, "running contrafold\n" );
	for( i=0; i<njob; i++ )
	{
		fprintf( stderr, "%d / %d\n", i+1, njob );
		commongappick_record( 1, nogap+i, gapmap[i] );
		infp = fopen( "_contrafoldin", "w" );
		fprintf( infp, ">in\n%s\n", nogap[i] );
		fclose( infp );
#if 0 // contrafold v1
		sprintf( com, "env PATH=%s contrafold predict _contrafoldin --posteriors 0.01 > _contrafoldout", whereiscontrafold );
#else // contrafold v2
		sprintf( com, "env PATH=%s contrafold predict _contrafoldin --posteriors 0.01   _contrafoldout", whereiscontrafold );
#endif
		res = system( com );
		if( res )
		{
			fprintf( stderr, "error in contrafold\n" );
			fprintf( stderr, "=================================================================\n" );
			fprintf( stderr, "=================================================================\n" );
			fprintf( stderr, "==\n" );
			fprintf( stderr, "== This version of MAFFT supports CONTRAfold v2.02.\n" );
			fprintf( stderr, "== If you have a lower version of CONTRAfold installed in the\n" );
			fprintf( stderr, "== %s directory,\n", whereiscontrafold );
			fprintf( stderr, "== please update it!\n" );
			fprintf( stderr, "==\n" );
			fprintf( stderr, "=================================================================\n" );
			fprintf( stderr, "=================================================================\n" );
			exit( 1 );
		}


		infp = fopen( "_contrafoldout", "r" );
		readcontrafold( infp, pairprob[i], nlenmax );
		fclose( infp );
		fprintf( stdout, ">%d\n", i );
		outcontrafold( stdout, pairprob[i], nlenmax );
	}

	for( i=0; i<njob; i++ )
	{
		for( j=0; j<nlen[i]; j++ ) for( pairprobpt=pairprob[i][j]; pairprobpt->bestpos!=-1; pairprobpt++ )
		{
			left = gapmap[i][j];
			right = gapmap[i][pairprobpt->bestpos];
			prob = pairprobpt->bestscore;

			for( pt=alnpairprob[left]; pt->bestpos!=-1; pt++ )
				if( pt->bestpos == right ) break;

			if( pt->bestpos == -1 )
			{
				alnpairprob[left] = (RNApair *)realloc( alnpairprob[left], (alnpairnum[left]+2) * sizeof( RNApair ) );
				adpos = alnpairnum[left];
				alnpairnum[left]++;
				alnpairprob[left][adpos].bestscore = 0.0;
				alnpairprob[left][adpos].bestpos = right;
				alnpairprob[left][adpos+1].bestscore = -1.0;
				alnpairprob[left][adpos+1].bestpos = -1;
				pt = alnpairprob[left]+adpos;
			}
			else
				adpos = pt-alnpairprob[left];

			pt->bestscore += prob;
			if( pt->bestpos != right )
			{
				fprintf( stderr, "okashii!\n" );
				exit( 1 );
			}
//			fprintf( stderr, "adding %d-%d, %f\n", left, right, prob );
		}
	}
	return( 0 );

#if 0
	fprintf( stdout, "result=\n" );

	for( i=0; i<nlenmax; i++ ) for( pairprobpt=alnpairprob[i]; pairprobpt->bestpos!=-1; pairprobpt++ )
	{
		pairprobpt->bestscore /= (double)njob;
		left = i;
		right = pairprobpt->bestpos;
		prob = pairprobpt->bestscore;
		fprintf( stdout, "%d-%d, %f\n", left, right, prob );
	}

	return( 0 );
#endif
}
