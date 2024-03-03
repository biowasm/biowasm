#include "mltaln.h"

#define DEBUG 0

char *alignmentfile;
int keeplength;
int mapout;

static void fillorichar( int nseq, int *oripos, char **a, char **s )
{
	int i;
	char *pta, *pts;
	for( i=0; i<nseq; i++ )
	{

		pta = a[i];
		pts = s[oripos[i]];
		while( *pta )
		{
			if( *pta != '-' ) *pta = *pts++;
			if( *pta++ == 0 )
			{
				fprintf( stderr, "ERROR!!\n" );
				fprintf( stderr, "alignment is inconsistent with the original sequences (1)\n" );
				exit( 1 );
			}
		}
		if( *pts != 0 )
		{
			fprintf( stderr, "ERROR!!\n" );
			fprintf( stderr, "alignment is inconsistent with the original sequences (2)\n" );
			exit( 1 );
		}
	}
}

void arguments( int argc, char *argv[] )
{
    int c;

	keeplength = 0;
	mapout = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					--argc;
					goto nextoption;
				case 'a':
					alignmentfile = *++argv;
					--argc;
					goto nextoption;
				case 'z': // add2ndharfarg wo tsukau tame.
					mapout = 2;
					break;
				case 'Z': // add2ndharfarg wo tsukau tame.
					mapout = 1;
					break;
				case 'p': // add2ndharfarg wo tsukau tame.
					break;
				case 'Y':
					keeplength = 1;
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
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
}



int main( int argc, char *argv[] )
{
	FILE *infp;
	FILE *alfp;
	FILE *dlfp;
	char **name;
	char **aname;
	char **oname;
	char **seq;
	char **aseq;
	int *nlen;
	int *oripos;
	char *npt, *npt0, *npt2, *pt, *pt2;
	int i, o, l, prelen;
	int nlenmin;
	int njobs, njoba;
//	int **dlist;
//	int *ndel;
	char *gett;
	char *insname;

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

	if( alignmentfile )
	{
		alfp = fopen( alignmentfile, "r" );
		if( !alfp )
		{
			fprintf( stderr, "Cannot open %s\n", alignmentfile );
			exit( 1 );
		}
	}
	else
	{
		fprintf( stderr, "No alignment is given.\n" );
		exit( 1 );
	}

	if( keeplength )
	{
		dlfp = fopen( "_deletelist", "r" );
		if( !dlfp )
		{
			fprintf( stderr, "Cannot open _deletefile\n" );
			exit( 1 );
		}
	}
	
	dorp = NOTSPECIFIED;
	getnumlen_casepreserve( infp, &nlenmin );
	njobs = njob;
//	fprintf( stderr, "in infp, %d x %d - %d %c\n", njob, nlenmin, nlenmax, dorp );

	seq = AllocateCharMtx( njob, nlenmax+1 );
	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	oripos = AllocateIntVec( njob );
	readData_pointer_casepreserve( infp, name, nlen, seq );
	fclose( infp );

	dorp = NOTSPECIFIED;
	getnumlen( alfp );
	njoba = njob;
//	fprintf( stderr, "in alfp, %d x %d %c\n", njob, nlenmax, dorp );
	aseq = AllocateCharMtx( njob, nlenmax+1 );
	aname = AllocateCharMtx( njob, B+1 );
	oname = AllocateCharMtx( njob, B+1 );
	insname = calloc( njob, sizeof( char ) );
	readData_pointer( alfp, aname, nlen, aseq );
	fclose( alfp );

	for( i=0; i<njob; i++ ) gappick_samestring( seq[i] );

//	reporterr( "seq[njob-1] = %s\n", seq[njob-1] );

	if( keeplength )
	{
		gett = calloc( 1000, sizeof( char ) );
//		ndel = (int *)calloc( njob, sizeof( int ) );
//		for( i=0; i<njob; i++ ) ndel[i] = 0;
//		dlist = (int **)calloc( njob+1, sizeof( int * ) );
//		for( i=0; i<njob; i++ ) dlist[i] = NULL;
//		dlist[njob] = NULL;

		while( 1 )
		{
			fgets( gett, 999, dlfp );
			if( feof( dlfp ) ) break;
//			sscanf( gett, "%d %d", &i, &o );
			sscanf( gett, "%d %d %d", &i, &o, &l );
//			reporterr( "%d, %d\n", i, o );
//			dlist[i] = realloc( dlist[i], sizeof( int ) * (ndel[i]+2) );
//			dlist[i][ndel[i]] = o;
//			ndel[i]++;

			while( l-- ) seq[i][o++] = '-';
//			seq[i][o] = '-';
			insname[i] = 1;

		}
		fclose( dlfp );
		free( gett );
	}

	for( i=0; i<njob; i++ )
	{
		gappick_samestring( seq[i] );
	}

	if( njoba != njobs )
	{
		fprintf( stderr, "ERROR!!\n" );
		fprintf( stderr, "In input file,\n" );
		fprintf( stderr, "njob = %d\n", njobs );
		fprintf( stderr, "but in alignment file,\n" );
		fprintf( stderr, "njob = %d\n", njoba );
		exit( 1 );
	}

	for( i=0; i<njob; i++ )
	{
#if 0
		if( strstr( aname[i], "_seed_" ) )
		{
			npt2 = aname[i] + 7;
			strcpy( oname[i], "=_seed_" );
		}
		else
		{
			npt2 = aname[i] + 1;
			strcpy( oname[i], "=" );
		}

		fprintf( stderr, "npt2 = %s\n", npt2 );

		o = oripos[i] = atoi( npt2 );
		npt = strstr( npt2, "_oe_" );
		if( npt == NULL )
		{
			fprintf( stderr, "Format error!\n" );
			exit( 1 );
		}
		npt += 4;
		strcat( oname[i], npt+1 );
#endif
		npt0 = strstr( aname[i], "_os_" );
		if( npt0 == NULL )
		{
			fprintf( stderr, "Format error!\n" );
			exit( 1 );
		}
		npt2 = npt0 + 4;
		o = oripos[i] = atoi( npt2 );

		npt = strstr( aname[i], "_oe_" );
		if( npt == NULL )
		{
			fprintf( stderr, "Format error!\n" );
			exit( 1 );
		}
		npt += 4;
		

		pt2 = npt;
		pt = npt2 - 4;
		while( *pt ) *pt++ = *pt2++; // okashii

		prelen = npt0-aname[i];
		strncpy( oname[i], aname[i], prelen ); oname[i][prelen] = 0;
		strcat( oname[i], npt0 );

		if( mapout == 2 && insname[i] )
		{
			npt0 = strstr( npt0, "|" );
			npt0 += 1;
		}
#if 0
		pt = strstr( aname[i], "_numo_e" );
		if( pt ) pt += 8;
		else pt = aname[i] + 1;

		if( strstr( pt, "_seed_" ) ) pt += 6;

		fprintf( stderr, "pt = :%s:\n", pt );

#endif
//		fprintf( stderr, "npt0 = :%s:\n", npt0 );
//
//		reporterr( "oname[i] = %s\n", oname[i] );
//		reporterr( "aname[i] = %s\n", aname[i] );
//		reporterr( " name[i] = %s\n",  name[i] );
//
//		fprintf( stderr, "aname[%d] = :%s:\n", i, aname[i] );
//		fprintf( stderr, "npt0 = :%s:\n", npt0 );
//		fprintf( stderr, "oname[%d] = :%s:\n", i, oname[i] );
//		fprintf( stderr, "name[%d] = :%s:\n", o, name[o] );

		if( strncmp( npt0, name[o]+1, 10 ) )
		{
			fprintf( stderr, "ERROR!!\n" );
			fprintf( stderr, "In input file,\n" );
			fprintf( stderr, "name[%d] = %s\n", o, name[o] );
			fprintf( stderr, "but in alignment file,\n" );
			fprintf( stderr, "oname[%d] = %s\n", i, oname[i] );
			fprintf( stderr, "npt0 = %s\n", npt0 );
			fprintf( stderr, "prelen = %d\n", prelen );
			fprintf( stderr, "name[%d] = %s\n", i, aname[i] );
			exit( 1 );
		}
#if 0
		else
		{
			fprintf( stderr, "OK!!\n" );
			fprintf( stderr, "In input file,\n" );
			fprintf( stderr, "name[%d] = %s\n", o, name[o] );
			fprintf( stderr, "and in alignment file,\n" );
			fprintf( stderr, "name[%d] = %s\n", i, aname[i] );
			fprintf( stderr, "\n" );
		}
#endif
	}
//	fprintf( stderr, "seq[0] = %s\n", seq[0] );
//	fprintf( stderr, "aseq[0] = %s\n", aseq[0] );

	fillorichar( njob, oripos, aseq, seq );


	writeData_pointer( stdout, njob, oname, nlen, aseq );

	FreeCharMtx( seq );
	FreeCharMtx( aseq );
	FreeCharMtx( name );
	FreeCharMtx( aname );
	FreeCharMtx( oname );
	free( nlen );
	free( oripos );
	free( insname );

	return( 0 );
}
