#include "mltaln.h"

#define DEBUG 0


static char *comment;
static char *orderfile;
static int format;
static int namelen;
static int excludedashseq;
static int extendedalphabet;

static void fillspace( char *seq, int lenmax )
{
	int len = strlen( seq );
	seq += len;
	lenmax -= len;
	while( lenmax-- ) *seq++ = ' ';
	*seq = 0;
}

void setmark_clustal( int nlen, int nseq, char **seq, char *mark )
{
	int i, j, k, nalpha;
	char firstletter;
	char *strong[9];
	char *weaker[11];
	int nstrong, nweaker;
	char s;

	if( dorp == 'd' ) 
	{
		strong[0] = "TU";
		nstrong = 1;
		weaker[0] = "AG";
		weaker[1] = "CT";
		weaker[2] = "CU";
		nweaker = 2;
		nalpha = 10;
	}
	else
	{
		strong[0] = "STA";
		strong[1] = "NEQK";
		strong[2] = "NHQK";
		strong[3] = "NDEQ";
		strong[4] = "QHRK";
		strong[5] = "MILV";
		strong[6] = "MILF";
		strong[7] = "HY";
		strong[8] = "FYW";
		nstrong = 9;
		weaker[0] = "CSA";
		weaker[1] = "ATV";
		weaker[2] = "SAG";
		weaker[3] = "STNK";
		weaker[4] = "STPA";
		weaker[5] = "SGND";
		weaker[6] = "SNDEQK";
		weaker[7] = "NDEQHK";
		weaker[8] = "NEQHRK";
		weaker[9] = "FVLIM";
		weaker[10] = "HFY";
		nweaker = 11;
		nalpha = 20;
	}

	for( i=0; i<nlen; i++ )
	{
		mark[i] = ' ';
		for( j=0; j<nseq; j++ )
		{
			s = seq[j][i];
			if( '-' == s || ' ' == s ) break;
		}
		if( j != nseq ) 
		{
			continue;
		}
		if( extendedalphabet )
		{
			firstletter = seq[0][i];
			if( amino_n[(unsigned char)firstletter] < 0 ) continue;
	
			for( j=0; j<nseq; j++ )
				if( seq[j][i] != firstletter ) break;
			if( j == nseq ) 
			{
				mark[i] = '*';
				continue;
			}
		}
		else 
		{
			firstletter = toupper( seq[0][i] );
			if( amino_n[(unsigned char)firstletter] >= nalpha || amino_n[(unsigned char)firstletter] < 0 ) continue;
	
			for( j=0; j<nseq; j++ )
				if( toupper( seq[j][i] ) != firstletter ) break;
			if( j == nseq ) 
			{
				mark[i] = '*';
				continue;
			}
			for( k=0; k<nstrong; k++ )
			{
				for( j=0; j<nseq; j++ )
				{
					if( !strchr( strong[k], toupper( seq[j][i] ) ) ) break;
				}
				if( j == nseq ) break;
			}
			if( k < nstrong )
			{
				mark[i] = ':';
				continue;
			}
			for( k=0; k<nweaker; k++ )
			{
				for( j=0; j<nseq; j++ )
				{
					if( !strchr( weaker[k], toupper( seq[j][i] ) ) ) break;
				}
				if( j == nseq ) break;
			}
			if( k < nweaker )
			{
				mark[i] = '.';
				continue;
			}
		}
	}
	mark[nlen] = 0;
}

void setmark( int nlen, int nseq, char **seq, char *mark )
{
	int i, j;

	for( i=0; i<nlen; i++ )
	{
		mark[i] = ' ';
		for( j=0; j<nseq; j++ )
			if( '-' == seq[j][i] ) break;
		if( j != nseq ) 
		{
			continue;
		}
		for( j=0; j<nseq; j++ )
			if( seq[0][i] != seq[j][i] ) break;
		if( j == nseq ) 
		{
			mark[i] = '*';
			continue;
		}
		for( j=0; j<nseq; j++ )
			if( amino_grp[(unsigned char)seq[0][i]] != amino_grp[(unsigned char)seq[j][i]] ) break;
		if( j == nseq ) 
		{
			mark[i] = '.';
			continue;
		}
	}
	mark[nlen] = 0;
}

void arguments( int argc, char *argv[] )
{
    int c;
	excludedashseq = 0;
	namelen = -1;
	scoremtx = 1;
	nblosum = 62;
	LineLengthInFASTA = -1;
	dorp = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	inputfile = NULL;
	comment = NULL;
	orderfile = NULL;
	format = 'c';
	extendedalphabet = 0;

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
				case 'c':
					comment = *++argv;
					fprintf( stderr, "comment = %s\n", comment );
					--argc;
					goto nextoption;
				case 'r':
					orderfile = *++argv;
					fprintf( stderr, "orderfile = %s\n", orderfile );
					--argc;
					goto nextoption;
				case 'n':
					namelen = myatoi( *++argv );
					fprintf( stderr, "namelen = %d\n", namelen );
					--argc;
					goto nextoption;
				case 'l':
					LineLengthInFASTA = myatoi( *++argv );
					if( LineLengthInFASTA == 0 )
					{
						reporterr( "LineLengthInFASTA = 0 ??\n" );
						exit( 1 );
					}
					--argc;
					goto nextoption;
				case 'f':
					format = 'f';
					break;
				case 'd':
					excludedashseq = 1;
					break;
				case 'y':
					format = 'y';
					break;
				case 'E':
					extendedalphabet = 1;
					nblosum = -2;
					break;
				case 'N':
					extendedalphabet = 0;
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
	static int  *nlen, *onlen;	
	static char **name, **oname, **seq, **oseq, *mark;
	static int *order, *oorder;
	int i, j;
	FILE *infp;
	FILE *orderfp;
	char gett[B];
	int nlenmin;
	int nout;

	arguments( argc, argv );


	if( inputfile )
	{
		infp = fopen( inputfile, "rb" );
		if( !infp )
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else
		infp = stdin;

	getnumlen_casepreserve( infp, &nlenmin );
	rewind( infp );

	seq = AllocateCharMtx( njob, nlenmax*2+1 );
	mark = AllocateCharVec( nlenmax*2+1 );
	order = AllocateIntVec( njob );
	name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob );

	if( orderfile )
	{
		orderfp = fopen( orderfile, "r" );
		if( !orderfp )
		{
			fprintf( stderr, "Cannot open %s\n", orderfile );
			exit( 1 );
		}
		for( i=0; i<njob; i++ )
		{
			fgets( gett, B-1, orderfp );
			order[i] = atoi( gett );
		}
		fclose( orderfp );
	}
	else
	{
		for( i=0; i<njob; i++ ) order[i] = i;
	}


	readData_pointer_casepreserve( infp, name, nlen, seq );
	fclose( infp );

	if( format == 'c' || format == 'y' ) for( i=0; i<njob; i++ ) fillspace( seq[i], nlenmax );
	constants( njob, seq );

//	initSignalSM();

//	initFiles();

	if( excludedashseq )
	{
		nout = 0;
		for( i=0; i<njob; i++ )
		{
			if( !strncmp( name[order[i]]+1, "DASH|", 5 ) || strstr( name[order[i]]+1, "_numo_e_DASH|" ) ) continue;
			nout++;
		}

		reporterr( "nout=%d\n", nout );

		oseq = calloc( sizeof( char * ), nout );
		oname = calloc( sizeof( char * ), nout );
   		onlen = AllocateIntVec( nout );
   		oorder = AllocateIntVec( nout );


		for( i=0,j=0; i<njob; i++ )
		{
			if( !strncmp( name[order[i]]+1, "DASH|", 5 ) || strstr( name[order[i]]+1, "_numo_e_DASH|" ) ) continue;
			oname[j] = name[order[i]];
			oseq[j] = seq[order[i]];
			onlen[j] = nlen[order[i]];
   			oorder[j] = j;
			j++;
		}
		commongappick( nout, oseq );
	}
	else
	{
		oseq = seq;
		oname = name;
   		onlen = nlen;
   		oorder = order;
		nout = njob;
	}

	reporterr( "seq = %p, oseq=%p\n", seq, oseq );



//	setmark( nlenmax, njob, seq, mark );
	setmark_clustal( nlenmax, nout, oseq, mark );

#if mingw
	setmode( fileno( stdout ), O_TEXT ); // windows deha saishuu tekina output nomi text mode
#endif

	if( format == 'f' )
		writeData_reorder_pointer( stdout, nout, oname, onlen, oseq, oorder );
	else if( format == 'c' )
		clustalout_pointer( stdout, nout, nlenmax, oseq, oname, mark, comment, oorder, namelen );
	else if( format == 'y' )
		phylipout_pointer( stdout, nout, nlenmax, oseq, oname, oorder, namelen );
	else
		fprintf( stderr, "Unknown format\n" );

	FreeCharMtx( seq ); seq = NULL;
	FreeCharMtx( name ); name = NULL;
	free( nlen ); nlen = NULL;
	free( order ); order = NULL;

	if( excludedashseq )
	{
		free( oseq ); oseq = NULL;
		free( oname ); oname = NULL;
		free( onlen ); onlen = NULL;
		free( oorder ); oorder = NULL;
	}

	free( mark ); mark = NULL;
	freeconstants();

//	SHOWVERSION;
	return( 0 );
}
