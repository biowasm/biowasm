#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

static int usecache;
static char *whereispairalign;
static char *odir;
static char *pdir;
static double scale;
static int *alreadyoutput;
static int equivthreshold;
static int equivwinsize;
static int equivshortestlen;

static void cutpath( char *s )
{
	char *pos;
	pos = s + strlen( s );

	while( --pos >= s )
	{
		if( *pos == '/' ) break;
	}

	strcpy( s, pos+1 );
}

static char getchainid( char *s )
{
	s += strlen( s ) - 2;
	if( isspace( s[0] ) && isalnum( s[1] ) )
		return( s[1] );
	else
		return( 'A' );
}

static void extractfirstword( char *s )
{
	while( *s )
	{
		if( isspace( *s ) ) break;
		s++;
	}
	*s = 0;
}

static char *strip( char *s )
{
	char *v;

	while( *s )
	{
		if( !isspace( *s ) ) break;
		s++;
	}
	v = s;

	s += strlen( v ) - 1;
	while( s>=v )
	{
		if( !isspace( *s ) ) 
		{
			*(s+1) = 0;
			break;
		}
		s--;
	}

	return( v );
}

#if 0
static void makeequivdouble( double *d, char *c )
{
	while( *c )
	{
		*d++ = (double)( *c++ - '0' );
	}
}

static void maskequiv( double *d, int n )
{
	int halfwin;
	int ok;
	int i, j;

	halfwin = (int)( equivwinsize / 2 );

	for( i=0; i<n; i++ )
	{
		ok = 1;
		for( j = i-halfwin; j<i+halfwin; j++ )
		{
			if( j<0 || n=<j ) continue;
			if( d[j] <= 0.0 )
			{
				ok = 0;
				break;
			}
		}
		if( ok == 0 ) d[i] = 0.0;
	}
}
#else
static void maskequiv( double *d, int n )
{
	int i, len;
	int count;
	len = 0;
	double *dbk, *dori, *dbkori;

	dbk = calloc( n, sizeof( double ) );

	dbkori = dbk;
	dori = d;
	count = n;
	while( count-- )
	{
		*dbk++ = *d++;
	}

	dbk = dbkori;
	d = dori;
	len = 0;


	for( i=0; i<n; i++ )
	{
		if( d[i] > 0.0 )
		{
			len += 1;
			d[i] = 0.0;
		}
		else
		{
			d[i] = 0.0;
			if( len >= equivshortestlen ) 
			{
				len++;
				while( len-- ) d[i-len] = dbk[i-len];
			}
			len = 0;
		}
	}

	if( len >= equivshortestlen )
	{
		len++;
		while( len-- ) d[n-len] = dbk[n-len];
	}

	free( dbk );
}
#endif

static void makeequivdouble_tmalign( double *d, char *c, int n )
{
	double tmpd;
	double *dbk;
	int tmpi;
	char s;
	dbk = d;
	while( *c )
	{
		if( ( s=*c++ ) == ':' )
			tmpi = 9;
		else if( s == '.' )
			tmpi = 4;
		else
			tmpi = 0;
//		tmpd = (double)( tmpi + 1 - equivthreshold ) / ( 10 - equivthreshold ) * 9.0;
//		if( tmpd < 0.0 ) tmpd = 0.0;
		tmpd = (double)( tmpi );
//		*d++ = (int)tmpd;
		*d++ = tmpd;
	}

	d = dbk;
//	maskequiv( d, n );
}

static void makeequivdouble_threshold( double *d, char *c, int n )
{
	double tmpd;
	double *dbk;
	int tmpi;
	dbk = d;
	while( *c )
	{
		tmpi = (int)( *c++ - '0' );
		tmpd = (double)( tmpi + 1 - equivthreshold ) / ( 10 - equivthreshold ) * 9.0;
		if( tmpd < 0.0 ) tmpd = 0.0;
//		*d++ = (int)tmpd;
		*d++ = tmpd;
	}

	d = dbk;
	maskequiv( d, n );
}

static void readtmalign( FILE *fp, char *seq1, char *seq2, double *equiv )
{
	static char *line = NULL;
	static char *equivchar = NULL;
	int n;

	
	if( equivchar == NULL )
	{
		equivchar = calloc( nlenmax * 2 + 1, sizeof( char ) );
		line = calloc( nlenmax * 2 + 1, sizeof( char ) );
	}
	seq1[0] = 0;
	seq2[0] = 0;
	equivchar[0] = 0;


//	system( "vi _tmalignout" );
	while( 1 )
	{
		if( feof( fp ) ) 
		{
			fprintf( stderr, "Error in TMalign\n" );
			exit( 1 );
		}
		fgets( line, 999, fp );
//		fprintf( stdout, "line = :%s:\n", line );
		if( !strncmp( line+5, "denotes the residue pairs", 20 ) ) break;
	}
	fgets( line, nlenmax*2, fp );
	strcat( seq1, strip( line ) );

	fgets( line, nlenmax*2, fp );
	strcat( equivchar, strip( line ) );

	fgets( line, nlenmax*2, fp );
	strcat( seq2, strip( line ) );

#if 0
	printf( "seq1=%s\n", seq1 );
	printf( "seq2=%s\n", seq2 );
	printf( "equi=%s\n", equivchar );
exit( 1 );
#endif
	n = strlen( seq1 );
	makeequivdouble_tmalign( equiv, equivchar, n );

#if 0
	fprintf( stdout, "\n" );
	for( i=0; i<n; i++ )
	{
		fprintf( stdout, "%1.0f", equiv[i] );
	}
	fprintf( stdout, "\n" );
	exit( 1 );
#endif
}
static void readrash( FILE *fp, char *seq1, char *seq2, double *equiv )
{
	char line[1000];
	static char *equivchar = NULL;
	int n;

	
	if( equivchar == NULL )
	{
		equivchar = calloc( nlenmax * 2, sizeof( char ) );
	}
	seq1[0] = 0;
	seq2[0] = 0;
	equivchar[0] = 0;

	while( 1 )
	{
		fgets( line, 999, fp );
//		fprintf( stdout, "line = :%s:\n", line );
		if( !strncmp( line, "Query ", 6 ) ) break;
		if( feof( fp ) ) break;
		if( !strncmp( line, "QUERY ", 6 ) )
		{
			strcat( seq1, strip( line+5 ) );
			fgets( line, 999, fp );
		}
		if( !strncmp( line, "TEMPL ", 6 ) )
		{
			strcat( seq2, strip( line+5 ) );
			fgets( line, 999, fp );
		}
		if( !strncmp( line, "Equiva", 6 ) )
		{
			strcat( equivchar, strip( line+11 ) );
			fgets( line, 999, fp );
		}
	}
#if 0
	printf( "seq1=:%s:\n", seq1 );
	printf( "seq2=:%s:\n", seq2 );
	printf( "equi=:%s:\n", equivchar );
exit( 1 );
#endif
	n = strlen( seq1 );
	makeequivdouble_threshold( equiv, equivchar, n );

#if 0
	fprintf( stdout, "\n" );
	for( i=0; i<n; i++ )
	{
		fprintf( stdout, "%1.0f", equiv[i] );
	}
	fprintf( stdout, "\n" );
#endif
}

static int checkcbeta( FILE *fp )
{
	char linec[1000];
	while( 1 )
	{
		fgets( linec, 999, fp );
		if( feof( fp ) ) break;
		if( !strncmp( "ATOM ", linec, 5 ) )
		{
			if( !strncmp( "CB ", linec+13, 3 ) ) return( 0 );
		}
	}
	return( 1 );
}


static double calltmalign( char **mseq1, char **mseq2, double *equiv, char *fname1, char *chain1, char *fname2, char *chain2, int alloclen )
{
	FILE *fp;
	int res;
	static char com[10000];
	double value;
	char cachedir[10000];
	char cachefile[10000];
	int runnow;


	if( usecache )
	{
		sprintf( cachedir, "%s/.tmalignoutcache", getenv( "HOME" ) );
		sprintf( com, "mkdir -p %s", cachedir );
		system( com );

		sprintf( cachefile, "%s/%s%s-%s%s", cachedir, fname1, chain1, fname2, chain2 );

		runnow = 0;
		fp = fopen( cachefile, "r" );
		if( fp == NULL ) runnow = 1;
		else
		{
			fgets( com, 100, fp );
			if( strncmp( com, "successful", 10 ) ) runnow = 1;
			fclose( fp );
		}
	}
	else
	{
		runnow = 1;
	}

	if( runnow )
	{
#if 0
		sprintf( com, "ln -s %s %s.pdb 2>_dum", fname1, fname1 );
		res = system( com );
		sprintf( com, "ln -s %s %s.pdb 2>_dum", fname2, fname2 );
		res = system( com );
#endif
		sprintf( com, "\"%s/TMalign\"  %s.pdb %s.pdb > _tmalignout 2>_dum", whereispairalign, fname1, fname2 );
		fprintf( stderr, "command = %s\n", com );
		res = system( com );
		if( res )
		{
			fprintf( stderr, "Error in TMalign\n" );
			exit( 1 );
		}
	
	}
	else
	{
		fprintf( stderr, "Cache is not supported!\n" );
		exit( 1 );
	}

	fp = fopen( "_tmalignout", "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _tmalignout\n" );
		exit( 1 );
	}

	readtmalign( fp, *mseq1, *mseq2, equiv );

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

	value = (double)naivepairscore11( *mseq1, *mseq2, penalty );

	return( value );
}

static double callrash( int mem1, int mem2, char **mseq1, char **mseq2, double *equiv, char *fname1, char *fname2, int alloclen )
{
	FILE *fp;
//	int res;
	static char com[10000];
	double value;
	char cachedir[10000];
	char cachefile[10000];
	int runnow;
	char pairid[1000];

	sprintf( pairid, "%d-%d", mem1, mem2 );
//	fprintf( stderr, "pairid = %s\n", pairid );

	if( usecache )
	{
//		sprintf( cachedir, "tmp" );
		sprintf( cachedir, "%s", pdir );

		sprintf( cachefile, "%s/%s.%s.rash", cachedir, fname1, fname2 );

		runnow = 0;
		fp = fopen( cachefile, "r" );
		if( fp == NULL )
		{
			fprintf( stderr, "Cannot open %s\n", cachefile );
			exit( 1 );
		}
		else
		{
			fclose( fp );
		}
	}
	else
	{
		fprintf( stderr, "Not supported!\n" );
		exit( 1 );
	}

#if 0
	if( 0 )
	{
#if 0
		sprintf( com, "ln -s %s %s.pdb 2>_dum", fname1, fname1 );
		res = system( com );
		sprintf( com, "ln -s %s %s.pdb 2>_dum", fname2, fname2 );
		res = system( com );
#endif
#if 0  // 091127, pdp nai!
		sprintf( com, "env PATH=%s PDP_ASH.pl --qf %s.pdb --qc %s --tf %s.pdb --tc %s > _rashout 2>_dum", whereispairalign, fname1, chain1, fname2, chain2 );
#else
		sprintf( com, "\"%s/rash\" --qf %s.pdb --qc %s --tf %s.pdb --tc %s --of %s.pdbpair > %s.rashout 2>%s.dum", whereispairalign, fname1, chain1, fname2, chain2, pairid, pairid, pairid );
#endif
		fprintf( stderr, "command = %s\n", com );
		res = system( com );
		if( res )
		{
			fprintf( stderr, "Error in structural alignment\n" );
			exit( 1 );
		}
		sprintf( com, "awk '/^REMARK/,/^TER/' %s.pdbpair > %s.%s-x-%s.%s.pdbpair", pairid, fname1, chain1, fname2, chain2 );
		res = system( com );

		sprintf( com, "awk '/^REMARK/,/^TER/{next} 1' %s.pdbpair > %s.%s-x-%s.%s.pdbpair", pairid, fname2, chain2, fname1, chain1 );
		res = system( com );

		sprintf( com, "rm %s.pdbpair", pairid );
		res = system( com );

	
	}
	else
#endif
	{
		fprintf( stderr, "Use cache! cachefile = %s\n", cachefile );
		sprintf( com, "cat %s > %s.rashout", cachefile, pairid );
		system( com );
	}

	if( usecache && runnow )
	{
		fprintf( stderr, "Okashii! usechache=%d, runnow=%d\n", usecache, runnow );
		exit( 1 );
	}

	sprintf( com, "%s.rashout", pairid );
	fp = fopen( com, "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open %s\n", com );
		exit( 1 );
	}

	readrash( fp, *mseq1, *mseq2, equiv );

	fclose( fp );

//	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
//	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );


	value = (double)naivepairscore11( *mseq1, *mseq2, penalty );

	return( value );
}

static void preparetmalign( FILE *fp, char ***strfiles, char ***chainids, char ***seqpt, char ***mseq1pt, char ***mseq2pt, double **equivpt, int *alloclenpt )
{
	int i, res;
	char *dumseq;
	char line[1000];
	char fname[1000];
	char command[1000];
	int linenum, istr, nstr;
	FILE *checkfp;
	char *sline; 
	int use[10000];
	linenum = 0;
	nstr = 0;
	while( 1 )
	{
		fgets( line, 999, fp );
		if( feof( fp ) ) break;
		sline = strip( line );
		use[linenum] = 1;
		if( sline[0] == '#' || strlen( sline ) < 2 )
		{
			use[linenum] = 0;
			linenum++;
			continue;
		}
		extractfirstword( sline );
		checkfp = fopen( sline, "r" );
		if( checkfp == NULL )
		{
			fprintf( stderr, "Cannot open %s.\n", sline );
			exit( 1 );
		}
#if 0
		fgets( linec, 999, checkfp );
		if( strncmp( "HEADER ", linec, 7 ) )
		{
			fprintf( stderr, "Check the format of %s.\n", sline );
			exit( 1 );
		}
#endif
		if( checkcbeta( checkfp ) ) 
		{
			fprintf( stderr, "%s has no C-beta atoms.\n", sline );
			exit( 1 );
		}
		else
			nstr++;
		fclose( checkfp );
		linenum++;
	}
	njob = nstr;
	fprintf( stderr, "nstr = %d\n", nstr );

	*strfiles = AllocateCharMtx( nstr, 1000 );
	*chainids = AllocateCharMtx( nstr, 2 );

	rewind( fp );
	istr = 0;
	linenum = 0;
	while( 1 )
	{
		fgets( line, 999, fp );
		if( feof( fp ) ) break;
		sline = strip( line );
		if( use[linenum++] ) 
		{
			(*chainids)[istr][0] = getchainid( sline );
			(*chainids)[istr][1] = 0;
			extractfirstword( sline );
			sprintf( fname, "%s", sline );
			cutpath( fname );
			sprintf( command, "cp %s %s.pdb", sline, fname );
			system( command );
			sprintf( command, "perl \"%s/clean.pl\" %s.pdb", whereispairalign, fname );
			res = system( command );
			if( res )
			{
				fprintf( stderr, "error: Install clean.pl\n" );
				exit( 1 );
			}
			strcpy( (*strfiles)[istr++], fname );
		}
	}

	*seqpt = AllocateCharMtx( njob, nlenmax*2+1 );
	*mseq1pt = AllocateCharMtx( njob, 0 );
	*mseq2pt = AllocateCharMtx( njob, 0 );
	*equivpt = AllocateDoubleVec( nlenmax*2+1 );
	*alloclenpt = nlenmax*2;
	dumseq = AllocateCharVec( nlenmax*2+1 );
	alreadyoutput = AllocateIntVec( njob );
	for( i=0; i<njob; i++ ) alreadyoutput[i] = 0;

	for( i=0; i<istr; i++ )
	{
		fprintf( stderr, "i=%d\n", i );
		(*seqpt)[i][0] = 0;

		(*mseq1pt)[0] = (*seqpt)[i];
		(*mseq2pt)[0] = dumseq;

		calltmalign( *mseq1pt, *mseq2pt, *equivpt, (*strfiles)[i], (*chainids)[i], (*strfiles)[i], (*chainids)[i], *alloclenpt );
		fprintf( stdout, ">%d_%s-%s\n%s\n", i+1, (*strfiles)[i], (*chainids)[i], (*seqpt)[i] );
		alreadyoutput[i] = 1;
	}
}

static void prepareash( FILE *fp, char *inputfile, char ***strfiles, char ***chainids, char ***seqpt, char ***mseq1pt, char ***mseq2pt, double **equivpt, int *alloclenpt )
{
	int i, res;
	char *dumseq;
	char line[1000];
	char fname[1000];
	char command[1000];
	int linenum, istr, nstr;
//	FILE *checkfp;
	char *sline; 
	int use[10000];
	linenum = 0;
	nstr = 0;

	fprintf( stderr, "inputfile = %s\n", inputfile );
	while( 1 )
	{
		fgets( line, 999, fp );
		if( feof( fp ) ) break;
		sline = strip( line );
		use[linenum] = 1;
		if( sline[0] == '#' || strlen( sline ) < 2 )
		{
			use[linenum] = 0;
			linenum++;
			continue;
		}
		extractfirstword( sline );
#if 0
		checkfp = fopen( sline, "r" );
		if( checkfp == NULL )
		{
			fprintf( stderr, "Cannot open %s.\n", sline );
			exit( 1 );
		}
		if( checkcbeta( checkfp ) ) 
		{
			fprintf( stderr, "%s has no C-beta atoms.\n", sline );
			exit( 1 );
		}
		else
			nstr++;
		fclose( checkfp );
#else
		nstr++;
#endif
		linenum++;
	}
	njob = nstr;
	fprintf( stderr, "nstr = %d\n", nstr );

	*strfiles = AllocateCharMtx( nstr, 1000 );
	*chainids = AllocateCharMtx( nstr, 2 );

	rewind( fp );
	istr = 0;
	linenum = 0;
	while( 1 )
	{
		fgets( line, 999, fp );
		if( feof( fp ) ) break;
		sline = strip( line );
		fprintf( stderr, "sline = %s\n", sline );
		if( use[linenum++] ) 
		{
			(*chainids)[istr][0] = getchainid( sline );
			(*chainids)[istr][1] = 0;
			extractfirstword( sline );
			sprintf( fname, "%s", sline );
			cutpath( fname );
#if 0
			sprintf( command, "cp %s %s.pdb", sline, fname );
			system( command );
			sprintf( command, "perl \"%s/clean.pl\" %s.pdb", whereispairalign, fname );
			res = system( command );
			if( res )
			{
				fprintf( stderr, "error: Install clean.pl\n" );
				exit( 1 );
			}
#endif
			strcpy( (*strfiles)[istr++], fname );
		}
	}

	*seqpt = AllocateCharMtx( njob, nlenmax*2+1 );
	*mseq1pt = AllocateCharMtx( njob, 0 );
	*mseq2pt = AllocateCharMtx( njob, 0 );
	*equivpt = AllocateDoubleVec( nlenmax*2+1 );
	*alloclenpt = nlenmax*2;
	dumseq = AllocateCharVec( nlenmax*2+1 );
	alreadyoutput = AllocateIntVec( njob );
	for( i=0; i<njob; i++ ) alreadyoutput[i] = 0;

	fprintf( stderr, "Running pdp_ash_batch.pl..\n" );
//	sprintf( command, "/opt/protein/share/domains/code/pdp_ash/pdp_ash_batch.pl -f %s -d tmp -i %d", inputfile, wheretooutput  );
	sprintf( command, "/opt/protein/share/mafftash/pdp_ash/pdp_ash_batch.pl -f %s -d %s -i %s", inputfile, pdir, odir );
	res = system( command );
	if( res )
	{
		fprintf( stderr, "Ask KM!\n" );
		exit( 1 );
	}
	fprintf( stderr, "done\n" );


	for( i=0; i<istr; i++ )
	{
		fprintf( stderr, "i=%d\n", i );
		(*seqpt)[i][0] = 0;

		(*mseq1pt)[0] = (*seqpt)[i];
		(*mseq2pt)[0] = dumseq;

		callrash( i, i, *mseq1pt, *mseq2pt, *equivpt, (*strfiles)[i], (*strfiles)[i], *alloclenpt );
		fprintf( stdout, ">%d_%s\n%s\n", i+1, (*strfiles)[i], (*seqpt)[i] );
		alreadyoutput[i] = 1;
	}
}

void arguments( int argc, char *argv[] )
{
    int c;

	usecache = 0;
	scale = 1.0;
	equivthreshold = 5;
	equivwinsize = 5;
	equivshortestlen = 1;
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
    alg = 'R';
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
	RNAppenalty = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	odir = "";
	pdir = "";

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
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
					fprintf( stderr, "jtt %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "TM %d\n", pamN );
					--argc;
					goto nextoption;
				case 'd':
					whereispairalign = *++argv;
					fprintf( stderr, "whereispairalign = %s\n", whereispairalign );
					--argc; 
					goto nextoption;
				case 'o':
					odir = *++argv;
					fprintf( stderr, "odir = %s\n", odir );
					--argc; 
					goto nextoption;
				case 'p':
					pdir = *++argv;
					fprintf( stderr, "pdir = %s\n", pdir );
					--argc; 
					goto nextoption;
				case 't':
					equivthreshold = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'w':
					equivwinsize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'l':
					equivshortestlen = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 's':
					scale = atof( *++argv );
					--argc;
					goto nextoption;
				case 'c':
					usecache = 1;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
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
#if 0
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
				case 'Q':
					calledByXced = 1;
					break;
				case 'x':
					disp = 1;
					break;
#if 0
				case 'a':
					alg = 'a';
					break;
#endif
				case 'S':
					alg = 'S';
					break;
				case 'L':
					alg = 'L';
					break;
				case 'B':
					alg = 'B';
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
				case 'N':
					alg = 'N';
					break;
				case 'K':
					alg = 'K';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'V':
					alg = 'V';
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
				case 'y':
					divpairscore = 1;
					break;
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
#if 0
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
#endif
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

static void pairalign( char **name, int nlen[M], char **seq, char **aseq, char **mseq1, char **mseq2, double *equiv, double *effarr, char **strfiles, char **chainids, int alloclen )
{
	int i, j, ilim;
	int clus1, clus2;
	int off1, off2;
	double pscore = 0.0; // by D.Mathog
	static char *indication1, *indication2;
	FILE *hat2p, *hat3p;
	static double **distancemtx;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	char *pt;
	char *hat2file = "hat2";
	LocalHom **localhomtable, *tmpptr;
	static char **pair;
//	int intdum;
	double bunbo;
	char **checkseq;


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
			localhomtable[i][j].nokori = 0;
		}
	}

	if( effarr1 == NULL ) 
	{
		distancemtx = AllocateDoubleMtx( njob, njob );
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		checkseq = AllocateCharMtx( njob, alloclen );
#if 0
#else
		pair = AllocateCharMtx( njob, njob );
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

	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) pair[i][j] = 0;
	for( i=0; i<njob; i++ ) pair[i][i] = 1;

	for( i=0; i<njob; i++ )
	{
		strcpy( checkseq[i], seq[i] );
//		fprintf( stderr, "checkseq[%d] = %s\n", i, checkseq[i] );
	}


	ilim = njob - 1;
	for( i=0; i<ilim; i++ ) 
	{
		fprintf( stderr, "% 5d / %d\r", i, njob );


		for( j=i+1; j<njob; j++ )
		{


#if 0
			if( strlen( seq[i] ) == 0 || strlen( seq[j] ) == 0 )
			{
				distancemtx[i][j] = pscore;
				continue;
			}
#endif

			strcpy( aseq[i], seq[i] );
			strcpy( aseq[j], seq[j] );
			clus1 = conjuctionfortbfast_old( pair, i, aseq, mseq1, effarr1, effarr, indication1 );
			clus2 = conjuctionfortbfast_old( pair, j, aseq, mseq2, effarr2, effarr, indication2 );
	//		fprintf( stderr, "mseq1 = %s\n", mseq1[0] );
	//		fprintf( stderr, "mseq2 = %s\n", mseq2[0] );
	
#if 0
			fprintf( stderr, "group1 = %.66s", indication1 );
			fprintf( stderr, "\n" );
			fprintf( stderr, "group2 = %.66s", indication2 );
			fprintf( stderr, "\n" );
#endif
//			for( l=0; l<clus1; l++ ) fprintf( stderr, "## STEP-eff for mseq1-%d %f\n", l, effarr1[l] );
	
#if 1
			{
				switch( alg )
				{
					case( 'T' ):
						fprintf( stderr, "  Calling tmalign %d-%d/%d    \r", i+1, j+1, njob );
						pscore = calltmalign( mseq1, mseq2, equiv, strfiles[i], chainids[i], strfiles[j], chainids[j], alloclen );
						off1 = off2 = 0;
						break;
					case( 'R' ):
						fprintf( stderr, "  Calling PDP_ASH.pl %d-%d/%d    \r", i+1, j+1, njob );
						pscore = callrash( i, j, mseq1, mseq2, equiv, strfiles[i], strfiles[j], alloclen );
						off1 = off2 = 0;
						break;
					ErrorExit( "ERROR IN SOURCE FILE" );
				}
			}
#endif
			distancemtx[i][j] = pscore;
#if SCOREOUT
			fprintf( stderr, "score = %10.2f (%d,%d)\n", pscore, i, j );
#endif

			putlocalhom_str( mseq1[0], mseq2[0], equiv, scale, localhomtable[i]+j, off1, off2, (int)pscore, strlen( mseq1[0] ), 'k' );
#if 1

			if( alreadyoutput[i] == 0 )
			{
				alreadyoutput[i] = 1;
				gappick0( seq[i], mseq1[0] );
				fprintf( stdout, ">%d_%s\n%s\n", i+1, strfiles[i], seq[i] );
				strcpy( checkseq[i], seq[i] );
			}
			else
			{
				gappick0( seq[i], mseq1[0] );
				fprintf( stderr, "checking seq%d\n", i );

//				fprintf( stderr, "     seq=%s\n", seq[i] );
//				fprintf( stderr, "checkseq=%s\n", checkseq[i] );

				if( strcmp( checkseq[i], seq[i] ) )
				{
					fprintf( stderr, "\n\nWARNING: Sequence changed!!\n" );
					fprintf( stderr, "i=%d\n", i );
					fprintf( stderr, "     seq=%s\n", seq[i] );
					fprintf( stderr, "checkseq=%s\n", checkseq[i] );
					exit( 1 );
				}
			}
			if( alreadyoutput[j] == 0 )
			{
				alreadyoutput[j] = 1;
				gappick0( seq[j], mseq2[0] );
				fprintf( stdout, ">%d_%s-%s\n%s\n", j+1, strfiles[j], chainids[j], seq[j] );
				strcpy( checkseq[j], seq[j] );
			}
			else
			{
				gappick0( seq[j], mseq2[0] );
				fprintf( stderr, "checking seq%d\n", j );
				if( strcmp( checkseq[j], seq[j] ) )
				{
					fprintf( stderr, "\n\nWARNING: Sequence changed!!\n" );
					fprintf( stderr, "j=%d\n", j );
					fprintf( stderr, "     seq=%s\n", seq[j] );
					fprintf( stderr, "checkseq=%s\n", checkseq[j] );
					exit( 1 );
				}
			}
#endif
		}
	}
	for( i=0; i<njob; i++ )
	{
		pscore = 0.0;
		for( pt=seq[i]; *pt; pt++ )
			pscore += amino_dis[(int)*pt][(int)*pt];
		distancemtx[i][i] = pscore;

	}

	ilim = njob-1;	
	for( i=0; i<ilim; i++ )
	{
		for( j=i+1; j<njob; j++ )
		{
			bunbo = MIN( distancemtx[i][i], distancemtx[j][j] );
			if( bunbo == 0.0 )
				distancemtx[i][j] = 2.0;
			else
				distancemtx[i][j] = ( 1.0 - distancemtx[i][j] / bunbo ) * 2.0;
		}
	}

	hat2p = fopen( hat2file, "w" );
	if( !hat2p ) ErrorExit( "Cannot open hat2." );
	WriteHat2_pointer( hat2p, njob, name, distancemtx );
	fclose( hat2p );

	fprintf( stderr, "##### writing hat3\n" );
	hat3p = fopen( "hat3", "w" );
	if( !hat3p ) ErrorExit( "Cannot open hat3." );
	ilim = njob-1;	
	for( i=0; i<ilim; i++ ) 
	{
		for( j=i+1; j<njob; j++ )
		{
			for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d k\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2 ); 
			}
		}
	}
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
	static char **mseq1, **mseq2;
	static char **aseq;
	static char **bseq;
	static double *eff;
	static double *equiv;
	char **strfiles;
	char **chainids;
	int i;
	FILE *infp;
	char c;
	int alloclen;

	arguments( argc, argv );

	if( equivthreshold < 1 || 9 < equivthreshold )
	{
		fprintf( stderr, "-t n, n must be 1..9\n" );
		exit( 1 );
	}

	if( ( equivwinsize + 1 ) % 2 != 0 )
	{
		fprintf( stderr, "equivwinsize = %d\n", equivwinsize );
		fprintf( stderr, "It must be an odd number.\n" );
		exit( 1 );
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

	nlenmax = 10000; // tekitou

	if( alg == 'R' )
		prepareash( infp, inputfile, &strfiles, &chainids, &seq, &mseq1, &mseq2, &equiv, &alloclen );
	else if( alg == 'T' )
		preparetmalign( infp, &strfiles, &chainids, &seq, &mseq1, &mseq2, &equiv, &alloclen );

	fclose( infp );

	name = AllocateCharMtx( njob, B+1 );
	aseq = AllocateCharMtx( njob, nlenmax*2+1 );
	bseq = AllocateCharMtx( njob, nlenmax*2+1 );
	eff = AllocateDoubleVec( njob );

	for( i=0; i<njob; i++ )
	{
		fprintf( stderr, "str%d = %s-%s\n", i, strfiles[i], chainids[i] );
	}

	if( njob < 1 )
	{
		fprintf( stderr, "No structure found.\n" ); 
		exit( 1 );
	}
	if( njob < 2 )
	{
		fprintf( stderr, "Only %d structure found.\n", njob ); 
		exit( 0 );
	}
	if( njob > M )
	{
		fprintf( stderr, "The number of structures must be < %d\n", M );
		fprintf( stderr, "Please try sequence-based methods for such large data.\n" );
		exit( 1 );
	}



#if 0
	readData( infp, name, nlen, seq );
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
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

//	writePre( njob, name, nlen, seq, 0 );

	for( i=0; i<njob; i++ ) eff[i] = 1.0;


	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );

	pairalign( name, nlen, bseq, aseq, mseq1, mseq2, equiv, eff, strfiles, chainids, alloclen );

	fprintf( trap_g, "done.\n" );
#if DEBUG
	fprintf( stderr, "closing trap_g\n" );
#endif
	fclose( trap_g );

//	writePre( njob, name, nlen, aseq, !contin );
#if 0
	writeData( stdout, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif
	SHOWVERSION;
	return( 0 );
}
