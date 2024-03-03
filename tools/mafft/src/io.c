#include "mltaln.h"

static int upperCase = 0;

#define DEBUG   0
#define IODEBUG 0

char creverse( char f )
{
	static TLS char *table = NULL;

	if( f == 0 )
	{
		free( table );
		table = NULL;
		return( 0 );
	}

	if( table == NULL )
	{
		int i;
		table = AllocateCharVec(0x80);
		for( i=0; i<0x80; i++ ) table[i] = i;
		table['A'] = 'T';
		table['C'] = 'G';
		table['G'] = 'C';
		table['T'] = 'A';
		table['U'] = 'A';
		table['M'] = 'K';
		table['R'] = 'Y';
		table['W'] = 'W';
		table['S'] = 'S';
		table['Y'] = 'R';
		table['K'] = 'M';
		table['V'] = 'B';
		table['H'] = 'D';
		table['D'] = 'H';
		table['B'] = 'V';
		table['N'] = 'N';
		table['a'] = 't';
		table['c'] = 'g';
		table['g'] = 'c';
		table['t'] = 'a';
		table['u'] = 'a';
		table['m'] = 'k';
		table['r'] = 'y';
		table['w'] = 'w';
		table['s'] = 's';
		table['y'] = 'r';
		table['k'] = 'm';
		table['v'] = 'b';
		table['h'] = 'd';
		table['d'] = 'h';
		table['b'] = 'v';
		table['n'] = 'n';
//		table['-'] = '-';
//		table['.'] = '.';
//		table['*'] = '*';
	}
	return( table[(int)f] );
}

static int countc( char *s, char q )
{
	int v = 0;
	while( *s )
		if( *s++ == q ) v++;
	return( v );
}

static void ttou( char *s )
{
	while( *s )
	{
		if( *s == 't' ) *s = 'u';
		else if( *s == 'T' ) *s = 'U';
		s++;
	}
}

void sreverse( char *r, char *s )
{
	int numt = countc( s, 't' ) + countc( s, 'T' );
	int numu = countc( s, 'u' ) + countc( s, 'U' );

//	reporterr( "numt=%d, numu=%d\n", numt, numu );
//	reporterr( "s=%s\n", s );

	r += strlen( s );
	*r-- = 0;
	while( *s )
		*r-- = creverse( *s++ );
//		*r-- = ( *s++ );
	if( numu > numt )
	{
//		reporterr( "RNA!\n" );
//		reporterr( "r before ttou =%s\n", r );
		ttou( r+1 );
//		reporterr( "r after ttou =%s\n", r );
	}
}

void gappick_samestring( char *seq )
{
	char *aseq = seq;

	for( ; *seq != 0; seq++ )
	{
		if( *seq != '-' )
			*aseq++ = *seq;
	}
	*aseq = 0;
}

#if 0

static int addlocalhom2( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, int skip )
{
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt;
	int st;
	int nlocalhom = 0;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	isumscore = 0;
	sumoverlap = 0;

#if 0
	fprintf( stderr, "nlocalhom = %d in addlocalhom\n", nlocalhom );
	fprintf( stderr, "al1 = %s, al2 = %s\n", al1, al2 );
	fprintf( stderr, "off1 = %d, off2 = %d\n", off1, off2 );
	fprintf( stderr, "localhopt = %p, skip = %d\n", localhompt, skip );
	fprintf( stderr, "pt1 = \n%s\n, pt2 = \n%s\n", pt1, pt2 );
#endif

	if( skip )
	{
		while( --skip > 0 ) localhompt = localhompt->next;
		localhompt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
		localhompt = localhompt->next;
//		fprintf( stderr, "tmppt = %p, localhompt = %p\n", tmppt, localhompt );
	}
	tmppt = localhompt;

	st = 0;
	iscore = 0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "In in while loop\n" );
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

#if 1
			isumscore += iscore;
			sumoverlap += end2-start2+1;
#else
			tmppt->overlapaa   = end2-start2+1;
			tmppt->opt = iscore * 5.8 / 600;
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "iscore (1)= %d\n", iscore );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			iscore = 0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			iscore += n_dis[(int)amino_n[(int)*pt1]][(int)amino_n[(int)*pt2]];
//			fprintf( stderr, "%c-%c, score(0) = %d\n", *pt1, *pt2, iscore );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}

	if( st )
	{
		if( nlocalhom++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}
		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;

#if 1
		isumscore += iscore;
		sumoverlap += end2-start2+1;
#else
		tmppt->overlapaa   = end2-start2+1;
		tmppt->opt = (double)iscore * 5.8 / 600;
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif
#if 0
		fprintf( stderr, "score (2)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
	{
		tmppt->overlapaa = sumoverlap;
		tmppt->opt = (double)sumscore * 5.8 / 600 / sumoverlap;
	}
	return( nlocalhom );
}

#endif



static int addlocalhom_r( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, int skip, char korh )
{
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	double score;
	double sumscore;
	int sumoverlap;
	LocalHom *tmppt = NULL; // by D.Mathog, a guess
	int st;
	int nlocalhom = 0;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	sumscore = 0.0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

#if 0
	fprintf( stderr, "nlocalhom = %d in addlocalhom\n", nlocalhom );
	fprintf( stderr, "al1 = %s, al2 = %s\n", al1, al2 );
	fprintf( stderr, "off1 = %d, off2 = %d\n", off1, off2 );
	fprintf( stderr, "localhopt = %p, skip = %d\n", localhompt, skip );
#endif
	fprintf( stderr, "pt1 = \n%s\n, pt2 = \n%s\n", pt1, pt2 );

	if( skip )
	{
		while( --skip > 0 ) localhompt = localhompt->next;
		localhompt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
		localhompt = localhompt->next;
		fprintf( stderr, "tmppt = %p, localhompt = %p\n", (void *)tmppt, (void *)localhompt );
	}
	tmppt = localhompt;

	st = 0;
	score = 0.0;
	while( *pt1 != 0 )
	{
		fprintf( stderr, "In in while loop\n" );
		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;
			tmppt->korh   = korh  ;

#if 1
			sumscore += score;
			sumoverlap += end2-start2+1;
#else
			tmppt->overlapaa   = end2-start2+1;
			tmppt->opt = score * 5.8 / 600;
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

			fprintf( stderr, "score (1)= %f\n", score );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
			score = 0.0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			score += (double)n_dis[(int)amino_n[(unsigned char)*pt1]][(int)amino_n[(unsigned char)*pt2]];
//			fprintf( stderr, "%c-%c, score(0) = %f\n", *pt1, *pt2, score );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( nlocalhom++ > 0 )
	{
//		fprintf( stderr, "reallocating ...\n" );
		tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//		fprintf( stderr, "done\n" );
		tmppt = tmppt->next;
		tmppt->next = NULL;
	}
	end1 = pos1 - 1;
	end2 = pos2 - 1;
	tmppt->start1 = start1;
	tmppt->start2 = start2;
	tmppt->end1   = end1  ;
	tmppt->end2   = end2  ;
	tmppt->korh   = korh  ;

#if 1
	sumscore += score;
	sumoverlap += end2-start2+1;
#else
	tmppt->overlapaa   = end2-start2+1;
	tmppt->opt = score * 5.8 / 600;
	tmppt->overlapaa   = overlapaa;
	tmppt->opt = (double)opt;
#endif

	fprintf( stderr, "score (2)= %f\n", score );
	fprintf( stderr, "al1: %d - %d\n", start1, end1 );
	fprintf( stderr, "al2: %d - %d\n", start2, end2 );

	for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
	{
		tmppt->overlapaa = sumoverlap;
		tmppt->opt = sumscore * 5.8 / 600 / sumoverlap;
	}
	return( nlocalhom );
}
void putlocalhom3( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, char korh )
{
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	double score;
	double sumscore;
	int sumoverlap;
	LocalHom *tmppt;
	LocalHom *subnosento;
	int st;
	int saisho;

	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	sumscore = 0.0;
	sumoverlap = 0;
	start1 = 0; // by Mathog, a guess
	start2 = 0; // by Mathog, a guess

	subnosento = localhompt;
	while( subnosento->next ) subnosento = subnosento->next;
	tmppt = subnosento;

	saisho = ( localhompt->nokori == 0 );

	fprintf( stderr, "localhompt = %p\n", (void *)localhompt );
	fprintf( stderr, "tmppt = %p\n", (void *)tmppt );
	fprintf( stderr, "subnosento = %p\n", (void *)subnosento );

	st = 0;
	score = 0.0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( localhompt->nokori++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;
			tmppt->korh   = korh  ;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				if( tmppt->overlapaa>0)
					tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
				else
					tmppt->opt = -1.0;
			}
			else
			{
				sumscore += score;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "score (1)= %f\n", score );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			score = 0.0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			score += (double)n_dis[(int)amino_n[(unsigned char)*pt1]][(int)amino_n[(unsigned char)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, score(0) = %f\n", *pt1, *pt2, score );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( *(pt1-1) != '-' && *(pt2-1) != '-'  )
	{
		if( localhompt->nokori++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}

		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;
		tmppt->korh   = korh  ;


#if 1
		if( divpairscore )
		{
			tmppt->overlapaa   = end2-start2+1;
			if( tmppt->overlapaa>0)
				tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
			else
				tmppt->opt = -1.0;
		}
		else
		{
			sumscore += score;
			sumoverlap += end2-start2+1;
		}
#else
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif

#if 0
		fprintf( stderr, "score (2)= %f\n", score );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	fprintf( stderr, "sumscore = %f\n", sumscore );
	if( !divpairscore )
	{

		if( !saisho ) subnosento = subnosento->next;
		for( tmppt=subnosento; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
			if( tmppt->overlapaa>0)
				tmppt->opt = sumscore * 5.8 / 600 / sumoverlap;
			else
				tmppt->opt = -1.0;
			fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
		}
	}
}
void putlocalhom_ext( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, char korh )
{
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;


	isumscore = 0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	st = 0;
	iscore = 0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;
			tmppt->korh   = korh  ;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				if( tmppt->overlapaa>0)
					tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
				else
					tmppt->opt = -1.0;
			}
			else
			{
				isumscore += iscore;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "iscore (1)= %d\n", iscore );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			iscore = 0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			iscore += n_dis[(int)amino_n[(unsigned char)*pt1]][(int)amino_n[(unsigned char)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, iscore(0) = %d\n", *pt1, *pt2, iscore );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( *(pt1-1) != '-' && *(pt2-1) != '-'  )
	{
		if( nlocalhom++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}
		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;
		tmppt->korh   = korh  ;
	
#if 1
		if( divpairscore )
		{
			tmppt->overlapaa   = end2-start2+1;
			if( tmppt->overlapaa>0)
				tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			else
				tmppt->opt = -1.0;
		}
		else
		{
			isumscore += iscore;
			sumoverlap += end2-start2+1;
		}
#else
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif

#if 0
		fprintf( stderr, "iscore (2)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	if( !divpairscore )
	{
		for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
//			tmppt->opt = (double)isumscore * 5.8 / ( 600 * sumoverlap );
			tmppt->opt = (double)600 * 5.8 / 600;
//			fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
		}
	}
}

void putlocalhom_str( char *al1, char *al2, double *equiv, double scale, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, char korh )
{
	int posinaln, pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
//	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	isumscore = 0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	posinaln = 0;
	while( *pt1 != 0 )
	{
		if( *pt1 != '-' && *pt2 != '-' && equiv[posinaln] > 0.0 )
		{
			start1 = end1 = pos1; start2 = end2 = pos2;
			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ... (posinaln=%d)\n", posinaln );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;
			tmppt->korh   = korh  ;

			tmppt->overlapaa   = 1;
//			tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			tmppt->opt = equiv[posinaln] * scale;
//			fprintf( stdout, "*pt1=%c, *pt2=%c, equiv=%f\n", *pt1, *pt2, equiv[posinaln] );

		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
		posinaln++;
	}
}



void putlocalhom2( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, char korh )
{
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;



	isumscore = 0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	st = 0;
	iscore = 0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;
			tmppt->korh   = korh  ;
			tmppt->nokori  += 1;
			localhompt->last  = tmppt;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				if(tmppt->overlapaa>0)
					tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
				else
					tmppt->opt = -1.0;
			}
			else
			{
				isumscore += iscore;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "iscore (1)= %d\n", iscore );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			iscore = 0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			iscore += n_dis[(int)amino_n[(unsigned char)*pt1]][(int)amino_n[(unsigned char)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, iscore(0) = %d\n", *pt1, *pt2, iscore );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( *(pt1-1) != '-' && *(pt2-1) != '-'  )
	{
		if( nlocalhom++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}
		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;
		tmppt->korh   = korh  ;
		tmppt->nokori  += 1;
		localhompt->last  = tmppt;
	
#if 1
		if( divpairscore )
		{
			tmppt->overlapaa   = end2-start2+1;
			if(tmppt->overlapaa>0)
				tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			else
				tmppt->opt = -1.0;
		}
		else
		{
			isumscore += iscore;
			sumoverlap += end2-start2+1;
		}
#else
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif

#if 0
		fprintf( stderr, "iscore (2)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	if( !divpairscore )
	{
		for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
			if(tmppt->overlapaa>0)
				tmppt->opt = (double)isumscore * 5.8 / ( 600 * sumoverlap );
			else
				tmppt->opt = -1.0;
//			fprintf( stderr, "tmpptr->opt = %f, sumoverlap=%d\n", tmppt->opt, sumoverlap );
		}
	}
}

#if 0
void putlocalhom( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, char korh )
{
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	double score;
	double sumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;


	sumscore = 0.0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	st = 0;
	score = 0.0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;
			tmppt->korh   = korh  ;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
			}
			else
			{
				sumscore += score;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "score (1)= %f\n", score );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			score = 0.0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			score += (double)n_dis[(int)amino_n[(unsigned char)*pt1]][(int)amino_n[(unsigned char)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, score(0) = %f\n", *pt1, *pt2, score );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( nlocalhom++ > 0 )
	{
//		fprintf( stderr, "reallocating ...\n" );
		tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//		fprintf( stderr, "done\n" );
		tmppt = tmppt->next;
		tmppt->next = NULL;
	}
	end1 = pos1 - 1;
	end2 = pos2 - 1;
	tmppt->start1 = start1;
	tmppt->start2 = start2;
	tmppt->end1   = end1  ;
	tmppt->end2   = end2  ;
	tmppt->korh   = korh  ;

#if 1
	if( divpairscore )
	{
		tmppt->overlapaa   = end2-start2+1;
		tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
	}
	else
	{
		sumscore += score;
		sumoverlap += end2-start2+1;
	}
#else
	tmppt->overlapaa   = overlapaa;
	tmppt->opt = (double)opt;
#endif

#if 0
	fprintf( stderr, "score (2)= %f\n", score );
	fprintf( stderr, "al1: %d - %d\n", start1, end1 );
	fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif

	if( !divpairscore )
	{
		for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
			tmppt->opt = sumscore * 5.8 / 600 / sumoverlap;
//			fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
		}
	}
}
#endif

char *cutal( char *al, int al_display_start, int start, int end )
{
	int pos;
	char *pt = al;
	char *val = NULL;

	pos = al_display_start;
	do
	{
		if( start == pos ) val = pt;
		if( end == pos ) break;
//		fprintf( stderr, "pos=%d, *pt=%c, val=%p\n", pos, *pt, val );
		if( *pt != '-' ) pos++;
	} while( *pt++ != 0 );
	*(pt+1) = 0;
	return( val );
}

void ErrorExit( char *message )
{
	fprintf( stderr, "%s\n", message );
	exit( 1 );
}

void strncpy_caseC( char *str1, char *str2, int len )
{
	if( dorp == 'd' && upperCase > 0 ) 
	{
		while( len-- )
			*str1++ = toupper( *str2++ );
	}
	else strncpy( str1, str2, len );
}
	
void seqUpper( int nseq, char **seq )
{
	int i, j, len;
	for( i=0; i<nseq; i++ ) 
	{
		len = strlen( seq[i] );
		for( j=0; j<len; j++ ) 
			seq[i][j] = toupper( seq[i][j] );
	}
}

void seqLower( int nseq, char **seq )
{
	int i, j, len;
	for( i=0; i<nseq; i++ ) 
	{
		len = strlen( seq[i] );
		for( j=0; j<len; j++ ) 
			seq[i][j] = tolower( seq[i][j] );
	}
}

int getaline_fp_eof( char *s, int l, FILE *fp )  /* end of file -> return 1 */
{
    int c, i = 0 ;
    int noteofflag = 0;
    for( i=0; i<l && ( noteofflag = ( (c=getc(fp)) != EOF ) ) && c != '\n'; i++ ) 
    	*s++ = c;
    *s = '\0' ;
     return( !noteofflag );
}

int getaline_fp_eof_new(s, l, fp)  /* end of file -> return 1 */
char    s[] ; int l ; FILE *fp ;
{
        int c = 0, i = 0 ;
		int noteofflag = 0;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( noteofflag = ( (c=getc(fp)) != EOF ) ) && c != '\n'; i++ ) 
        { *s++ = c ; }
        *s = '\0' ;
		if( c != '\n' && c != EOF ) while( getc(fp) != '\n' )
			;
		return( !noteofflag );
}

int myfgets(s, l, fp)  /* l以上は、行末まで読み飛ばす */
char    s[] ; int l ; FILE *fp ;
{
        int     c = 0, i = 0 ;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( c=getc( fp ) ) != '\n'; i++ ) 
        	*s++ = c;
        *s = '\0' ;
		if( c != '\n' ) 
			while( getc(fp) != '\n' )
				;
		return( 0 );
}

double input_new( FILE *fp, int d )
{
	char mojiretsu[10];
	int i, c;

	c = getc( fp );
	if( c != '\n' )
		ungetc( c, fp );

	for( i=0; i<d; i++ )
		mojiretsu[i] = getc( fp );
	mojiretsu[i] = 0;

	return( atof( mojiretsu ) );
}


void PreRead( FILE *fp, int *locnjob, int *locnlenmax )
{
	int i, nleni;
	char b[B];

	fgets( b, B-1, fp ); *locnjob = atoi( b );
	*locnlenmax = 0;
	i=0; 
	while( i<*locnjob )
	{
		fgets( b, B-1, fp );
		if( !strncmp( b, "=", 1 ) )
		{
			i++;
			fgets( b, B-1, fp ); nleni = atoi( b );
			if( nleni > *locnlenmax ) *locnlenmax = nleni;
		}
	}
	if( *locnlenmax > N )
	{
		fprintf( stderr, "TOO LONG SEQUENCE!\n" );
		exit( 1 );
	}
	if( njob > M  ) 
	{
		fprintf( stderr, "TOO MANY SEQUENCE!\n" );
		fprintf( stderr, "%d > %d\n", njob, M );
		exit( 1 );
	}
}	

int allSpace( char *str )
{
	int value = 1;
	while( *str ) value *= ( !isdigit( *str++ ) );
	return( value );
}
	
void Read( char name[M][B], int nlen[M], char **seq )
{
	extern void FRead( FILE *x, char y[M][B], int z[M], char **w );
	FRead( stdin, name, nlen, seq );
}


void FRead( FILE *fp, char name[][B], int nlen[], char **seq )
{
	int i, j; 
	char b[B];

	fgets( b, B-1, fp );
#if DEBUG
	fprintf( stderr, "b = %s\n", b );
#endif

    if( strstr( b, "onnet" ) ) scoremtx = 1;
    else if( strstr( b, "DnA" ) ) 
	{
		scoremtx = -1; 
		upperCase = -1;
	}
    else if( strstr( b, "dna" ) ) 
	{
		scoremtx = -1; 
		upperCase = 0;
	}
	else if( strstr( b, "DNA" ) )
	{
		scoremtx = -1; 
		upperCase = 1;
	}
    else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2; 
    else scoremtx = 0;
#if DEBUG
	fprintf( stderr, " %s->scoremtx = %d\n", b, scoremtx );
#endif

	geta2 = GETA2;

#if 0
	if( strlen( b ) >=25 )
	{
		b[25] = 0;
	#if DEBUG
		fprintf( stderr, "kimuraR = %s\n", b+20 );
	#endif
		kimuraR = atoi( b+20 );

		if( kimuraR < 0 || 20 < kimuraR ) ErrorExit( "Illeagal kimuraR value.\n" );
		if( allSpace( b+20 ) ) kimuraR = NOTSPECIFIED;
	}
	else kimuraR = NOTSPECIFIED;
	#if DEBUG
	fprintf( stderr, "kimuraR = %d\n", kimuraR );
	#endif

	if( strlen( b ) >=20 )
	{
		b[20] = 0;
	#if DEBUG
		fprintf( stderr, "pamN = %s\n", b+15 );
	#endif
		pamN = atoi( b+15 );
		if( pamN < 0 || 400 < pamN ) ErrorExit( "Illeagal pam value.\n" );
		if( allSpace( b+15 ) ) pamN = NOTSPECIFIED;
	}
	else pamN = NOTSPECIFIED;

	if( strlen( b ) >= 15 )
	{
		b[15] = 0;
	#if DEBUG
		fprintf( stderr, "poffset = %s\n", b+10 );
	#endif
		poffset = atoi( b+10 );
		if( poffset > 500 ) ErrorExit( "Illegal extending gap ppenalty\n" );
		if( allSpace( b+10 ) ) poffset = NOTSPECIFIED;
	}
	else poffset = NOTSPECIFIED;

	if( strlen( b ) >= 10 )
	{
		b[10] = 0;
	#if DEBUG
		fprintf( stderr, "ppenalty = %s\n", b+5 );
	#endif
		ppenalty = atoi( b+5 );
		if( ppenalty > 0 ) ErrorExit( "Illegal opening gap ppenalty\n" );
		if( allSpace( b+5 ) ) ppenalty = NOTSPECIFIED;
	}
	else ppenalty = NOTSPECIFIED;
#endif

	for( i=0; i<njob; i++ )
	{
		getaline_fp_eof_new( b, B-1, fp );
		strcpy( name[i], b );
#if DEBUG
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		fgets( b, B-1, fp ); nlen[i] = atoi( b );      /* seq i no nagasa */
		seq[i][0] = 0;
		if( nlen[i] ) for( j=0; j <= (nlen[i]-1)/C; j++ )
		{
			getaline_fp_eof_new( b, B-1, fp );
			/*	b[C] = 0;  */
			strcat( seq[i], b );
		} 
		seq[i][nlen[i]] = 0;
	}
	if( scoremtx == -1 && upperCase != -1 ) seqLower( njob, seq );
}


static int countKUorWA( FILE *fp )
{
	int value;
	int c, b;

	value= 0;
	b = '\n';
	while( ( c = getc( fp ) ) != EOF )
	{
		if( b == '\n' && ( c == '>' ) )
			value++;
		b = c;
	}
	rewind( fp );
	return( value );
}

void searchKUorWA( FILE *fp )
{
	int c, b;
	b = '\n';
	while( !( ( ( c = getc( fp ) ) == '>' || c == EOF ) && b == '\n' ) )
		b = c;
	ungetc( c, fp );
}

#if 0
static int onlyGraph( char *str )
{
	char tmp;
	char *res = str;
	char *bk = str;

//	while( (tmp=*str++) ) if( isgraph( tmp ) ) *res++ = tmp;
	while( (tmp=*str++) ) 
	{
		if( 0x20 < tmp && tmp < 0x7f ) *res++ = tmp;
		if( tmp == '>' || tmp == '(' )
		{
			fprintf( stderr, "========================================================\n" );
			fprintf( stderr, "========================================================\n" );
			fprintf( stderr, "=== \n" );
			fprintf( stderr, "=== ERROR!! \n" );
//			fprintf( stderr, "=== In the '--anysymbol' and '--preservecase' modes, \n" );
			fprintf( stderr, "=== '>' and '(' are acceptable only in title lines.\n" );
			fprintf( stderr, "=== \n" );
			fprintf( stderr, "========================================================\n" );
			fprintf( stderr, "========================================================\n" );
			exit( 1 );
		}
	}
	*res = 0;
	return( res - bk );
}
#endif

static int charfilter( unsigned char *str )
{
	unsigned char tmp;
	unsigned char *res = str;
	unsigned char *bk = str;

	while( (tmp=*str++) )
	{
//		if( tmp == '=' || tmp == '*' || tmp == '<' || tmp == '>' || tmp == '(' || tmp == ')' )
		if( tmp == '=' || tmp == '<' || tmp == '>' )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "Characters '= < >' can be used only in the title lines in the --anysymbol or --text mode.\n" );
			fprintf( stderr, "\n" );
			exit( 1 );
		}
//		if( 0x20 < tmp && tmp < 0x7f )
//		if( 0x0 <=tmp && tmp < 0x100 && 
		if( tmp != 0x0a && tmp != 0x20 && tmp != 0x0d )
//		if( tmp != '\n' && tmp != ' ' && tmp != '\t' ) // unprintable characters mo ok.
		{
			*res++ = tmp;
//			reporterr( "tmp=%d (%c)\n", tmp, tmp );
		}
	}
	*res = 0;
	return( res - bk );
}


static int onlyAlpha_lower( char *str )
{
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
		if( isalpha( tmp ) || tmp == '-' || tmp == '*' || tmp == '.' )
			*res++ = tolower( tmp );
	*res = 0;
	return( res - bk );
}
static int onlyAlpha_upper( char *str )
{
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
		if( isalpha( tmp ) || tmp == '-' || tmp == '*' || tmp == '.' )
			*res++ = toupper( tmp );
	*res = 0;
	return( res - bk );
}

void kake2hiku( char *str )
{
	do
		if( *str == '*' ) *str = '-';
	while( *str++ );
}

char *load1SeqWithoutName_realloc_casepreserve( FILE *fpp )
{
	int c, b;
	char *cbuf;
	int size = N;
	char *val;

	val = malloc( (size+1) * sizeof( char ) );
	cbuf = val;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&           
          !( ( c == '>' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;  /* 長すぎてもしらない */
		if( cbuf - val == size )
		{
			size += N;
			fprintf( stderr, "reallocating...\n" );
			val = (char *)realloc( val, (size+1) * sizeof( char ) );
			if( !val )
			{
				fprintf( stderr, "Allocation error in load1SeqWithoutName_realloc \n" );
				exit( 1 );
			}
			fprintf( stderr, "done.\n" );
			cbuf = val + size-N;
		}
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;
//	onlyGraph( val );
	charfilter( (unsigned char *) val );
//	kake2hiku( val );
	return( val );
}

char *load1SeqWithoutName_realloc( FILE *fpp )
{
	int c, b;
	char *cbuf;
	int size = N;
	char *val;

	val = malloc( (size+1) * sizeof( char ) );
	cbuf = val;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&           
          !( ( c == '>' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;  /* 長すぎてもしらない */
		if( cbuf - val == size )
		{
			size += N;
			fprintf( stderr, "reallocating...\n" );
			val = (char *)realloc( val, (size+1) * sizeof( char ) );
			if( !val )
			{
				fprintf( stderr, "Allocation error in load1SeqWithoutName_realloc \n" );
				exit( 1 );
			}
			fprintf( stderr, "done.\n" );
			cbuf = val + size-N;
		}
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;

	if( nblosum == -2 )
	{
		charfilter( (unsigned char *) val );
	}
	else
	{
		if( dorp == 'd' )
			onlyAlpha_lower( val );
		else
			onlyAlpha_upper( val );
		kake2hiku( val );
	}
	return( val );
}

int load1SeqWithoutName_new( FILE *fpp, char *cbuf )
{
	int c, b;
	char *bk = cbuf;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&                    /* by T. Nishiyama */
          !( ( c == '>' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;  /* 長すぎてもしらない */
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;
	if( dorp == 'd' )
		onlyAlpha_lower( bk );
	else
		onlyAlpha_upper( bk );
	kake2hiku( bk );
	return( 0 );
}


void readDataforgaln( FILE *fp, char **name, int *nlen, char **seq )
{
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		nlen[i] = strlen( seq[i] );
		free( tmpseq );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}

void readData_varlen( FILE *fp, char **name, int *nlen, char **seq )
{
	int i; 
	static char *tmpseq = NULL;

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		nlen[i] = strlen( tmpseq );
//		fprintf( stderr, "nlen[%d] = %d\n", i+1, nlen[i] );
		seq[i] = calloc( nlen[i]+1, sizeof( char ) );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}

void readData_pointer2( FILE *fp, int nseq, char **name, int *nlen, char **seq )
{
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<nseq; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
		nlen[i] = strlen( seq[i] );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( nseq, seq );
#if 0
	free( tmpseq );
#endif
	if( outnumber )
	{
		char *namebuf;
		char *cptr;
		namebuf = calloc( B+100, sizeof( char ) );
		for( i=0; i<nseq; i++ )
		{
			namebuf[0] = '=';
			cptr = strstr( name[i], "_numo_e_" );
			if( cptr )
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, cptr+8 );
			else
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, name[i]+1 );
			strncpy( name[i], namebuf, B );
			name[i][B-1] = 0;
		}
		free( namebuf );
//		exit( 1 );
	}
}


void readData_pointer_casepreserve( FILE *fp, char **name, int *nlen, char **seq )
{
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
		nlen[i] = strlen( seq[i] );
	}
}


int copydatafromgui( char **namegui, char **seqgui, char **name, int *nlen, char **seq )
{
	int i; 


	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; 
		strncpy( name[i]+1, namegui[i], B-2 );
		name[i][B-1] = 0;

		strcpy( seq[i], seqgui[i] );
		nlen[i] = strlen( seq[i] );
	}
	if( dorp == 'd' ) 
		seqLower( njob, seq );
	else if( dorp == 'p' ) 
		seqUpper( njob, seq );
	else
	{
		reporterr( "DNA or Protein?\n" );
		return( 1 );
	}
#if 0
	free( tmpseq );
#endif
	if( outnumber )
	{
		char *namebuf;
		char *cptr;
		namebuf = calloc( B+100, sizeof( char ) );
		for( i=0; i<njob; i++ )
		{
			namebuf[0] = '=';
			cptr = strstr( name[i], "_numo_e_" );
			if( cptr )
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, cptr+8 );
			else
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, name[i]+1 );
			strncpy( name[i], namebuf, B );
			name[i][B-1] = 0;
		}
		free( namebuf );
	}
	return( 0 );
}

void readData_pointer( FILE *fp, char **name, int *nlen, char **seq )
{
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
//		reporterr( "B=%d, name[i]+1=%s\n", B, name[i]+1 );
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
		nlen[i] = strlen( seq[i] );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
	if( outnumber )
	{
		char *namebuf;
		char *cptr;
		namebuf = calloc( B+100, sizeof( char ) );
		for( i=0; i<njob; i++ )
		{
			namebuf[0] = '=';
			cptr = strstr( name[i], "_numo_e_" );
			if( cptr )
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, cptr+8 );
			else
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, name[i]+1 );
			strncpy( name[i], namebuf, B );
			name[i][B-1] = 0;
		}
		free( namebuf );
//		exit( 1 );
	}
}

void readData( FILE *fp, char name[][B], int nlen[], char **seq )
{
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		nlen[i] = strlen( seq[i] );
		free( tmpseq );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}

void cutAlignment( FILE *fp, int **regtable, char **revtable, int *outtable, char **name, char **outseq )
{
	int i, j; 
	int outlen;
	static char *tmpseq = NULL;
	static char *dumname = NULL;
	char *fs, *rs;
	int npos, lpos;
	int startpos, endpos, seqlen;

	if( dumname == NULL )
	{
		dumname = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );


	npos = 0;
	for( i=0; i<njob; i++ )
	{
		dumname[0] = '>'; getc( fp ); 
		myfgets( dumname+1, B-1, fp ); 
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );

		if( outtable[i] )
		{
//			putc( '>', stdout );
//			puts( dumname+1 );

	
			strncat( name[npos], dumname, B-1 );
			name[npos][B-1] = 0;
	
			if( dorp == 'd' && upperCase != -1 ) seqLower( 1, &tmpseq );
			seqlen = strlen( tmpseq );
			lpos = 0;
			for( j=0; j<5; j++ )
			{
				if( regtable[0][j*2] == -1 && regtable[0][j*2+1] == -1 ) continue;

				startpos = regtable[0][j*2];
				endpos   = regtable[0][j*2+1];
				if( startpos > endpos )
				{
					endpos   = regtable[0][j*2];
					startpos = regtable[0][j*2+1];
				}

				if( startpos < 0 ) startpos = 0;
				if( endpos   < 0 ) endpos   = 0;
				if( endpos   >= seqlen ) endpos   = seqlen-1;
				if( startpos >= seqlen ) startpos = seqlen-1;

//				fprintf( stderr, "startpos = %d, endpos = %d\n", startpos, endpos );

				outlen = endpos - startpos+1;
				if( revtable[0][j] == 'f' )
				{
//					fprintf( stderr, "regtable[%d][st] = %d\n", i, regtable[0][j*2+0] );
//					fprintf( stderr, "regtable[%d][en] = %d\n", i, regtable[0][j*2+1] );
//					fprintf( stderr, "outlen = %d\n", outlen );
//					fprintf( stdout, "%.*s\n", outlen, tmpseq+regtable[0][j*2] );
					strncpy( outseq[npos] + lpos, tmpseq+startpos, outlen );
					lpos += outlen;
				}
				else
				{
					fs = AllocateCharVec( outlen+1 );
					rs = AllocateCharVec( outlen+1 );

					fs[outlen] = 0;
					strncpy( fs, tmpseq+startpos, outlen );
					sreverse( rs, fs );
//					fprintf( stdout, "%s\n", rs );
					strncpy( outseq[npos] + lpos, rs, outlen );
					lpos += outlen;
					free( fs );
					free( rs );
				}
				outseq[npos][lpos] = 0;
			}
			npos++;
		}
		free( tmpseq );
	}
}

void cutData( FILE *fp, int **regtable, char **revtable, int *outtable )
{
	int i, j; 
	int outlen, seqlen, startpos, endpos;
	static char *tmpseq = NULL;
	static char *dumname = NULL;
	char *fs, *rs;

	if( dumname == NULL )
	{
		dumname = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		dumname[0] = '='; getc( fp ); 
		myfgets( dumname+1, B-2, fp ); 
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );

		if( outtable[i] )
		{
			gappick_samestring( tmpseq );
			putc( '>', stdout );
			puts( dumname+1 );

			seqlen = strlen( tmpseq );

			if( dorp == 'd' && upperCase != -1 ) seqLower( 1, &tmpseq );
			if( outtable[i] == 2 )
			{
				startpos = 0;
				endpos   = seqlen-1;
				outlen = endpos - startpos + 1;
				fprintf( stdout, "%.*s\n", outlen, tmpseq+startpos );
			}
			else
			{
				for( j=0; j<5; j++ )
				{
					if( regtable[i][j*2] == -1 && regtable[i][j*2+1] == -1 ) continue;
	
					startpos = regtable[i][j*2];
					endpos   = regtable[i][j*2+1];
	
					if( startpos > endpos )
					{
						endpos   = regtable[i][j*2];
						startpos = regtable[i][j*2+1];
					}
	
					if( startpos < 0 ) startpos = 0;
					if( endpos   < 0 ) endpos   = 0;
					if( endpos   >= seqlen ) endpos   = seqlen-1;
					if( startpos >= seqlen ) startpos = seqlen-1;
	
					outlen = endpos - startpos + 1;
					if( revtable[i][j] == 'f' )
					{
						fprintf( stderr, "startpos = %d\n", startpos );
						fprintf( stderr, "endpos   = %d\n", endpos );
						fprintf( stderr, "outlen = %d\n", outlen );
						fprintf( stdout, "%.*s\n", outlen, tmpseq+startpos );
					}
					else
					{
						fs = AllocateCharVec( outlen+1 );
						rs = AllocateCharVec( outlen+1 );
	
						fs[outlen] = 0;
						strncpy( fs, tmpseq+startpos, outlen );
						sreverse( rs, fs );
						fprintf( stdout, "%s\n", rs );
						free( fs );
						free( rs );
					}
				}
			}
		}
		free( tmpseq );
	}
}

void catData( FILE *fp )
{
	int i; 
	static char *tmpseq = NULL;
	static char *dumname = NULL;
//	char *cptr;

	if( dumname == NULL )
	{
		dumname = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		dumname[0] = '='; getc( fp ); 
		myfgets( dumname+1, B-2, fp ); 
		if( outnumber )
		{
			fprintf( stdout, ">_numo_s_%08d_numo_e_", i+1 );
		}
		else
		{
			putc( '>', stdout );
		}
		puts( dumname+1 );
		tmpseq = load1SeqWithoutName_realloc( fp );
		if( dorp == 'd' && upperCase != -1 ) seqLower( 1, &tmpseq );
		puts( tmpseq );
		free( tmpseq );
	}
}

int countATGCandN( char *s, int *countN, int *total )
{
	int nATGC;
	int nChar;
	int nN;
	char c;
	nN = nATGC = nChar = 0;

	if( *s == 0 ) 
	{
		*total = 0;
		return( 0 );
	}

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
			if( c == 'n' )
				nN++;
		}
	}
	while( *++s );

//	reporterr( "nN = %d", nN );

	*total = nChar;
	*countN = nN;
	return( nATGC );
}

int countATGC( char *s, int *total )
{
	int nATGC;
	int nChar;
	char c;
	nATGC = nChar = 0;

	if( *s == 0 ) 
	{
		*total = 0;
		return( 0 );
	}

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
		}
	}
	while( *++s );

	*total = nChar;
	return( nATGC );
}

double countATGCbk( char *s )
{
	int nATGC;
	int nChar;
	char c;
	nATGC = nChar = 0;

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
		}
	}
	while( *++s );
	return( (double)nATGC / nChar );
}


int countnogaplen( char *seq )
{
	int val = 0;
	while( *seq )
		if( *seq++ != '-' ) val++;
	return( val );
}

int countnormalletters( char *seq, char *ref )
{
	int val = 0;
	while( *seq )
		if( strchr( ref, *seq++ ) ) val++;
	return( val );
}

void getnumlen_casepreserve( FILE *fp, int *nlenminpt )
{
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq, *tmpname;
	double atgcfreq;

#if mingw
	setmode( fileno( fp ), O_BINARY );
	setmode( fileno( stdout ), O_BINARY );
#endif

	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		tmp = strlen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		if( total < 1000000 )
		{
			atgcnum += countATGC( tmpseq, &nsite );
			total += nsite;
		}
		free( tmpseq );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
//	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
}

void getnumlen_nogap_countn( FILE *fp, int *nlenminpt, double *nfreq )
{
	int total;
	int nsite = 0;
	int atgcnum, nnum, nN;
	int i, tmp;
	char *tmpseq, *tmpname;
	double atgcfreq;
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	nnum = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpseq = load1SeqWithoutName_realloc( fp );
		tmp = countnogaplen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		if( total < 100000 )
		{
			atgcnum += countATGCandN( tmpseq, &nN, &nsite );
			total += nsite;
		}
		nnum += nN;
		free( tmpseq );
//		if( i % 10000 == 0 ) reporterr( "atgcnum=%d, total=%d\n", atgcnum, total );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
	*nfreq = (double)nnum / atgcnum;
//	fprintf( stderr, "##### nnum = %d\n", nnum );
//	fprintf( stderr, "##### atgcfreq = %f, *nfreq = %f\n", atgcfreq, *nfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
}
void getnumlen_nogap( FILE *fp, int *nlenminpt )
{
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq, *tmpname;
	double atgcfreq;
#if mingw // web nomi de shiyou suru node nakutemo ii
	setmode( fileno( fp ), O_BINARY );
	setmode( fileno( stdout ), O_BINARY );
#endif
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpseq = load1SeqWithoutName_realloc( fp );
		tmp = countnogaplen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		if( total < 100000 )
		{
			atgcnum += countATGC( tmpseq, &nsite );
			total += nsite;
		}
		free( tmpseq );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
//	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
}


void getnumlen_nogap_outallreg( FILE *fp, int *nlenminpt )
{
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq, *tmpname;
	double atgcfreq;
#if mingw // web nomi de shiyou suru node nakutemo ii
	setmode( fileno( fp ), O_BINARY );
	setmode( fileno( stdout ), O_BINARY );
#endif
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		fprintf( stdout, "%s\n", tmpname );
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		tmp = countnogaplen( tmpseq );
		fprintf( stdout, "%d\n", tmp );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		if( total < 100000 )
		{
			atgcnum += countATGC( tmpseq, &nsite );
			total += nsite;
		}
		free( tmpseq );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
//	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
}

static void escapehtml( char *res, char *ori, int maxlen )
{
	char *res0 = res;
	while( *ori )
	{
		if( *ori == '<' ) 
		{
			strcpy( res, "&lt;" );
			res += 3;
		}
		else if( *ori == '>' ) 
		{
			strcpy( res, "&gt;" );
			res += 3;
		}
		else if( *ori == '&' ) 
		{
			strcpy( res, "&amp;" );
			res += 4;
		}
		else if( *ori == '"' ) 
		{
			strcpy( res, "&quot;" );
			res += 5;
		}
		else if( *ori == ' ' ) 
		{
			strcpy( res, "&nbsp;" );
			res += 5;
		}
		else
		{
			*res = *ori;
		}
		res++;
		ori++;

		if( res - res0 -10 > N ) break;
	}
	*res = 0;
}

void getnumlen_nogap_outallreg_web( FILE *fp, FILE *ofp, int *nlenminpt, int *isalignedpt )
{
	int total;
	int nsite = 0;
	int atgcnum;
	int alnlen = 0, alnlen_prev;
	int i, tmp, lennormalchar;
	char *tmpseq, *tmpname, *tmpname2;
	double atgcfreq;
#if mingw // web nomi de shiyou suru node nakutemo ii
	setmode( fileno( fp ), O_BINARY );
	setmode( fileno( stdout ), O_BINARY );
#endif
	tmpname = AllocateCharVec( N );
	tmpname2 = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	alnlen_prev = -1;
	*isalignedpt = 1;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpname2[0] = tmpname[0];
		escapehtml( tmpname2+1, tmpname+1, N );
//		fprintf( stdout, "%s\n", tmpname );
//		fprintf( stdout, "%s\n", tmpname2 );
//		exit(1);
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		tmp = countnogaplen( tmpseq );
//		fprintf( stdout, "%d\n", tmp );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		atgcnum += countATGC( tmpseq, &nsite );
		total += nsite;

		alnlen = strlen( tmpseq );
//		fprintf( stdout, "##### alnlen, alnlen_prev = %d, %d\n", alnlen, alnlen_prev );
		if( i>0 && alnlen_prev != alnlen ) *isalignedpt = 0;
		alnlen_prev = alnlen;

		atgcfreq = (double)atgcnum / total;
//		fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
//		if( dorp == NOTSPECIFIED ) // you kentou
		{
			if( atgcfreq > 0.75 ) 	
			{
				dorp = 'd';
				upperCase = -1;
			}
			else                  
			{
				dorp = 'p';
				upperCase = 0;
			}
		}
	
		if( dorp == 'd' ) lennormalchar = countnormalletters( tmpseq, "atgcuATGCU" );
		else              lennormalchar = countnormalletters( tmpseq, "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv" );
		free( tmpseq );

		fprintf( ofp, " <label for='s%d'><span id='ss%d'><input type='checkbox' id='s%d' name='s%d' checked></span> <input type='text' class='ll' id='ll%d' style='display:none' size='6' value='%d' readonly='readonly'>%s</label>\n", i, i, i, i, i, lennormalchar, tmpname2 );
		fprintf( ofp, "<span id='t%d-0' style='display:none'>", i );
		fprintf( ofp, " <a href='javascript:void(0)' onclick='ddcycle(this.form,\"t%d\")'>+reg</a>", i );
		fprintf( ofp, " Begin:<input type='text' name='b%d-0' size='8' value='1' class='ie'> End:<input type='text' name='e%d-0' size='8' value='%d' class='ie'>", i, i, tmp );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-0'><input type='checkbox' name='r%d-0' id='r%d-0'>Reverse</label>", i, i, i );
//		fprintf( ofp, "  Sequence Length:<input type='text' name='l%d' size='8' value='%d' readonly='readonly'>", i, tmp );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-1' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-1' size='8' value='' class='ie'> End:<input type='text' name='e%d-1' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-1'><input type='checkbox' name='r%d-1' id='r%d-1'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-2' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-2' size='8' value='' class='ie'> End:<input type='text' name='e%d-2' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-2'><input type='checkbox' name='r%d-2' id='r%d-2'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-3' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-3' size='8' value='' class='ie'> End:<input type='text' name='e%d-3' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-3'><input type='checkbox' name='r%d-3' id='r%d-3'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-4' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-4' size='8' value='' class='ie'> End:<input type='text' name='e%d-4' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-4'><input type='checkbox' name='r%d-4' id='r%d-4'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
	}
	free( tmpname );
	free( tmpname2 );
	atgcfreq = (double)atgcnum / total;
	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
//	if( dorp == NOTSPECIFIED ) // you kentou
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
	fprintf( ofp, "\n" );
	if( *isalignedpt )
	{
		fprintf( ofp, "<span id='tall-0' style='display:none'>" );
		fprintf( ofp, "Cut the alignment\n" );
		fprintf( ofp, " <a href='javascript:void(0)' onclick='ddcycle(this.form,\"tall\")'>+reg</a>" );
		fprintf( ofp, " Begin:<input type='text' name='ball-0' size='8' value='1'> End:<input type='text' name='eall-0' size='8' value='%d'>", alnlen );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-0'><input type='checkbox' name='rall-0' id='rall-0'>Reverse</label>" );
		fprintf( ofp, "  Alignment length:<input type='text' name='lall' size='8' value='%d' readonly='readonly'>", alnlen );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-1' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-1' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-1'><input type='checkbox' name='rall-1' id='rall-1'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-2' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-2' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-2'><input type='checkbox' name='rall-2' id='rall-2'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-3' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-3' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-3'><input type='checkbox' name='rall-3' id='rall-3'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-4' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-4' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-4'><input type='checkbox' name='rall-4' id='rall-4'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
	}

}

void getnumlen( FILE *fp )
{
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq;
	char *tmpname;
	double atgcfreq;

#if mingw
	setmode( fileno( fp ), O_BINARY );
	setmode( fileno( stdout ), O_BINARY );
#endif

	
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpseq = load1SeqWithoutName_realloc( fp );
		tmp = strlen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( total < 1000000 )
		{
			atgcnum += countATGC( tmpseq, &nsite );
			total += nsite;
		}
//		fprintf( stderr, "##### total = %d\n", total );
		free( tmpseq );
	}


	atgcfreq = (double)atgcnum / total;
//	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
	free( tmpname );
}
	


void WriteGapFill( FILE *fp, int locnjob, char name[][B], int nlen[M], char **aseq )
{
	static char b[N];
	int i, j;
	int nalen[M];
	static char gap[N];
	static char buff[N];

#if IODEBUG
	fprintf( stderr, "IMAKARA KAKU\n" );
#endif
	nlenmax = 0;
	for( i=0; i<locnjob; i++ )
	{
		int len = strlen( aseq[i] );
		if( nlenmax < len ) nlenmax = len;
	}

	for( i=0; i<nlenmax; i++ ) gap[i] = '-';
	gap[nlenmax] = 0;

	fprintf( fp, "%5d", locnjob );
	fprintf( fp, "\n" );

	for( i=0; i<locnjob; i++ )
	{
		strcpy( buff, aseq[i] );
		strncat( buff, gap, nlenmax-strlen( aseq[i] ) );
		buff[nlenmax] = 0;
		nalen[i] = strlen( buff );
		fprintf( fp, "%s\n", name[i] );
		fprintf( fp, "%5d\n", nalen[i] );
		for( j=0; j<nalen[i]; j=j+C )
		{
			strncpy_caseC( b, buff+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
		}
	}
#if DEBUG
	fprintf( stderr, "nalen[0] = %d\n", nalen[0] );
#endif
#if IODEBUG
	fprintf( stderr, "KAKIOWATTA\n" );
#endif
}

void writeDataforgaln( FILE *fp, int locnjob, char **name, int *nlen, char **aseq )
{
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		nalen = strlen( aseq[i] );
		fprintf( fp, ">%s\n", name[i]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[i]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[i]+j );
#endif
		}
	}
}

void writeData_pointer( FILE *fp, int locnjob, char **name, int *nlen, char **aseq )
{
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
#if DEBUG
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		fprintf( fp, ">%s\n", name[i]+1 );
#if 1
		if( LineLengthInFASTA < 0 )
			fprintf( fp, "%s\n", aseq[i] );
		else // V7.510 deha tsukawanai
		{
			nalen = strlen( aseq[i] );
			for( j=0; j<nalen; j=j+LineLengthInFASTA ) fprintf( fp, "%.*s\n", LineLengthInFASTA, aseq[i]+j );
		}
#else
		fprintf( fp, "%s\n", aseq[i] );
#endif
	}
}

void writeData( FILE *fp, int locnjob, char name[][B], int nlen[], char **aseq )
{
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
#if DEBUG
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		nalen = strlen( aseq[i] );
		fprintf( fp, ">%s\n", name[i]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[i]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[i]+j );
#endif
		}
	}
}


void write1seq( FILE *fp, char *aseq )
{
	int j;
	int nalen;

	nalen = strlen( aseq );
	for( j=0; j<nalen; j=j+C )
		fprintf( fp, "%.*s\n", C, aseq+j );
}

void readhat2_doublehalf_part_pointer( FILE *fp, int nseq, int nadd, char **name, double **mtx )
{
    int i, j, nseq0, norg;
    char b[B];

    fgets( b, B, fp );
//  fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) 
    fgets( b, B, fp ); nseq0 = atoi( b );  // 2020/Oct/23
	if( nseq != nseq0 ) 
	{
		fprintf( stderr, "%d != %d\n", nseq, nseq0 );
		ErrorExit( "hat2 is wrong." );
	}
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
	norg = nseq-nadd;
    for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ )
    {
        mtx[i][j] = ( input_new( fp, D ) );
    }
}

void readhat2_doublehalf_pointer( FILE *fp, int nseq, char **name, double **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) 
	{
		fprintf( stderr, "%d != %d\n", nseq, nseq0 );
		ErrorExit( "hat2 is wrong." );
	}
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j-i] = ( input_new( fp, D ) );
    }
}
void readhat2_doublehalf( FILE *fp, int nseq, char name[M][B], double **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j-i] = ( input_new( fp, D ) );
    }
}
void readhat2_double( FILE *fp, int nseq, char name[M][B], double **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = ( input_new( fp, D ) );
    }
}
void readhat2_int( FILE *fp, int nseq, char name[M][B], int **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = (int)( input_new( fp, D ) * INTMTXSCALE + 0.5 );
    }
}

void readhat2_pointer( FILE *fp, int nseq, char **name, double **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = (double)input_new( fp, D);
    }
}
void readhat2( FILE *fp, int nseq, char name[M][B], double **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = (double)input_new( fp, D);
    }
}

void WriteFloatHat2_pointer_halfmtx( FILE *hat2p, int locnjob, char **name, double **mtx )
{
	int i, j, ijsa;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=1; j<locnjob-i; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob; i++ )
	{
		for( j=i+1; j<njob; j++ )
		{
			fprintf( hat2p, DFORMAT, mtx[i][j-i] );
			ijsa = j-i;
			if( ijsa % 12 == 0 || ijsa == locnjob-i-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteFloatHat2_pointer( FILE *hat2p, int locnjob, char **name, double **mtx )
{
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=1; j<locnjob-i; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob; i++ )
	{
		for( j=1; j<locnjob-i; j++ ) 
		{
			fprintf( hat2p, DFORMAT, mtx[i][j] );
			if( j % 12 == 0 || j == locnjob-i-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteFloatHat2( FILE *hat2p, int locnjob, char name[M][B], double **mtx )
{
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=1; j<locnjob-i; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob; i++ )
	{
		for( j=1; j<locnjob-i; j++ ) 
		{
			fprintf( hat2p, DFORMAT, mtx[i][j] );
			if( j % 12 == 0 || j == locnjob-i-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2_int( FILE *hat2p, int locnjob, char name[M][B], int **mtx )
{
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];
	max /= INTMTXSCALE;

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, DFORMAT, (double)mtx[i][j] / INTMTXSCALE );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2_part_pointer( FILE *hat2p, int locnjob, int nadd, char **name, double **mtx )
{
	int i, j;
	int norg = locnjob-nadd;
	double max = 0.0;
//	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<norg; i++ )
	{
		for( j=0; j<nadd; j++ ) 
		{
			fprintf( hat2p, DFORMAT, mtx[i][j] );
			if( (j+1) % 12 == 0 || j == nadd-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2_pointer( FILE *hat2p, int locnjob, char **name, double **mtx )
{
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, DFORMAT, mtx[i][j] );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2( FILE *hat2p, int locnjob, char name[M][B], double **mtx )
{
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, DFORMAT, mtx[i][j] );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

#if 0
void WriteHat2plain( FILE *hat2p, int locnjob, double **mtx )
{
	int i, j, ilim;

	ilim = locnjob-1;
	for( i=0; i<ilim; i++ )
	{
		fprintf( hat2p, "%d-%d d=%.3f\n", i+1, i+1, 0.0 );
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, "%d-%d d=%.3f\n", i+1, j+1, mtx[i][j] );
		}
	}
}
#endif

int ReadFasta_sub( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
            count++;
        }
    }

	for( i=0; i<nseq; i++ ) dis[i] = 0.0;
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
		if( fgets( b, B-1, fp ) ) break;
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            fgets( b, B-1, fp );
            dis[junban[count]] = atof( b );
            count++;
        }
    }
    return 0;
}


int ReadSsearch( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];
	int opt;

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
			sscanf( b+75, "%d", &opt ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }

/*
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            dis[junban[count]] = atof( b+65 );
            count++;
        }
    }
*/
    return 0;
}

int ReadBlastm7_avscore( FILE *fp, double *dis, int nin )
{
    int count=0;
    char b[B];
	char *pt;
    int *junban;
	double score, sumscore;
	double len, sumlen;
	int qstart, qend, tstart, tend;
	double scorepersite;
	static char qal[N], tal[N], al[N];
	int nlocalhom;

	junban = calloc( nin, sizeof( int ) );

	count = 0;
	sumscore = 0.0;
	sumlen = 0.0;
	score = 0.0;
	len = 0.0;
	scorepersite = 0.0; // by D.Mathog, a guess
    while( 1 )
	{

		if( feof( fp ) ) break;

		while( fgets( b, B-1, fp ) )
		{
			if( !strncmp( "          <Hit_def>", b, 19 ) || !strncmp( "              <Hsp_num>", b, 23 ) ) break;
		}

		if( !strncmp( "          <Hit_def>", b, 19 ) )
		{
			junban[count] = atoi( b+31 );
			nlocalhom = 0;
		}


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_score>", b, 25 ) ) break;
		pt = b + 25;
		score = atof( pt );
		sumscore += score;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-from>", b, 30 ) ) break;
		pt = b + 30;
		qstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-to>", b, 28 ) ) break;
		pt = b + 28;
		qend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-from>", b, 28 ) ) break;
		pt = b + 28;
		tstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-to>", b, 26 ) ) break;
		pt = b + 26;
		tend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_align-len>", b, 29 ) ) break;
		pt = b + 29;
		len = atoi( pt );
		sumlen += len;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_qseq>", al, 24 ) ) break;

		strcpy( qal, al+24 );
		pt = qal;
		while( *++pt != '<' )
			;
		*pt = 0;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_hseq>", al, 24 ) ) break;

		strcpy( tal, al+24 );
		pt = tal;
		while( *++pt != '<' )
			;
		*pt = 0;


//		fprintf( stderr, "t=%d, score = %f, qstart=%d, qend=%d, tstart=%d, tend=%d, overlapaa=%d\n", junban[count], score, qstart, qend, tstart, tend, overlapaa );


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "            </Hsp>:", b, 18 ) ) break;


		fgets( b, B-1, fp );


		if( !strncmp( "          </Hit_hsps>", b, 21 ) )
		{
			dis[junban[count++]] = sumscore;
			sumscore = 0.0;
			fgets( b, B-1, fp );
			fgets( b, B-1, fp );
			scorepersite = sumscore / sumlen;
			if( scorepersite != (int)scorepersite )
			{
				fprintf( stderr, "ERROR! sumscore=%f, sumlen=%f, and scorepersite=%f\n", sumscore, sumlen, scorepersite );
				exit( 1 );
			}

			if( !strncmp( "      </Iteration_hits>", b, 23 ) ) break;
		}
	}

	free( junban );

    return (int)scorepersite;
}
int ReadBlastm7_scoreonly( FILE *fp, double *dis, int nin )
{
    int count=0;
    char b[B];
	char *pt;
    int *junban;
	int overlapaa;
	double score, sumscore;
	int qstart, qend, tstart, tend;
	static char qal[N], tal[N], al[N];
	int nlocalhom;

	junban = calloc( nin, sizeof( int ) );

	count = 0;
	sumscore = 0.0;
	score = 0.0;
    while( 1 )
	{

		if( feof( fp ) ) break;

		while( fgets( b, B-1, fp ) )
		{
			if( !strncmp( "          <Hit_def>", b, 19 ) || !strncmp( "              <Hsp_num>", b, 23 ) ) break;
		}

		if( !strncmp( "          <Hit_def>", b, 19 ) )
		{
			junban[count] = atoi( b+31 );
			nlocalhom = 0;
		}


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_score>", b, 25 ) ) break;
		pt = b + 25;
		score = atof( pt );
		sumscore += score;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-from>", b, 30 ) ) break;
		pt = b + 30;
		qstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-to>", b, 28 ) ) break;
		pt = b + 28;
		qend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-from>", b, 28 ) ) break;
		pt = b + 28;
		tstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-to>", b, 26 ) ) break;
		pt = b + 26;
		tend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_align-len>", b, 29 ) ) break;
		pt = b + 29;
		overlapaa = atoi( pt );


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_qseq>", al, 24 ) ) break;

		strcpy( qal, al+24 );
		pt = qal;
		while( *++pt != '<' )
			;
		*pt = 0;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_hseq>", al, 24 ) ) break;

		strcpy( tal, al+24 );
		pt = tal;
		while( *++pt != '<' )
			;
		*pt = 0;


//		fprintf( stderr, "t=%d, score = %f, qstart=%d, qend=%d, tstart=%d, tend=%d, overlapaa=%d\n", junban[count], score, qstart, qend, tstart, tend, overlapaa );

//		nlocalhom += addlocalhom_r( qal, tal, localhomlist+junban[count], qstart, tstart, score, overlapaa, nlocalhom );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "            </Hsp>:", b, 18 ) ) break;


		fgets( b, B-1, fp );


		if( !strncmp( "          </Hit_hsps>", b, 21 ) )
		{
			dis[junban[count++]] = sumscore;
			sumscore = 0.0;
			fgets( b, B-1, fp );
			fgets( b, B-1, fp );
			if( !strncmp( "      </Iteration_hits>", b, 23 ) ) break;
		}
	}

	free( junban );

    return count;
}

int ReadBlastm7( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	double score, sumscore;
	int qstart, qend, tstart, tend;
	static char qal[N], tal[N], al[N];
	int nlocalhom;



	count = 0;
	sumscore = 0.0;
	score = 0.0;
	nlocalhom = 0;
    while( 1 )
	{

		if( feof( fp ) ) break;

		while( fgets( b, B-1, fp ) )
		{
			if( !strncmp( "          <Hit_def>", b, 19 ) || !strncmp( "              <Hsp_num>", b, 23 ) ) break;
		}

		if( !strncmp( "          <Hit_def>", b, 19 ) )
		{
			junban[count] = atoi( b+31 );
			nlocalhom = 0;
		}


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_score>", b, 25 ) ) break;
		pt = b + 25;
		score = atof( pt );
		sumscore += score;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-from>", b, 30 ) ) break;
		pt = b + 30;
		qstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-to>", b, 28 ) ) break;
		pt = b + 28;
		qend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-from>", b, 28 ) ) break;
		pt = b + 28;
		tstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-to>", b, 26 ) ) break;
		pt = b + 26;
		tend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_align-len>", b, 29 ) ) break;
		pt = b + 29;
		overlapaa = atoi( pt );


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_qseq>", al, 24 ) ) break;

		strcpy( qal, al+24 );
		pt = qal;
		while( *++pt != '<' )
			;
		*pt = 0;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_hseq>", al, 24 ) ) break;

		strcpy( tal, al+24 );
		pt = tal;
		while( *++pt != '<' )
			;
		*pt = 0;


//		fprintf( stderr, "t=%d, score = %f, qstart=%d, qend=%d, tstart=%d, tend=%d, overlapaa=%d\n", junban[count], score, qstart, qend, tstart, tend, overlapaa );

		nlocalhom += addlocalhom_r( qal, tal, localhomlist+junban[count], qstart, tstart, score, overlapaa, nlocalhom, 'h' );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "            </Hsp>:", b, 18 ) ) break;


		fgets( b, B-1, fp );


		if( !strncmp( "          </Hit_hsps>", b, 21 ) )
		{
			dis[junban[count++]] = sumscore;
			sumscore = 0.0;
			fgets( b, B-1, fp );
			fgets( b, B-1, fp );
			if( !strncmp( "      </Iteration_hits>", b, 23 ) ) break;
		}
	}
    return count;
}

int ReadFasta34noalign( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int opt;
	double z, bits;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
    }

    return count;
}
int ReadFasta34m10_nuc( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	int opt, qstart, qend, tstart, tend;
	double z, bits;
	int qal_display_start, tal_display_start;
	static char qal[N], tal[N];
	char *qal2, *tal2;
	int c;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			if( strchr( b, 'r' ) ) continue;

			pt = strchr( b, ']' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
		else if( 0 == strncmp( ">>+==========+", b, 14 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

	count = 0;
    while( 1 )
	{
		if( strncmp( ">>+==========+", b, 14 ) )
		{
			fgets( b, B-1, fp );
			if( feof( fp ) ) break;
			continue;
		}
		junban[count++] = atoi( b+14 );
//		fprintf( stderr, "t = %d\n", atoi( b+14 ) );
		while( fgets( b, B-1, fp ) )
			if( !strncmp( "; fa_opt:", b, 9 ) || !strncmp( "; sw_s-w opt:", b, 13 ) ) break;
		pt = strstr( b, ":" ) +1;
		opt = atoi( pt );


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_overlap:", b+4, 9 ) ) break;
		pt = strstr( b, ":" ) +1;
		overlapaa = atoi( pt );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) +1;
		qstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) +1;
		qend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) +1;
		qal_display_start = atoi( pt ) - 1;

		pt = qal;
		while( (c = fgetc( fp )) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tal_display_start = atoi( pt ) - 1;

		pt = tal;
		while( ( c = fgetc( fp ) ) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

//		fprintf( stderr, "(%d-%d:%d-%d)\n", qstart, qend, tstart, tend );
//		fprintf( stderr, "qal_display_start = %d, tal_display_start = %d\n", qal_display_start, tal_display_start );

//		fprintf( stderr, "qal = %s\n", qal );
//		fprintf( stderr, "tal = %s\n", tal );

		qal2 = cutal( qal, qal_display_start, qstart, qend );
		tal2 = cutal( tal, tal_display_start, tstart, tend );

//		fprintf( stderr, "qal2 = %s\n", qal2 );
//		fprintf( stderr, "tal2 = %s\n", tal2 );

//		fprintf( stderr, "putting   %d - %d, opt = %d\n", qmem, junban[count-1], opt );
		putlocalhom2( qal2, tal2, localhomlist+junban[count-1], qstart, tstart, opt, overlapaa, 'h' );
	}
//	fprintf( stderr, "count = %d\n", count );
    return count;
}
int ReadFasta34m10( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	int opt, qstart, qend, tstart, tend;
	double z, bits;
	int qal_display_start, tal_display_start;
	static char qal[N], tal[N];
	char *qal2, *tal2;
	int c;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
		else if( 0 == strncmp( ">>+==========+", b, 14 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

	count = 0;
    while( 1 )
	{
		if( strncmp( ">>+==========+", b, 14 ) )
		{
			fgets( b, B-1, fp );
			if( feof( fp ) ) break;
			continue;
		}
		junban[count++] = atoi( b+14 );
//		fprintf( stderr, "t = %d\n", atoi( b+14 ) );
		while( fgets( b, B-1, fp ) )
			if( !strncmp( "; fa_opt:", b, 9 ) || !strncmp( "; sw_s-w opt:", b, 13 ) ) break;
		pt = strstr( b, ":" ) +1;
		opt = atoi( pt );


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_overlap:", b+4, 9 ) ) break;
		pt = strstr( b, ":" ) +1;
		overlapaa = atoi( pt );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) +1;
		qstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) +1;
		qend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) +1;
		qal_display_start = atoi( pt ) - 1;

		pt = qal;
		while( (c = fgetc( fp )) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tal_display_start = atoi( pt ) - 1;

		pt = tal;
		while( ( c = fgetc( fp ) ) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

//		fprintf( stderr, "(%d-%d:%d-%d)\n", qstart, qend, tstart, tend );
//		fprintf( stderr, "qal_display_start = %d, tal_display_start = %d\n", qal_display_start, tal_display_start );

//		fprintf( stderr, "qal = %s\n", qal );
//		fprintf( stderr, "tal = %s\n", tal );

		qal2 = cutal( qal, qal_display_start, qstart, qend );
		tal2 = cutal( tal, tal_display_start, tstart, tend );

//		fprintf( stderr, "qal2 = %s\n", qal2 );
//		fprintf( stderr, "tal2 = %s\n", tal2 );

//		fprintf( stderr, "putting   %d - %d, opt = %d\n", qmem, junban[count-1], opt );
		putlocalhom2( qal2, tal2, localhomlist+junban[count-1], qstart, tstart, opt, overlapaa, 'h' );
	}
//	fprintf( stderr, "count = %d\n", count );
    return count;
}
int ReadFasta34m10_scoreonly_nucbk( FILE *fp, double *dis, int nin )
{
    int count=0;
    char b[B];
	char *pt;
    int pos;
	int opt;
	double z, bits;

    count = 0;
    while( !feof( fp ) )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+===========+", b, 13 ) )
        {
            pos = atoi( b+13 );

			if( strchr( b, 'r' ) ) continue;

//			pt = strchr( b, ')' ) + 1;
			pt = strchr( b, ']' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[pos] += (double)opt;
            count++;
#if 0
			fprintf( stderr, "b=%s\n", b );
			fprintf( stderr, "opt=%d\n", opt );
			fprintf( stderr, "pos=%d\n", pos );
			fprintf( stderr, "dis[pos]=%f\n", dis[pos] );
#endif

        }
		else if( 0 == strncmp( ">>><<<", b, 6 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

    return count;
}

int ReadFasta34m10_scoreonly_nuc( FILE *fp, double *dis, int nin )
{
    int count=0;
    char b[B];
	char *pt;
    int pos;
	int opt;
	double z, bits;
	int c;
	int *yonda;


	yonda = AllocateIntVec( nin );
	for( c=0; c<nin; c++ ) yonda[c] = 0;
	for( c=0; c<nin; c++ ) dis[c] = 0.0;

    count = 0;
    while( !feof( fp ) )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+===========+", b, 13 ) )
        {
            pos = atoi( b+13 );

			if( strchr( b, 'r' ) ) continue;

//			pt = strchr( b, ')' ) + 1;
			pt = strchr( b, ']' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
			if( yonda[pos] == 0 )
			{
	            dis[pos] += (double)opt;
				yonda[pos] = 1;
			}
            count++;
#if 0
			fprintf( stderr, "b=%s\n", b );
			fprintf( stderr, "opt=%d\n", opt );
			fprintf( stderr, "pos=%d\n", pos );
			fprintf( stderr, "dis[pos]=%f\n", dis[pos] );
#endif

        }
        else if( !strncmp( ">>>", b, 3 ) )
		{
			for( c=0; c<nin; c++ ) yonda[c] = 0;
		}
		else if( 0 == strncmp( ">>><<<", b, 6 ) )
		{
			break;
		}

    }

	free( yonda );

	if( !count ) return -1;

    return count;
}

int ReadFasta34m10_scoreonly( FILE *fp, double *dis, int nin )
{
    int count=0;
    char b[B];
	char *pt;
    int pos;
	int opt;
	double z, bits;
	int c;
	int *yonda;


	yonda = AllocateIntVec( nin );
	for( c=0; c<nin; c++ ) yonda[c] = 0;
	for( c=0; c<nin; c++ ) dis[c] = 0.0;

    count = 0;
    while( !feof( fp ) )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+===========+", b, 13 ) )
        {
            pos = atoi( b+13 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
			if( yonda[pos] == 0 )
			{
	            dis[pos] += (double)opt;
				yonda[pos] = 1;
			}
            count++;
#if 0
			fprintf( stderr, "b=%s\n", b );
			fprintf( stderr, "opt=%d\n", opt );
			fprintf( stderr, "pos=%d\n", pos );
			fprintf( stderr, "dis[pos]=%f\n", dis[pos] );
#endif

        }
        else if( !strncmp( ">>>", b, 3 ) )
		{
			for( c=0; c<nin; c++ ) yonda[c] = 0;
		}
		else if( 0 == strncmp( ">>><<<", b, 6 ) )
		{
			break;
		}

    }

	free( yonda );

	if( !count ) return -1;

    return count;
}
int ReadFasta34( FILE *fp, double *dis, int nseq, char name[M][B], LocalHom *localhomlist )
{
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	int opt, qstart, qend, tstart, tend;
	double z, bits;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
		else if( 0 == strncmp( ">>+==========+", b, 14 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

	count = 0;
    while( !feof( fp ) )
	{
		if( !strncmp(">>+==========+", b, 14 ) )
		{
            junban[count] = atoi( b+14 );
            count++;
        	fgets( b, B-1, fp ); // initn:
			pt = strstr( b, "opt: " ) + 5;
			localhomlist[junban[count-1]].opt = atof( pt );
        	fgets( b, B-1, fp ); // Smith-Waterman score
			pt = strstr( b, "ungapped) in " ) + 13;
			sscanf( pt, "%d", &overlapaa ); 
			fprintf( stderr, "pt = %s, overlapaa = %d\n", pt, overlapaa );
			pt = strstr( b, "overlap (" ) + 8;
			sscanf( pt, "(%d-%d:%d-%d)", &qstart, &qend, &tstart, &tend ); 
			localhomlist[junban[count-1]].overlapaa = overlapaa;
			localhomlist[junban[count-1]].start1 = qstart-1;
			localhomlist[junban[count-1]].end1   = qend-1;
			localhomlist[junban[count-1]].start2 = tstart-1;
			localhomlist[junban[count-1]].end2   = tend-1;
		}
        fgets( b, B-1, fp );
	}
	fprintf( stderr, "count = %d\n", count );
    return count;
}

int ReadFasta3( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int count=0;
    char b[B];
	char *pt;
    int junban[M];
	int initn, init1, opt;
	double z;

    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %d %d %lf", &initn, &init1, &opt, &z ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }
    return 0;
}

int ReadFasta( FILE *fp, double *dis, int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];
	int initn, init1, opt;

    count = 0;
	for( i=0; i<nseq; i++ ) dis[i] = 0.0;
    for( i=0; !feof( fp ) && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			sscanf( b+50, "%d %d %d", &initn, &init1, &opt ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }

/*
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            dis[junban[count]] = atof( b+65 );
            count++;
        }
    }
*/
    return 0;
}


int ReadOpt( FILE *fp, int opt[M], int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];
	int optt, initn, init1;

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
			sscanf( b+50, "%d %d %d", &initn, &init1, &optt ); 
            opt[junban[count]] = (double)optt;
            count++;
        }
    }
    return 0;
}

int ReadOpt2( FILE *fp, int opt[M], int nseq, char name[M][B] )
{
    int i, count=0;
    char b[B];
    int junban[M];

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
            opt[junban[count]] = atoi( b+65 );
            count++;
        }
    }
    return 0;
}



int writePre( int nseq, char **name, int nlen[M], char **aseq, int force )
{
#if USE_XCED
	int i, value;
	if( !signalSM )
	{
		if( force ) 
		{
			rewind( prep_g );
			if( devide ) dvWrite( prep_g, nseq, name, nlen, aseq );
#if 0
			else    WriteGapFill( prep_g, nseq, name, nlen, aseq );
#else
			else    writeData( prep_g, nseq, name, nlen, aseq );
#endif
		}
		return( 0 );
	}
	for( i=0; i<10; i++ )
	{
#if IODEBUG
		fprintf( stderr, "SEMAPHORE = %d\n", signalSM[SEMAPHORE] );
#endif
		if( signalSM[SEMAPHORE]-- > 0 )
		{
#if 0 /* /tmp/pre の関係ではずした */
			if( ferror( prep_g ) ) prep_g = fopen( "pre", "w" );
			if( !prep_g ) ErrorExit( "Cannot re-open pre." ); 
#endif
			rewind( prep_g );
			signalSM[STATUS] = IMA_KAITERU;
#if IODEBUG
			if( force ) fprintf( stderr, "FINAL " );
#endif
			if( devide ) dvWrite( prep_g, nseq, name, nlen, aseq );
			else    WriteGapFill( prep_g, nseq, name, nlen, aseq );
			/*
			fprintf( prep_g, '\EOF' );
			*/
			fflush( prep_g );
			if( force ) signalSM[STATUS] = OSHIMAI;
			else        signalSM[STATUS] = KAKIOWATTA;
			value = 1;
			signalSM[SEMAPHORE]++;
#if IODEBUG
			fprintf( stderr, "signalSM[STATUS] = %c\n", signalSM[STATUS] );
#endif
			break;
		}
		else
		{
#if IODEBUG
			fprintf( stderr, "YONDERUKARA_AKIRAMERU\n" );
#endif
			value = 0;
			signalSM[SEMAPHORE]++;
			if( !force ) break;
#if IODEBUG
			fprintf( stderr, "MATSU\n" );
#endif
			sleep( 1 );
		}
	}
	if( force && !value ) ErrorExit( "xced ga pre wo hanasanai \n" );
	return( value );
#else
	if( force ) 
	{
		rewind( prep_g );
		writeData_pointer( prep_g, nseq, name, nlen, aseq );
	}
#endif
	return( 0 );
}


void readOtherOptions( int *ppidptr, int *fftThresholdptr, int *fftWinSizeptr )
{
	if( calledByXced )
	{
		FILE *fp = fopen( "pre", "r" );
		char b[B];
		if( !fp ) ErrorExit( "Cannot open pre.\n" );
		fgets( b, B-1, fp );
		sscanf( b, "%d %d %d", ppidptr, fftThresholdptr, fftWinSizeptr );
		fclose( fp );
#if IODEBUG
	fprintf( stderr, "b = %s\n", b );
	fprintf( stderr, "ppid = %d\n", ppid );
	fprintf( stderr, "fftThreshold = %d\n", fftThreshold );
	fprintf( stderr, "fftWinSize = %d\n", fftWinSize );
#endif
	}
	else
	{
		*ppidptr = 0;
		*fftThresholdptr = FFT_THRESHOLD;
		if( dorp == 'd' )
			*fftWinSizeptr = FFT_WINSIZE_D;
		else
			*fftWinSizeptr = FFT_WINSIZE_P;
	}
#if 0
	fprintf( stderr, "fftThresholdptr=%d\n", *fftThresholdptr );
	fprintf( stderr, "fftWinSizeptr=%d\n", *fftWinSizeptr );
#endif
}

void initSignalSM( void )
{
//	int signalsmid;

#if IODEBUG
	if( ppid ) fprintf( stderr, "PID of xced = %d\n", ppid );
#endif
	if( !ppid )
	{
		signalSM = NULL;
		return;
	}

#if 0
	signalsmid = shmget( (key_t)ppid, 3, IPC_ALLOC | 0666 );
	if( signalsmid == -1 ) ErrorExit( "Cannot get Shared memory for signal.\n" );
	signalSM = shmat( signalsmid, 0, 0 );
	if( (int)signalSM == -1 ) ErrorExit( "Cannot attatch Shared Memory for signal!\n" );
	signalSM[STATUS] = IMA_KAITERU;
	signalSM[SEMAPHORE] = 1;
#endif
}

void initFiles( void )
{
	char pname[100];
	if( ppid )
		sprintf( pname, "/tmp/pre.%d", ppid );
	else
		sprintf( pname, "pre" );
	prep_g = fopen( pname, "w" );
	if( !prep_g ) ErrorExit( "Cannot open pre" );

#if mingw
	setmode( fileno( prep_g ), O_BINARY );
#endif
	trap_g = fopen( "trace", "w" );
	if( !trap_g ) ErrorExit( "cannot open trace" );
	fprintf( trap_g, "PID = %d\n", getpid() );
	fflush( trap_g );
}

void closeFiles( void )
{
	fclose( prep_g );
	fclose( trap_g );
}


void WriteForFasta( FILE *fp, int locnjob, char **name, int nlen[M], char **aseq )
{
    static char b[N];
    int i, j;
    int nalen[M];

    for( i=0; i<locnjob; i++ )
    {
        nalen[i] = strlen( aseq[i] );
        fprintf( fp, ">%s\n", name[i] );
        for( j=0; j<nalen[i]; j=j+C ) 
        {
            strncpy( b, aseq[i]+j, C ); b[C] = 0;
            fprintf( fp, "%s\n",b );
        }
    }
}

void readlocalhomtable2_target( FILE*fp, int njob, LocalHom **localhomtable, char *kozoarivec, int *targetmap )
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2, it, jt;
	LocalHom *tmpptr1, *tmpptr2;

//	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif
		it = targetmap[i];
		if( it == -1 )
		{
			reporterr( "hat3 ga okashii.  _target_ \n" );
			exit( 1 );
		}
		jt = targetmap[j];



//		if( i < j )
		{
			if( localhomtable[it][j].nokori++ > 0 )
			{
				tmpptr1 = localhomtable[it][j].last;
//				fprintf( stderr, "reallocating, localhomtable[%d][%d].nokori = %d\n", i, j, localhomtable[i][j].nokori );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->extended = -1;
				tmpptr1->next = NULL;
				localhomtable[it][j].last = tmpptr1;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j].nokori );
			}
			else
			{
				tmpptr1 = localhomtable[it]+j;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j].nokori );
			}
	
			tmpptr1->start1 = start1;
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1;
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
		}
//		else
		if( jt != -1 )
		{
			if( localhomtable[jt][i].nokori++ > 0 )
			{
				tmpptr2 = localhomtable[jt][i].last;
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->extended = -1;
				tmpptr2->next = NULL;
				localhomtable[jt][i].last = tmpptr2;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", j, i, localhomtable[j][i].nokori );
			}
			else
			{
				tmpptr2 = localhomtable[jt]+i;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", j, i, localhomtable[j][i].nokori );
			}
	
			tmpptr2->start2 = start1;
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1;
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, st1=%d, en1=%d, opt = %f\n", i, j, tmpptr1->start1, tmpptr1->end1, opt );
		}

	}
}

void readlocalhomtable2_half( FILE*fp, int njob, LocalHom **localhomtable, char *kozoarivec )
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	LocalHom *tmpptr1;

//	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif

		if( j <= i || i >= njob || i >= njob )
		{
			reporterr( "Check hat3.  The first sequence must be younger than the second one.\n" );
			exit( 1 );
		}
		{
			if( localhomtable[i][j-i].nokori++ > 0 )
			{
				tmpptr1 = localhomtable[i][j-i].last;
//				fprintf( stderr, "reallocating, localhomtable[%d][%d].nokori = %d\n", i, j, localhomtable[i][j].nokori );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->extended = -1;
				tmpptr1->next = NULL;
				localhomtable[i][j-i].last = tmpptr1;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j-i].nokori );
			}
			else
			{
				tmpptr1 = localhomtable[i]+j-i;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j-i].nokori );
			}
	
			tmpptr1->start1 = start1;
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1;
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
		}
	}
}

void readlocalhomtable2( FILE*fp, int njob, LocalHom **localhomtable, char *kozoarivec )
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	LocalHom *tmpptr1, *tmpptr2;

//	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif

//		if( i < j )
		{
			if( localhomtable[i][j].nokori++ > 0 )
			{
				tmpptr1 = localhomtable[i][j].last;
//				fprintf( stderr, "reallocating, localhomtable[%d][%d].nokori = %d\n", i, j, localhomtable[i][j].nokori );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->extended = -1;
				tmpptr1->next = NULL;
				localhomtable[i][j].last = tmpptr1;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j].nokori );
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j].nokori );
			}
	
			tmpptr1->start1 = start1;
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1;
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
		}
//		else
		{
			if( localhomtable[j][i].nokori++ > 0 )
			{
				tmpptr2 = localhomtable[j][i].last;
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->extended = -1;
				tmpptr2->next = NULL;
				localhomtable[j][i].last = tmpptr2;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", j, i, localhomtable[j][i].nokori );
			}
			else
			{
				tmpptr2 = localhomtable[j]+i;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", j, i, localhomtable[j][i].nokori );
			}
	
			tmpptr2->start2 = start1;
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1;
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, st1=%d, en1=%d, opt = %f\n", i, j, tmpptr1->start1, tmpptr1->end1, opt );
		}

	}
}

#if 0
void readlocalhomtable_target( FILE*fp, int ntarget, int njob, LocalHom **localhomtable, char *kozoarivec, int *targetmap )
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2, it, jt;
	int **nlocalhom = NULL;
	LocalHom *tmpptr1=NULL, *tmpptr2=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( njob, njob );
	for( i=0; i<ntarget; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif

		printf( "reading %d-%d\n", i, j );

		it = targetmap[i];
		if( it == -1 )
		{
			reporterr( "hat3 ga okashii.  _target_ \n" );
			exit( 1 );
		}
		jt = targetmap[j];

//		if( i < j )
		{
			if( nlocalhom[it][j]++ > 0 )
			{
				printf( "extending %d-%d, ->%d\n", i, j, nlocalhom[it][j] );
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[it]+j;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );


		}
//		else

		if( jt != -1 )
		{
			if( nlocalhom[jt][i]++ > 0 )
			{
				printf( "extending %d-%d, ->%d\n", i, j, nlocalhom[jt][i] );
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->next = NULL;
			}
			else
				tmpptr2 = localhomtable[jt]+i;
	
			tmpptr2->start2 = start1; // CHUUI!!!!
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1; // CHUUI!!!!
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;

//			fprintf( stderr, "j=%d, i=%d, opt = %f\n", j, i, opt );
		}

	}
	LocalHom *tmpptr;
	for( tmpptr = localhomtable[1]+11; tmpptr; tmpptr=tmpptr->next )
		fprintf( stdout, "reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f, next=%p\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt / 600 * 5.8, tmpptr->next );
	FreeIntMtx( nlocalhom );
}

void readlocalhomtable_half( FILE*fp, int njob, LocalHom **localhomtable, char *kozoarivec )
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	int **nlocalhom = NULL;
	LocalHom *tmpptr1=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( njob, njob );
	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif


		if( j <= i )
		{
			reporterr( "Check hat3.  The first sequence must be younger than the second one.\n" );
			exit( 1 );
		}
		{
			if( nlocalhom[i][j]++ > 0 )
			{
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[i]+j-i;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );
		}
	}
	FreeIntMtx( nlocalhom );
}
#endif

void readlocalhomtable( FILE*fp, int njob, LocalHom **localhomtable, char *kozoarivec )
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	int **nlocalhom = NULL;
	LocalHom *tmpptr1=NULL, *tmpptr2=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( njob, njob );
	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif


		if( j <= i )
		{
			reporterr( "Check hat3.  The first sequence must be younger than the second one.\n" );
			exit( 1 );
		}
		{
			if( nlocalhom[i][j]++ > 0 )
			{
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );
		}
//		else
		{
			if( nlocalhom[j][i]++ > 0 )
			{
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->next = NULL;
			}
			else
				tmpptr2 = localhomtable[j]+i;
	
			tmpptr2->start2 = start1; // CHUUI!!!!
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1; // CHUUI!!!!
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;

//			fprintf( stderr, "j=%d, i=%d, opt = %f\n", j, i, opt );
		}

	}
	FreeIntMtx( nlocalhom );
}


void readlocalhomtable_two( FILE*fp, int norg, int nadd, LocalHom **localhomtable, LocalHom **localhomtablex, char *kozoarivec ) // for test only
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	int **nlocalhom = NULL;
	int **nlocalhomx = NULL;
	LocalHom *tmpptr1=NULL, *tmpptr2=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( norg, nadd );
	for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ ) nlocalhom[i][j] = 0;
	nlocalhomx = AllocateIntMtx( nadd, norg );
	for( i=0; i<nadd; i++ ) for( j=0; j<norg; j++ ) nlocalhomx[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) 
		{
			fprintf( stderr, "Not supported!\n" );
			exit( 1 );
		}
		j -= norg;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif


//		if( i < j )
		{
			if( nlocalhom[i][j]++ > 0 )
			{
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );
		}

		{
			if( nlocalhomx[j][i]++ > 0 )
			{
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->next = NULL;
			}
			else
				tmpptr2 = localhomtablex[j]+i;
	
			tmpptr2->start2 = start1+1; // CHUUI!!!!
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1+1; // CHUUI!!!!
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;

//			fprintf( stderr, "j=%d, i=%d, opt = %f\n", j, i, opt );
		}

	}
	FreeIntMtx( nlocalhom );
	FreeIntMtx( nlocalhomx );
}

void readlocalhomtable_one( FILE*fp, int norg, int nadd, LocalHom **localhomtable, char *kozoarivec ) // for test only
{
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	int **nlocalhom = NULL;
	LocalHom *tmpptr1=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( norg, nadd );
	for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) 
		{
			fprintf( stderr, "Not supported!\n" );
			exit( 1 );
		}
		j -= norg;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif


//		if( i < j )
		{
			if( nlocalhom[i][j]++ > 0 )
			{
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );
		}

	}
	FreeIntMtx( nlocalhom );
}

void outlocalhom_part( LocalHom **localhom, int norg, int nadd )
{
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ )
	{
		tmpptr = localhom[i]+j;
		fprintf( stdout, "%d-%d\n", i, j+norg );
		do
		{
			fprintf( stdout, "reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt / 600 * 5.8 );
		}
		while( (tmpptr=tmpptr->next) );
	}
}

void outlocalhom_target( LocalHom **localhom, int norg, int nadd )
{
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ )
	{
		tmpptr = localhom[i]+j;
		fprintf( stdout, "%d-%d\n", i, j );
		for( ; tmpptr; tmpptr=tmpptr->next )
		{
			fprintf( stdout, "reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f, next=%p\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt / 600 * 5.8, (void *)tmpptr->next );
		}
//		while( (tmpptr=tmpptr->next) );
	}
}

void outlocalhom_half( LocalHom **localhom, int nseq )
{
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j-i;
		fprintf( stdout, "%d-%d\n", i, j );
		do
		{
			fprintf( stdout, "reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f, next=%p\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt / 600 * 5.8, (void *)tmpptr->next );
		}
		while( (tmpptr=tmpptr->next) );
	}
}

void outlocalhom( LocalHom **localhom, int nseq )
{
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j;
		fprintf( stderr, "%d-%d\n", i, j );
		do
		{
			fprintf( stderr, "reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt );
		}
		while( (tmpptr=tmpptr->next) );
	}
}

void outlocalhompt( LocalHom ***localhom, int n1, int n2 )
{
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		tmpptr = localhom[i][j];
//		fprintf( stdout, "%d-%d\n", i, j );
		do
		{
			fprintf( stdout, "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f\n", i, j, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt );
		}
		while( (tmpptr=tmpptr->next) );
	}
}

void initlocalhom1( LocalHom *lh )
{
	lh->start1 = -1;
	lh->end1 = -1;
	lh->start2 = -1; 
	lh->end2 = -1; 
	lh->opt = -1.0;
	lh->next = NULL;
	lh->nokori = 0;
	lh->extended = -1;
	lh->last = lh;
	lh->korh = 'h';
}

void freelocalhom1( LocalHom *lh )
{
	if( lh == NULL ) return;
	LocalHom *tmpptr = lh;
	LocalHom *ppp;
	for( ; tmpptr; tmpptr=ppp )
	{
		ppp = tmpptr->next;
		if( tmpptr!=lh ) 
		{
			free( tmpptr );
			continue;
		}
		tmpptr->start1 = -1;
		tmpptr->end1 = -1;
		tmpptr->start2 = -1; 
		tmpptr->end2 = -1; 
		tmpptr->opt = -1.0;
		tmpptr->next = NULL;
		tmpptr->nokori = 0;
		tmpptr->extended = -1;
		tmpptr->last = tmpptr;
		tmpptr->korh = 'h';
	}

}

void FreeLocalHomTable_part( LocalHom **localhomtable, int n, int m ) 
{
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<m; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}
#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

void FreeLocalHomTable_two( LocalHom **localhomtable, int n, int m ) 
{
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<m; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}

	for( i=n; i<n+m; i++ ) 
	{
		for( j=0; j<n; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}
#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

void FreeLocalHomTable_one( LocalHom **localhomtable, int n, int m ) 
{
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<m; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}

#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

void FreeLocalHomTable_half( LocalHom **localhomtable, int n ) 
{
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<n-i; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}
#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}
void FreeLocalHomTable( LocalHom **localhomtable, int n ) 
{
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<n; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}
#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

char *progName( char *str )
{
    char *value; 
    if( ( value = strrchr( str, '/' ) ) != NULL )
        return( value+1 );
    else    
        return( str );
}

static void tabtospace( char *str )
{
	char *p;
//	fprintf( stderr, "before = %s\n", str );
	while( NULL != ( p = strchr( str , '\t' ) ) )
	{
		*p = ' ';
	}
//	fprintf( stderr, "after = %s\n", str );
}

static char *extractfirstword( char *str )
{
	char *val = str;

	tabtospace( str );
	while( *str )
	{
		if( val == str && *str == ' ' )
		{
			val++; str++;
		}
		else if( *str != ' ' )
		{
			str++;
		}
		else if( *str == ' ' )
		{
			*str = 0;
		}
	}
	return( val );
}

void phylipout_pointer( FILE *fp, int nseq, int maxlen, char **seq, char **name, int *order, int namelen )
{
	int pos, pos2, j;
	if( namelen == -1 ) namelen = 10;
	pos = 0;

	fprintf( fp, " %d %d\n", nseq, maxlen );
	
	while( pos < maxlen )
	{
		for( j=0; j<nseq; j++ )
		{
			if( pos == 0 )
				fprintf( fp, "%-*.*s", namelen, namelen, extractfirstword( name[order[j]]+1 ) );
			else
				fprintf( fp, "%-*.*s", namelen, namelen, "" );

			pos2 = pos;
			while( pos2 < pos+41 && pos2 < maxlen )
			{
				fprintf( fp, " %.10s", seq[order[j]]+pos2 );
				pos2 += 10;
			}
			fprintf( fp, "\n" );
		}
		fprintf( fp, "\n" );
		pos += 50;
	}
}

void clustalout_pointer( FILE *fp, int nseq, int maxlen, char **seq, char **name, char *mark, char *comment, int *order, int namelen )
{
	int pos, j;
	if( namelen == -1 ) namelen = 15;
	pos = 0;
	if( comment == NULL )
		fprintf( fp, "CLUSTAL format alignment by MAFFT (v%s)\n\n", VERSION );
	else
		fprintf( fp, "CLUSTAL format alignment by MAFFT %s (v%s)\n\n", comment, VERSION );
	
	while( pos < maxlen )
	{
		fprintf( fp, "\n" );
		for( j=0; j<nseq; j++ )
		{
			fprintf( fp, "%-*.*s ", namelen, namelen, extractfirstword( name[order[j]]+1 ) );
			fprintf( fp, "%.60s\n", seq[order[j]]+pos ); // 長さが違うとだめ。
		}
		if( mark )
		{
			fprintf( fp, "%-*.*s ", namelen, namelen, "" );
			fprintf( fp, "%.60s\n", mark + pos ); // 長さが違うとだめ。
		}
		pos += 60;
	}
}


void writeData_reorder_pointer( FILE *fp, int locnjob, char **name, int *nlen, char **aseq, int *order )
{
	int i, j, k;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		k = order[i];
#if DEBUG
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		fprintf( fp, ">%s\n", name[k]+1 );
#if 1
		if( LineLengthInFASTA < 0 )
			fprintf( fp, "%s\n", aseq[k] );
		else
		{
			nalen = strlen( aseq[k] );
			for( j=0; j<nalen; j=j+LineLengthInFASTA ) fprintf( fp, "%.*s\n", LineLengthInFASTA, aseq[k]+j );
		}
#else
		fprintf( fp, "%s\n", aseq[k] );
#endif
	}
}
void writeData_reorder( FILE *fp, int locnjob, char name[][B], int nlen[], char **aseq, int *order )
{
	int i, j, k;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		k = order[i];
#if DEBUG
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		nalen = strlen( aseq[k] );
		fprintf( fp, ">%s\n", name[k]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[k]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[k]+j );
#endif
		}
	}
}
static void showaamtxexample()
{
	fprintf( stderr, "Format error in aa matrix\n" );
	fprintf( stderr, "# Example:\n" );
	fprintf( stderr, "# comment\n" );
	fprintf( stderr, "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n" );
	fprintf( stderr, "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0\n" );
	fprintf( stderr, "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3\n" );
	fprintf( stderr, "...\n" );
	fprintf( stderr, "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4\n" );
	fprintf( stderr, "frequency 0.07 0.05 0.04 0.05 0.02 .. \n" );
	fprintf( stderr, "# Example end\n" );
	fprintf( stderr, "Only the lower half is loaded\n" );
	fprintf( stderr, "The last line (frequency) is optional.\n" );
	exit( 1 );
}

double *loadaamtx( int *rescalept )
{
	int i, j, k, ii, jj;
	double *val;
	double **raw;
	int *map;
	char *aaorder = "ARNDCQEGHILKMFPSTWYV";
	char *inorder;
	char *line;
	char *ptr1;
	char *ptr2;
	char *mtxfname = "_aamtx";
	FILE *mf;
	char key[1000];


	raw = AllocateDoubleMtx( 21, 20 );
	val = AllocateDoubleVec( 420 );
	map = AllocateIntVec( 20 );

	if( dorp != 'p' )
	{
		fprintf( stderr, "User-defined matrix is not supported for DNA\n" );
		exit( 1 );
	}

	mf = fopen( mtxfname, "r" );
	if( mf == NULL )
	{
		fprintf( stderr, "Cannot open the _aamtx file\n" );
		exit( 1 );
	}

	inorder = calloc( 1000, sizeof( char ) );
	line = calloc( 1000, sizeof( char ) );
	

	while( !feof( mf ) )
	{
		fgets( inorder, 999, mf );
		if( inorder[0] != '#' ) break;
	}
	ptr1 = ptr2 = inorder;
	while( *ptr2 )
	{
		if( isalpha( *ptr2 ) )
		{
			*ptr1 = toupper( *ptr2 );
			ptr1++;
		}
		ptr2++;
	}
	inorder[20] = 0;

	for( i=0; i<20; i++ )
	{
		ptr2 = strchr( inorder, aaorder[i] );
		if( ptr2 == NULL )
		{
			fprintf( stderr, "%c: not found in the first 20 letters.\n", aaorder[i] );
			showaamtxexample();
		}
		else
		{
			map[i] = ptr2 - inorder;
		}
	}

	i = 0;
	while( !feof( mf ) )
	{
		fgets( line, 999, mf );
//		fprintf( stderr, "line = %s\n", line );
		if( line[0] == '#' ) continue;
		ptr1 = line;
//		fprintf( stderr, "line = %s\n", line );
		for( j=0; j<=i; j++ )
		{
			while( !isdigit( *ptr1 ) && *ptr1 != '-' && *ptr1 != '.' )
				ptr1++;

			raw[i][j] = atof( ptr1 );
//			fprintf( stderr, "raw[][]=%f, %c-%c %d-%d\n", raw[i][j], inorder[i], inorder[j], i, j );
			ptr1 = strchr( ptr1, ' ' );
			if( ptr1 == NULL && j<i) showaamtxexample();
		}
		i++;
		if( i > 19 ) break;
	}

	*rescalept = 1;
	for( i=0; i<20; i++ ) raw[20][i] = -1.0;
	while( !feof( mf ) )
	{
		fgets( line, 999, mf );

		sscanf( line, "%s", key );
		
		if( !strcmp( key, "norescale" ) )
		{
			reporterr( "no rescale\n" );
			*rescalept = 0;
			break;
		}
//		else if( line[0] == 'f' )
		else if( !strcmp( key, "frequency" ) )
		{
//			fprintf( stderr, "found! line = %s\n", line );
			ptr1 = line;
			for( j=0; j<20; j++ )
			{
				while( !isdigit( *ptr1 ) && *ptr1 != '-' && *ptr1 != '.' )
					ptr1++;
	
				raw[20][j] = atof( ptr1 );
//				fprintf( stderr, "raw[20][]=%f, %c %d\n", raw[20][j], inorder[i], j );
				ptr1 = strchr( ptr1, ' ' );
				if( ptr1 == NULL && j<19) showaamtxexample();
			}
			break;
		}
	}

	k = 0;
	for( i=0; i<20; i++ )
	{
		for( j=0; j<=i; j++ )
		{
			if( i != j )
			{
				ii = MAX( map[i], map[j] );
				jj = MIN( map[i], map[j] );
			}
			else ii = jj = map[i];
			val[k++] = raw[ii][jj];
//			fprintf( stderr, "%c-%c, %f\n", aaorder[i], aaorder[j], val[k-1] );
		}
	}
	for( i=0; i<20; i++ ) val[400+i] = raw[20][map[i]];

	fprintf( stderr, "inorder = %s\n", inorder );
	fclose( mf );
	free( inorder );
	free( line );
	FreeDoubleMtx( raw );
	free( map );
	return( val );
}

static void tab2space( char *s ) // nen no tame
{
	while( *s )
	{
		if( *s == '\t' ) *s = ' ';
		s++;
	}
}

static int readasubalignment( char *s, int *t, int *preservegaps )
{
	int v = 0;
	char status = 's';
	char *pt = s;
	*preservegaps = 0;
	tab2space( s );
	while( *pt )
	{
		if( *pt == ' ' )
		{
			status = 's';
		}
		else
		{
			if( status == 's' ) 
			{
				if( *pt == '\n' || *pt == '#' ) break;
				status = 'n';
				t[v] = atoi( pt );
				if( t[v] == 0 )
				{
					fprintf( stderr, "Format error? Sequences must be specified as 1, 2, 3...\n" );
					exit( 1 );
				}
				if( t[v] < 0 ) *preservegaps = 1;
				t[v] = abs( t[v] );
				t[v] -= 1;
				v++;
			}
		}
		pt++;
	}
	t[v] = -1;
	return( v );
}

static int countspace( char *s )
{
	int v = 0;
	char status = 's';
	char *pt = s;
	tab2space( s );
	while( *pt )
	{
		if( *pt == ' ' )
		{
			status = 's';
		}
		else
		{
			if( status == 's' ) 
			{
				if( *pt == '\n' || *pt == '#' ) break;
				v++;
				status = 'n';
				if( atoi( pt ) == 0 )
				{
					fprintf( stderr, "Format error? Sequences should be specified as 1, 2, 3...\n" );
					exit( 1 );
				}
			}
		}
		pt++;
	}
	return( v );
}


void readsubalignmentstable( int nseq, int **table, int *preservegaps, int *nsubpt, int *maxmempt )
{
	FILE *fp;
	char *line;
	int linelen = 1000000;
	int nmem;
	int lpos;
	int i, p;
	int *tab01;

	line = calloc( linelen, sizeof( char ) );
	fp = fopen( "_subalignmentstable", "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _subalignmentstable\n" );
		exit( 1 );
	}
	if( table == NULL )
	{
		*nsubpt = 0;
		*maxmempt = 0;
		while( 1 )
		{
			fgets( line, linelen-1, fp );
			if( feof( fp ) ) break;
			if( line[strlen(line)-1] != '\n' )
			{
				fprintf( stderr, "too long line? \n" );
				exit( 1 );
			}
			if( line[0] == '#' ) continue;
			if( atoi( line ) == 0 ) continue;
			nmem = countspace( line );
			if( nmem > *maxmempt ) *maxmempt = nmem;
			(*nsubpt)++;
		}
	}
	else
	{
		tab01 = calloc( nseq, sizeof( int ) );
		for( i=0; i<nseq; i++ ) tab01[i] = 0;
		lpos = 0;
		while( 1 )
		{
			fgets( line, linelen-1, fp );
			if( feof( fp ) ) break;
			if( line[strlen(line)-1] != '\n' )
			{
				fprintf( stderr, "too long line? \n" );
				exit( 1 );
			}
			if( line[0] == '#' ) continue;
			if( atoi( line ) == 0 ) continue;
			readasubalignment( line, table[lpos], preservegaps+lpos );
			for( i=0; (p=table[lpos][i])!=-1; i++ )
			{
				if( tab01[p] )
				{
					fprintf( stderr, "\nSequence %d appears in different groups.\n", p+1 );
					fprintf( stderr, "Hierarchical grouping is not supported.\n\n" );
					exit( 1 );
				}
				tab01[p] = 1;
				if( p > nseq-1 )
				{
					fprintf( stderr, "Sequence %d does not exist in the input sequence file.\n", p+1 );
					exit( 1 );
				}
			}
			lpos++;
		}
		free( tab01 );
	}
	fclose( fp );
	free( line );
}


void readmccaskill( FILE *fp, RNApair **pairprob, int length )
{
	char gett[1000];
	int *pairnum;
	int i;
	int left, right;
	double prob;
	int c;

	pairnum = (int *)calloc( length, sizeof( int ) );
	for( i=0; i<length; i++ ) pairnum[i] = 0;

	c = getc( fp );
	{
		if( c != '>' )
		{
			fprintf( stderr, "format error in hat4 - 1\n" );
			exit( 1 );
		}
	}
	fgets( gett, 999, fp );
	while( 1 )
	{
		if( feof( fp ) ) break;
		c = getc( fp );
		ungetc( c, fp );
		if( c == '>' || c == EOF )
		{
			break;
		}
		fgets( gett, 999, fp );
//		fprintf( stderr, "gett = %s\n", gett );
		sscanf( gett, "%d %d %lf", &left, &right, &prob );

		if( left >= length || right >= length )
		{
			fprintf( stderr, "format error in hat4 - 2\n" );
			exit( 1 );
		}

		if( prob < 0.01 ) continue; // 080607, mafft ni dake eikyou

		if( left != right && prob > 0.0 )
		{
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

void readpairfoldalign( FILE *fp, char *s1, char *s2, char *aln1, char *aln2, int q1, int q2, int *of1, int *of2, int sumlen )
{
	char gett[1000];
	int *maptoseq1;
	int *maptoseq2;
	char dumc;
	int dumi;
	char sinseq[100], sinaln[100];
	int posinseq, posinaln;
	int alnlen;
	int i;
	int pos1, pos2;
	char *pa1, *pa2;
	char qstr[1000];

	*of1 = -1;
	*of2 = -1;

	maptoseq1 = AllocateIntVec( sumlen+1 );
	maptoseq2 = AllocateIntVec( sumlen+1 );

	posinaln = 0; // foldalign ga alingment wo kaesanaitok no tame.

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; ALIGNING", 10 ) ) break;
	}
	sprintf( qstr, "; ALIGNING            %d against %d\n", q1+1, q2+1 );
	if( strcmp( gett, qstr ) )
	{
		fprintf( stderr, "Error in FOLDALIGN\n" );
		fprintf( stderr, "qstr = %s, but gett = %s\n", qstr, gett );
		exit( 1 );
	}

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; --------", 10 ) ) break;
	}


	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; ********", 10 ) ) break;
//		fprintf( stderr, "gett = %s\n", gett );
		sscanf( gett, "%c %c %s %s %d %d", &dumc, &dumc, sinseq, sinaln, &dumi, &dumi );
		posinaln = atoi( sinaln );
		posinseq = atoi( sinseq );
//		fprintf( stderr, "posinseq = %d\n", posinseq );
//		fprintf( stderr, "posinaln = %d\n", posinaln );
		maptoseq1[posinaln-1] = posinseq-1;
	}
	alnlen = posinaln;

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; --------", 10 ) ) break;
	}

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; ********", 10 ) ) break;
//		fprintf( stderr, "gett = %s\n", gett );
		sscanf( gett, "%c %c %s %s %d %d", &dumc, &dumc, sinseq, sinaln, &dumi, &dumi );
		posinaln = atof( sinaln );
		posinseq = atof( sinseq );
//		fprintf( stderr, "posinseq = %d\n", posinseq );
//		fprintf( stderr, "posinaln = %d\n", posinaln );
		maptoseq2[posinaln-1] = posinseq-1;
	}
	if( alnlen != posinaln )
	{
		fprintf( stderr, "Error in foldalign?\n" );
		exit( 1 );
	}

	pa1 = aln1;
	pa2 = aln2;
	for( i=0; i<alnlen; i++ )
	{
		pos1 = maptoseq1[i];
		pos2 = maptoseq2[i];

		if( pos1 > -1 )
			*pa1++ = s1[pos1];
		else
			*pa1++ = '-';

		if( pos2 > -1 )
			*pa2++ = s2[pos2];
		else
			*pa2++ = '-';
	}
	*pa1 = 0;
	*pa2 = 0;

	*of1 = 0;
	for( i=0; i<alnlen; i++ )
	{
		*of1 = maptoseq1[i];
		if( *of1 > -1 ) break;
	}
	*of2 = 0;
	for( i=0; i<alnlen; i++ )
	{
		*of2 = maptoseq2[i];
		if( *of2 > -1 ) break;
	}

//	fprintf( stderr, "*of1=%d, aln1 = :%s:\n", *of1, aln1 );
//	fprintf( stderr, "*of2=%d, aln2 = :%s:\n", *of2, aln2 );

	free( maptoseq1 );
	free( maptoseq2 );
}

int myatoi( char *in )
{
	if( in == NULL )
	{
		fprintf( stderr, "Error in myatoi()\n" );
		exit( 1 );
	}
	return( atoi( in ) );
}

unsigned long long myatoll( char *in )
{
	if( in == NULL )
	{
		fprintf( stderr, "Error in myatoi()\n" );
		exit( 1 );
	}

	unsigned long long tanni;

	if( strchr( in, 'G' ) ) tanni = 1000 * 1000 * 1000;
	else if( strchr( in, 'M' ) ) tanni = 1000 * 1000;
	else if( strchr( in, 'k' ) ) tanni = 1000;
	else tanni = 1;

	return( tanni * atoi( in ) );
}


double myatof( char *in )
{
	if( in == NULL )
	{
		fprintf( stderr, "Error in myatof()\n" );
		exit( 1 );
	}
	return( atof( in ) );
}

void reporterr( const char *str, ... )
{
//	static int loglen = 0;
	va_list args;

	if( gmsg )
	{
# if 1  // ato de sakujo
		static FILE *errtmpfp = NULL;
		if( errtmpfp == NULL )
			errtmpfp = fopen( "maffterr", "w" );
		else
			errtmpfp = fopen( "maffterr", "a" );
		va_start( args, str );
		vfprintf( errtmpfp, str, args );
		va_end( args );
		fclose( errtmpfp );
#endif

#if 0
		char *tmpptr;
		tmpptr = (char *)realloc( *gmsg, (loglen+10000) * sizeof( char ) );
		if( tmpptr == NULL )
		{
			fprintf( stderr, "Cannot relloc *gmsg\n" );
			exit( 1 );
		}
		*gmsg = tmpptr;
		va_start( args, str );
		loglen += vsprintf( *gmsg + loglen, str, args );
		va_end( args );


		va_start( args, str );
		loglen += vsprintf( *gmsg + loglen, str, args );
		va_end( args );
		*(*gmsg + loglen) = 0;
		if( loglen > gmsglen - 100 ) loglen = 0; // tekitou
#endif

	}
	else
	{
		va_start( args, str );
		vfprintf( stderr, str, args );
		va_end( args );
//		fflush( stderr ); // iru?
	}
	return;
}


#if !defined(mingw) && !defined(_MSC_VER)
void setstacksize(rlim_t kStackSize )
{
//	const rlim_t kStackSize = 100 * 1024 * 1024;   // min stack size = 10MB
	struct rlimit rl;
	int result;
	rlim_t originalsize;

	result = getrlimit(RLIMIT_STACK, &rl);
	if (result == 0)
	{
		originalsize = rl.rlim_cur;
		if (rl.rlim_cur < kStackSize)
		{
			rl.rlim_cur = kStackSize;
			reporterr( "stacksize: %d kb->%d kb\n", originalsize/1024, rl.rlim_cur/1024 );
			result = setrlimit(RLIMIT_STACK, &rl);
			if (result != 0)
			{
				reporterr( "Warning: Failed to extend stack size. It's ok in most cases but there may be problems in --pileup and --chainedtree.\n" );
			}
		}
		else
			reporterr( "stacksize: %d kb\n", rl.rlim_cur / 1024 );
	}
	else
		reporterr( "Warning: Cannot check stack size.\n" );
}
#endif

void treeout_bin( FILE *fp, int n, int ***topol, double **len, Treedep *dep, int *nfilesfornode )
// dep ha nakutemo topol kara saigen dekiru.
{

	int i;
	char c = '\n';
	for( i=0; i<n-1; i++ )
	{
		if( fwrite( topol[i][0], sizeof( int ), 1, fp ) != 1 || 
		    fwrite( topol[i][1], sizeof( int ), 1, fp ) != 1 )
		{
			reporterr( "write error in treeout_bin(), topol, i=%d\n", i );
			exit( 1 );
		}
		if( fwrite( len[i]+0, sizeof( double ), 1, fp ) != 1 || 
		    fwrite( len[i]+1, sizeof( double ), 1, fp ) != 1 )
		{
			reporterr( "write error in treeout_bin(), len, i=%d\n", i );
			exit( 1 );
		}
		if( fwrite( &(dep[i].child0), sizeof( int ), 1, fp ) != 1 || 
		    fwrite( &(dep[i].child1), sizeof( int ), 1, fp ) != 1 ||
		    fwrite( &(nfilesfornode[i]), sizeof( int ), 1, fp ) != 1 ||
		    fwrite( &(dep[i].distfromtip), sizeof( double ), 1, fp ) != 1 )
		{
			reporterr( "write error in treeout_bin(), dep, i=%d\n", i );
			exit( 1 );
		}
		if( fwrite( &c, sizeof( char ), 1, fp ) != 1 )
		{
			reporterr( "write error in treeout_bin(), c, i=%d\n", i );
			exit( 1 );
		}
	}
}

void treein_bin( FILE *fp, int n, int ***topol, double **len, Treedep *dep, int *nfilesfornode )
// dep ha nakutemo topol kara saigen dekiru.
{
	int i;
	char c;
	for( i=0; i<n-1; i++ )
	{

		topol[i][0] = calloc( sizeof( int ), 2 );
		topol[i][1] = calloc( sizeof( int ), 2 );

		topol[i][0][1] = -1;
		topol[i][1][1] = -1;

		if( fread( topol[i][0], sizeof( int ), 1, fp ) != 1 || 
		    fread( topol[i][1], sizeof( int ), 1, fp ) != 1 )
		{
			reporterr( "read error in treein_bin(), topol, i=%d\n", i );
			exit( 1 );
		}
		if( fread( len[i]+0, sizeof( double ), 1, fp ) != 1 || 
		    fread( len[i]+1, sizeof( double ), 1, fp ) != 1 )
		{
			reporterr( "read error in treein_bin(), len, i=%d\n", i );
			exit( 1 );
		}
		if( fread( &(dep[i].child0), sizeof( int ), 1, fp ) != 1 || 
		    fread( &(dep[i].child1), sizeof( int ), 1, fp ) != 1 ||
		    fread( &(nfilesfornode[i]), sizeof( int ), 1, fp ) != 1 ||
		    fread( &(dep[i].distfromtip), sizeof( double ), 1, fp ) != 1 )
		{
			reporterr( "read error in treein_bin(), dep, i=%d\n", i );
			exit( 1 );
		}
		if( fread( &c, sizeof( char ), 1, fp ) != 1 )
		{
			reporterr( "read error in treein_bin(), c, i=%d\n", i );
			exit( 1 );
		}

		if( c != '\n' )
		{
			reporterr( "Error in tree file\n" );
			exit( 1 );
		}
	}
}

void uselhout( FILE *fp, int n, int *uselh )
{

	if( fwrite( uselh, sizeof( int ), n, fp ) != n ) 
	{
		reporterr( "write error in uselhout()\n" );
		exit( 1 );
	}
}

int uselhin( FILE *fp, int n, int *uselh )
{
	if( fread( uselh, sizeof( int ), n, fp ) != n ) 
	{
		reporterr( "read error in uselhout()\n" );
		exit( 1 );
	}

	while( n-- ) if( *uselh++ == 0 ) return( 0 );
	return( 1 );
}

void getweightfromname( int n, double *w, char **name )
{
	int i;
	int res;
	char tmp[3][1000];
	for( i=0; i<n; i++ )
	{
		res = sscanf( name[i]+1, "n0=%s n1=%s n2=%s", tmp[0], tmp[1], tmp[2] );

		if( res == EOF )
		{
			reporterr( "error in reading \">%s\"\n", name[i]+1 );
			reporterr( "Format has to be \">n0=1000 n1=51 n2=2 sequencename\"\n" );
			exit( 1 );
		}

		//reporterr( "tmp[0] = %s\n", tmp[0] );
		//reporterr( "tmp[1] = %s\n", tmp[1] );
		//reporterr( "tmp[2] = %s\n", tmp[2] );

		w[i] = atof( tmp[0] );

		if( w[i] == 0.0 )
		{
			reporterr( "Error in reading \">%s\". n0=0?\n", name[i]+1 );
			reporterr( "Format has to be \">n0=1000 n1=51 n2=2 sequencename\"\n" );
			exit( 1 );
		}

		w[i] = 1.0 / atof( tmp[0] );

		if( w[i] == 0.0 || w[i] > 1.0 )
		{
			reporterr( "Warning: weight for \">%s\" is %f\n", name[i]+1, w[i] );
		}
		//reporterr( "w[%d] = %f\n", i, w[i] );
	}
}

#if 0
#include <sys/time.h>
#include <sys/resource.h>

void use_getrusage(void)
{
	struct rusage r;
	if (getrusage(RUSAGE_SELF, &r) != 0) {
		/*Failure*/
	}
	fprintf(stderr, "\nmaxrss = %ld MB\n", r.ru_maxrss/1000);
}

#endif

void commongappick( int nseq, char **seq )
{
	int i, j, count;
	int len = strlen( seq[0] );
#if 1

	int *mapfromnewtoold;


	mapfromnewtoold = calloc( len+1, sizeof( int ) );

	for( i=0, count=0; i<=len; i++ ) 
	{
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			mapfromnewtoold[count++] = i;
	 	}
	}
//	mapfromnewtoold[count] = -1; // iranai
	for( j=0; j<nseq; j++ )
	{
		for( i=0; i<count; i++ )
		{
			seq[j][i] = seq[j][mapfromnewtoold[i]];
		}
	}
	free( mapfromnewtoold );
#else
--------------------------

	int *mapfromoldtonew;
	int pos;

	mapfromoldtonew = calloc( len+1, sizeof( int ) );
	for( i=0; i<=len; i++ ) mapfromoldtonew[i] = -1;

	for( i=0, count=0; i<=len; i++ ) 
	{
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			mapfromoldtonew[i] = count;
			count++;
	 	}
	}
	for( j=0; j<nseq; j++ )
	{
		for( i=0; i<=len; i++ ) 
		{
			if( (pos=mapfromoldtonew[i]) != -1 )
				seq[j][pos] = seq[j][i];
		}
	}
	free( mapfromoldtonew );
--------------------------

	for( i=0, count=0; i<=len; i++ ) 
	{
	/*
		allgap = 1;
		for( j=0; j<nseq; j++ ) 
			allgap *= ( seq[j][i] == '-' );
		if( !allgap )
	*/
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			for( j=0; j<nseq; j++ )
			{
				seq[j][count] = seq[j][i];
			}
			count++;
	 	}
	}

#endif
}

void readexternalanchors( ExtAnch **extanch, int nseq, int *nogaplen )
{
	FILE *fp;
	int size, lineno;
	char buf[10000];
	fp = fopen( "_externalanchors", "r" );

	if( fp == NULL )
	{
		reporterr( "Cannot open _externalanchors\n" );
		exit( 1 );
	}

	size = 0;
	lineno = 0;
	while( 1 )
	{
		lineno++;
//		reporterr( "size = %d\n", size );
		fgets( buf, 9999, fp );
		if( feof( fp ) ) break;

		if( buf[0] == '#' ) continue;

//		reporterr( "buf=%s\n", buf );
		*extanch = realloc( *extanch, sizeof( ExtAnch ) * (size+2) );
		if( *extanch == NULL )
		{
			reporterr( "Cannot realloc *extanch\n" );
			exit( 1 );
		}

		sscanf( buf, "%d %d %d %d %d %d %d", &(((*extanch)+size)->i), &(((*extanch)+size)->j), &(((*extanch)+size)->starti), &(((*extanch)+size)->endi), &(((*extanch)+size)->startj), &(((*extanch)+size)->endj), &(((*extanch)+size)->score) );
//		reporterr( "i=%d, j=%d, %d-%d, %d-%d, score=%d\n", (*extanch)[size].i, (*extanch)[size].j, (*extanch)[size].starti, (*extanch)[size].endi, (*extanch)[size].startj, (*extanch)[size].endj, (*extanch)[size].score );

		((*extanch)+size)->i -= 1; // 1-origin -> 0-origin
		((*extanch)+size)->j -= 1; // 1-origin -> 0-origin
		((*extanch)+size)->starti -= 1; 
		((*extanch)+size)->startj -= 1;
		((*extanch)+size)->endi -= 1;
		((*extanch)+size)->endj -= 1;

		if( (*extanch)[size].i >= nseq || (*extanch)[size].j >= nseq )
		{
			reporterr( "\nOut of range?  The input file has %d sequences but pair %d-%d was specified in line %d.\nNote that sequence IDs are counted from 1.\n", nseq, (*extanch)[size].i+1, (*extanch)[size].j+1, lineno );
			exit( 1 );
		}
		if( (*extanch)[size].i >= (*extanch)[size].j )
		{
			reporterr( "\nFormat problem?  \"%d %d\" in line %d.\nThe sequence id of the first column must be less than the second.\n", (*extanch)[size].i+1, (*extanch)[size].j+1, lineno );
			exit( 1 );
		}
		if( (*extanch)[size].starti > nogaplen[(*extanch)[size].i] )
		{
			reporterr( "\nOut of range?  len(seq%d)=%d, but anchor=%d in line %d.\nNote that position is counted from 1.\n", (*extanch)[size].i+1, nogaplen[(*extanch)[size].i], (*extanch)[size].starti+1, lineno );
			exit( 1 );
		}
		if( (*extanch)[size].startj > nogaplen[(*extanch)[size].j] )
		{
			reporterr( "\nOut of range?  len(seq%d)=%d, but anchor=%d in line %d.\nNote that position is counted from 1.\n", (*extanch)[size].j, nogaplen[(*extanch)[size].j]+1, (*extanch)[size].startj+1, lineno );
			exit( 1 );
		}

		size++;
		(*extanch)[size].i = (*extanch)[size].j = -1;

	}
	fclose( fp );
}

static char id2nuc( int id ) // just for error message
{
	if( id == 0 ) return 't';
	if( id == 1 ) return 'c';
	if( id == 2 ) return 'a';
	if( id == 3 ) return 'g';
	return -1;
}
static void id2codon( int id, char *codon ) // just for error message
{
	codon[0] = id2nuc( id/16 );
	id = id - 16 * id/16;
	codon[1] = id2nuc( id/4 );
	id = id - 4 * id/4;
	codon[2] = id2nuc( id );
	return;
}

static int nuc2id( char nuc )
{
	if( nuc == 't' ) return 0;
	if( nuc == 'c' ) return 1;
	if( nuc == 'a' ) return 2;
	if( nuc == 'g' ) return 3;
	return -1;
}

int codon2id( char *codon )
{
	int id0, id1, id2, id;
	id0 = nuc2id( codon[0] );
	id1 = nuc2id( codon[1] );
	id2 = nuc2id( codon[2] );

	if( id0 < 0 || id1 < 0 | id2 < 0 ) return( -1 );

	return id0 * 16 + id1 * 4 + id2;
}

double **loadcodonscore( FILE *fp, double **scoremtx )
{
	int i, j;
	char *buf = calloc( sizeof( char ), 1000 );
	char codon1[3], codon2[3];
	char aa1[1000], aa2[1000];
	int id1, id2;
	double score;
	for( i=0; i<64; i++ ) for( j=0; j<64; j++ ) scoremtx[i][j] = -99999;
	i = 0;
	while( fgets( buf, 1000, fp ) )
	{
//		reporterr( "%s", buf );
		if( buf[0] == '#' ) continue;
		if( buf[strlen(buf)-1] != '\n' )
		{
			reporterr( "%s: too long in codonscore file.\n", buf );
			exit( 1 );
		}
		sscanf( buf, "%3s %s %3s %s %lf", codon1, aa1, codon2, aa2, &score );
//		reporterr( "codon1=%s, id=%d\n", codon1, codon2id( codon1 ) );
//		reporterr( "codon2=%s, id=%d\n", codon2, codon2id( codon2 ) );
//		reporterr( "score=%f\n", score );
		id1 = codon2id( codon1 );
		id2 = codon2id( codon2 );
		if( id1 < 0 || id2 < 0 )
		{
			reporterr( "Cannot use codon pair %s - %s: Use small letter, a, c, g, t (instead of u)\n", codon1, codon2 );
			exit( 1 );
		}
		scoremtx[id1][id2] = scoremtx[id2][id1] = score;
		i++;
	}
	free( buf );

	for( i=0; i<64; i++ ) for( j=0; j<64; j++ )
	{
		char codon[3];
		if( scoremtx[i][j] == -99999 )
		{
			id2codon( i, codon );
			reporterr( "\nCodon score for %s", codon );
			id2codon( j, codon );
			reporterr( "-%s (id%d-id%d) is not given.\n", codon, i, j );
			exit( 1 );
		}
//		reporterr( "id%d-id%d = %f\n", i, j, scoremtx[i][j] );
	}
	return( scoremtx );
}
