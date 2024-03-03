#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

int main( int ac, char **av )
{
	unsigned int c;
	FILE *fp;
	unsigned char buf[10000];
	unsigned char *bpt;
	char format;

	if( ac == 1 || ( ac == 2 && av[1][0] == '-' ) )
	{
		fp = stdin;
	}
	else if( ac == 2 )
	{
		fp = fopen( av[1], "rb" );
		if( fp == NULL )
		{
			fprintf( stderr, "%s: Cannot open %s.\n", av[0], av[1] );
			exit( 1 );
		}
	}
	else
	{
		fprintf( stderr, "Usage %s input > output\n", av[0] );
		exit( 1 );
	}

#if mingw
	setmode( fileno( fp ), O_BINARY );
#endif

#if 1
	format = 'f';
	c = fgetc( fp );
	if( c == 'C' ) format = 'c';
	ungetc( c, fp );
#else
	fgets( (char *)buf, 999, fp );
	rewind( fp );
	if( !strncmp( (char *)buf, "CLUSTAL", 7 ) ) format = 'c';
	else format = 'f';
#endif

	if( format == 'c' ) // clustal
	{
		int ln = 0;
		int titlelen = -1;
		while( 1 )
		{
			fgets( (char *)buf, 999, fp );
			if( feof( fp ) ) break;
			if( ln == 0 ) 
			{
				ln = 1;
				printf( "%s", buf );
				continue;
			}
	
			bpt = (unsigned char *)strchr( (char *)buf, ' ' );
			if( bpt == NULL ) 
			{
				printf( "\n" );
				continue;
			}
			if( titlelen == -1 ) 
			{
				while( *++bpt == ' ' ) 
					;
				titlelen = bpt - buf;
			}
			else
			{
				bpt = buf + titlelen;
			}
			*(bpt-1) = 0;
			printf( "%s ", buf );
	
			while( (c=(unsigned int)*bpt++)!='\n' ) 
			{
				if( c == '-' ) printf( "-- " );
				else if( c == '=' ) printf( "== " );
				else if( c == '*' ) printf( "** " );
				else if( c == ' ' ) printf( "   " );
				else printf( "%02x ", c );
			}
			printf( "\n" );
		}
	}
	else // fasta
	{
		while( 1 )
		{
			c = fgetc( fp );
			if( c == EOF ) break;
			else if( c == '\n' ) printf( "\n" );
			else if( c == '\r' ) printf( "\r" ); // nai
			else if( c == '-' ) printf( "-- " );
			else if( c == '=' ) printf( "== " ); // nai
			else if( c == ' ' ) printf( "   " ); // nai
			else if( c == '>' || c == '<' )
			{
				printf( "%c", c );
				while( 1 )
				{
					c = fgetc( fp );
					printf( "%c", c );
					if( c == '\n' ) break;
				}
			}
			else printf( "%02x ", c );
		}
	}
	fclose( fp );
	return( 0 );
}
