#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void fgetstilspace( unsigned char *b, FILE *fp )
{
	unsigned int c;
	int alreadyread = 0;
	while( 1 )
	{
		c = fgetc( fp );
		if( c == ' ' && alreadyread == 0 ) continue;
		alreadyread = 1;
		if( c == ' ' || c == '\n' || c == '\r' || c == EOF ) 
		{
			ungetc( c, fp );
			break;
		}
		*b++ = (unsigned char)c;
	}
	*b = 0;
}

int main( int ac, char **av )
{
	unsigned int c;
	unsigned char buf[100];
	FILE *fp;
	int res;

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
	{
		while( 1 )
		{
			c = fgetc( fp );
			if( c == EOF ) break;
			else if( c == '\n' ) printf( "\n" );
			else if( c == '\r' ) printf( "\r" );
#if 0
			else if( c == '-' ) 
			{
				printf( "-" );
				fprintf( stderr, "Warning: '-' will be removed\n" );
			}
			else if( c == '=' )
			{
				printf( "=" );
				fprintf( stderr, "Warning: '=' will be removed\n" );
			}
//			else if( c == ' ' ) printf( " " ); // nai
#endif
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
			else
			{
				buf[0] = (unsigned char)c;
				fgetstilspace( buf+1, fp );
				//fprintf( stderr, "buf=%s\n", buf );
				if( strchr( (const char *)buf, '-' ) ) // added cast, 2019/Jan/25
				{
					printf( "-" );
					continue;
				}
				//res = sscanf( buf, " %x ", &c );
				res = sscanf( (const char *)buf, " %x ", &c ); // added cast, 2019/Jan/25
				if( res == EOF )
				{
					//fprintf( stderr, "%s was ignored.\n", buf );
				}
				else if( res != 1 )
				{
					fprintf( stderr, "Error in reading %s\n", buf );
					exit( 1 );
				}
				else if( c <= 0 || c > 0xff )
				{
					fprintf( stderr, "Out of range: 0x%x\n", c );
					//exit( 1 );
				}
				else if( c == 0x0d || c == 0x0a )
				{
					fprintf( stderr, "Warning: skipped 0x%x (CR or LF) that cannot be used in mafft --text.\n", c );
					//printf( "%c", c );
				}
				else if( c == 0x20 || c == 0x3E || c == 0x3C || c == 0x3D )
				{
					fprintf( stderr, "Warning: skipped 0x%x (%c) that cannot be used in mafft --text.\n", c, c );
					//printf( "%c", c );
				}
				else if( c == 0x2D )
				{
					fprintf( stderr, "Warning: put 0x%x (%c) that is interpreted as gap in mafft --text.\n", c, c );
					printf( "%c", c );
				}
				else
					printf( "%c", c );
			}
		}
	}
	fclose( fp );
	return( 0 );
}
