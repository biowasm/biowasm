 /* iteration  ( algorithm C ) */
#include "mltaln.h"

#define DEBUG 0

static void Writeoptions( FILE *fp )
{
    if( scoremtx == 1 )
        fprintf( fp, "Dayhoff( ... )\n" );
    else if( scoremtx == -1 )
        fprintf( fp, "DNA\n" );
    else if( scoremtx == 2 )
        fprintf( fp, "Miyata-Yasunaga\n" );
	else
		fprintf( fp, "JTT %dPAM\n", pamN );

	if( scoremtx == 0 )
	    fprintf( fp, "Gap Penalty = %+d, %+d\n", penalty, offset );
	else
	    fprintf( fp, "Gap Penalty = %+d\n", penalty );

    fprintf( fp, "marginal score to search : best - %f\n", cut );
    if( scmtd == 3 )
        fprintf( fp, "score of rnd or sco\n" );
    else if( scmtd == 4 )
        fprintf( fp, "score = sigma( score for a pair of homologous amino acids ) / ( number of amino acids pairs )\n" );
    else if( scmtd == 5 )
        fprintf( fp, "score : SP\n" );
    if( mix )
        fprintf( fp, "?\n" );
    else
    { 
        if( weight == 2 )
            fprintf( fp, "weighted,  geta2 = %f\n", geta2 );
        else if( weight == 3 )
        {
            if( scmtd == 4 )
                fprintf( fp, "reversely weighted in function 'align', unweighted in function 'score_calc'\n" );
            else
                fprintf( fp, "weighted like ClustalW," );
        }
        else
            fprintf( fp, "unweighted\n" );
    }
    if( weight && utree )
    {
        fprintf( fp, "using tree defined by the file hat2 with simplified UPG method\n" );
    }
    if( weight && !utree )
        fprintf( fp, "using temporary tree by simplified UPG method\n" );
    fprintf( fp, "Algorithm %c\n", alg );
}




char **align0( double *wm, char **aseq, char *seq, double effarr[M], int icyc, int ex )
{
    char **result;

    if( alg == 'B' )
    {
		ErrorExit( "Sorry!" );
	/*
        if( outgap == 0 )
        {
            result = alignm1_o( wm, aseq, seq, scmx, effarr, icyc, ex );
        }
        if( outgap == 1 )
        {
            result = alignm1( wm, aseq, seq, scmx, effarr, icyc, ex );
        }
	*/
    }
    else if( alg == 'C' )
    {
        result = Calignm1( wm, aseq, seq, effarr, icyc, ex );
    }
    return( result );
}
    

double score_m_1_0( char **aseq, int locnjob, int s, double **eff, double effarr[M] )
{
    double x;

    if( alg == 'B' )
    {
		ErrorExit( "Sorry!" );
    }
    if( alg == 'C' )
    {
        x = Cscore_m_1( aseq, locnjob, s, eff );
    }
    fprintf( stderr, "in score_m_1_0 %f\n", x );
    return( x );
}

int iteration( int locnjob, char name[M][B], int nlen[M], char **aseq, char **bseq, int ***topol, double **len, double **eff ) 
{
    double tscore, mscore;
    int identity;
    static char *mseq1, **mseq2 = NULL;
	static char **result;
	int i, l;
	static double effarr[M];
	int s;
	int sss[2];
	char ou;
	int alloclen; 
	int resultlen;
	int nlenmax0 = nlenmax;
	FILE *prep;
	char sai[M];
	char sai1[M];
	char sai2[M];
#if 0
	double his[2][M][MAXITERATION/locnjob+1];
#else
	double ***his;
#endif
	int cyc[2];
	char shindou = 0;
	double wm;
	int returnvalue;

    for( i=0; i<locnjob; i++ ) 
    {
		sai[i] = 1;
        sai1[i] = 1;
        sai2[i] = 2;
    }
    sai[locnjob] = sai1[locnjob] = sai2[locnjob] = 0;


	Writeoptions( trap_g );

	his = AllocateDoubleCub( 2, M, MAXITERATION/locnjob+1 );

	if( mseq2 == NULL )
	{
    	alloclen = nlenmax * 2.0;
   		AllocateTmpSeqs( &mseq2, &mseq1, alloclen );
	}

	if( !tbitr && !tbweight )
	{
		writePre( locnjob, name, nlen, aseq, 0 );

#if 0
		prep = fopen( "best", "w" );
		Write( prep, locnjob, name, nlen, aseq );
		fclose( prep );
#endif
	}
	



	treeconstruction( aseq, locnjob, topol, len, eff );
	tscore = score_calc0( aseq, locnjob, eff, 0 );

#if DEBUG
    fprintf( stderr, "eff mtx in iteration\n" );
    for( i=0; i<locnjob; i++ )
    {
        for( j=0; j<locnjob; j++ ) 
        {
            fprintf( stderr, "%5.3f ", eff[i][j] );
        }
        fprintf( stderr, "\n" );
    }
#endif

    fprintf( stderr, "\n" );
	if( disp )
	{
    	fprintf( stderr, "aseq(initial)\n" );
		display( aseq, locnjob );
	}
	fprintf( stderr, "initial score = %f     \n", tscore );
	fprintf( stderr, "\n" );
	for( i=0; i<locnjob; i++ ) strcpy( bseq[i], aseq[i] );
	mscore = tscore;
    srand( time(NULL) );

	sss[1] = 0;
	sss[0] = locnjob-1;
/*
	sss[0] = (int)( (double)locnjob/2.0 );
*/
	ou = 1;
	cyc[0] = 0; cyc[1] = 0;

    for( s=-1, l=0; l<MAXITERATION; l++ )
    {
        int ss;
        double tmpscore, tmpscore1;

		if( strlen( aseq[0] ) > nlenmax )
			nlenmax = strlen( aseq[0] );

		/*
        s = ( int )( rnd() * locnjob );
		s++; 
		if( s == locnjob ) s = 0;
		ou = 0;
		*/
		if( ou == 0 )
		{
			ou = 1;
			s = sss[0];
			/*
			sss[0]++;
			if( sss[0] == locnjob )
			{
				sss[0] = 0;
				cyc[0]++;
			}
			*/
			sss[0]--;
			if( sss[0] == -1 )
			{
				sss[0] = locnjob-1;
				cyc[0]++;
			}
		}
		else
		{
			ou = 0;
			s = sss[1];
			sss[1]++;
			if( sss[1] == locnjob ) 
			{
				sss[1] = 0;
				cyc[1]++;
			}
		}
		fprintf( trap_g, "%d  ", weight );

/*
        for( i=0, count=0; i<strlen( aseq[s] ); i++ ) 
        {
            if( aseq[s][i] != '-' )
            {
                mseq1[count] = aseq[s][i];
                count++;
            }
        }
        mseq1[count] = 0;
*/
		gappick0( mseq1, aseq[s] );

		if( checkC )
			tmpscore = score_m_1_0( aseq, locnjob, s, eff, effarr );

		gappick( locnjob, s, aseq, mseq2, eff, effarr );

        result = align0( &wm, mseq2, mseq1, effarr, locnjob-2, s );
		resultlen = strlen( result[0] );
		if( resultlen > alloclen )
		{
			if( resultlen > nlenmax0*3 || resultlen > N )
			{
				fprintf(stderr, "Error in main1\n");
				exit( 1 );
			}
			FreeTmpSeqs( mseq2, mseq1 );
			alloclen = strlen( result[0] ) * 2.0;
			fprintf( stderr, "\n\ntrying to allocate TmpSeqs\n\n" );
			AllocateTmpSeqs( &mseq2, &mseq1, alloclen );
		}
		for( i=0; i<locnjob; i++ ) strcpy( mseq2[i], result[i] ); 

		if( checkC  )
			fprintf( stderr, "wm in iteration == %f\n", wm );

		strcpy( mseq1, mseq2[locnjob-1] );
/*
		Write( stdout, locnjob, name, nlen, mseq2 );
*/
        for( i=locnjob-2; i>=s; i-- ) strcpy( mseq2[i+1], mseq2[i] );
        strcpy( mseq2[s], mseq1 );
		if( checkC )
		{
			tmpscore1= score_m_1_0( mseq2, locnjob, s, eff, effarr );
			fprintf( stderr, "pick up %d, before ALIGNM1 score_m_1_0 = %f\n", s+1, tmpscore );
			fprintf( stderr, "pick up %d,  after ALIGNM1 score_m_1_0 = %f\n", s+1, tmpscore1 );
			if( tmpscore1 < tmpscore ) 
			{
				fprintf( stderr, "\7" );
				fprintf( trap_g, ">>>>>>>n\n" );
			}
			if( fabs( wm - tmpscore1 ) / wm  > 0.001 ) 
			{
				fprintf( stderr, "\7sorry\n" );
				exit( 1 );
			}
		}

        identity = !strcmp( mseq2[s], aseq[s] );
        if( s == locnjob - 1 ) ss = 0; else ss=s+1;

        identity *= !strcmp( mseq2[ss], aseq[ss] );

   	    if( !identity ) 
		{
			tmpscore = score_calc0( mseq2, locnjob, eff, s );
		}
		else tmpscore = tscore;

		if( disp )
		{
       		fprintf( stderr, "% 3d    % 3d / the rest   \n", l+1, s+1 );
       		display( mseq2, locnjob );
		}
       	fprintf( stderr, "% 3d    % 3d / the rest   \n", l+1, s+1 );
       	fprintf( stderr, "score = %f    mscore =  %f  ", tmpscore, mscore );

       	fprintf( trap_g, "%#4d    %#4d / the rest     ", l+1, s+1 );
       	fprintf( trap_g, "score = %f    mscore =  %f  ", tmpscore, mscore );

		if( identity ) 
		{
			fprintf( stderr, "( identical )\n" );
			fprintf( trap_g, "( identical )\n" );
			sai[s] = 2;
		}

        else if( tmpscore > mscore - cut )
        {
            fprintf( stderr, "accepted\n" );
            fprintf( trap_g, "accepted\n" );
            for( i=0; i<locnjob; i++ ) strcpy( aseq[i], mseq2[i] );
			strcpy( sai, sai1 );   /* kokoka ? */
			if( !tbitr && !tbweight )
			{
				writePre( locnjob, name, nlen, aseq, 0 );
			}
			strcpy( sai, sai1 );
			tscore = tmpscore;
			/*
			tscore = tmpscore = score_calc0( aseq, locnjob, eff, s );   * ? *
			*/
    		if( tmpscore > mscore ) 
			{
            	for( i=0; i<locnjob; i++ ) strcpy( bseq[i], mseq2[i] );
				treeconstruction( bseq, locnjob, topol, len, eff );
				tscore = mscore = score_calc0( bseq, locnjob, eff, s );
				fprintf( trap_g, "                                    -> %f\n", mscore );
				strcpy( sai, sai1 );   /* kokoka ? */
#if 0
				if( !tbitr && !tbweight )
				{	prep = fopen( "best", "w" );
					Write( prep, locnjob, name, nlen, bseq );
					fclose( prep );
				}
#endif
			}
        }

		else
		{
			if( tmpscore == tscore )
			{
				fprintf( stderr, "occational coincidence \n" );
				fprintf( trap_g, "occational coincidence\n" );
			}
			else
			{
				fprintf( stderr, "rejected\n" );
           	    fprintf( trap_g, "rejected\n" );
			}
			for( i=0; i<locnjob; i++ ) strcpy( aseq[i], bseq[i] );
			tscore = mscore;
			sai[s] = 2;
		}

/*
		prep = fopen( "cur", "w" );
		Write( prep, locnjob, name, nlen, mseq2 );
		fclose( prep );
*/

		his[ou][s][cyc[ou]] = tmpscore;
		if( !strcmp( sai, sai2 ) )
		{
			returnvalue = 0;
			fprintf( trap_g, "converged\n" );
			break;
		}
		for( i=cyc[ou]-1; i>0; i-- ) 
		{
			if( tmpscore == his[ou][s][i] ) 
			{
				shindou = 1;
				break;
			}
		}
		fprintf( stderr, "\n" );
		if( shindou == 1 )
		{
			returnvalue = -1;
			fprintf( trap_g, "oscillating\n" );
			break;
		}
	}
	if( l == MAXITERATION ) returnvalue = -2;
	FreeDoubleCub( his );
	return( returnvalue );
}

