#include "mltaln.h"
#include "dp.h"

TLS int commonAlloc1 = 0;
TLS int commonAlloc2 = 0;
TLS int **commonIP = NULL;
TLS int **commonJP = NULL;
int nthread = 1;
int nthreadpair = 1;
int randomseed = 0;
int parallelizationstrategy = BAATARI1;


char modelname[500];
int njob, nlenmax;
int amino_n[0x100];
char amino_grp[0x100];
//int amino_dis[0x100][0x100];
int **amino_dis = NULL;
double **n_disLN = NULL;
//double amino_dis_consweight_multi[0x100][0x100];
double **amino_dis_consweight_multi = NULL;
int **n_dis = NULL;
int **n_disFFT = NULL;
double **n_dis_consweight_multi = NULL;
unsigned char amino[0x100];
double polarity[0x100];
double volume[0x100];
int ribosumdis[37][37];

int ppid;
double thrinter;
double fastathreshold;
int pslocal, ppslocal;
int constraint;
int divpairscore;
int fmodel; // 1-> fmodel 0->default -1->raw
int nblosum; // 45, 50, 62, 80
int kobetsubunkatsu;
int bunkatsu;
int dorp = NOTSPECIFIED; // arguments de shitei suruto, tbfast -> pairlocalalign no yobidashi de futsugou
int niter;
int contin;
int calledByXced;
int devide;
int scmtd;
int weight;
int utree;
int tbutree;
int refine;
int check;
double cut;
int cooling;
int trywarp = 0;
int penalty, ppenalty, penaltyLN;
int penalty_dist, ppenalty_dist;
int RNApenalty, RNAppenalty;
int RNApenalty_ex, RNAppenalty_ex;
int penalty_ex, ppenalty_ex, penalty_exLN;
int penalty_EX, ppenalty_EX;
int penalty_OP, ppenalty_OP;
int penalty_shift, ppenalty_shift;
double penalty_shift_factor = 100.0;
int RNAthr, RNApthr;
int offset, poffset, offsetLN, offsetFFT;
int scoremtx;
int TMorJTT;
char use_fft;
char force_fft;
int nevermemsave;
int fftscore;
int fftWinSize;
int fftThreshold;
int fftRepeatStop;
int fftNoAnchStop;
int divWinSize;
int divThreshold;
int disp;
int outgap = 1;
char alg;
int cnst;
int mix;
int tbitr;
int tbweight;
int tbrweight;
int disopt;
int pamN;
int checkC;
double geta2;
int treemethod;
int kimuraR;
char *swopt;
int fftkeika;
int score_check;
int makedistmtx;
char *inputfile;
char *addfile;
int addprofile = 1;
int rnakozo;
char rnaprediction;
int scoreout = 0;
int spscoreout = 0;
int outnumber = 0;
int legacygapcost = 0;
double minimumweight = 0.0005;
int nwildcard = 0;

char *signalSM;
FILE *prep_g;
FILE *trap_g;
char **seq_g;
char **res_g;

double consweight_multi = 1.0;
double consweight_rna = 0.0;
char RNAscoremtx = 'n';

TLS char *newgapstr = "-";

int nalphabets = 26;
int nscoredalphabets = 20;

double specificityconsideration = 0.0;
int ndistclass = 10;
int maxdistclass = -1;

int gmsg = 0;

double sueff_global = SUEFF;

double lenfaca, lenfacb, lenfacc, lenfacd;
int maxl, tsize;

char codonpos = 0;
char codonscore = 0;


void initglobalvariables()
{
	commonAlloc1 = 0;
	commonAlloc2 = 0;
	commonIP = NULL;
	commonJP = NULL;
	nthread = 1;
	randomseed = 0;
	parallelizationstrategy = BAATARI1;

	trywarp = 0;
	penalty_shift_factor = 100.0;
	outgap = 1;
	addprofile = 1;
	scoreout = 0;
	outnumber = 0;
	legacygapcost = 0;
	consweight_multi = 1.0;
	consweight_rna = 0.0;
	RNAscoremtx = 'n';

	newgapstr = "-";

	nalphabets = 26;
	nscoredalphabets = 20;

	specificityconsideration = 0.0;
	ndistclass = 10;
	maxdistclass = -1;
	
	gmsg = 0;
}

// for usetmpfile
int compacttree = 0;
int lhlimit = INT_MAX;
int specifictarget = 0;
int nadd = 0; // <- static in tbfast.c, pairlocalalign.c
int usenaivescoreinsteadofalignmentscore = 0;
int nthreadreadlh = 1;
int LineLengthInFASTA = -1;
