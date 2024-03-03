void MtxuntDouble( double **, int );
void MtxmltDouble( double **, double **, int );

char *AllocateCharVec( int );
void FreeCharVec( char * );

char **AllocateCharMtx( int, int);
void ReallocateCharMtx( char **, int, int);
void FreeCharMtx( char ** );

double *AllocateFloatVec( int );
void FreeFloatVec( double * );

double **AllocateFloatHalfMtx( int );
double **AllocateFloatMtx( int, int );
void FreeFloatHalfMtx( double **, int );
void FreeFloatMtx( double ** );

double **AlocateFloatTri( int );
void FreeFloatTri( double ** );

int *AllocateIntVec( int );
int *AllocateIntVecLarge( unsigned long long );
void FreeIntVec( int * );

int **AllocateIntMtx( int, int );
int **AllocateIntMtxLarge( unsigned long long, unsigned long long );
void FreeIntMtx( int ** );

char ***AllocateCharCub( int, int, int );
void FreeCharCub( char *** );

int ***AllocateIntCub( int, int, int );
void FreeIntCub( int *** );

double *AllocateDoubleVec( int );
void FreeDoubleVec( double * );

double **AllocateDoubleHalfMtx( int );
double **AllocateDoubleMtx( int, int );
void FreeDoubleHalfMtx( double **, int );
void FreeDoubleMtx( double ** );

double ***AllocateDoubleCub( int, int, int );
void FreeDoubleCub( double *** );

double ***AllocateFloatCub( int, int, int );
void FreeFloatCub( double *** );

short *AllocateShortVec( int );
void FreeShortVec( short * );

short **AllocateShortMtx( int, int );
void FreeShortMtx( short ** );

void freeintmtx( int **, int );
