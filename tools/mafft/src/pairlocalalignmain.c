#include "mltaln.h"
int main( int argc, char **argv )
{
	int res = pairlocalalign( 0, 0, NULL, NULL, NULL, NULL, argc, argv, NULL );
	if( res == GUI_CANCEL ) res = 0; // treeout de goto chudan wo riyousuru
	return res;
}
