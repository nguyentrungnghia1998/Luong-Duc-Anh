#include "lbm_D2Q9.h"


void D2Q9_streaming_step(double ***f, double ***fpc, double **type)
{
	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		for(int j=1; j<jm+1; j++){

			if(type[i][j]==FIELD) // 1~8
			{
				f[i][j][0]=fpc[i][j][0];
				if(type[i-1][j-1]==FIELD){ f[i][j][1]=fpc[i-1][j-1][ 1]; }else if(type[i-1][j-1]==BOUND){ f[i][j][1]=fpc[i][j][8]; }
				if(type[i-1][j  ]==FIELD){ f[i][j][2]=fpc[i-1][j  ][ 2]; }else if(type[i-1][j  ]==BOUND){ f[i][j][2]=fpc[i][j][7]; }
				if(type[i-1][j+1]==FIELD){ f[i][j][3]=fpc[i-1][j+1][ 3]; }else if(type[i-1][j+1]==BOUND){ f[i][j][3]=fpc[i][j][6]; }
				if(type[i  ][j-1]==FIELD){ f[i][j][4]=fpc[i  ][j-1][ 4]; }else if(type[i  ][j-1]==BOUND){ f[i][j][4]=fpc[i][j][5]; }
				if(type[i  ][j+1]==FIELD){ f[i][j][5]=fpc[i  ][j+1][ 5]; }else if(type[i  ][j+1]==BOUND){ f[i][j][5]=fpc[i][j][4]; }
				if(type[i+1][j-1]==FIELD){ f[i][j][6]=fpc[i+1][j-1][ 6]; }else if(type[i+1][j-1]==BOUND){ f[i][j][6]=fpc[i][j][3]; }
				if(type[i+1][j  ]==FIELD){ f[i][j][7]=fpc[i+1][j  ][ 7]; }else if(type[i+1][j  ]==BOUND){ f[i][j][7]=fpc[i][j][2]; }
				if(type[i+1][j+1]==FIELD){ f[i][j][8]=fpc[i+1][j+1][ 8]; }else if(type[i+1][j+1]==BOUND){ f[i][j][8]=fpc[i][j][1]; }
			}
		}
	}
}
