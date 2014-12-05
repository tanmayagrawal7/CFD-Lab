#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define grid 128

__global__ void conduction (float *dev_a) 
{
	
	int x = blockIdx.y * blockDim.y + threadIdx.y;
	int y = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(x==0||x==grid-1||y==0||y==grid-1)
	{
	dev_a[x*grid+y]=dev_a[x*grid+y];
	}
	
	else
	{
	float Value;
	Value = 0.25f*( dev_a[x*grid + y + 1] + dev_a[x*grid + y-1] + dev_a[grid + x*grid + y] + dev_a[x*grid + y - grid] );
	dev_a[x*grid+y]=Value;
	}
}
	

int main()
{
	float T[grid][grid];
	float temp[grid*grid];
	for(int i=0;i<grid;i++)
		{
		for(int j=0;j<grid;j++)
			{
			T[i][j]=0.0;
			}
		}
	for ( int e = 0 ; e < grid ; e++ )
		{
		T[0][e] = 1.0;
		}

	float *dev_a,*output;

	float size = grid * grid * sizeof(float);

	cudaMalloc( (void**)&dev_a, size );
	
	FILE *fp;                                  
	fp = fopen("124.dat","w");
	
	cudaMemcpy(dev_a,T,size,cudaMemcpyHostToDevice);
	
	int step=1;
	while(step<10000)
	{
	dim3 dimBlock(128,1);
    dim3 dimGrid(1,128);
	
	//Kernel launch
	conduction<<<dimGrid, dimBlock>>>(dev_a);
	
	cudaMemcpy(output,dev_a,size,cudaMemcpyDeviceToDevice);
	cudaFree(dev_a);
	cudaMalloc( (void**)&dev_a, size );
	cudaMemcpy(dev_a,output,size,cudaMemcpyDeviceToDevice);
	
	step++;
	
	}
	
	cudaMemcpy(temp,dev_a,size,cudaMemcpyDeviceToHost);
	
		fprintf(fp, "VARIABLES=\"X\",\"Y\",\"Texact\"\nZONE F=POINT\nI=%d,J=%d\n",grid,grid);
	
		for (int i = 0; i < grid ; i++)
		{
			for( int j =0; j < grid ; j++ )
			{

				fprintf(fp,"%-6d%-6d%-15.10f\n",i,j,temp[i*grid+j]);
			}
		}

		printf("writing\n\n\n");

}