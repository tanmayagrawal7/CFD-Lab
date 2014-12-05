// Code written by Tanmay Agrawal for simulation of two dimensional heat conduction problem with second order finite difference scheme.
// Left, right and bottom walls are at a temperature of 20 units while the top wall has a sinusoidal temperature distribution.
// Simulated with CUDA-C

#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define grid 512

#define THREADX 128
#define THREADY 1


void init( double *T, double *Texact )
{
	int i, j, k;
	double pi;
	pi=acos(-1);
	for (i=0; i<grid; i++)
	{
		for (j=0; j<grid; j++)
		{
			k=j*grid+i;
			if (i==0 && j>=0 && j<grid)
			{
			T[k]=20.;
			}
			
			else if (i==(grid-1) && j>=0 && j<grid)
			{
			T[k]=20.;
			}
			
			else if (j==0 && i>=0 && i<grid)
			{
			T[k]=20.;
			}
			
			else if (j==(grid-1) && i>0 && i<(grid-1))
			{
			T[k]=20.+80.*sin((pi*i)/grid);
			}
			else
			{
			T[k]=20.;
			}
			
			Texact[k]=20. + (80.*sin((pi*i)/grid)*sinhf((pi*j)/grid))/(sinhf(pi));
		}
	}
}

__global__ void calculation1 (double *T, double *Tnew)
{
	int i, j, k;
	i = blockDim.x * blockIdx.x + threadIdx.x;
	j = blockDim.y * blockIdx.y + threadIdx.y;
	//printf("Inside the kernel\n");
	k=j*grid+i;
	
	if (i>0 && i<(grid-1) && j>0 && j<(grid-1) )
	{
	Tnew[k] = 0.25*(T[k-1] + T[k+1] + T[k-grid] + T[k+grid]);
	//printf("gets into loop");
	}
	else
	{
	Tnew[k] = T[k];
	}
}


__global__ void calculation2 (double *T, double *Tnew)
{
	int i, j, k;
	i = blockDim.x * blockIdx.x + threadIdx.x;
	j = blockDim.y * blockIdx.y + threadIdx.y;
	k=j*grid+i;
	T[k]=Tnew[k];
}


int main (void)
{
	double *T_h, *Tnew_h, *Texact;
	double *T_d, *Tnew_d, terror;
	int i, j, k, timeStep, maxStep;
	
	size_t NG = grid*grid*sizeof(double);

	dim3 dimGrid( (grid/THREADX), (grid/THREADY));
	dim3 dimBlock( THREADX, THREADY );
	
	
	// LOCATE MEMORY--HOST
	T_h = (double*)malloc(NG);
	Tnew_h = (double*)malloc(NG);
	Texact = (double*)malloc(NG);
	
	
	// LOCATE MEMORY--DEVICE
	cudaMalloc((void**)&T_d, NG);
	cudaMalloc((void**)&Tnew_d, NG);
	
	// INITIALIZATION
	init( T_h, Texact );
	
	// File writing after initialisation
	FILE *fout1;
	fout1 = fopen("Initialisation.dat","w+t");
	if ( fout1 == NULL )
	{
    printf("\nERROR when opening file\n");
    fclose( fout1 );
	}

  else
	{
	fprintf( fout1, "VARIABLES=\"T\"\n");
	fprintf( fout1, "I=%d, J=%d\n", grid, grid );

	for ( j = 0 ; j < grid ; j++ )
	{
    for ( i = 0 ; i < grid ; i++ )
    {
	k = j*grid + i;
		fprintf( fout1, "%5.8lf\n", T_h[k] );
    }
	}
	}
	
	fclose( fout1 );
	
	
	// COPY MEMORY FROM HOST TO DEVICE
	cudaMemcpy( T_d, T_h, NG, cudaMemcpyHostToDevice );
	cudaMemcpy( Tnew_d, T_h, NG, cudaMemcpyHostToDevice );
	
	printf("ENTER MAXIMUM SIMULATION STEPS\n");
	scanf("%d", &maxStep);
	
	
	
	// START ITERATION
  
	for ( timeStep = 1 ; timeStep <= maxStep ; timeStep++ )
	{

	
	//Kernel launch
	calculation1<<<dimGrid, dimBlock>>>( T_d, Tnew_d   );
	
	cudaThreadSynchronize();
	calculation2<<<dimGrid, dimBlock>>>( T_d, Tnew_d );
	
	cudaThreadSynchronize();

	calculation1<<<dimGrid, dimBlock>>>( T_d , Tnew_d  );
	
	cudaThreadSynchronize();
	
	
	
	
	if (timeStep % 10000==3)
        {
			//calculate error
			cudaMemcpy(Tnew_h, T_d    ,NG , cudaMemcpyDeviceToHost);
			cudaMemcpy(T_h   , Tnew_d ,NG , cudaMemcpyDeviceToHost);	

			
			
			
			terror = 0;
			for ( j = 0 ; j < grid  ; j++ )
			{
				for ( i = 0 ; i < grid  ; i++ )
				{
				k = j * grid + i;
				terror += fabs( T_h[k]-Tnew_h[k] );
				}
			}
	        printf("timestep = %d \n", timeStep);
		    printf("T_error  = %e \n", terror);
			
			if ( (terror < 1E-6) ) goto END;	
        }
	
	}
	END:
	cudaThreadSynchronize();
  
	// COPY MEMORY FROM DEVICE TO HOST
	cudaMemcpy( Tnew_h, T_d, NG, cudaMemcpyDeviceToHost );
	
	// OUTPUT DATA
	FILE *fout2, *fout3, *fout4;
	fout2 = fopen("TemperatureDistribution.dat","w+t");
	fout3 = fopen("CentreLineTemperature.dat","w+t");
	fout4 = fopen("CentreLineExact.dat","w+t");
	
	if ( fout2 == NULL )
	{
    printf("\nERROR when opening file\n");
    fclose( fout2 );
	}

  else
	{
	fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"T\"\n");
	fprintf( fout2, "ZONE  F=POINT\n");
	fprintf( fout2, "I=%d, J=%d\n", grid, grid );

	for ( j = 0 ; j < grid ; j++ )
	{
    for ( i = 0 ; i < grid ; i++ )
    {
		k = j*grid + i;
		double dx, dy, xpos, ypos;
		dx=1.0/grid;
		dy=dx;
		xpos = i*dx;
		ypos = j*dy;
		
		fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, Tnew_h[k] );
    }
	}
	}
	
	fclose( fout2 );
	
	
	
	
	// CENTRAL TEMP--T
  fprintf(fout3, "VARIABLES=\"X / L\",\"T\"\n");
  fprintf(fout3, "ZONE F=POINT\n");
  fprintf(fout3, "I=%d\n", grid );

  for ( i = 0 ; i < grid ; i++ )
  {
    int k1 = (grid/2)*grid + i;
    int k2 = (grid/2-1)*grid + i;
	double dx, xpos;
	dx=1.0/grid;
    xpos = (double) i*dx;

    fprintf( fout3, "%5.8lf\t%5.8lf\n", xpos, (Tnew_h[k1] + Tnew_h[k2])/(2.0) );
  }
  
  
  
  // EXACT SOLUTION CENTERLINE
  fprintf(fout4, "VARIABLES=\"X / L\",\"Texact\"\n");
  fprintf(fout4, "ZONE F=POINT\n");
  fprintf(fout4, "I=%d\n", grid );
  
  for ( i = 0 ; i < grid ; i++ )
  {
    int k1 = (grid/2)*grid + i;
    int k2 = (grid/2-1)*grid + i;
	double dx, xpos;
	dx=1.0/grid;
    xpos = (double) i*dx;
	
    fprintf( fout4, "%5.8lf\t%5.8lf\n", xpos, (Texact[k1] + Texact[k2])/(2.0) );
  }
  
  
  
	// FREE MEMORY
	
	cudaFree( T_d );
  
  
}