#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define grid 128

int main (void)
{
	double u[grid][grid+1], um[grid][grid+1], uold[grid][grid+1], uolder[grid][grid+1], uwrite[grid][grid+1];
	double v[grid+1][grid], vm[grid+1][grid], vold[grid+1][grid], volder[grid+1][grid], vwrite[grid+1][grid];
	double p[grid+1][grid+1], pm[grid+1][grid+1], pold[grid+1][grid+1], polder[grid+1][grid+1], pwrite[grid+1][grid+1];
	double m[grid+1][grid+1], m2[grid+1][grid+1];
	int i, j, step1, step2;
	double dx, dy, dt, tau, delta, error1, error2, Re;
	dx = 1.0/(grid-1);
	dy = 1.0/(grid-1);
	dt = 0.0005;
	tau = 0.0001;
	delta = 15.5;
	Re = 1000.0;
	
	// Initializing u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				u[i][j] = 0.0;
				u[i][grid] = 1.0;
				u[i][grid-1] = 1.0;
				um[i][j] = 0.0;
				um[i][grid] = 1.0;
				um[i][grid-1] = 1.0;
				uold[i][j] = 0.0;
				uold[i][grid] = 1.0;
				uold[i][grid-1] = 1.0;
				uolder[i][j] = 0.0;
				uolder[i][grid] = 1.0;
				uolder[i][grid-1] = 1.0;
			}
		}
		
	// Initializing v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{
				v[i][j] = 0.0;
				vm[i][j] = 0.0;
				vold[i][j] = 0.0;
				volder[i][j] = 0.0;
			}
		}
		
	// Initializing p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				p[i][j] = 1.0;
				pm[i][j] = 1.0;
				pold[i][j] = 1.0;
				polder[i][j] = 1.0;
			}
		}
	
	error2 = 1.0;
	step2 =1;
	while (step2 < 500000)
	{
		step1=1;
		error1 = 1.0;
		while (error1 > 0.0000001)
		{
		// Solve u-momentum
		for (i=1; i<=(grid-2); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				u[i][j] = (2*dt*um[i][j] + 4*tau*uold[i][j] - tau*uolder[i][j]
							- (2*dt*tau)*((um[i+1][j]*um[i+1][j]-um[i-1][j]*um[i-1][j])/2.0/dx 
							+0.25*( (um[i][j]+um[i][j+1])*(vm[i][j]+vm[i+1][j])-(um[i][j]+um[i][j-1])*(vm[i+1][j-1]+vm[i][j-1]) )/dy  )
								- (2*dt*tau)/dx*(pm[i+1][j]-pm[i][j]) 
									+ (2*dt*tau)*1.0/Re*( (um[i+1][j]-2.0*um[i][j]+um[i-1][j])/dx/dx +(um[i][j+1]-2.0*um[i][j]+um[i][j-1])/dy/dy ))/(2*dt + 3*tau);
			}
		}
		
		for (j=1; j<=(grid-1); j++)
		{
			u[0][j] = 0.0;
			u[grid-1][j] = 0.0;
		}
		
		for (i=0; i<=(grid-1); i++)
		{
			u[i][0] = -u[i][1];
			u[i][grid] = 2 - u[i][grid-1];
		}
		
		
		// Solves v-momentum
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-2); j++)
			{
				v[i][j] = (2*dt*vm[i][j] + 4*tau*vold[i][j] - tau*volder[i][j]
							- 2*dt*tau*(0.25*( (um[i][j]+um[i][j+1])*(vm[i][j]+vm[i+1][j])-(um[i-1][j]+um[i-1][j+1])*(vm[i][j]+vm[i-1][j]) )/dx 
							+(vm[i][j+1]*vm[i][j+1]-vm[i][j-1]*vm[i][j-1])/2.0/dy ) 
								- (2*dt*tau)/dy*(pm[i][j+1]-pm[i][j]) 
									+ (2*dt*tau)*1.0/Re*( (vm[i+1][j]-2.0*vm[i][j]+vm[i-1][j])/dx/dx+(vm[i][j+1]-2.0*vm[i][j]+vm[i][j-1])/dy/dy ))/(2*dt + 3*tau);
			}
		}
		
		for (j=1; j<=(grid-2); j++)
		{
			v[0][j] = -v[1][j];
			v[grid][j] = -v[grid-1][j];
		}		

		for (i=0; i<=(grid); i++)
		{
			v[i][0] = 0.0;
			v[i][grid-1] = 0.0;
		}		
	
		// Solves continuity equation
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				p[i][j] = pm[i][j]-tau*delta*(  ( u[i][j]-u[i-1][j] )/dx + ( v[i][j]-v[i][j-1] ) /dy  );
			}
		}
		
		for (i=1; i<=(grid-1); i++)
		{
			p[i][0] = p[i][1];
			p[i][grid] = p[i][grid-1];
		}
		
		for (j=0; j<=(grid); j++)
		{
			p[0][j] = p[1][j];
			p[grid][j] = p[grid-1][j];
		}		
		
		// Displaying error
		error1 = 0.0;
		
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				m[i][j] = dy*(  ( u[i][j]-u[i-1][j] ) + ( v[i][j]-v[i][j-1] )  );
				error1 = error1 + fabs(m[i][j]);
			}
		}
		
		if (step1%10000 ==1)
		{
	    printf("Error is %5.8lf for the step %d\n", error1, step1);
		}
		
		
		// Iterating u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				um[i][j] = u[i][j];
			}
		}
		
		// Iterating v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{
				vm[i][j] = v[i][j];
			}
		}
		
		// Iterating p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				pm[i][j] = p[i][j];
			}
		}	
		
		step1 = step1 + 1;	
		}
		
		
		// Iterating u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				uolder[i][j] = uold[i][j];
				uold[i][j] = u[i][j];
				um[i][j] = uold[i][j];
			}
		}
		
		// Iterating v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{
				volder[i][j] = vold[i][j];
				vold[i][j] = v[i][j];
				vm[i][j] = vold[i][j];
			}
		}
		
		// Iterating p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				polder[i][j] = pold[i][j];
				pold[i][j] = p[i][j];
				pm[i][j] = pold[i][j];
			}
		}
		
		error2=0.0;
		// Displaying error
		
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				m[i][j] = dy*(  ( uold[i][j]-uold[i-1][j] ) + ( vold[i][j]-vold[i][j-1] )  );
				error2 = error2 + fabs(m[i][j]);
			}
		}
		
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{	
			    uwrite[i][j] = 0.5*(u[i][j]+u[i][j+1]);
                vwrite[i][j] = 0.5*(v[i][j]+v[i+1][j]);
                pwrite[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
			}
		}
		
		
		if (step2%100 ==1)
		{
	    //printf("Error in second loop is %5.8lf for the step %d\n", error2, step2);
		
		char filename[64], filename2[64];
		sprintf(filename, "UVP%d.dat", step2);
		sprintf(filename2, "CentralU%d.dat", step2);
		
		FILE *file, *file2;
		file= fopen(filename, "wb");
		file2= fopen(filename2, "wb");
		
		// Write UVP
		{
		fprintf( file, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
		fprintf( file, "ZONE  F=POINT\n");
		fprintf( file, "I=%d, J=%d\n", grid, grid );

		for ( j = 0 ; j < (grid) ; j++ )
		{
		for ( i = 0 ; i < (grid) ; i++ )
		{
		double xpos, ypos;
		xpos = i*dx;
		ypos = j*dy;

		fprintf( file, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, uwrite[i][j], vwrite[i][j], pwrite[i][j] );
		}
		}
		}
		fclose(file);	
		
		// Write Central U velocity
		{
		fprintf(file2, "VARIABLES=\"U\",\"Y\"\n");
		fprintf(file2, "ZONE F=POINT\n");
		fprintf(file2, "I=%d\n", grid );

		for ( j = 0 ; j < grid ; j++ )
		{
		double ypos;
		ypos = (double) j*dy;

		fprintf( file2, "%5.8lf\t%5.8lf\n", (uwrite[grid/2][j] + uwrite[(grid/2)+1][j])/(2.), ypos );
		}
		}
		fclose(file2);	
		
		}	
		step2 = step2 + 1;
		}

	
	
	for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{	
			    u[i][j] = 0.5*(um[i][j]+um[i][j+1]);
                v[i][j] = 0.5*(vm[i][j]+vm[i+1][j]);
                p[i][j] = 0.25*(pm[i][j]+pm[i+1][j]+pm[i][j+1]+pm[i+1][j+1]);
			}
		}
	
	
	
	// OUTPUT DATA
	FILE *fout2, *fout3;
	fout2 = fopen("UVP.plt","w+t");
	fout3 = fopen("Central_U.plt","w+t");

	if ( fout2 == NULL )
	{
    printf("\nERROR when opening file\n");
    fclose( fout2 );
	}

  else
	{
	fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
	fprintf( fout2, "ZONE  F=POINT\n");
	fprintf( fout2, "I=%d, J=%d\n", grid, grid );

	for ( j = 0 ; j < (grid) ; j++ )
	{
    for ( i = 0 ; i < (grid) ; i++ )
    {
		double xpos, ypos;
		xpos = i*dx;
		ypos = j*dy;

		fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, u[i][j], v[i][j], p[i][j] );
    }
	}
	}

	fclose( fout2 );
	
	// CENTRAL --U
  fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
  fprintf(fout3, "ZONE F=POINT\n");
  fprintf(fout3, "I=%d\n", grid );

  for ( j = 0 ; j < grid ; j++ )
  {
	double ypos;
    ypos = (double) j*dy;

    fprintf( fout3, "%5.8lf\t%5.8lf\n", (u[grid/2][j] + u[(grid/2)+1][j])/(2.), ypos );
  }

}