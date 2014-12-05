#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define grid 128

int main (void)
{
	double *u, *u1, *u2, *umid;
	double *v, *v1, *v2, *vmid;
	double *p, *p1, *pc;
	double *m, *realu, *realv, *apu, *apv;
	int i, j, step;
	double Re, nu, omega, dx, dy, De, Dw, Dn, Ds, AE, AW, AN, AS, Ce, Cw, Cn, Cs, Ce1, Cw1, Cn1, Cs1, Aep, Awp, Anp, Asp, App, error, alpha, k1 ;
	size_t NG = grid*grid*sizeof(double);
	
// Memory allocation

	u = (double*)malloc(NG);
	v = (double*)malloc(NG);
	p = (double*)malloc(NG);
	u1 = (double*)malloc(NG);
	v1 = (double*)malloc(NG);
	p1 = (double*)malloc(NG);
	umid = (double*)malloc(NG);
	vmid = (double*)malloc(NG);
	pc = (double*)malloc(NG);
	realu = (double*)malloc(NG);
	realv = (double*)malloc(NG);
	m = (double*)malloc(NG);
	apu = (double*)malloc(NG);
	apv = (double*)malloc(NG);
	
//Initialisation

	for (i=0; i<grid; i++)
	{
		for (j=0; j<grid; j++)
		{
			int k=j*grid+i;
			if ( j==(grid-1) && i>0 && i<(grid-1) )
			{
				u[k]=1.;
				realu[k]=1.;
				v[k]=0.;
				p[k]=1.;
			}
			else
			{
				u[k]=0.;
				u1[k]=0.;
				realu[k]=0.;
				v[k]=0.;
				v1[k]=0.;
				realv[k]=0.;
				p[k]=1.;
			}
		}
	}
	

		
	printf("ENTER Reynolds NUMBER\n");
	scanf("%lf", &Re);
	omega=0.999;
	error=1.;
	step=1;
	
	while (error > 0.000000001)
	{
		alpha=0.001;
		dx=1./grid;
		dy=dx;
		k1=0.8;
		
	// u-velocity
		
		for (i=2; i<(grid-1); i++)
		{
			for (j=1; j<(grid-1); j++)
			{
				int k=j*grid+i;
				De=(dy/dx)/Re;
				Dw=(dy/dx)/Re;
				Dn=(dx/dy)/Re;
				Ds=(dx/dy)/Re;
		
				Ce= dy*( u[k+1] + u[k]    )/2.0;
				Cw= dy*( u[k]   + u[k-1]  )/2.0;
				Cn= dx*( v[k+grid] + v[k+grid-1])/2.0;
				Cs= dx*( v[k]   + v[k-1]  )/2.0;
    
				if(j == 1)
				{
                   Ds=2.0*Ds;
				}

				if(j == (grid-2))
				{
                   Dn=2.0*Dn;
				}
    
				AE=De-Ce/2.0;
				AW=Dw+Cw/2.0;
				AN=Dn-Cn/2.0;
				AS=De+Cs/2.0;
				apu[k]=AE + AW + AN + AS;
    
				// procedure1 solve momentum eq. !
     
				u1[k]=alpha*( -(dy)*(p[k]-p[k-1]) + AE*u[k+1] + AW*u[k-1] + AN*u[k+grid] + AS*u[k-grid] )/apu[k]+(1-alpha)*u[k];
			}
		}
	
	// v-velocity
    
        for (i=1; i<(grid-1); i++)
		{
			for (j=2; j<(grid-1); j++)
			{
				int k=j*grid+i;
				De=(dy/dx)/Re;
				Dw=(dy/dx)/Re;
				Dn=(dx/dy)/Re;
				Ds=(dx/dy)/Re;
		
				Ce1=  dy*( u[k+1]+ u[k+1-grid])/2.0;
				Cw1=  dy*( u[k]  + u[k-grid]  )/2.0;
				Cn1=  dx*( v[k+grid]+ v[k]    )/2.0;
				Cs1=  dx*( v[k]  + v[k-grid]  )/2.0;
    
				if(i == 2)
				{
                   Dw=2.0*Dw; 
				}

				if(i == (grid-2))
				{
                   De=2.0*De;
				}
    
				AE=De-(Ce1/2.0);
				AW=Dw+(Cw1/2.0);
				AN=Dn-(Cn1/2.0);
				AS=Ds+(Cs1/2.0);
				apv[k]=AE + AW + AN + AS;
				
	 // procedure1 solve momentum eq !

				v1[k]=alpha*(-(dy)*(p[k]-p[k-grid])+ AE*v[k+1]+AW*v[k-1]+AN*v[k+grid]+AS*v[k-grid])/apv[k]+(1-alpha)*v[k];      	 
			} 
		}

		
	error=0.0;   
				
				
      // procedure3 solve pressure correction!
        
		for (i=0; i<grid; i++)
		{
			for (j=0; j<grid; j++)
			{
					int k=j*grid+i;
					pc[k]=0.;
			}
		}
		
					
		for (i=1; i<(grid-1); i++)
		{
			for (j=1; j<(grid-1); j++)
				{		   
					int k=j*grid+i;
					Aep=(dy*dy)/apu[k+1];   
					Awp=(dy*dy)/apu[k];   
					Anp=(dx*dx)/apv[k+grid];  
					Asp=(dx*dx)/apv[k];
					
					if (i==(grid-2) )
					{
						Aep=0.;
					}
					
					if(i==1)
					{
						Awp=0.; 
					}

					if(j==(grid-2))
					{
						Anp=0.;
					}

					if(j==1)
					{
						Asp=0.;  
					}
            
					
					App=Aep+Awp+Anp+Asp;
					//printf("App is %lf for i %d and j %d\n", App, i, j);
					
					m[k] = (u1[k]-u1[k+1])*dy + (v1[k]-v1[k+grid])*dx;
					//printf("m is %lf\n", m[k]);
					pc[k]=(Aep*pc[k+1]+Awp*pc[k-1]+Anp*pc[k+grid]+Asp*pc[k-grid]+m[k])/App;
					
					error=error+fabs(m[k]);
					
				}
        }
             
		if (step%1000 ==1)
		{
	    printf("Error is %5.8lf for the step %d\n", error, step);
		}
		
		for (i=2; i<(grid-1); i++)
		{
			for (j=1; j<(grid-1); j++)
			{
				int k=j*grid+i;
				umid[k]=(dy/apu[k])*(pc[k-1]-pc[k]);
				u[k]=u1[k]+umid[k];
			}
		}
  
 
		for (i=1; i<(grid-1); i++)
		{
			for (j=2; j<(grid-1); j++)
			{
				int k=j*grid+i;
				vmid[k]=dx*(pc[k-grid]-pc[k])/apv[k];
				v[k]=v1[k]+vmid[k];
			}
		}
	
	// procedure4 correct pressure !
	
		for (i=1; i<grid; i++)
		{
			for (j=1; j<grid; j++)
			{	
				int k=j*grid+i;
				p[k]=p[k]+k1*pc[k];
			}
		}
		
		step=step+1;
	}
	
	
	for (i=1; i<grid; i++)
	{
		for (j=1; j<grid; j++)
		{
			int k=j*grid+i;
			if (j== (grid-1))
			{
			realu[k]=1.0;
			}
		}
	}
	
	for (i=1; i<grid; i++)
	{
		for (j=1; j<grid; j++)
		{
			int k=j*grid+i;
			realu[k]=(u[k]+u[k+1])/2.;
			realv[k]=(v[k]+v[k+grid])/2.;
		}
    }
	
	// OUTPUT DATA
	FILE *fout2, *fout3, *fout4;
	fout2 = fopen("UVP.dat","w+t");
	fout3 = fopen("Central U.dat","w+t");
	fout4 = fopen("Central V.dat","w+t");
	
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

	for ( j = 0 ; j < grid ; j++ )
	{
    for ( i = 0 ; i < grid ; i++ )
    {
		int k = j*grid + i;
		double dx, dy, xpos, ypos;
		dx=1./grid;
		dy=dx;
		xpos = i*dx;
		ypos = j*dy;
		
		fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, realu[k], realv[k], p[k] );
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
    int k1 = j*grid + grid/2;
    int k2 = j*grid + grid/2 - 1;
	double dy, ypos;
	dy=1./grid;
    ypos = (double) j*dy;

    fprintf( fout3, "%5.8lf\t%5.8lf\n", (realu[k1] + realu[k2])/(2.), ypos );
  }
  
  
  // CENTRAL --V
  fprintf(fout4, "VARIABLES=\"Y / L\",\"V\"\n");
  fprintf(fout4, "ZONE F=POINT\n");
  fprintf(fout4, "I=%d\n", grid );

  for ( i = 0 ; i < grid ; i++ )
  {
    int k1 = (grid/2)*grid + i;
    int k2 = (grid/2-1)*grid + i;
	double dx, xpos;
	dx=1./grid;
    xpos = (double) i*dx;

    fprintf( fout4, "%5.8lf\t%5.8lf\n", xpos, (realv[k1] + realv[k2])/(2.) );
  }
  
}