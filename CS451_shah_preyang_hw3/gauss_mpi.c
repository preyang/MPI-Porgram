#include <mpi.h>
#include <sys/times.h>
#include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

void MatrixIntialization();
void backSubstitution();
void DisplayMatrixValue();
void gaussian_mpi(int N);
int proc,id,N;
double X[2000][2000],y[2000],Z[2000];
void MatrixIntialization()
{
	int i,j;
	for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++)
		{
						
			X[i][j]=rand()/1000.0;
			
	//	scanf("%lf",&X[i][j]);
				
		}
	y[i] = rand()/32768.0;

	}
}

void DisplayMatrixValue()
{	
	int i,j;
	printf("Displaying Initial Matrix.\n");
	for (i=0;i<N;i++)
	{
		printf("| ");
		for(j=0;j<N; j++)
		{	
			printf("%lf ",X[i][j]);
		}
		printf("| | %lf |\n",y[i]);
	}
}

void backSubstitution()
{
	int i,j;
	for (i=N-1;i>=0;i--)
	{
		int count=0;
		Z[i] = y[i];
		for (j=i+1;j<N;j++)
		{
			Z[i]-=X[i][j]*Z[j];
			
		}
		Z[i] = Z[i]/X[i][i];
	}
}

void gaussian_mpi(int N)
{	
	double wp_time,wa_time=0;
	MPI_Request request;
	MPI_Status status;
	int p,k,i,j;
	float mp;	

	MPI_Barrier(MPI_COMM_WORLD);	
	if(id==0)
	{
		wa_time = MPI_Wtime();
	}	
	for (k=0;k<N-1;k++)
 	{	
		MPI_Bcast(&X[k][0],N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&y[k],1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(id==0)
		{
		for (p=1;p<proc;p++)
			{
		  		for (i=k+1+p;i<N;i+=proc)
		  		{
				   MPI_Isend(&X[i], N, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request);
				   MPI_Wait(&request, &status);
				   MPI_Isend(&y[i], 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request);
				   MPI_Wait(&request, &status);
		  		}
			}
			for (i=k+1 ; i<N ; i += proc)
			{
	  			mp = X[i][k] / X[k][k];
	  			for (j = k; j < N; j++)
	 			{
	   				X[i][j] -= X[k][j] * mp;
	 			}
	   			y[i] -= y[k] * mp;
			}
			/* MPI_Barrier(MPI_COMM_WORLD);*/
			for (p = 1; p < proc; p++)
			{
			  for (i = k + 1 + p; i < N; i += proc)
			  {
			    MPI_Recv(&X[i], N, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &status);
			    MPI_Recv(&y[i], 1, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, &status);
			  }
			}
			
			if (k == N - 2)
			{
 				wp_time = MPI_Wtime();
 				printf("elapsed time = %f\n", wp_time-wa_time);
			}
		}
		else
		{
		
			for (i = k + 1 + id; i < N; i += proc)
			{
				MPI_Recv(&X[i], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);	
				MPI_Recv(&y[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);				
				mp = X[i][k] / X[k][k];
				for (j = k; j < N; j++)
				{
				    X[i][j] -= X[k][j] * mp;
				}
				y[i] -= y[k] * mp;
				MPI_Isend(&X[i], N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);						    
				MPI_Wait(&request, &status);		
				MPI_Isend(&y[i], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
				MPI_Wait(&request, &status);
			}
		}
		 MPI_Barrier(MPI_COMM_WORLD);
	}
}
int main(int argc,char *argv[])
{
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&proc);	
	MPI_Status status;
	if (argc >= 2) 
	{
	    N = atoi(argv[1]);
   	}
	struct timeval etstart, etstop;  
 	struct timezone tzdummy;
 	clock_t StartElapsedTime, StopElapsedTime;
  	unsigned long long usecstart, usecstop;
  	struct tms StartCPUTime, StopCPUTime;


	if(id==0)
	{
		printf("\n matrix Initialized working on send and recieve operation\n");
		MatrixIntialization();
 		
		/*printf("\n matrix Display\n");
		DisplayMatrixValue();*/
		
  		gettimeofday(&etstart, &tzdummy);
  		StartElapsedTime = times(&StartCPUTime);					
	}
	
	gaussian_mpi(N);
	
	if(id==0)
	{
		backSubstitution();
  		gettimeofday(&etstop, &tzdummy);
  		StopElapsedTime = times(&StopCPUTime);
  		printf("Stopped clock.\n");
  		usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
  		usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;
	

  		printf("\nElapsed time = %g ms.\n",
	 	(float)(usecstop - usecstart)/(float)1000);

  		printf("(CPU times are accurate to the nearest %g ms)\n",
	 	1.0/(float)CLOCKS_PER_SEC * 1000.0);
	
  		printf("My total CPU time for parent = %g ms.\n",
	 	(float)( (StopCPUTime.tms_utime + StopCPUTime.tms_stime) -
		  (StartCPUTime.tms_utime + StartCPUTime.tms_stime) ) /
	 	(float)CLOCKS_PER_SEC * 1000);
  	
		printf("My system CPU time for parent = %g ms.\n",
	 	(float)(StopCPUTime.tms_stime - StartCPUTime.tms_stime) /
	 	(float)CLOCKS_PER_SEC * 1000);
  	
		printf("My total CPU time for child processes = %g ms.\n",
	 	(float)( (StopCPUTime.tms_cutime + StopCPUTime.tms_cstime) -
		  (StartCPUTime.tms_cutime + StartCPUTime.tms_cstime) ) /
	 	(float)CLOCKS_PER_SEC * 1000);
     	
		  /* Contrary to the man pages, this appears not to include the parent */
		 printf("--------------------------------------------\n");

	}
	MPI_Finalize(); 
  	return 0;
}
