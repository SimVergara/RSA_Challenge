#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

int nextprime(long int n);

int main(int argc, char** argv)
{
	int 		my_rank,
				procs;

	MPI_Status	status;

	int 		q,
				p,
				sqn;

	int			*k;

	double 		*timearray;

	long int 	mag;
	long int 	n = atoi(argv[1]);
	int 		done = 0,
				found = 0;

	clock_t		begin, end;
	double 		time_spent;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);


begin = clock();

	sqn = sqrt(n);
	p = 1;
	q = nextprime(sqn);
	long int gap = (n/2)-sqn;
	long int sdsd = log10(gap) -1;

	mag = (int) pow((double)10,(double)sdsd);

	k = (int*)malloc(sizeof(int)*mag);

	for (int i=0;i<mag;i++)
	{
		k[i]=sqn + (i* (int)(gap/mag));
	}
	k[mag] = n/2;


//printf("log-1=%d,  	mag=%d\n",sdsd,mag);

//for (int i=0;i<=mag;i++)
//	printf("k[%d]=%d\n", i,k[i]);



	for (int i=my_rank; i<mag; i=i+procs)
	{
//printf("P%d: i=%d\n", my_rank, i);
		q=nextprime(k[i]-1);

		while ((q <= k[i+1])&&(!done)&&(!found))
		{//finding the prime numbers
			p = nextprime(p);

			if (p*q == n)
			{
				found = 1;
			}
			else if ((p*q > n)||(p>q))
			{
				p=1;
				q=nextprime(q);
			}
		}//done finding primes

		MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		if (found || done) break;
	}

	if (found)
	{
		printf("P%d: Finished: p=%d q=%d\n", my_rank,p,q);
	}
	/*else
		printf("P%d: Numbers not found\n",my_rank);*/
	end = clock();

	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	if (my_rank!=0)
	{
		MPI_Send(&time_spent, sizeof(double),MPI_INT,0,0,MPI_COMM_WORLD);
	}
	else
	{
		timearray = (double*)malloc(sizeof(double)*procs);
		timearray[0] = time_spent;

		int j; 
		for (j = 1; j<procs;j++)
		{
			MPI_Recv(&timearray[j], sizeof(double),MPI_INT,j,0,MPI_COMM_WORLD,&status);
//			printf("Printing to file %dtime:%lf\n",j,timearray[j]);
		}

		
		writeTime(n,timearray,procs);
	}

	MPI_Finalize();
	return 0;
}



void writeTime(int n, double *timearray, int j)
{
  int i;

// open file for writing
  char filename[100];
  sprintf(filename, "time_%d",n);

  FILE *fd;
  fd = fopen(filename, "w");

// write the image
  for(i = 0; i < j; i++)
  {
    fprintf(fd, "%d\t%lf\n", i,timearray[i]);
  }
  fclose(fd);
}



int nextprime(long int n)
{
	mpz_t i;
	mpz_init(i);
	mpz_set_d(i,n);
	mpz_nextprime(i, i);
	long int k = mpz_get_d(i);
	return k;
}