#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>


void writeTime(char *n, double *timearray, int j);

int main(int argc, char** argv)
{
	int 		my_rank,
				procs;

	MPI_Status	status;

	mpz_t 		q,
				p,
				pq,
				n,
				gap,
				sqn;

	mpz_t			*k;

	double 		*timearray;
	double		begin, end;
	double 		time_spent;
	
	int 		done = 0,
				found = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);


	begin = MPI_Wtime();

	mpz_init(n);
	mpz_init(q);
	mpz_init(p);
	mpz_init(pq);
	mpz_init(sqn);
	mpz_init(gap);

	mpz_set_str(n,argv[1],10);
	mpz_sqrt(sqn,n);
	mpz_set(gap, sqn);

	mpz_set_ui(p,1);
	mpz_nextprime(q, p);


	size_t mag = mpz_sizeinbase (gap, 10);
	mag=mag*procs;

	k = (mpz_t*)malloc(sizeof(mpz_t)*(mag+1));
	
	mpz_t temp;
	mpz_init(temp);
	mpz_tdiv_q_ui(temp,gap,mag);	

	for (int i=0;i<=mag;i++)
	{
		mpz_init(k[i]);
		mpz_mul_ui(k[i],temp,i);
	}
	mpz_set(k[mag],sqn);





	int counter=0;


	for (int i=my_rank; i<mag; i=i+procs)
	{
		mpz_set(q,k[i]);

		mpz_sub_ui(q,q,1);
		mpz_nextprime(q,q);


		while (( mpz_cmp(q,k[i+1]) <= 0 )&&(!done)&&(!found)&&(mpz_cmp(q,sqn)<=0))
		{//finding the prime numbers
			counter++;
			if (counter%5000==0)MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			if (mpz_divisible_p(n,q)==0) 
			{//if n is not divisible by q
				mpz_nextprime(q,q);
				continue;
			}

			//since it is divisible, try n/q and see if result is prime
			mpz_divexact(p,n,q);
			int reps;
			if (mpz_probab_prime_p(p,reps)!=0)
			{
				found = 1;
				done = 1;
			}
			else
			{
				mpz_nextprime(q,q);
			}
		}//done finding primes

		
		if (found || done) break;
		
	}
	
	end = MPI_Wtime();



	if (found)
	{
		MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		gmp_printf("*************************\nP%d: Finished\np*q=n\np=%Zd q=%Zd n=%Zd\n*************************\n", my_rank,p,q,n);
	}
	else if (!done)
	{
		while (!done&&!found)
			MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	}



	time_spent = (double)(end - begin);





	if (my_rank!=0)
	{
		MPI_Send(&time_spent, sizeof(double),MPI_CHAR,0,0,MPI_COMM_WORLD);
	}
	else
	{
		timearray = (double*)malloc(sizeof(double)*procs);
		
		
		timearray[0] = time_spent;

		int j; 
		for (j = 1; j<procs;j++)
		{
			double *ptr = timearray+j;
			MPI_Recv(ptr, sizeof(double),MPI_CHAR,j,0,MPI_COMM_WORLD,&status);
		}


		writeTime(argv[1],timearray,procs);
		free(timearray);

	}


	
	free(k);

	MPI_Finalize();
	return 0;
}



void writeTime(char *n, double *timearray, int j)
{
  int i;

// open file for writing
  char filename[100];
  sprintf(filename, "time_%s",n);

  FILE *fd;
  fd = fopen(filename, "w");

// write the image
  for(i = 0; i < j; i++)
  {
    fprintf(fd, "%d\t%le\n", i,timearray[i]);
  }
  fclose(fd);
}

