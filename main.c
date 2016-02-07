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
				sqn,
				kiplusone;

	mpz_t			*k;

	double 		*timearray;
	double		begin, end;
	double 		time_spent;
	
	char 		*ptr;
	long int 	N = strtol(argv[1],&ptr,10);
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
	mpz_init(kiplusone);

if(!my_rank) printf("N=%ld\n", N);
	mpz_set_str(n,argv[1],10);
if(!my_rank) gmp_printf("n=%Zd\n",n);
	mpz_sqrt(sqn,n);
if(!my_rank) gmp_printf("sqn:%Zd\n", sqn);
	mpz_set(gap, sqn);
if(!my_rank) gmp_printf("gap:%Zd\n", gap);

	mpz_set_ui(p,1);
	mpz_nextprime(q, p);


	size_t mag = mpz_sizeinbase (gap, 10);
	mag=mag*procs;

if(!my_rank) printf("mag:%zu\n", mag);





	k = (mpz_t*)malloc(sizeof(mpz_t)*(mag+1));




	
	mpz_t temp;
	mpz_init(temp);
	mpz_tdiv_q_ui(temp,gap,mag);	

//unsigned long int spacing = mpz_get_ui(temp);

	for (int i=0;i<=mag;i++)
	{
		mpz_init(k[i]);
		mpz_mul_ui(k[i],temp,i);
		//k[i] = spacing*i;
	}
	//k[mag] = mpz_get_ui(sqn);;
	mpz_set(k[mag],sqn);


if (!my_rank){
	for (int i =0;i<mag;i++)
		gmp_printf("k[%d]=%Zd\n", i,k[i+1]);
}



	int counter=0;

int tem=0;

	for (int i=my_rank; i<mag; i=i+procs)
	{
		//mpz_sub_ui(kiplusone,k[i+1],1);

		mpz_set(q,k[i]);

		mpz_sub_ui(q,q,1);
		mpz_nextprime(q,q);
gmp_printf("P%d n%Zd q%Zd\n", my_rank, n,q);

		while (( mpz_cmp(q,k[i+1]) <= 0 )&&(!done)&&(!found)&&(mpz_cmp(q,sqn)<=0))
		{//finding the prime numbers
			counter++;
			if (counter%500==0)MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
if (counter%500==0)gmp_printf("q:%Zd k[%d+1]%Zd done%d, found%d\n", q,i,k[i+1],done,found);

			if (mpz_divisible_p(n,q)==0) 
			{//if n is not divisible by q
//gmp_printf("P%d found that %Zd is NOT divisible by %Zd\n", my_rank, n,q);
				mpz_nextprime(q,q);
				continue;
			}
if (tem==0)
{
gmp_printf("P%d found that %Zd is divisible by %Zd\n", my_rank, n,q);
tem=1;
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
tem=0;
			}
			///end try
/*
			mpz_nextprime(p,p);

			mpz_mul(pq,p,q);

			if (mpz_cmp(pq,n)==0)
			{
				found = 1;
				done=1;
			}
			else if ( mpz_cmp(pq,n)>0)
			{
				mpz_set(p,q);
				mpz_nextprime(q,q);
tem=0;
			}*/
		}//done finding primes

		
		if (found || done) break;
		
	}
	
	end = MPI_Wtime();



	if (found)
	{
		MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		gmp_printf("P%d: Finished: p=%Zd q=%Zd*************************\n", my_rank,p,q);
	}
	else if (!done)
	{
		while (!done&&!found)
			MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	}



	time_spent = (double)(end - begin);





	if (my_rank!=0)
	{
printf("P%d sending\n", my_rank);
		MPI_Send(&time_spent, sizeof(double),MPI_CHAR,0,0,MPI_COMM_WORLD);
printf("P%d sent\n", my_rank);
	}
	else
	{
		timearray = (double*)malloc(sizeof(double)*procs);
		
		
		timearray[0] = time_spent;


printf("timearray[0]=%le\n", timearray[0]);

		int j; 
		for (j = 1; j<procs;j++)
		{
			double *ptr = timearray+j;
			MPI_Recv(ptr, sizeof(double),MPI_CHAR,j,0,MPI_COMM_WORLD,&status);
printf("timearray[%d]=%le\n", j,timearray[j]);
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

