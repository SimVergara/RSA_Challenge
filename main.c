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
				halfn,
				kiplusone;

	mpz_t			*k;

	double 		*timearray;

	
	char 		*ptr;
	long int 	N = strtol(argv[1],&ptr,10);
	int 		done = 0,
				found = 0;

	double		begin, end;
	double 		time_spent;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);


	begin = MPI_Wtime();

	mpz_init(n);
	mpz_init(q);
	mpz_init(p);
	mpz_init(pq);
	mpz_init(sqn);
	mpz_init(halfn);
	mpz_init(gap);
	mpz_init(kiplusone);

if(!my_rank) printf("N=%ld\n", N);
	mpz_set_d(n,N);
if(!my_rank) gmp_printf("n=%Zd\n",n);
	mpz_set_ui(halfn,2);
	mpz_sqrt(sqn,n);
if(!my_rank) gmp_printf("sqn:%Zd\n", sqn);
	mpz_tdiv_q(halfn,n,halfn);
	mpz_set(gap, sqn);
if(!my_rank) gmp_printf("gap:%Zd\n", gap);

	mpz_set_ui(p,1);
	mpz_nextprime(q, p);


	size_t mag = mpz_sizeinbase (gap, 10);
	mag=mag*procs;

if(!my_rank) printf("mag:%zu\n", mag);



	k = (mpz_t*)malloc(sizeof(mpz_t)*(mag+1));


	

	for (int i=0;i<=mag;i++)
	{
		mpz_init(k[i]);
		mpz_t temp;
		mpz_init(temp);
		mpz_tdiv_q_ui(temp,gap,mag);
		mpz_mul_ui(k[i],temp,i);
		//mpz_add(k[i],sqn,temp);
		//k[i]=sqn + (i* (int)(gap/mag));
	}
	//	k[mag] = n/2;
	mpz_set(k[mag],sqn);

/*if (!my_rank){
for (int i = 0; i<=mag; i++)
	gmp_printf("k[%d]=%Zd\n",i,k[i]);
}*/

int counter=0;
//printf("mag%d\n",mag);


	for (int i=my_rank; i<=mag; i=i+procs)
	{
		mpz_sub_ui(kiplusone,k[i+1],1);

		mpz_set(q,k[i]);
		mpz_sub_ui(q,q,1);
		mpz_nextprime(q,q);//		q=nextprime(k[i]-1);

//gmp_printf("i:%d p%Zd q%Zd\n", i,p,q);


		while (( mpz_cmp(q,kiplusone) <= 0/*q <= k[i+1]*/  )&&(!done)&&(!found))
		{//finding the prime numbers
/*if (i==0)
	gmp_printf("p%Zd q%Zd\n",p,q );
*/
			MPI_Allreduce(&found, &done, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			mpz_t temp;
			mpz_mod(temp,n,q);
			if (mpz_cmp_ui(temp,0)!=0) 
			{
				//gmp_printf("%Zd NOT DIVISIBLE BY %Zd\n", n,q);
				mpz_nextprime(q,q);
				continue;
			}
			mpz_nextprime(p,p);

			mpz_mul(pq,p,q);

			if (mpz_cmp(pq,n)==0 /*p*q == n*/)
			{
				found = 1;
				done=1;
			}
			else if ( mpz_cmp(pq,n)>0/*(p*q > n)*/)// || mpz_cmp(p,q)>0/*(p>q)*/)
			{
				mpz_set_ui(p,1);
				mpz_nextprime(q,q);
			}

			counter++;
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
	/*else
		printf("P%d: Numbers not found\n",my_rank);*/
	


printf("P%d finished with done=%d found=%d.\n", my_rank, done, found);

	time_spent = (double)(end - begin);

	if (my_rank!=0)
	{
		MPI_Send(&time_spent, sizeof(double),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
//printf("P%d sent their time %le\n", my_rank, time_spent);
	}
	else
	{
		timearray = (double*)malloc(sizeof(double)*procs);
		timearray[0] = time_spent;

		int j; 
		for (j = 1; j<procs;j++)
		{
			MPI_Recv(&timearray[j], sizeof(double),MPI_DOUBLE,j,0,MPI_COMM_WORLD,&status);
//			printf("Printing to file %dtime:%lf\n",j,timearray[j]);
		}


		writeTime(argv[1],timearray,procs);
	}

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

