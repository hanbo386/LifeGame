/*
 *   MyMPI.c -- A library of matrix/vector
 *   input/output/redistribution functions
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 4 September 2002
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MyMPI.h"


/***************** MISCELLANEOUS FUNCTIONS *****************/

/*
 *   Given MPI_Datatype 't', function 'get_size' returns the
 *   size of a single datum of that data type.
 */

int get_size (MPI_Datatype t) {
   if (t == MPI_BYTE) return sizeof(char);
   if (t == MPI_DOUBLE) return sizeof(double);
   if (t == MPI_FLOAT) return sizeof(float);
   if (t == MPI_INT) return sizeof(int);
   printf ("Error: Unrecognized argument to 'get_size'\n");
   fflush (stdout);
   MPI_Abort (MPI_COMM_WORLD, TYPE_ERROR);
}


/*
 *   Function 'my_malloc' is called when a process wants
 *   to allocate some space from the heap. If the memory
 *   allocation fails, the process prints an error message
 *   and then aborts execution of the program.
 */

void *my_malloc (
   int id,     /* IN - Process rank */
   int bytes)  /* IN - Bytes to allocate */
{
   void *buffer;
   if ((buffer = malloc ((size_t) bytes)) == NULL) {
      printf ("Error: Malloc failed for process %d\n", id);
      fflush (stdout);
      MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
   }
   return buffer;
}


/*
 *   Function 'terminate' is called when the program should
 *   not continue execution, due to an error condition that
 *   all of the processes are aware of. Process 0 prints the
 *   error message passed as an argument to the function.
 *
 *   All processes must invoke this function together!
 */

void terminate (
   int   id,            /* IN - Process rank */
   char *error_message) /* IN - Message to print */
{
   if (!id) {
      printf ("Error: %s\n", error_message);
      fflush (stdout);
   }
   MPI_Finalize();
   exit (-1);
}

/******************** OUTPUT FUNCTIONS ********************/

/*
 *   Print elements of a doubly-subscripted array.
 */

void print_submatrix (
   void       **a,       /* OUT - Doubly-subscripted array */
   MPI_Datatype dtype,   /* OUT - Type of array elements */
   int          rows,    /* OUT - Matrix rows */
   int          cols)    /* OUT - Matrix cols */
{
   int i, j;

   for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
         if (dtype == MPI_DOUBLE)
            printf ("%6.3f ", ((double **)a)[i][j]);
         else {
            if (dtype == MPI_FLOAT)
               printf ("%6.3f ", ((float **)a)[i][j]);
            else if (dtype == MPI_INT)
               printf ("%6d ", ((int **)a)[i][j]);
         }
      }
      putchar ('\n');
   }
}


void read_row_striped_matrix (
   char        *s,        /* IN - File name */
   void      ***subs,     /* OUT - 2D submatrix indices */
   int       **storage,  /* OUT - Submatrix stored here */
   MPI_Datatype dtype,    /* IN - Matrix element type */
   int         *m,        /* OUT - Matrix rows */
   int         *n,        /* OUT - Matrix cols */
   MPI_Comm     comm)     /* IN - Communicator */
{
   int          datum_size;   /* Size of matrix element */
   int          i,j,k;
   int          id;           /* Process rank */
   FILE        *infileptr;    /* Input file pointer */
   int          local_rows;   /* Rows on this proc */
   void       **lptr;         /* Pointer into 'subs' */
   int          p;            /* Number of processes */
   void        *rptr;         /* Pointer into 'storage' */
   MPI_Status   status;       /* Result of receive */
   int          x;            /* Result of read */
   *m = 0;
   *n = 0;   //printf("%d",*m);
   MPI_Comm_size (comm, &p);
   MPI_Comm_rank (comm, &id);
   datum_size = get_size (dtype);

   /* Process p-1 opens file, reads size of matrix,
      and broadcasts matrix dimensions to other procs */

	if(id == (p-1))
	{
		infileptr = fopen(s,"r");
		if(infileptr == NULL) *m = 0;
		else
		{
			while(!feof(infileptr))			//读文件，读到文件末尾结束
			{
				int ch;
				ch = fgetc(infileptr);
				if(ch != 10)
					*n = *n + 1;		//记录文件的列数
				else	
					*m = *m + 1;		//记录文件的行数
			}
			*n = (*n-1)/(*m);			//会把最后的文件结束符读进来；//计算出实际的列数
		}
		rewind(infileptr);				//这里文件指针已到末尾，要移回文件开头
	}


  /* if (id == (p-1)) {
      infileptr = fopen (s, "r");		//第二个参数是文件的打开方式，r为只读方式打开，该文件必须存在 //打开文件后，指向该流的文件指针被返回
      if (infileptr == NULL) *m = 0;
      else {
         fread (m, sizeof(int), 1, infileptr);		//从数据流中读数据：用于接受数据地址；单个元素的大小；元素个数；提供数据的文件指针
         fread (n, sizeof(int), 1, infileptr);		//获取矩阵行数和列数
      }      						//二进制流读文件，不懂，因此替换为以字符读文件
   }*/
//	printf("%d  ",*m);
//	printf("%d  ",*n);
   MPI_Bcast (m, 1, MPI_INT, p-1, comm);
   if (!(m)) MPI_Abort (MPI_COMM_WORLD, OPEN_FILE_ERROR);
   MPI_Bcast (n, 1, MPI_INT, p-1, comm);
   local_rows = BLOCK_SIZE(id,p,*m);;
   /* Dynamically allocate matrix. Allow double subscripting
      through 'a'. */
/*--------------问题在这儿-------------------*/		//这里不清楚，先跳过了
   *storage = (int *) malloc(local_rows * *n * datum_size);
  //ch = fgetc(infileptr);storage = (void *) my_malloc (id,		//原来的my_malloc函数会报错，原因不明，待求证，此处改为malloc
  //     local_rows * *n * datum_size);

   //subs = (void **) my_malloc (id, local_rows * PTR_SIZE);
   *subs = (void **) malloc(local_rows * PTR_SIZE);
   lptr = (void *) &(*subs[0]);
   rptr = (void *) *storage;
   for (i = 0; i < local_rows; i++) {
      *(lptr++)= (void *) rptr;
      rptr += *n * datum_size;
   }		
/*--------------------------------------------------------------------	//谁能告诉我这段代码怎么理解？？？？？？先不管了
   /* Process p-1 reads blocks of rows from file and
      sends each block to the correct destination process.
      The last block it keeps. */
//	printf("%d  ",row);
/*--------------------------------------这里也有问题---------------------------------------------------------------*/
   if (id == (p-1)) {						//最后一个进程负责按行读取矩阵并将读出的数据发送到每个进程
      for (i = 0; i <= p-1; i++) {				//因为最后要给自己留数据，所以用=号
        /* x = fread (*storage, datum_size,
            BLOCK_SIZE(i,p,*m) * *n, infileptr);*/		//fread函数什么的都替换掉
	int *temp = *storage;
	for(j = 0; j < BLOCK_SIZE(i,p,*m) * *n; j++)		//这里要格外注意，控制循环的应该是每一个进程各自的大小，不能用local_rows
	{
		x = fgetc(infileptr);
		if(x != 10)					//换行符的值是10
			temp[j] = x-48;				//记得把ASCII码转换成整数
		else
			j--;					//遇到换行符之后记得把变量回退1
	}
	x = fgetc(infileptr);					//最后一个换行符记得吃掉
	if(i != p-1)
        	MPI_Send (*storage, BLOCK_SIZE(i,p,*m) * *n, dtype, i, DATA_MSG, comm);		//把属于各个进程的数据发送出去，自己除外
      }
	
    // x = fread (*storage, datum_size, local_rows * *n,		//问题在这句话 //最后剩下的数据，由最后一个进程负责（原来这里不晓得为什么有问题）
     //    infileptr);
      fclose (infileptr);
   } else
      MPI_Recv (*storage, local_rows * *n, dtype, p-1,		//其它进程负责接收自己的数据就好了
         DATA_MSG, comm, &status);
   
}
/*-------------------------------------------还没看这个函数怎么实现的，但是真的很好用---------------------------------------------------------*/
/*
 *   Print a matrix that is distributed in row-striped
 *   fashion among the processes in a communicator.
 */

void print_row_striped_matrix (
   void **a,            /* IN - 2D array */
   MPI_Datatype dtype,  /* IN - Matrix element type */
   int m,               /* IN - Matrix rows */
   int n,               /* IN - Matrix cols */
   MPI_Comm comm)       /* IN - Communicator */
{
   MPI_Status  status;          /* Result of receive */
   void       *bstorage;        /* Elements received from
                                   another process */
   void      **b;               /* 2D array indexing into
                                   'bstorage' */
   int         datum_size;      /* Bytes per element */
   int         i;
   int         id;              /* Process rank */
   int         local_rows;      /* This proc's rows */
   int         max_block_size;  /* Most matrix rows held by
                                   any process */
   int         prompt;          /* Dummy variable */
   int         p;               /* Number of processes */

   MPI_Comm_rank (comm, &id);
   MPI_Comm_size (comm, &p);
   local_rows = BLOCK_SIZE(id,p,m);
   if (!id) {
      print_submatrix (a, dtype, local_rows, n);
      if (p > 1) {
         datum_size = get_size (dtype);
         max_block_size = BLOCK_SIZE(p-1,p,m);
	bstorage = malloc(max_block_size * n * datum_size);
	b = (void **)malloc(max_block_size * datum_size);
         /*bstorage = my_malloc (id,
            max_block_size * n * datum_size);			//这里记得替换
         b = (void **) my_malloc (id,
            max_block_size * datum_size);*/
         b[0] = bstorage;
         for (i = 1; i < max_block_size; i++) {
            b[i] = b[i-1] + n * datum_size;
         }
         for (i = 1; i < p; i++) {
            MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG,
               MPI_COMM_WORLD);
            MPI_Recv (bstorage, BLOCK_SIZE(i,p,m)*n, dtype,
               i, RESPONSE_MSG, MPI_COMM_WORLD, &status);
            print_submatrix (b, dtype, BLOCK_SIZE(i,p,m), n);
         }
         free (b);
         free (bstorage);
      }
      putchar ('\n');
   } else {
      MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG,
         MPI_COMM_WORLD, &status);
      MPI_Send (*a, local_rows * n, dtype, 0, RESPONSE_MSG,
         MPI_COMM_WORLD);
   }
}


