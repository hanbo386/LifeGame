#include <stdio.h>
#include <mpi.h>
#include "MyMPI.h"
#include "NewMPI.c"

typedef int dtype;
#define MPI_TYPE MPI_INT
void DoLife(int id, int p, dtype **a, int m,int n,dtype *send_1,dtype *send_2,dtype *recv_1,dtype *recv_2); 		//执行一次迭代函数
void DataExchange(int id, int p, int m, int n, MPI_Status status,dtype **a,dtype *send_temp_1,dtype *send_temp_2,dtype *recv_temp_1,dtype *recv_temp_2);
//执行一次进程间数据通信

int main(int argc, char *argv[])
{
	int j,k;					//循环控制变量，j是迭代次数，每k次迭代打印一次结果
	int i;
	int m;						//读入矩阵的行数
	int n;  					//读入矩阵的列数
	int** a;					//每个子进程子矩阵指针
	int* storage;					//每个子进程数据域？？？？
	int id;						
	int p;		
	int *send_temp_1,*send_temp_2;			//用于临时保存发送的数据
	int *recv_temp_1,*recv_temp_2;			//用于临时保存接收的数据
	send_temp_1 = (dtype*)malloc(n * sizeof(dtype));//分配内存
	send_temp_2 = (dtype*)malloc(n * sizeof(dtype));
	recv_temp_1 = (dtype*)malloc(n * sizeof(dtype));
	recv_temp_2 = (dtype*)malloc(n * sizeof(dtype));
	j = *argv[2]-48;					//这里如果命令行缺少参数会出错 还缺少参数判断
	k = *argv[3]-48;					
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);					//获得进程总数
	MPI_Comm_rank(MPI_COMM_WORLD, &id);					//获得每个进程编号
	
	//memset();
	read_row_striped_matrix(argv[1],(void*) &a,(void*) &storage, MPI_TYPE, &m, &n, MPI_COMM_WORLD);
	//int row = BLOCK_SIZE(id,p,m);						//暂时还用不到，先留着吧，省得注销了
	for(i = 0; i < j; i++)					//这里不知道为什么
	{
		DataExchange(id,p,m,n,status,(dtype**) a,send_temp_1,send_temp_2,recv_temp_1,recv_temp_2);
		DoLife(id,p,(dtype**) a,m,n,send_temp_1,send_temp_2,recv_temp_1,recv_temp_2);
		if(i%k == 0) 
			print_row_striped_matrix((void**) a, MPI_TYPE, m, n, MPI_COMM_WORLD);		
	}
	MPI_Finalize();
	
}
/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void DataExchange(int id, int p, int m, int n,MPI_Status status,dtype **a,dtype *send_temp_1,dtype *send_temp_2,dtype *recv_temp_1,dtype *recv_temp_2)
{
	int i,j,k;
	int root;		//决定哪些变量用来广播
	int offset;		//进程内局部下标
	int row;		//保存块的大小，即每个进程分配的行数
	row = BLOCK_SIZE(id,p,m);//进程内行数
	for(k = 0; k < m; k++)
	{
		root = BLOCK_OWNER(k,p,m);
		if(root == id)							//遍历行索引，过滤不属于本进程的行
		{
			offset = k - BLOCK_LOW(id,p,m);				//offset是每个进程内部的指针偏移量
			if(id == 0 && offset == row - 1 && id != p-1)		//是0号进程，并且是当前进程最后一行，则向下一进程发送本行数据		
			{
				for(j = 0; j < n; j++)
				{
					//printf("%d ",a[offset][j]);
					send_temp_2[j] = a[offset][j];
				}
				//printf("\n");
				MPI_Send(send_temp_2,n,MPI_INT,id+1,1,MPI_COMM_WORLD);	  //进程间通信
				MPI_Recv(recv_temp_2,n,MPI_INT,id+1,1,MPI_COMM_WORLD,&status);	//注意：收发进程tag保持一致；传递变量个数正确指定					
			}				
			
			if((id == p - 1) && (offset == 0)&&(id != 0))		//是最后一个进程，并且是当前进程第一行，则向上一进程发送本行数据
			{								
				for(j = 0; j < n; j++)
					send_temp_1[j] = a[offset][j];
				MPI_Send(send_temp_1,n,MPI_INT,id-1,1,MPI_COMM_WORLD);
				MPI_Recv(recv_temp_1,n,MPI_INT,id-1,1,MPI_COMM_WORLD,&status);
				/*for(j = 0; j < n; j++) 
				{
					if(id == 1)
						printf("%d ",recv_temp_1[j]);
					//B[i][j] = recv_temp_1[j-1];
				}*/
				//printf("\n");
				
			}
			else if(id != 0 && id != p-1)				//其它中间进程收发本进程第一行以及最后一行数据
			{
				if(offset == 0)	//是首行数据
				{
					for(j = 0; j < n; j++)
						send_temp_1[j] = a[offset][j];
					MPI_Send(send_temp_1,n,MPI_INT,id-1,1,MPI_COMM_WORLD);
					MPI_Recv(recv_temp_1,n,MPI_INT,id-1,1,MPI_COMM_WORLD,&status);
				}
				if(offset == row-1)//是末行数据
				{
					for(j = 0; j < n; j++)
						send_temp_2[j] = a[offset][j];
					MPI_Send(send_temp_2,n,MPI_INT,id+1,1,MPI_COMM_WORLD);
					MPI_Recv(recv_temp_2,n,MPI_INT,id+1,1,MPI_COMM_WORLD,&status);
				}
			}
		}
	}
}

/*------------------------------------------------------------------一次迭代的实现－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－－*/

void DoLife(int id, int p,dtype **a, int m,int n,dtype *send_temp_1,dtype *send_temp_2,dtype *recv_temp_1,dtype *recv_temp_2)						
{
	int i,j,k;			//循环控制变量													
	int row;			//保存块的大小，即每个进程分配的行数										
	int **B, *Bstorage;		//动态分配内存	
	row = BLOCK_SIZE(id,p,m);	//进程内行数
	int row_new, col_new;		//中间变量，方便为中间矩阵分配内存
	row_new = row+2;
	col_new = n+2; 
	Bstorage = (dtype *)malloc(row_new * col_new * sizeof(dtype));	//动态分配一片内存，用作矩阵存档
	B = (dtype **)malloc(row_new * sizeof(dtype *));
	for(i = 0;i < row_new; i++)
	{
		B[i] = &Bstorage[i*col_new];				
	}
	for(j = 0; j < row_new*col_new; j++)			//初始化存档矩阵
	{
		Bstorage[j] = 0;
	}
	//printf("ddd ");			
	//MPI_Bcast(temp,m,MPI_TYPE,root,MPI_COMM_WORLD);
		
	//下面是具体的迭代实现
	int temp;							//用于保存邻居点加和之后的值
	i = 0;
	for(j = 1; j <= n; j++)  B[i][j] = recv_temp_1[j-1];		//把从其它进程收到的数据保存到存档矩阵第一行
	for(i = 1; i <= row; i++)					//把本进程子矩阵顺延保存至存档矩阵	{
		for(j = 1; j <=n; j++) B[i][j] = a[i-1][j-1];
	}	
	for(j = 1; j <=n; j++) B[i][j] = recv_temp_2[j-1];		//把从其它进程收到的数据保存到存档矩阵末行
	/*if(id == 1)
	for(i = 0; i < row+2;i++)
	{
		for(j = 0; j<n+2;j++)
			printf("%d ",B[i][j]);
		printf("\n");
	}*/
	for(i = 1; i <= row; i++)					//迭代的过程
	{
		for(j = 1; j <= n; j++)
		{
			temp = B[i][j] + B[i-1][j-1] + B[i-1][j] + B[i-1][j+1] + B[i][j-1] + B[i][j+1] + B[i+1][j-1] + B[i+1][j] + B[i+1][j+1];
			if(B[i][j] == 1)				//当前细胞为活细胞
			{
				if(temp == 3 || temp == 4)
					a[i-1][j-1] = 1;
				else
					a[i-1][j-1] = 0;
			}
			else						//当前细胞为死细胞
			{
				if(temp == 3)
					a[i-1][j-1] = 1;
				else
					a[i-1][j-1] = 0;
			}
		}
		//printf("\n");
	}	
}
