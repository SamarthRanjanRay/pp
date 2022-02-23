#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

void merge(int a[], int x1, int y1, int x2, int y2) {
    int p = x1, q = x2, temp[100], i = 0, j = 0, k = 0;
    while (p <= y1 && q <= y2) {
        if (a[p] < a[q]) {
            temp[k++] = a[p++];
        }
        else {
            temp[k++] = a[q++];
        }
    }

    while (p <= y1) {
        temp[k++] = a[p++];
    }
    while (q <= y2) {
        temp[k++] = a[q++];
    }
    for(i = x1, j = 0; i <= y2; i++, j++) {
        a[i] = temp[j];
    }
}

void mergesort(int a[], int low, int high) {
    if (low < high) {
        int mid = (low+high)/2;
        #pragma omp parallel sections {
            #pragma omp section {
                mergesort(a, low, mid);
            }
            #pragma omp section {
                mergesort(a, mid+1, high);
            }
        }
        merge(a,low,mid,mid+1,high);
    }
}

void main() {
    int n = 100;
    scanf("%d",&n);
    int arr[n];
    for(int i = 0; i < n; i++) {
        arr[i] = rand()%100;
    }
    
    printf("Unsorted :\n");
    for(int i = 0; i < n; i++)
        printf("%d ", a[i]);
    print("\n\n\n");
        
    mergesort(arr,0,n-1);
    
    printf("Sorted :\n");
    for(int i = 0; i < n; i++)
        printf("%d ", a[i]);
    print("\n\n\n");
}

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

void main() {
    int rows, cols, size;
    scanf("%d%d%d",&rows, &cols, &size);
    if (rows <= 0 || cols <= 0 || size <= 0 || cols != size) {
        printf("Invalid sizes");
        return;
    }

    int matrix[rows][cols], vector[size], result[rows];

    printf("Matrix is :\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = i+j;
            printf("%d ", (i+j));
        }
        printf("\n");
        result[i] = 0;
    }

    printf("Vector is :\n");
    for (int i = 0; i < size; i++) {
        vector[i] = i;
        printf("%d ", i);
    }

    #pragma omp parallel for private(j) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i] += matrix[i][j]*vector[j];
            }
        }
    }

    printf("\nResult is :\n");
    for (int i = 0; i < rows; i++) {
        printf("%d ", result[i]);
    }
        
}

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<math.>

int primes[1000];
double sines[1000];

void findPrime(int max) {
    int i = 2, j = 2, flag = 1, k=0;
    while (k < max) {
        for(i = 2; i*i < j; i++) {
            flag = 1;
            if (j%i == 0) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            primes[k++] = j;
        }
        j++;
    }
}

void findSines(int max) {
    for(int i = 0; i < max; i++) {
        sines[i] = (double) sin(i*3.14/180);
    }
}

void main() {
    int max = 0;
    scanf("%d",&max);
    #pragma omp parallel sections {
        #pragma omp section {
            findPrime(max);
            for(int i = 0; i < max; i++) {
                printf("Prime no at index %d is %lf", i, primes[i]);
            }
             
        }
        #pragma omp section {
            findSines(max);
            for(int i = 0; i < max; i++) {
                printf("Sine no at index %d is %lf", i, sines[i]);
            }
        }
    }
    
}

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

void main() {
    int a, i, fact1 = 1, fact2 = 1;
    scanf("%d",&a);

    //wihtout firstprivate
    #pragma omp parallel for firstprivate(a) {
        for (i = 2 ; i <= a; i++) {
            fact1 *= i;
        }
    }
    printf("\nFactorial without firstprivate %d\n",fact1);
    
    //wihtout firstprivate
    #pragma omp parallel for firstprivate(a,fact2) {
        for (i = 2 ; i <= a; i++) {
            fact2 *= i;
        }
    }
    printf("\nFactorial with firstprivate %d\n",fact2)
}


#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<math.h>
#include<time.h>

int clusterCenter[4][2]={{25,25},{25,75},{75,75},{75,25}}, clusterCount[4], points[100][2];

void fill(int max) {
    #pragma omp parallel {
        for (int i = 0; i < max; i++) {
            points[i][0] = rand()%100;
            points[i][1] = rand()%100;
        }
    }
    for (int i = 0; i < 4; i++) {
        clusterCount[i] = 0;
    }
}

double dis(int x1, int y1, int x2, int y2) {
    int dx = x1-x2;
    int dy = y1-y2;
    return (double)sqrt(dx*dx+dy*dy)
}

void classify(int max) {
    int i = 0, j = 0, clusterIndex = 0;
    double curDist = 0, minDist = 9999;
    #pragma omp parallel for private(i,j,clusterIndex,curDist,minDist) {
        for(i = 0; i < max; i++) {
            minDist = 9999;
            clusterIndex = -1;
            for (j = 0; j < 4; j++) {
                cur_dist = dis(points[i][0],points[i][1],clusterCenter[j][0],clusterCenter[j][1]);
                if (curDist < minDist) {
                    minDist = curDist;
                    clusterIndex = j;
                }
            }
            clusterCount[clusterIndex]++;
        }
    } 
}

void main() {
    int max;
    scanf("%d",&max);
    fill(max);
    double t1 = omp_get_wtime();
    classify(max);
    double t2 = omp_get_wtime();
    printf("\nCluster Counts\n");
    for(int i = 0; i < max; i++)
        printf("%d index has %d points\n", i, clusterCount[i]);
}

==============================

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>

#define BUFFER 32

void main(int argc, char* argv[]) {
    int rank, size, temp = 0, root = 3;
    char* msg[BUFFER];
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status;

    if (rank == 3) {
        strcpy(msg, "Hello");
        for (temp = 0, temp < size; temp++) {
            if (temp != 3) {
                MPI_Send(msg, BUFFER, MPI_CHAR, temp, 0, MPI_COMM_WORLD);
            }
        }
    }
    else {
        MPI_Recv(msg, BUFFER, MPI_CHAR, root, 0, MPI_COMM_WORLD, &status);
        printf("\nmsg %s at processor %d from processor %d\n", msg, rank, root);
    }
    MPI_Finalize();
}

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

void main(int argc, char* argv[]) {
    int i, rank, size;
    double a[100], b[100];
    double mysum, allsum;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("\nstarting with %d proc\n",size);
    }

    for(i = 0; i < 100; i++) {
        a[i] = 1.0;
        b[i] = a[i];
    }

    mysum = 0.0;
    
    for(int i = 0; i < 100; i++) {
        mysum += a[i]*b[i];
    }
    printf("\nprocessor %d's sum is %d\n",rank,mysum);
    MPI_Reduce(&mysum,&allsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if (rank == 0) {
        printf("\nAllsum is %d\n", allsum);
    }
    MPI_Finalize();
}

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>

#define SEED 12345678
#define N 100000

void main(int argc, char*argv[]) {
    int rank, size;
    double count = 0, final = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    srand(SEED+rank);
    for(int i = 0; i < N; i++) {
        double x = (double)rand()/RAND_MAX;
        double y = (double)rand()/RAND_MAX;
        z = x*x+y*y;
        if (z < 1) {
            count++;
        }
    }
    MPI_Reduce(&count,&final,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if (rank == 0) {
        printf("\nPI value is %ld\n", final/N);
    }
    MPI_Finalize();
}

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

void main(int argc, char* argv[]) {
    double evenp, eventotal, oddp, oddtotal;
    int rank, size;
    MPI_Comm even_comm, odd_comm;
    MPI_Group world_group, even_group, odd_group;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size)
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);

    int even_size = size/2;
    int even_p = {0,2};

    MPI_Group_incl(world_group,even_size,even_p,&even_group);
    MPI_Comm_create(MPI_COMM_WORLD,even_group,&even_comm);

    int odd_size = size/2;
    int odd_p = {1,3};

    MPI_Group_incl(world_group,even_size,even_p,&even_group);
    MPI_Comm_create(MPI_COMM_WORLD,even_group,&even_comm);

    if (rank == even_p[0] || rank == even_p[1]) {
        evenp = rank;
        MPI_Reduce(&evenp,&eventotal,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }
    else {
        oddp = rank;
        MPI_Reduce(&oddp,&oddtotal,1,MPI_DOUBLE,MPI_SUM,1,MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("\nEven total is %d\n",eventotal);
    }
    else if (rank == 1) {
        printf("\nOdd total is %d\n",oddtotal);
    }
    MPI_Finalize();
}

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

#define BUFFER sizeof(int)*4

int a[2][2] = {{1,2},{3,4}};
int b[2][2] = {{1,0}, {1,0}};
int c[2][2] = {{0,0},{0,0}};

void main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if (rank == 0) {
        printf("Matrix mul at proc %d", rank);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                c[i][j] = 0;
                for (int k = 0; k < 2; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        MPI_Send(&c,BUFFER,MPI_INT,1,0,MPI_COMM_WORLD);
    }
    else if (rank == 1) {
        MPI_Recv(&c,BUFFER,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        printf("\nResult is : \n");
        for(int i = 0; i < 2; i++) {
            for(int j = 0; j < 2; j++) {
                printf("%d ", c[i][j]);
            }
            printf("\n");
        }
    }
    MPI_Finalize();
}
