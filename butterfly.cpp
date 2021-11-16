#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <mpi.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>
using namespace std;

const int numall = 8; //要求和的数的数目

int ceil_log2(int n) {
    return ceil(log(n) / log(2));
}

int main(int argc, char **argv) {
    //生成numall个数字
    int data[numall] = {0,1,2,3,4,5,6,7};
    //for (int i = 0; i < numall; i++) {
    //  data[i] = rand() % 100;
    //}

    //MPI multi thread
    int comm_sz;
    int my_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //串行求和，用于验证是否正确
    if (my_rank == 0) {
        int test_sum = 0;
        for (int i = 0; i < numall; i++) {
            test_sum += data[i];
        }
        cout << "求得的和应为: " << test_sum << endl;
    }


    //分派数字
    int local_n = numall / comm_sz;  //每个线程被分配到的数字的数目
    if (numall %comm_sz > my_rank) local_n += 1;
    int local_init_index = my_rank * local_n;
    if (my_rank == numall % comm_sz) local_init_index = my_rank * (local_n + 1);
    if (numall%comm_sz < my_rank) local_init_index = (numall / comm_sz + 1)*(numall%comm_sz) + (my_rank - numall % comm_sz)*local_n;

    //对每个线程分配的数进行求和
    int local_sum = 0;
    for (int i = 0; i < local_n; i++) {
        local_sum += data[local_init_index + i];
    }

    //cout << my_rank << ": " << local_sum << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    //将每个线程求得的局部和进行求和，蝶形
    int k = ceil_log2(comm_sz);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < comm_sz; j += int(pow(2, i + 1))) {
            for (int k = 0; k < (1 << i); k++) {
                int px = j+k;
                int py = (1 << i) + px;
                if (px >= comm_sz || py >= comm_sz)break;
                if (px == my_rank) {
                    int recv_sum = 0;
                    MPI_Send(&local_sum, 1, MPI_INT, py, 0, MPI_COMM_WORLD);
                    MPI_Recv(&recv_sum, 1, MPI_INT, py, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    local_sum += recv_sum;
                }
                else if (py == my_rank) {
                    int recv_sum = 0;
                    MPI_Send(&local_sum, 1, MPI_INT, px, 0, MPI_COMM_WORLD);
                    MPI_Recv(&recv_sum, 1, MPI_INT, px, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    local_sum += recv_sum;
                }
            }
        }
        cout << my_rank << ": " << local_sum << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    //cout << my_rank << ": " << local_sum << endl;
    if (my_rank == 0) cout << local_sum << endl;
    MPI_Finalize();
    return 0;
}
