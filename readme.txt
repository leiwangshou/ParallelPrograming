---------------------------------------------ReadMe------------------------------------------
scattergather_mpi.c: all processors get same amount data via MPI_Scatter and return data to processor 0 via MPI_Gather. These two methods handle position for each part data.

distributed_mpi.c: In manager-worker mode, workers get data almost equally and encrypt it, which means processor 0 doesn’t encrypt it. However, processor 0 distributes data to others based on their id. Thus, when this manager collects, he knows position of each part data.

Fox.c: Uses Fox's algorithm to multiply two square matrices

cyclic.c: Cyclic Reduction to solve heat equation

recursive.c: Recursive Doubling to solve heat equation

prefix.cu: prefix sum