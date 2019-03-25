
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void encryptText(char* ori, char* enTxt, int enkey, int size);
void decryptText(char* enTxt, char* deTxt, int enkey, int size);

int main(int argc, char* argv[]) {

	int myId, numProc;
	int ttSize;
	int chunkSize;
	int enkey = 3;
	double start_time, end_time, total_time;

	char* origTxt = NULL;
	char* enTxt = NULL;
	char* deTxt = NULL;
	char* chunkTxt = NULL;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	start_time = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);

	if (myId == 0) {
		printf("\n-----------------------------------------\n");
		printf("size of test file is 1kb, the number of processors is %d\n", numProc);
		printf("\n-----------------------------------------\n");
		FILE* fp;
		long fsize;
		fp = fopen("./test_1kb.txt", "r");
		if (fp == NULL) {
			printf("File not found! \n");
		}
		else {
			fseek(fp, 0, SEEK_END);
			fsize = ftell(fp);
			fseek(fp, 0, SEEK_SET);
			ttSize = fsize + 1;
			origTxt = (char*)malloc(ttSize * sizeof(char));
			memset(origTxt, 0, ttSize);
			fread(origTxt, fsize, 1, fp);
			origTxt[fsize] = '\0';
		}
		fclose(fp);

		enTxt = (char*)malloc(ttSize * sizeof(char));
		deTxt = (char*)malloc(ttSize * sizeof(char));
		memset(enTxt, 0, ttSize);
		memset(deTxt, 0, ttSize);

		if (numProc == 1) {
			//encrypt text
			encryptText(origTxt, enTxt, enkey, ttSize);
			//decrypt text
			decryptText(enTxt, deTxt, enkey, ttSize);
/*			printf("\n--------------------------------------\n");
			printf("encrypted text \n");
			printf("--------------------------------------\n");
			printf("%s", enTxt);

			printf("\n--------------------------------------\n");
			printf("decrypted text \n");
			printf("--------------------------------------\n");
			printf("%s\n", deTxt);*/

			printf("\n--------------------------------------\n");

			if (strcmp(origTxt, deTxt) == 0) {
				printf("Encryption and decryption are both correct.\n");
			}
			else {
				printf("There are something wrong in encryption and decrytion!\n");
			}

			printf("\n--------------------------------------\n");
			printf("\n--------------------------------------\n");
			end_time = MPI_Wtime();
			total_time = end_time - start_time;
			printf("total time is %f s", total_time);
			printf("\n---------------------------------------\n");
		}
		else {
			chunkSize = ((ttSize - 1) / numProc) + 1;
			chunkTxt = (char*)malloc(chunkSize * sizeof(char));
			memset(chunkTxt, 0, chunkSize);
			int i;
			for (i = 1; i < numProc; ++i)
				MPI_Send(&ttSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		

	}
	else {
		MPI_Recv(&ttSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		chunkSize = ((ttSize - 1) / numProc) + 1;
		chunkTxt = (char*)malloc(chunkSize * sizeof(char));
		memset(chunkTxt, 0, chunkSize);
		enTxt = (char*)malloc(ttSize * sizeof(char));
		memset(enTxt, 0, ttSize);

	}

	if (numProc > 1) {
		if (myId == 0) {
			MPI_Scatter(origTxt, (chunkSize - 1), MPI_CHAR, chunkTxt,
				(chunkSize - 1), MPI_CHAR, 0, MPI_COMM_WORLD);
			chunkTxt[chunkSize - 1] = '\0';
		}
		else {
			MPI_Scatter(origTxt, (chunkSize - 1), MPI_CHAR, chunkTxt,
				(chunkSize - 1), MPI_CHAR, 0, MPI_COMM_WORLD);
			chunkTxt[chunkSize - 1] = '\0';
		}

		//encrypt text
		encryptText(chunkTxt, chunkTxt, enkey, chunkSize); 

		if (myId == 0) {
			MPI_Gather(chunkTxt, chunkSize - 1, MPI_CHAR, enTxt, chunkSize - 1, MPI_CHAR,
				0, MPI_COMM_WORLD);
		}
		else {
			MPI_Gather(chunkTxt, chunkSize - 1, MPI_CHAR, enTxt, chunkSize - 1, MPI_CHAR,
				0, MPI_COMM_WORLD);
		}


		if (myId == 0) {
			enTxt[ttSize - 1] = '\0';
/*			printf("\n--------------------------------------\n");
			printf("encrypted text \n");
			printf("--------------------------------------\n");
			printf("%s\n", enTxt);*/

			//decrypt text
			decryptText(enTxt, deTxt, enkey, ttSize);
			/*
			printf("\n--------------------------------------\n");
			printf("decrypted text \n");
			printf("----------------------------------------\n");
			printf("%s\n", deTxt);                   */
			printf("\n--------------------------------------\n");
			if (strcmp(origTxt, deTxt) == 0) {
				printf("Encryption and decryption are both correct.\n");
			}
			else {
				printf("There are something wrong in encryption and decrytion!\n");
			}

			printf("\n--------------------------------------\n");
			printf("\n--------------------------------------\n");
			end_time = MPI_Wtime();
			total_time = end_time - start_time;
			printf("total time is %f s\n", total_time);
			printf("\n--------------------------------------\n");
		}  

		if (myId == 0) {
			memset(origTxt, 0, ttSize);
			memset(enTxt, 0, ttSize);
			memset(deTxt, 0, ttSize);
			memset(chunkTxt, 0, chunkSize);
			free(chunkTxt);
			free(enTxt);
			free(deTxt);
			free(origTxt);
		}
		else {
			memset(enTxt, 0, ttSize);
			memset(chunkTxt, 0, chunkSize);
			free(chunkTxt);
			free(enTxt);
		}
	}  
	else {
		memset(origTxt, 0, ttSize);
		memset(enTxt, 0, ttSize);
		memset(deTxt, 0, ttSize);
		free(enTxt);
		free(deTxt);
		free(origTxt);
	}

	MPI_Finalize();

	return 0;
}

void encryptText(char* ori, char* enTxt, int key, int size) {
	int i;
	for (i = 0; i < (size - 1); ++i) {
		enTxt[i] = (char)(((int)(ori[i]) + key) % 256);
	}
	enTxt[size - 1] = '\0';

}

void decryptText(char* enTxt, char* deTxt, int key, int size) {
	int i;
	for (i = 0; i < size - 1; ++i) {
		deTxt[i] = (char)(((int)(enTxt[i]) - key) % 256);
	}
	deTxt[size - 1] = '\0';
}
