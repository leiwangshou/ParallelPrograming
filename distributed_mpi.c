#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void encryptText(char* ori, char* enTxt, int enkey, int size);
void decryptText(char* enTxt, char* deTxt, int enkey, int size);

int main(int argc, char* argv[]) {
	
	int myId, numProc;
	int enkey = 3;
	long ttSize = 0;
	long chunkSize = 0;
	long lastSize = 0;
	char* origTxt = NULL; //orignal text
	char* enTxt = NULL;   //encrypted text
	char* deTxt = NULL;   //decrypted text
	char* chunkTxt = NULL; //chunk text
	char* enChunk = NULL; //chunk encrypted text
	char* lastTxt = NULL;
	double start_time, end_time, total_time;

	MPI_Init(&argc, &argv);
	start_time = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);
	MPI_Status status;
	
	if (myId == 0) {
		printf("\n-----------------------------------------------\n");
		printf("size of test file is 1kb, and number of process is %d\n", numProc);
		printf("-------------------------------------------------\n");
		FILE* fp;
		long fsize;
		fp = fopen("./test_1kb.txt", "r");
		if(fp == NULL) {
			printf("File not found! \n");
		}
		else {
			fseek(fp, 0, SEEK_END);
			fsize = ftell(fp);
			fseek(fp, 0, SEEK_SET);
			ttSize = fsize + 1;
			//Read text to origTxt
			origTxt = (char*) malloc(ttSize * sizeof(char));
			memset(origTxt, 0, ttSize);
			fread(origTxt, fsize, 1, fp);
			origTxt[fsize] = '\0';
			enTxt = (char*)malloc(ttSize * sizeof(char));
			deTxt = (char*)malloc(ttSize * sizeof(char));
			memset(enTxt, 0, ttSize);
			memset(deTxt, 0, ttSize);
		}	
		fclose(fp);
	/*	printf("\n--------------------------------------\n");
		printf("original text \n");
		printf("--------------------------------------\n");
		printf("%s\n", origTxt);
		*/
		if (numProc == 1) {
			//encrypt text
			encryptText(origTxt, enTxt, enkey, ttSize);
			//decrypt text
			decryptText(enTxt, deTxt, enkey, ttSize);
		/*	printf("\n--------------------------------------\n");
			printf("encrypted text \n");
			printf("--------------------------------------\n");
			printf("%s", enTxt);

			printf("\n--------------------------------------\n");
			printf("decrypted text \n");
			printf("--------------------------------------\n");
			printf("%s\n", deTxt);

			printf("\n---------------------------------------\n");
                        */
			if (strcmp(origTxt, deTxt) == 0) {
				printf("Encryption and decryption are both correct.\n");
			}
			else {
				printf("There are something wrong in encryption and decrytion!\n");
			}
			printf("\n---------------------------------------------\n");
			end_time = MPI_Wtime();
			total_time = end_time - start_time;
			printf("\n----------------------------------------------\n");
			printf("total time is %f s\n", total_time);
		}	
		else {
			chunkSize = (fsize / (numProc - 1)) + 1;
			int i;
			for (i = 0; i < (numProc - 2); ++i) {
				chunkTxt = (char*)malloc(chunkSize * sizeof(char));
				memset(chunkTxt, 0, chunkSize);
				strncpy(chunkTxt, origTxt + i * (chunkSize - 1), chunkSize - 1);
				chunkTxt[chunkSize - 1] = '\0';
				MPI_Send(&chunkSize, 1, MPI_LONG, (i + 1), 0, MPI_COMM_WORLD);
				MPI_Send(chunkTxt, chunkSize, MPI_CHAR, (i + 1), 1, MPI_COMM_WORLD);
			}
			lastSize = fsize - (chunkSize - 1) * (numProc - 2) + 1;
			lastTxt = (char*)malloc(lastSize * sizeof(char));
			memset(lastTxt, 0, lastSize);
			strncpy(lastTxt, origTxt + (numProc - 2) * (chunkSize - 1), lastSize - 1);
			lastTxt[lastSize - 1] = '\0';
			MPI_Send(&lastSize, 1, MPI_LONG, (numProc - 1), 0, MPI_COMM_WORLD);
			MPI_Send(lastTxt, lastSize, MPI_CHAR, (numProc - 1), 1, MPI_COMM_WORLD);
		}
		
	}
	else {
		if (myId != numProc - 1) {
			MPI_Recv(&chunkSize, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
			chunkTxt = (char*) malloc(chunkSize * sizeof(char));
			memset(chunkTxt, 0, chunkSize);

			MPI_Recv(chunkTxt, chunkSize, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status);
			enChunk = (char*) malloc(chunkSize * sizeof(char));
			memset(enChunk, 0, chunkSize);
			encryptText(chunkTxt, enChunk, enkey, chunkSize);
			MPI_Send(enChunk, chunkSize, MPI_CHAR, 0, 2, MPI_COMM_WORLD);
		}
		else {
			MPI_Recv(&lastSize, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
			lastTxt = (char*)malloc(lastSize * sizeof(char));
			memset(lastTxt, 0, lastSize);
			MPI_Recv(lastTxt, lastSize, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status);
			enChunk = (char*)malloc(lastSize * sizeof(char));
			memset(enChunk, 0, lastSize);
			encryptText(lastTxt, enChunk, enkey, lastSize);
			MPI_Send(enChunk, lastSize, MPI_CHAR, 0, 2, MPI_COMM_WORLD);
		}
	}
	
	if ((myId == 0) && (numProc > 1)) {
		enChunk = (char*)malloc(chunkSize * sizeof(char));
		memset(enChunk, 0, chunkSize);
		int idx = 0;
		int j;
		for (j = 1; j < numProc - 1; ++j) {
			MPI_Recv(enChunk, chunkSize, MPI_CHAR, j, 2, MPI_COMM_WORLD, &status);
			idx = (j - 1) * (chunkSize - 1);
			int t;
			for (t = 0; t < chunkSize - 1; ++t) {
				enTxt[idx + t] = enChunk[t];
			}
		}
		MPI_Recv(enChunk, lastSize, MPI_CHAR, numProc - 1, 2, MPI_COMM_WORLD, &status);
		idx = (numProc - 2) * (chunkSize - 1);
		int m;
		for (m = 0; m < lastSize - 1; ++m) {
			enTxt[idx + m] = enChunk[m];
		}
		enTxt[idx + lastSize - 1] = '\0';
		
/*		printf("\n--------------------------------------\n");
		printf("encrypted text \n");
		printf("--------------------------------------\n");
		printf("%s", enTxt);*/
		//decrypt text
		decryptText(enTxt, deTxt, enkey, ttSize);
/*		printf("\n--------------------------------------\n");
		printf("decrypted text \n");
		printf("--------------------------------------\n");
		printf("%s\n", deTxt);   */
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
		printf("total_time is %f \n", total_time);

	}
	
	if (numProc > 1){
		if (myId == 0) {
			free(deTxt);
			free(enTxt);
			if(chunkTxt != NULL) {
				free(chunkTxt);
			}
			free(lastTxt);
			free(origTxt);
			memset(enChunk, '\0', strlen(enChunk));
		}
		else if(myId == numProc - 1) { 
			free(enChunk);
			free(lastTxt);
		}
		else {
			free(enChunk);
			free(chunkTxt);
		}
	}
	else {
		free(deTxt);
		free(enTxt);
		free(origTxt);
	}
	
	MPI_Finalize();
	
	
	return 0;
}

void encryptText(char* ori, char* enTxt, int key, int size) {
	int i;	
	for (i = 0; i < size - 1; ++i) {
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
