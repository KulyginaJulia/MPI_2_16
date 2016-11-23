#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include <random>


double* CreateAndFillMatrix(int n, int m) {
	int i,j;
	double* matrix = (double*)malloc(n * m * sizeof(double));
	for (i = 0; i < n*m; i++) {
		matrix[i] = rand() / 10000.0;
	}
	return matrix;
}

void PrintMatrix(double* matrix, int n, int m) {
	int i, j;
	if ((n <= 10) && (m <= 10)) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++)
				printf_s("%f ", matrix[i*m + j]);
			printf_s("\n");
		}
	}
}
void PrintVector(double* vector, int m) {
	int i;
	if (m<=10)
		for (i = 0; i < m; i++)
			printf_s("%f ", vector[i]);
}

double* CreateAndFillVector(int m){
	int i;
	double *vector = (double*)malloc(m * sizeof(double));
	for (i = 0; i < m; i++)
		vector[i] = rand() / 10000.0;
	return vector;
}
void DeleteMatrix(double* matrix) {
	free(matrix);
}
double* Consistent(double* matrix, double* vector, double* result, int m, int n) {
	int i, j;
	double *tmp = (double*)malloc(n*m*sizeof(double));
	for (j = 0; j < m; j++) {
		for (i = 0; i < n; i++)
			tmp[i*m + j] = matrix[i*m + j] * vector[j];
		}
	printf_s("\n");
	printf_s("tmp: \n");
	PrintMatrix(tmp,n,m);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++)
			result[i] += tmp[i*m + j];
	}
	printf_s("\n");
	printf_s("Result: \n");
	PrintVector(result, n);
	free(tmp);
	return result;
}

bool CheckResult(double* A, double* B, int m) {
	int i;
	for (i = 0; i < m; i++){
		if (A[i] != B[i]){
			return false;
		}
	}
	return true;
}

using namespace std;
int main(int argc, char* argv[]) {
	double time1, time2, delta_time_1, delta_time_2;

	int i, n, j, m;
	double *vector = NULL;
	int ProcNum, ProcRank;
	int dataSize, deltaSize;
	double *matrix = NULL;
	double *result1 = NULL;
	double *result2 = NULL;
	
	double *mbuff = NULL;
	double *vbuff = NULL;
	double *resultbuff = NULL;
	
	if (argc >= 3) {
		n = atoi(argv[1]);
		m = atoi(argv[2]);
	}
	else {
		printf_s("Error with argv");
	}
	MPI_Status status;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	dataSize = (int)(m / (ProcNum - 1));
	deltaSize = m % (ProcNum - 1);

	result1 = (double *)malloc(n * sizeof(double));
	result2 = (double *)malloc(n * sizeof(double));
	for (i = 0; i < n; i++) {
		result1[i] = 0;
		result2[i] = 0;
	}
	mbuff = (double*)malloc(n*dataSize* sizeof(double));
	vbuff = (double*)malloc(m * sizeof(double));
	resultbuff = (double*)malloc(n * sizeof(double));
	
	if (ProcRank == 0) {
		matrix = CreateAndFillMatrix(n, m);
		vector = CreateAndFillVector(m);
		if ((n <= 10) && (m <= 10)) {
			printf_s("Matrix: \n");
			PrintMatrix(matrix, n, m);
			printf("\n");
			printf_s("Vector: \n");
			PrintVector(vector, m);
			printf_s("\n");
		}
	}

	if (ProcRank == 0) {
		time1 = MPI_Wtime();
		result1 = Consistent(matrix, vector, result1, m, n);
		time2 = MPI_Wtime();
		printf_s("Consistent: ");
		PrintVector(result1,n);
		delta_time_1 = time2 - time1;
		printf_s("Time: %f \n", delta_time_1);
	}

	MPI_Datatype ColomnMatrix;
	int blocklength = 1;	//число элементов базового типа в каждом блоке
	int stride = m;			//шаг между началами соседних блоков
	int count = m;			//число блоков
	MPI_Datatype oldtype = MPI_DOUBLE;
	MPI_Type_vector(count, blocklength, stride, oldtype, &ColomnMatrix);
	MPI_Type_commit(&ColomnMatrix);

	if (ProcRank == 0) {
		//printf_s("RankProc = %d", ProcRank);
		printf_s("\n");
			time1 = MPI_Wtime();

			double *temp_vector = vector;
			for (i = 1; i < ProcNum; i++) {
				MPI_Send(temp_vector, dataSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				temp_vector = temp_vector + dataSize;
			}
			for(i = 0; i < m; i++)
				for (j = 1; j < ProcNum; j++)
					MPI_Send(matrix + i, 1, ColomnMatrix, j, 0, MPI_COMM_WORLD);

	
			//Приняли обратно результат
			double* temp_result2 = result2;
			for (i = 1; i < ProcNum; i++) {
					MPI_Recv(temp_result2, dataSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
					temp_result2 = temp_result2 + dataSize;
			}

			//Обработка конца матрицы
			if (deltaSize != 0) {
				double *tmp = (double*)malloc(n*dataSize * sizeof(double));
				for (j = ProcNum; j < m; j++)
					for (i = n - deltaSize; i < n; i++) {
						tmp[i*m+j] = matrix[i*m+j] * vector[j];
					}
				for (i = n - deltaSize; i < n; i++)
					for (j = ProcNum; j < m; j++)
						temp_result2[i] += tmp[i*m + j];
				free(tmp);
				}
			time2 = MPI_Wtime();
			delta_time_2 = time2 - time1;
	}
	if (ProcRank != 0) {
		//printf_s("RankProc = %d", ProcRank);
		MPI_Recv(vbuff, dataSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		double *tmp = (double*)malloc(n*dataSize * sizeof(double));
		for (i = 0; i < dataSize; i++) {
			MPI_Recv(mbuff, dataSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			for (j = 0; j < n; j++)
				tmp[j*dataSize +i] = mbuff[j] * vbuff[i];
		}
		for (i = 0; i < n; i++)
			for (j = 0; j < dataSize; j++)
				resultbuff[i] += tmp[i*dataSize + j];
		free(tmp);
		MPI_Send(resultbuff, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	if (ProcRank == 0){	
		printf_s("MPI: ");
		delta_time_2 = time2 - time1;
		printf_s("Time: %f \n", delta_time_2);
		//printf_s("Acceleration(parallel):  %f\n", (delta_time_1 / delta_time_2));
		
		bool flag = true;
		flag = CheckResult(result1, result2, n);
		if (flag == true)
			printf_s("Consistent = MPI");
		else printf_s("Consistent != MPI");
		}
	MPI_Type_free(&ColomnMatrix);

	DeleteMatrix(matrix);
	free(vector);
	free(result1);
	free(result2);
	free(mbuff);
	free(vbuff);
	free(resultbuff);
	
	MPI_Finalize();
	return 0;
}