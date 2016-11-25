/*							Exercício 1
Faça um programa usando MPI que implemente o padrão mes-
tre/escravo. Inicialmente, o processo mestre deverá enviar uma mensagem
para cada um dos escravos. Ao receberem a mensagem, os processos es-
cravos deverão respondê-la, enviando uma mensagem com o seu rank. Ao
receber uma mensagem de um processo escravo, o processo mestre deverá
imprimir o rank do processo escravo.								 */

#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
	
	int size, rank;
	MPI_Status st;

	MPI_Init(&argc, &argv); 						// inicializa o ambiente MPI

	MPI_Comm_size(MPI_COMM_WORLD, &size);		// pega o numero de processo	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);		// pega o id de cada processo

	int message[size];

	if(rank == 0) {
		for(int i = 1; i < size; i++) {
			MPI_Send(&message, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		MPI_Recv(&message, size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
		for(int i = 1; i < size; i++) {
			printf("Rank do processo: %d\n", message[i]);
		}
	} else {
		MPI_Recv(&message, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
		//int MPI_Iprobe(0, i, MPI_COMM_WORLD, int *flag, *st);
		// if(flag == 1)
		if(rank != 0)
		MPI_Send(&message, size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
	}
	MPI_Finalize();									// Finaliza o ambiente MPI
	return 0;										// Retorna 0 quando o programa exeutou sem problemas
}
