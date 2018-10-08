#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <time.h> 

#include "hungarianAlg.h"

int main(int argc, char** argv)
{
	const size_t& nOfRows = 4; const size_t& nOfColumns = 4;
	int size_t = nOfRows*nOfColumns;

	//srand((unsigned)time(NULL));
	distMatrix_t distMatrix_Random(size_t);
	assignments_t assignment;
	for (int i = 0; i < size_t;i++)
	{
		distMatrix_Random[i] = rand() / double(RAND_MAX)*100.0f;
	}

	for (int j = 0; j < nOfColumns; j++)
		std::cout << std::setw(10) << j<< "   ";
	std::cout << std::setw(10) << std::endl;
	AssignmentProblemSolver APS;
	for (int i = 0; i < nOfRows; i++)
	{
		for (int j = 0; j < nOfColumns; j++)
			std::cout << std::setw(10) << distMatrix_Random[i*nOfColumns + j] << "   ";
		std::cout<< std::setw(10) << std::endl;
	}

	APS.Solve(distMatrix_Random, nOfRows, nOfColumns, assignment, AssignmentProblemSolver::TMethod::optimal);

	for (int i = 0; i < nOfRows; i++)
	{
		std::cout << std::setw(10) << assignment[i] << std::endl;
	}

	return 0;
}
