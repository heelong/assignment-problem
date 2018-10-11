#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <time.h> 
#include <windows.h>
#include "hungarianAlg.h"

int main(int argc, char** argv)
{
	clock_t time1, time2;
	srand((unsigned)time(NULL));
	int a = 1, b = 30;
	while (1)
	{
		const size_t& nOfRows = 30;// rand() % (b - a);
		const size_t& nOfColumns = 30;// rand() % (b - a);
		int size_t = nOfRows*nOfColumns;
		//float data[] = { 0.125, 56.000, 19.000, 80.000, 0.135, 152.000, 36.000, 33.000,0,0,
		//	58.000, 0.100, 35.000, 89.000, 6.500, 4.000, 254.000, 323.000, 0, 0,
		//	82.000, 74.000, 17.000, 85.000, 1.200, 152.000, 23.000, 15.000, 0, 0,
		//	71.000, 51.000, 2.000, 1.400, 24.000, 112.000, 23.000, 0.340, 0, 0,
		//	9.000, 36.000, 14.000, 16.000, 36.000, 151.000, 25.000, 158.000, 0, 0,
		//	0.100, 11.000, 13.000, 0.500, 15.000, 36.000, 33.000, 15.000, 0, 0,
		//	136.000, 152.000, 12.000, 0.100, 15.000, 30.000, 365.000, 156.000, 0, 0,
		//	101.000, 11.000, 0.100, 56.000, 15.000, 1456.000, 562.000, 10.100, 0, 0,
		//	112.000, 6896.000, 156.000, 2336.000, 235.000, 1.000, 11.000, 0.010, 0, 0,
		//	150.000, 11.000, 135.000, 158.000, 0.100, 0.001, 1.000, 2.000,0, 0 };
		float data[] = { 0.125, 56.000, 19.000, 80.000, 0.135, 152.000, 36.000, 33.000, 
			58.000, 0.100, 35.000, 89.000, 6.500, 4.000, 254.000, 323.000,
			82.000, 74.000, 17.000, 85.000, 1.200, 152.000, 23.000, 15.000, 
			71.000, 51.000, 2.000, 1.400, 24.000, 112.000, 23.000, 0.340,
			9.000, 36.000, 14.000, 16.000, 36.000, 151.000, 25.000, 158.000, 
			0.100, 11.000, 13.000, 0.500, 15.000, 36.000, 33.000, 15.000,
			136.000, 152.000, 12.000, 0.100, 15.000, 30.000, 365.000, 156.000,
			101.000, 11.000, 0.100, 56.000, 15.000, 1456.000, 562.000, 10.100,
			112.000, 6896.000, 156.000, 2336.000, 235.000, 1.000, 11.000, 0.010, 
			150.000, 11.000, 135.000, 158.000, 0.100, 0.001, 1.000, 2.000 };
		//srand((unsigned)time(NULL));
		distMatrix_t distMatrix_Random(size_t);
		assignments_t assignment;
		for (int i = 0; i < size_t; i++)
		{
			distMatrix_Random[i] = rand() / double(RAND_MAX)*100.0f;
		}

		//for (int j = 0; j < nOfColumns; j++)
		//	std::cout << std::setw(10) << j << "   ";
		//std::cout << std::setw(10) << std::endl;

		//for (int i = 0; i < nOfRows; i++)
		//{
		//	for (int j = 0; j < nOfColumns; j++)
		//		std::cout << std::setw(10) << distMatrix_Random[i*nOfColumns + j] << "   ";
		//	std::cout << std::setw(10) << std::endl;
		//}
		time1 = clock();
		AssignmentProblemSolver APS;
		float cost = APS.Solve(distMatrix_Random, nOfRows, nOfColumns, assignment, AssignmentProblemSolver::TMethod::optimal);
		std::cout << cost << std::endl << std::endl;
		time2 = clock();
		std::cout << "Time 1 :" << (double)(time2 - time1) / CLOCKS_PER_SEC*1000.0 << "ms" << std::endl;


		//if (cost<0.0001)
		//{
		//for (int i = 0; i < nOfRows; i++)
		//{
		//	for (int j = 0; j < nOfColumns; j++)
		//		std::cout << std::setw(10) << distMatrix_Random[i*nOfColumns + j] << "   ";
		//	std::cout << std::setw(10) << std::endl;
		//}
		for (int i = 0; i < nOfRows; i++)
		{
			std::cout << "第" << std::setw(4) << i << "号轨迹匹配第" << std::setw(10) << assignment[i] << "对象" << std::endl;
		}
		Sleep(1000);
		//	Sleep(10000);
		//}

	}


	return 0;
}
