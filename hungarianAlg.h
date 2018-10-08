/************************************************************************
*
*  hungarianAlg.h
author  Andrey Smorodov
modified by Andrea Pennisi

header file for Hungarian algorithm
*
**************************************************************************/

#ifndef _HUNGARIAN_ALG_H_
#define _HUNGARIAN_ALG_H_

#include <vector>
#include <iostream>
#include <iomanip> 
#include <limits>
#include <time.h>
#include <assert.h>

typedef unsigned int uint;
typedef std::vector<uint> UIntVec;
typedef std::vector<UIntVec> UIntMat;
typedef std::vector<bool> BoolVec;
typedef float track_t;
typedef std::vector<int> assignments_t;//配对
typedef std::vector<track_t> distMatrix_t;//距离矩阵


class AssignmentProblemSolver
{
private:
	typedef std::vector<uint> BoolVec;
private:
	// --------------------------------------------------------------------------
	// 通过匈牙利算法计算最小代价匹配.
	// --------------------------------------------------------------------------
	void assignmentoptimal(assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows, const size_t& nOfColumns);
	void buildassignmentvector(assignments_t& assignment, BoolVec& starMatrix, const size_t& nOfRows, const size_t& nOfColumns);
	void computeassignmentcost(const assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows);
	void step2a(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t&  nOfColumns, const size_t& minDim);
	void step2b(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim);
	void step3(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim);
	void step4(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim, const size_t& row, const size_t& col);
	void step5(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim);

public:
	enum TMethod
	{
		optimal,
		many_forbidden_assignments,
		without_forbidden_assignments
	};

public:
	AssignmentProblemSolver() { ; }
	/**
	 * \brief
	 * \param distMatrixIn
	 * \param nOfRows 轨迹的数目
	 * \param nOfColumns 跟踪目标的数目
	 * \param assignment
	 * \param Method
	 * \return
	 */
	track_t Solve(const distMatrix_t& distMatrixIn, const size_t& nOfRows, const size_t& nOfColumns, assignments_t& assignment, const TMethod& Method = optimal);
};


#endif