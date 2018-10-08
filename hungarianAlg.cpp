#include "hungarianAlg.h"

track_t AssignmentProblemSolver::Solve(const distMatrix_t& distMatrixIn, const size_t& nOfRows,
	const size_t& nOfColumns, std::vector<int>& assignment, const TMethod& Method)
{
	assignment.resize(nOfRows, -1);

	track_t cost = 0;

	switch (Method)
	{
	case optimal:
		assignmentoptimal(assignment, cost, distMatrixIn, nOfRows, nOfColumns);
		break;
	}

	return cost;
}
// --------------------------------------------------------------------------
// 通过匈牙利算法计算最小代价匹配
// --------------------------------------------------------------------------
void AssignmentProblemSolver::assignmentoptimal(assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows, const size_t& nOfColumns)
{

//	if (nOfRows != nOfColumns)
	std::cout << "原始矩阵，轨迹（行）：" << nOfRows << "    跟踪目标（列）：" << nOfColumns << std::endl;

	const size_t& nOfElements = nOfRows * nOfColumns;
	distMatrix_t distMatrix(nOfElements);
	// 指向最后一个元素
	track_t* distMatrixEnd = distMatrix.data() + nOfElements;

	//生成距离矩阵，并检查元素
	track_t value;
	for (size_t row = 0; row < nOfElements; ++row)
	{
		value = distMatrixIn[row];
		assert(value >= 0);
		distMatrix[row] = value;
	}

	for (int i = 0; i < nOfRows; i++)
	{
		for (int j = 0; j < nOfColumns; j++)
			std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
		std::cout << std::endl;
	}

	// Memory allocation
	BoolVec coveredColumns(nOfColumns, 0);//标记对象被分配了（列中有0）
	BoolVec coveredRows(nOfRows, 0);//标记轨迹已经有对象匹配了
	BoolVec starMatrix(nOfElements, 0);//标记 对象与轨迹进行配对
	BoolVec primeMatrix(nOfElements, 0);
	BoolVec newStarMatrix(nOfElements, 0); /* used in step4 */

	std::cout << " 1、当行小于列时，找出每一行中值最小的元素，然后把该行所有元素都减去这一最小值" << std::endl;
	std::cout << " 2、当行大于列时，找出每一列中值最小的元素，然后把该列所有元素都减去这一最小值" << std::endl;
	if (nOfRows <= nOfColumns)//行小于列
	{
		std::cout << "跟踪目标多于轨迹" << std::endl;
		track_t  minValue;
		track_t* distMatrixTemp;
		track_t* rowEnd;
		track_t* distMatrix_ptr = distMatrix.data();
		for (size_t row = 0; row < nOfRows; row++)
		{
			/* 找出每一行中值最小的元素 */
			distMatrixTemp = distMatrix_ptr + row*nOfColumns;
			rowEnd = distMatrixTemp + nOfColumns-1;//得到行的末尾
			minValue = *distMatrixTemp;
			distMatrixTemp += 1;
			while (distMatrixTemp < rowEnd)
			{
				track_t value = *distMatrixTemp;
				if (value < minValue)
				{
					minValue = value;
				}
				distMatrixTemp += 1;
			}
			/* 然后把该行所有元素都减去这一最小值 */
			distMatrixTemp = distMatrix_ptr + row*nOfColumns;
			while (distMatrixTemp < rowEnd)
			{
				*distMatrixTemp -= minValue;
				distMatrixTemp += 1;
			}
		}
		/* Steps 1 and 2a，将所有包含0的列进行标记 */
		for (size_t row = 0; row < nOfRows; row++)
		{
			for (size_t col = 0; col < nOfColumns; col++)
			{
				if (distMatrix[row*nOfColumns +col] == 0)
				{
					if (!coveredColumns[col])
					{
						starMatrix[row*nOfColumns + col] = true;
						coveredColumns[col] = true;
						coveredRows[row] = true;//标记所在行
						break;
					}
				}
			}
		}
	}
	else /* if(nOfRows > nOfColumns) 行较多，以行为基准进行遍历*/
	{
		std::cout << "跟踪目标少于轨迹" << std::endl;
		track_t* distMatrixTemp;
		track_t* columnEnd;
		track_t  minValue;
		track_t* distMatrix_ptr = distMatrix.data();
		int detal_size = (nOfRows - 1)*nOfColumns;//列末尾
		for (size_t col = 0; col < nOfColumns; col++)
		{
			/*找出每一列中值最小的元素 */
			distMatrixTemp = distMatrix_ptr + col;
			columnEnd = distMatrixTemp + detal_size;
			minValue = *distMatrixTemp;
			distMatrixTemp += nOfColumns;//转到下一行
			while (distMatrixTemp < columnEnd)
			{
				track_t value = *distMatrixTemp;
				if (value < minValue)
				{
					minValue = value;
				}
				distMatrixTemp += nOfColumns;//转到下一行
			}

			/* 把该列中的每一个元素都减去该最小值 */
			distMatrixTemp = distMatrix_ptr + col;
			while (distMatrixTemp < columnEnd)
			{
				*distMatrixTemp -= minValue;
				distMatrixTemp += nOfColumns;//转到下一行
			}
		}
		/* Steps 1 and 2a */
		for (size_t col = 0; col < nOfColumns; col++)
		{
			for (size_t row = 0; row < nOfRows; row++)
			{
				if (distMatrix[row*nOfColumns + col] == 0)
				{
					if (!coveredRows[row])//判断轨迹是否已有对象进行匹配
					{
						starMatrix[row*nOfColumns + col] = true;
						coveredColumns[col] = true;//标记所在列
						coveredRows[row] = true;//标记所在行
						break;
					}//如果当前行已经被覆盖过，当前行的其他0元素不进行列覆盖
				}
			}
		}

		//for (size_t row = 0; row < nOfRows; row++)
		//{
		//	coveredRows[row] = false;
		//}
	}

	for (int i = 0; i < nOfRows; i++)
	{
		for (int j = 0; j < nOfColumns; j++)
			std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
		std::cout << std::endl;
	}
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, (nOfRows <= nOfColumns) ? nOfRows : nOfColumns);
	/* 计算分配成本*/
	computeassignmentcost(assignment, cost, distMatrixIn, nOfRows);

	return;
}
// --------------------------------------------------------------------------
// 生成轨迹-》障碍物对应关系
// --------------------------------------------------------------------------
void AssignmentProblemSolver::buildassignmentvector(assignments_t& assignment, BoolVec& starMatrix, const size_t& nOfRows, const size_t& nOfColumns)
{
	for (size_t row = 0; row < nOfRows; row++)
	{
		for (size_t col = 0; col < nOfColumns; col++)
		{
			if (starMatrix[row*nOfColumns + col])
			{
				assignment[row] = static_cast<int>(col);//轨迹对应-》障碍物
				break;
			}
		}
	}
}
// --------------------------------------------------------------------------
// 生成轨迹-》障碍物对应关系，计算分配成本
// --------------------------------------------------------------------------
void AssignmentProblemSolver::computeassignmentcost(const assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows)
{
	for (size_t row = 0; row < nOfRows; row++)
	{
		const int col = assignment[row];
		if (col >= 0)
		{
			cost += distMatrixIn[row + nOfRows * col];
		}
	}
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step2a(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t&  nOfColumns, const size_t& minDim)
{
	uint *starMatrixTemp, *columnEnd;
	/* 覆盖所有包含起始0的列 */
	uint* starMatrix_ptr = starMatrix.data();
	int detal_size = (nOfRows - 1)*nOfColumns;//列末尾
	for (size_t col = 0; col < nOfColumns; col++)
	{
		starMatrixTemp = starMatrix_ptr + col;//转到下一列
		columnEnd = starMatrixTemp + detal_size;//列结尾
		while (starMatrixTemp < columnEnd)
		{
			if (*starMatrixTemp)
			{
				coveredColumns[col] = true;
				break;
			}
			starMatrixTemp += nOfColumns;//转到下一行
		}
	}
	/* move to step 3 */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

// --------------------------------------------------------------------------
// 判断是否完成分配，没有完成分配则进入步骤3（用尽量少的横线、竖线覆盖矩阵中的所有0）
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step2b(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	/* 统计被覆盖的列 */
	size_t nOfCoveredColumns = 0;
	for (size_t col = 0; col < nOfColumns; col++)
	{
		if (coveredColumns[col])
		{
			nOfCoveredColumns++;
		}
	}
	if (nOfCoveredColumns == minDim)//包含0的列数等于最小维数，则完成分配
	{
		/*完成分配*/
		buildassignmentvector(assignment, starMatrix, nOfRows, nOfColumns);
		std::cout << "完成分配" << std::endl;
	}
	else//存在对象没有完成分类
	{
		step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}
}

// --------------------------------------------------------------------------
// 用尽量少的横线、竖线覆盖矩阵中的所有0
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step3(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	bool zerosFound = true;
	while (zerosFound)
	{
		zerosFound = false;
		for (size_t col = 0; col < nOfColumns; col++)//遍历列
		{
			if (!coveredColumns[col])//当前对象还没有被分配
			{
				for (size_t row = 0; row < nOfRows; row++)
				{
					if ((!coveredRows[row]) && (distMatrix[row*nOfColumns + col] == 0))//当前轨迹还没有被分配，且存在对象与其匹配
					{
						/* prime zero */
						primeMatrix[row*nOfColumns + col] = true;
						/* 找到当前行中的第一个0 */
						size_t starCol = 0;
						for (; starCol < nOfColumns; starCol++)
						{
							if (starMatrix[row*nOfColumns + starCol])
							{
								break;
							}
						}
						if (starCol == nOfColumns) /* 当前轨迹没有对象与其进行配对*/
						{
							/* move to step 4 */
							step4(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
							return;
						}
						else
						{
							coveredRows[row] = true;//标记当前轨迹有配对
							coveredColumns[starCol] = false;//标记当前对象没有配对
							zerosFound = true;//找到配对
							break;
						}
					}
				}
			}
		}
	}
	/* move to step 5 */
	step5(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step4(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim, const size_t& row, const size_t& col)
{
	const size_t& nOfElements = nOfRows * nOfColumns;
	/* 对 starMatrix 进行拷贝*/
	for (size_t n = 0; n < nOfElements; n++)
	{
		newStarMatrix[n] = starMatrix[n];
	}
	/* star current zero */
	newStarMatrix[row*nOfColumns + col] = true;
	/* 找到当前列的起始0*/
	size_t starCol = col;
	size_t starRow = 0;
	for (; starRow < nOfRows; starRow++)
	{
		if (starMatrix[starRow * nOfColumns + starCol])
		{
			break;
		}
	}
	while (starRow < nOfRows)
	{
		/* 对标记的0取消标记 */
		newStarMatrix[starRow*nOfColumns + starCol] = false;
		/* 找到当前行中的第一个0 */
		size_t primeRow = starRow;
		size_t primeCol = 0;
		for (; primeCol < nOfColumns; primeCol++)
		{
			if (primeMatrix[primeRow*nOfColumns + primeCol])
			{
				break;
			}
		}
		/* star the primed zero */
		newStarMatrix[primeRow*nOfColumns + primeCol] = true;
		/* 找到当前列中的第一个0 */
		starCol = primeCol;
		for (starRow = 0; starRow < nOfRows; starRow++)
		{
			if (starMatrix[starRow*nOfColumns + starCol])
			{
				break;
			}
		}
	}
	/* use temporary copy as new starMatrix */
	/* delete all primes, uncover all rows */
	for (size_t n = 0; n < nOfElements; n++)
	{
		primeMatrix[n] = false;
		starMatrix[n] = newStarMatrix[n];
	}
	//释放所有行
	for (size_t n = 0; n < nOfRows; n++)
	{
		coveredRows[n] = false;
	}
	/* move to step 2a */
	step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

// --------------------------------------------------------------------------
// 从上一步中未被覆盖的元素中找到最小值，然后把这些元素都减去最这一小值、给覆盖交叉点的元素加上这一最小值
// 重复步骤3、4
// 被覆盖元素中的最小值实际上是完成所有任务过程中不可避免的开销。 
// 这一步的作用是增加开销矩阵中0的个数，使得任务更易分配
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step5(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	/* 从上一步中未被覆盖的元素中找到最小值 h */
	float h = std::numeric_limits<track_t>::max();
	for (size_t row = 0; row < nOfRows; row++)
	{
		if (!coveredRows[row])
		{
			for (size_t col = 0; col < nOfColumns; col++)
			{
				if (!coveredColumns[col])
				{
					const float& value = distMatrix[row*nOfColumns + col];
					if (value < h)
					{
						h = value;
					}
				}
			}
		}
	}
	//给没被覆盖的元素统一减去最小值，给被十字交叉覆盖的元素加上最小值

	//给被覆盖的行上所有元素加上这一最小值
	for (size_t row = 0; row < nOfRows; row++)
	{
		if (coveredRows[row])
		{
			for (size_t col = 0; col < nOfColumns; col++)
			{
				distMatrix[row + nOfRows*col] += h;
			}
		}
	}
	//给没被覆盖的列上的所有元素都减去最这一小值
	for (size_t col = 0; col < nOfColumns; col++)
	{
		if (!coveredColumns[col])
		{
			for (size_t row = 0; row < nOfRows; row++)
			{
				distMatrix[row + nOfRows*col] -= h;
			}
		}
	}
	std::cout << "从上一步中未被覆盖的元素中找到最小值，然后把这些元素都减去最这一小值、给覆盖交叉点的元素加上这一最小值" << std::endl;
	for (int i = 0; i < nOfRows; i++)
	{
		for (int j = 0; j < nOfColumns; j++)
			std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
		std::cout << std::endl;
	}
	step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}
