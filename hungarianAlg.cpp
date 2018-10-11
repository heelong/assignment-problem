#include "hungarianAlg.h"

track_t AssignmentProblemSolver::Solve(const distMatrix_t& distMatrixIn, const size_t& nOfRows,
	const size_t& nOfColumns, std::vector<int>& assignment, const TMethod& Method)
{
	size_t maxDim = (nOfRows >= nOfColumns) ? nOfRows : nOfColumns;

	const size_t& nOfRows_new = maxDim; const size_t& nOfColumns_new = maxDim;
	const size_t& nOfElements = nOfRows_new * nOfColumns_new;
	if (nOfElements == 0)
	{
		if (nOfRows == 0)
		{
			std::cout << "无轨迹进行匹配" << std::endl;
		}
		if (nOfColumns == 0)
		{
			std::cout << "无对象进行匹配" << std::endl;
		}
		return 0;
	}
	distMatrix_t distMatrix_Square(nOfElements,0);
	std::vector<int> assignment_new(maxDim, -1);
	//生成距离矩阵，并检查元素
	track_t value;
	for (size_t row = 0; row < nOfRows; row++)
	{
		for (size_t col = 0; col < nOfColumns; col++)
		{
			value = distMatrixIn[row*nOfColumns + col];
			if (value < 0.0)
				std::cout << "距离矩阵计算有误" << std::endl;
			distMatrix_Square[row*nOfColumns_new + col] = value;
		}
	}
	//std::cout << std::setw(10) << "   ";
	//for (int j = 0; j < nOfColumns_new; j++)
	//	std::cout << std::setw(10) << j << "   ";
	//std::cout << std::setw(10) << std::endl;
	//for (int i = 0; i < nOfRows_new; i++)
	//{
	//	std::cout << std::setw(10) << i << "   ";
	//	for (int j = 0; j < nOfColumns_new; j++)
	//		std::cout << std::setw(10) << distMatrix_Square[i*nOfColumns_new + j] << "   ";
	//	std::cout << std::endl;
	//}
	track_t cost = 0;
	switch (Method)
	{
	case optimal:
		assignmentoptimal(assignment_new, cost, distMatrix_Square, nOfRows_new, nOfColumns_new);
		break;
	default:
		break;
	}
	assignment.resize(nOfRows, -1);
	for (size_t row = 0; row < nOfRows; row++)
	{
		assignment[row] = (assignment_new[row] >= nOfColumns) ? -1 : assignment_new[row];//轨迹对应-》障碍物
	}
	return cost;
}
// --------------------------------------------------------------------------
// 通过匈牙利算法计算最小代价匹配
// --------------------------------------------------------------------------
void AssignmentProblemSolver::assignmentoptimal(assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows, const size_t& nOfColumns)
{

	//std::cout << "原始矩阵，轨迹（行）：" << nOfRows << "    跟踪目标（列）：" << nOfColumns << std::endl;

	const size_t& nOfElements = nOfRows * nOfColumns;

	distMatrix_t distMatrix(nOfElements);
	// 指向最后一个元素
	track_t* distMatrixEnd = distMatrix.data() + nOfElements;

	//生成距离矩阵，并检查元素
	track_t value;
	for (size_t row = 0; row < nOfElements; ++row)
	{
		value = distMatrixIn[row];
		if (value < 0.0)
			std::cout << "距离矩阵计算有误" << std::endl;
		distMatrix[row] = value;
	}

	//for (int i = 0; i < nOfRows; i++)
	//{
	//	for (int j = 0; j < nOfColumns; j++)
	//		std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
	//	std::cout << std::endl;
	//}

	BoolVec coveredColumns(nOfColumns, 0);//标记对象被分配了（列中有0）
	BoolVec coveredRows(nOfRows, 0);//标记轨迹已经有对象匹配了
	BoolVec starMatrix(nOfElements, 0);//标记 对象与轨迹进行配对
	BoolVec primeMatrix(nOfElements, 0);
	BoolVec newStarMatrix(nOfElements, 0); /* used in step4 */

	 markRows.resize(nOfRows, 0);//标记仍未分配对象的轨迹
	 markColumns.resize(nOfColumns, 0);//标记所有 **刚被标记过的行中** 0所在的列
	 CircleZero.resize(nOfElements, 0);//圈住的0元素
	 SlashZero.resize(nOfElements, 0); //划掉的0元素

	//std::cout << " 1、当行小于列时，找出每一行中值最小的元素，然后把该行所有元素都减去这一最小值" << std::endl;
	//std::cout << " 2、当行大于列时，找出每一列中值最小的元素，然后把该列所有元素都减去这一最小值" << std::endl;
	if (nOfRows <= nOfColumns)//行小于列
	{
		//std::cout << "跟踪目标多于轨迹" << std::endl;
		track_t  minValue;
		track_t* distMatrixTemp;
		track_t* rowEnd;
		track_t* distMatrix_ptr = distMatrix.data();
		for (size_t row = 0; row < nOfRows; row++)
		{
			/* 找出每一行中值最小的元素 */
			distMatrixTemp = distMatrix_ptr + row*nOfColumns;
			rowEnd = distMatrixTemp + nOfColumns;//得到行的末尾
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
		for (size_t row = 0; row < nOfRows; row++)
		{
			for (size_t col = 0; col < nOfColumns; col++)
			{
				if (distMatrix[row*nOfColumns +col] < 0.00001)
				{
					if (!coveredColumns[col])
					{
						starMatrix[row*nOfColumns + col] = true;
						coveredColumns[col] = true;
						//coveredRows[row] = true;//标记所在行
						break;
					}
				}
			}
		}
	}
	else /* 行较多，以行为基准进行遍历*/
	{
		track_t* distMatrixTemp;
		track_t* columnEnd;
		track_t  minValue;
		track_t* distMatrix_ptr = distMatrix.data();
		int detal_size = (nOfRows - 1)*nOfColumns+1;//列末尾
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

		for (size_t col = 0; col < nOfColumns; col++)
		{
			for (size_t row = 0; row < nOfRows; row++)
			{
				if (distMatrix[row*nOfColumns + col] < 0.00001)
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
		for (size_t row = 0; row < nOfRows; row++)
		{
			coveredRows[row] = false;
		}
	}
	//std::cout << "减去每行/列最小值"<<std::endl;
	//for (int i = 0; i < nOfRows; i++)
	//{
	//	for (int j = 0; j < nOfColumns; j++)
	//		std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
	//	std::cout << std::endl;
	//}
	step2b(assignment,distMatrixIn ,distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, (nOfRows <= nOfColumns) ? nOfRows : nOfColumns);
	/* 计算分配成本*/
	computeassignmentcost(assignment, cost, distMatrixIn, nOfRows, nOfColumns);

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
void AssignmentProblemSolver::computeassignmentcost(const assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows, const size_t& nOfColumns)
{
	for (size_t row = 0; row < nOfRows; row++)
	{
		const int col = assignment[row];
		if (col >= 0)
		{
			cost += distMatrixIn[row*nOfColumns +  col];
		}
	}
}

// --------------------------------------------------------------------------
// 可以完成分配进入2a
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step2a(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t&  nOfColumns, const size_t& minDim)
{
	uint *starMatrixTemp, *columnEnd;
	/* 检查每一列是否包含独立0元素*/
	uint* starMatrix_ptr = starMatrix.data();
	int detal_size = (nOfRows - 1)*nOfColumns+1;//列末尾
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
	step2b(assignment, distMatrixIn, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

// --------------------------------------------------------------------------
// 判断是否完成分配，没有完成分配则进入步骤3（用尽量少的横线、竖线覆盖矩阵中的所有0）
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step2b(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
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
		//std::cout << "完成分配" << std::endl;
	}
	else//存在对象没有完成分类
	{
		step3(assignment,distMatrixIn, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}
}

// --------------------------------------------------------------------------
// 用尽量少的横线、竖线覆盖矩阵中的所有0
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step3(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	//用尽量少的横线竖线去覆盖矩阵，当覆盖列数等于矩阵最小纬度时，完成匹配
	const size_t& nOfElements = nOfRows * nOfColumns;
	while (1)
	{
		for (size_t row = 0; row < nOfRows; row++)
		{
			markRows[row] = false;
			for (size_t col = 0; col < nOfColumns; col++)//遍历row_O行的所有列,找到原始0
			{
				CircleZero[row*nOfColumns + col] = 0;
				SlashZero[row*nOfColumns + col] = 0;
				markColumns[col] = false;
			}
		}
		size_t CircleZeroNum = 0;          //圈住的0元素的数目
		//找出独立0元素
		step3b(assignment, distMatrixIn, distMatrix, starMatrix, CircleZero, SlashZero, CircleZeroNum, nOfRows, nOfColumns, minDim);

		if (CircleZeroNum == minDim)//判断独立0元素的个数，如果与距离矩阵最小维度相等，则
		{
			step2a(assignment, distMatrixIn, distMatrix, CircleZero, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
			return;
		}

		//不存在未被标记过的零元素，但圈0的个数<n
		//作最少直线覆盖当前所有零元素，便于下步增加独立零元素的个数。

		//2、打钩
		//2.1、无圈0的行上打钩（已在上一步完成）
		for (size_t row = 0; row < nOfRows; row++)
		{
			markRows[row] = true;
			for (size_t col = 0; col < nOfColumns; col++)//遍历row_O行的所有列,找到原始0
				if (CircleZero[row*nOfColumns + col])
				{
					markRows[row] = false;
					break;
				}
		}
		//2.2、打钩的行上被划掉的0所在列打钩，进入下一步
		//2.3、打钩的列上有圈0的行打钩，进入上一步
		bool zeorFrond = true;
		while (zeorFrond)
		{
			zeorFrond = false;
			//2.2、打钩的行上被划掉的0所在列打钩----------------
			for (size_t row = 0; row < nOfRows; row++)
			{
				if (markRows[row])
					for (size_t col = 0; col < nOfColumns; col++)//遍历列
					{
						if (!markColumns[col] && (distMatrix[row*nOfColumns + col] < 0.00001) && SlashZero[row*nOfColumns + col])
						{
							markColumns[col] = true;
							zeorFrond = true;
						}
					}
			}
			//2.3、打钩的列上有圈0的行打钩-----------------------
			for (size_t col = 0; col < nOfColumns; col++)//遍历列
			{
				if (markColumns[col])
					for (size_t row = 0; row < nOfRows; row++)
					{
						if (!markRows[row] && (distMatrix[row*nOfColumns + col] < 0.00001) && CircleZero[row*nOfColumns + col])
						{
							markRows[row] = true;
							zeorFrond = true;
						}
					}
			}
		}

		size_t nOfCoveredColumns = 0;
		size_t nOfCoveredRow = 0;
		//4、在所有标记过的列和未标记的行上画线 
		for (size_t col = 0; col < nOfColumns; col++)//遍历列
		{
			coveredColumns[col] = false;
			if (markColumns[col])
			{
				coveredColumns[col] = true;
				nOfCoveredColumns++;
			}
		}
		for (size_t row = 0; row < nOfRows; row++)
		{
			coveredRows[row] = false;
			if (!markRows[row])
			{
				coveredRows[row] = true;
				nOfCoveredRow++;
			}
		}

		//step4(assignment,distMatrixIn, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
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
		//给没被覆盖的元素统一减去最小值，给被十字交叉覆盖的元素加上最小值-----------------------

		//给被覆盖的行上所有元素加上这一最小值
		for (size_t row = 0; row < nOfRows; row++)
		{
			if (coveredRows[row])
			{
				for (size_t col = 0; col < nOfColumns; col++)
				{
					distMatrix[row*nOfColumns + col] += h;
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
					distMatrix[row*nOfColumns + col] -= h;
				}
			}
		}
	}

}

// --------------------------------------------------------------------------
// 进行试分配，判断是否存在n个独立零元素
// 尝试对所有零元素做标记，确定独立零元素。
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step3b(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& CircleZero, BoolVec& SlashZero, size_t& CircleZeroNum, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	size_t notMatchedZeroNum = 0;//没有被标记的0的数目
	for (size_t row_O = 0; row_O < nOfRows; row_O++)
	{
		for (size_t col = 0; col < nOfColumns; col++)//遍历row_O行的所有列,找到原始0
			if ((distMatrix[row_O*nOfColumns + col] < 0.00001) && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
				notMatchedZeroNum++;
	}
	//用尽量少的横线竖线去覆盖矩阵，当覆盖列数等于矩阵最小维度时，完成匹配
	do
	{
		bool zeorFrond = true;
		//1.1、逐行检查、找到每行 >>>只有一个<<< 没有被圈住或划掉的0元素，将其圈住，并将其所在行与列的其他未标记的0元素划掉---------------------
		while (zeorFrond&&notMatchedZeroNum>0)
		{
			zeorFrond = false;
			for (size_t row_O = 0; row_O < nOfRows; row_O++)
			{
				size_t col_O = 0;
				size_t notMatchedZeroNumInRow_O = 0;
				for (size_t col = 0; col < nOfColumns; col++)//遍历row_O行的所有列,找到原始0
					if ((distMatrix[row_O*nOfColumns + col] < 0.00001) && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
					{
						notMatchedZeroNumInRow_O++;
						col_O = col;
					}
				if (notMatchedZeroNumInRow_O == 0 || notMatchedZeroNumInRow_O > 1)
					continue;
				else
				{
					CircleZero[row_O*nOfColumns + col_O] = 1;//找到圈住的0
					notMatchedZeroNum--;
					starMatrix[row_O*nOfColumns + col_O] = 1;
					CircleZeroNum++;
					zeorFrond = true;
				}
				//将被圈起的零元素所在行与列的其他未标记的零元素用记号 / 划去
				//for (size_t col = 0; col < nOfColumns; col++)////遍历row_O行的所有列,将未标记的0划掉
				//{
				//	if (distMatrix[row_O*nOfColumns + col] < 0.00001 && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
				//		SlashZero[row_O*nOfColumns + col] = 1;
				//}
				for (size_t row = 0; row < nOfRows; row++)////遍历col_O列的所有行,将未标记的0划掉
				{
					if (distMatrix[row*nOfColumns + col_O] < 0.00001 && (!CircleZero[row*nOfColumns + col_O] && !SlashZero[row*nOfColumns + col_O]))
					{
						SlashZero[row*nOfColumns + col_O] = 1;
						notMatchedZeroNum--;
					}
				}
			}
		}
		//1.2、逐列检查、找到每列 >>>只有一个<<< 没有被圈住或划掉的0元素，将其圈住，并将其所在行与列的其他未标记的0元素划掉---------------------
		zeorFrond = true;
		while (zeorFrond&&notMatchedZeroNum>0)
		{
			zeorFrond = false;
			for (size_t col_O = 0; col_O < nOfColumns; col_O++)
			{
				size_t row_O = 0;
				size_t notMatchedZeroNumInCol_O = 0;//没有标记的0
				for (size_t row = 0; row < nOfRows; row++)//遍历row_O行的所有列,找到原始0
					if ((distMatrix[row*nOfColumns + col_O] < 0.00001) && (!CircleZero[row*nOfColumns + col_O] && !SlashZero[row*nOfColumns + col_O]))
					{
						notMatchedZeroNumInCol_O++;
						row_O = row;
					}
				if (notMatchedZeroNumInCol_O == 0 || notMatchedZeroNumInCol_O > 1)
					continue;
				else
				{
					CircleZero[row_O*nOfColumns + col_O] = 1;//找到圈住的0
					notMatchedZeroNum--;
					starMatrix[row_O*nOfColumns + col_O] = 1;
					CircleZeroNum++;
					zeorFrond = true;
				}
				//将被圈起的零元素所在行与列的其他未标记的零元素用记号 / 划去
				for (size_t col = 0; col < nOfColumns; col++)////遍历row_O行的所有列,将未标记的0划掉
				{
					if (distMatrix[row_O*nOfColumns + col] < 0.00001 && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
					{
						SlashZero[row_O*nOfColumns + col] = 1;
						notMatchedZeroNum--;
					}
				}
				//for (size_t row = 0; row < nOfRows; row++)////遍历col_O列的所有行,将未标记的0划掉
				//{
				//	if (distMatrix[row*nOfColumns + col_O] < 0.00001 && (!CircleZero[row*nOfColumns + col_O] && !SlashZero[row*nOfColumns + col_O]))
				//		SlashZero[row*nOfColumns + col_O] = 1;
				//}
			}
		}
		//if ((CircleZeroNum < minDim))//每一行均有圈0出现，圈0的个数恰好等于n
		//	break;
		if (notMatchedZeroNum>0)
		{
			size_t row_O = -1;
			size_t col_O = -1;
			float minCost = (std::numeric_limits<float>::max)();//遍历、找到所有
			for (size_t row=0; row < nOfRows; row++)
			{
				for (size_t col=0; col < nOfColumns; col++)//遍历row_O行的所有列,找到原始0
					if ((distMatrix[row*nOfColumns + col] < 0.00001) && (!CircleZero[row*nOfColumns + col] && !SlashZero[row*nOfColumns + col]))
					{
						if (minCost>distMatrixIn[row*nOfColumns + col])
						{
							minCost = distMatrixIn[row*nOfColumns + col];
							row_O = row;
							col_O = col;
						}
					}
			}
			if ((row_O < nOfRows&&row_O >= 0)&&
				(col_O < nOfColumns&&col_O >= 0))
			{
				CircleZero[row_O*nOfColumns + col_O] = 1;
				CircleZeroNum++;
				notMatchedZeroNum--;
				//将被圈起的零元素所在行与列的其他未标记的零元素用记号 / 划去
				for (size_t col = 0; col < nOfColumns; col++)////遍历row_O行的所有列,将未标记的0划掉
				{
					if (distMatrix[row_O*nOfColumns + col] < 0.00001 && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
					{
						notMatchedZeroNum--;
						SlashZero[row_O*nOfColumns + col] = 1;
					}
				}
				for (size_t row = 0; row < nOfRows; row++)////遍历col_O列的所有行,将未标记的0划掉
				{
					if (distMatrix[row*nOfColumns + col_O] < 0.00001 && (!CircleZero[row*nOfColumns + col_O] && !SlashZero[row*nOfColumns + col_O]))
					{
						notMatchedZeroNum--;
						SlashZero[row*nOfColumns + col_O] = 1;
					}
				}
			}
		}
	} while ((CircleZeroNum < minDim) && notMatchedZeroNum>0);

	return ;
}

// --------------------------------------------------------------------------
// 从上一步中未被覆盖的元素中找到最小值，然后把这些元素都减去最这一小值、给覆盖交叉点的元素加上这一最小值
// 重复步骤3、4
// 被覆盖元素中的最小值实际上是完成所有任务过程中不可避免的开销。 
// 这一步的作用是增加开销矩阵中0的个数，使得任务更易分配
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step4(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
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
	//给没被覆盖的元素统一减去最小值，给被十字交叉覆盖的元素加上最小值-----------------------

	//给被覆盖的行上所有元素加上这一最小值
	for (size_t row = 0; row < nOfRows; row++)
	{
		if (coveredRows[row])
		{
			for (size_t col = 0; col < nOfColumns; col++)
			{
				distMatrix[row*nOfColumns + col] += h;
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
				distMatrix[row*nOfColumns + col] -= h;
			}
		}
	}
	//std::cout << "从上一步中未被覆盖的元素中找到最小值，然后把这些元素都减去最这一小值、给覆盖交叉点的元素加上这一最小值" << std::endl;
	//for (int i = 0; i < nOfRows; i++)
	//{
	//	for (int j = 0; j < nOfColumns; j++)
	//		std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
	//	std::cout << std::endl;
	//}

	step3(assignment,distMatrixIn, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}
