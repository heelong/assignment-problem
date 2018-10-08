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
// ͨ���������㷨������С����ƥ��
// --------------------------------------------------------------------------
void AssignmentProblemSolver::assignmentoptimal(assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows, const size_t& nOfColumns)
{

//	if (nOfRows != nOfColumns)
	std::cout << "ԭʼ���󣬹켣���У���" << nOfRows << "    ����Ŀ�꣨�У���" << nOfColumns << std::endl;

	const size_t& nOfElements = nOfRows * nOfColumns;
	distMatrix_t distMatrix(nOfElements);
	// ָ�����һ��Ԫ��
	track_t* distMatrixEnd = distMatrix.data() + nOfElements;

	//���ɾ�����󣬲����Ԫ��
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
	BoolVec coveredColumns(nOfColumns, 0);//��Ƕ��󱻷����ˣ�������0��
	BoolVec coveredRows(nOfRows, 0);//��ǹ켣�Ѿ��ж���ƥ����
	BoolVec starMatrix(nOfElements, 0);//��� ������켣�������
	BoolVec primeMatrix(nOfElements, 0);
	BoolVec newStarMatrix(nOfElements, 0); /* used in step4 */

	std::cout << " 1������С����ʱ���ҳ�ÿһ����ֵ��С��Ԫ�أ�Ȼ��Ѹ�������Ԫ�ض���ȥ��һ��Сֵ" << std::endl;
	std::cout << " 2�����д�����ʱ���ҳ�ÿһ����ֵ��С��Ԫ�أ�Ȼ��Ѹ�������Ԫ�ض���ȥ��һ��Сֵ" << std::endl;
	if (nOfRows <= nOfColumns)//��С����
	{
		std::cout << "����Ŀ����ڹ켣" << std::endl;
		track_t  minValue;
		track_t* distMatrixTemp;
		track_t* rowEnd;
		track_t* distMatrix_ptr = distMatrix.data();
		for (size_t row = 0; row < nOfRows; row++)
		{
			/* �ҳ�ÿһ����ֵ��С��Ԫ�� */
			distMatrixTemp = distMatrix_ptr + row*nOfColumns;
			rowEnd = distMatrixTemp + nOfColumns-1;//�õ��е�ĩβ
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
			/* Ȼ��Ѹ�������Ԫ�ض���ȥ��һ��Сֵ */
			distMatrixTemp = distMatrix_ptr + row*nOfColumns;
			while (distMatrixTemp < rowEnd)
			{
				*distMatrixTemp -= minValue;
				distMatrixTemp += 1;
			}
		}
		/* Steps 1 and 2a�������а���0���н��б�� */
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
						coveredRows[row] = true;//���������
						break;
					}
				}
			}
		}
	}
	else /* if(nOfRows > nOfColumns) �н϶࣬����Ϊ��׼���б���*/
	{
		std::cout << "����Ŀ�����ڹ켣" << std::endl;
		track_t* distMatrixTemp;
		track_t* columnEnd;
		track_t  minValue;
		track_t* distMatrix_ptr = distMatrix.data();
		int detal_size = (nOfRows - 1)*nOfColumns;//��ĩβ
		for (size_t col = 0; col < nOfColumns; col++)
		{
			/*�ҳ�ÿһ����ֵ��С��Ԫ�� */
			distMatrixTemp = distMatrix_ptr + col;
			columnEnd = distMatrixTemp + detal_size;
			minValue = *distMatrixTemp;
			distMatrixTemp += nOfColumns;//ת����һ��
			while (distMatrixTemp < columnEnd)
			{
				track_t value = *distMatrixTemp;
				if (value < minValue)
				{
					minValue = value;
				}
				distMatrixTemp += nOfColumns;//ת����һ��
			}

			/* �Ѹ����е�ÿһ��Ԫ�ض���ȥ����Сֵ */
			distMatrixTemp = distMatrix_ptr + col;
			while (distMatrixTemp < columnEnd)
			{
				*distMatrixTemp -= minValue;
				distMatrixTemp += nOfColumns;//ת����һ��
			}
		}
		/* Steps 1 and 2a */
		for (size_t col = 0; col < nOfColumns; col++)
		{
			for (size_t row = 0; row < nOfRows; row++)
			{
				if (distMatrix[row*nOfColumns + col] == 0)
				{
					if (!coveredRows[row])//�жϹ켣�Ƿ����ж������ƥ��
					{
						starMatrix[row*nOfColumns + col] = true;
						coveredColumns[col] = true;//���������
						coveredRows[row] = true;//���������
						break;
					}//�����ǰ���Ѿ������ǹ�����ǰ�е�����0Ԫ�ز������и���
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
	/* �������ɱ�*/
	computeassignmentcost(assignment, cost, distMatrixIn, nOfRows);

	return;
}
// --------------------------------------------------------------------------
// ���ɹ켣-���ϰ����Ӧ��ϵ
// --------------------------------------------------------------------------
void AssignmentProblemSolver::buildassignmentvector(assignments_t& assignment, BoolVec& starMatrix, const size_t& nOfRows, const size_t& nOfColumns)
{
	for (size_t row = 0; row < nOfRows; row++)
	{
		for (size_t col = 0; col < nOfColumns; col++)
		{
			if (starMatrix[row*nOfColumns + col])
			{
				assignment[row] = static_cast<int>(col);//�켣��Ӧ-���ϰ���
				break;
			}
		}
	}
}
// --------------------------------------------------------------------------
// ���ɹ켣-���ϰ����Ӧ��ϵ���������ɱ�
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
	/* �������а�����ʼ0���� */
	uint* starMatrix_ptr = starMatrix.data();
	int detal_size = (nOfRows - 1)*nOfColumns;//��ĩβ
	for (size_t col = 0; col < nOfColumns; col++)
	{
		starMatrixTemp = starMatrix_ptr + col;//ת����һ��
		columnEnd = starMatrixTemp + detal_size;//�н�β
		while (starMatrixTemp < columnEnd)
		{
			if (*starMatrixTemp)
			{
				coveredColumns[col] = true;
				break;
			}
			starMatrixTemp += nOfColumns;//ת����һ��
		}
	}
	/* move to step 3 */
	step2b(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

// --------------------------------------------------------------------------
// �ж��Ƿ���ɷ��䣬û����ɷ�������벽��3���þ����ٵĺ��ߡ����߸��Ǿ����е�����0��
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step2b(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	/* ͳ�Ʊ����ǵ��� */
	size_t nOfCoveredColumns = 0;
	for (size_t col = 0; col < nOfColumns; col++)
	{
		if (coveredColumns[col])
		{
			nOfCoveredColumns++;
		}
	}
	if (nOfCoveredColumns == minDim)//����0������������Сά��������ɷ���
	{
		/*��ɷ���*/
		buildassignmentvector(assignment, starMatrix, nOfRows, nOfColumns);
		std::cout << "��ɷ���" << std::endl;
	}
	else//���ڶ���û����ɷ���
	{
		step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}
}

// --------------------------------------------------------------------------
// �þ����ٵĺ��ߡ����߸��Ǿ����е�����0
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step3(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	//�þ����ٵĺ�������ȥ���Ǿ��󣬵������������ھ�����Сγ��ʱ�����ƥ��
	bool zerosFound = true;
	while (zerosFound)
	{
		zerosFound = false;
		for (size_t col = 0; col < nOfColumns; col++)//������
		{
			if (!coveredColumns[col])//��ǰ����û�б�����
			{
				for (size_t row = 0; row < nOfRows; row++)
				{
					if ((!coveredRows[row]) && (distMatrix[row*nOfColumns + col] == 0))//��ǰ�켣��û�б����䣬�Ҵ��ڶ�������ƥ��
					{
						/* prime zero */
						primeMatrix[row*nOfColumns + col] = true;
						/* �ҵ���ǰ���еĵ�һ��0 */
						size_t starCol = 0;
						for (; starCol < nOfColumns; starCol++)
						{
							if (starMatrix[row*nOfColumns + starCol])
							{
								break;
							}
						}
						if (starCol == nOfColumns) /* ��ǰ�켣û�ж�������������*/
						{
							/* move to step 4 */
							step4(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim, row, col);
							return;
						}
						else
						{
							coveredRows[row] = true;//��ǵ�ǰ�켣�����
							coveredColumns[starCol] = false;//��ǵ�ǰ����û�����
							zerosFound = true;//�ҵ����
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
	/* �� starMatrix ���п���*/
	for (size_t n = 0; n < nOfElements; n++)
	{
		newStarMatrix[n] = starMatrix[n];
	}
	/* star current zero */
	newStarMatrix[row*nOfColumns + col] = true;
	/* �ҵ���ǰ�е���ʼ0*/
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
		/* �Ա�ǵ�0ȡ����� */
		newStarMatrix[starRow*nOfColumns + starCol] = false;
		/* �ҵ���ǰ���еĵ�һ��0 */
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
		/* �ҵ���ǰ���еĵ�һ��0 */
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
	//�ͷ�������
	for (size_t n = 0; n < nOfRows; n++)
	{
		coveredRows[n] = false;
	}
	/* move to step 2a */
	step2a(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

// --------------------------------------------------------------------------
// ����һ����δ�����ǵ�Ԫ�����ҵ���Сֵ��Ȼ�����ЩԪ�ض���ȥ����һСֵ�������ǽ�����Ԫ�ؼ�����һ��Сֵ
// �ظ�����3��4
// ������Ԫ���е���Сֵʵ���������������������в��ɱ���Ŀ����� 
// ��һ�������������ӿ���������0�ĸ�����ʹ��������׷���
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step5(assignments_t& assignment, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	/* ����һ����δ�����ǵ�Ԫ�����ҵ���Сֵ h */
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
	//��û�����ǵ�Ԫ��ͳһ��ȥ��Сֵ������ʮ�ֽ��渲�ǵ�Ԫ�ؼ�����Сֵ

	//�������ǵ���������Ԫ�ؼ�����һ��Сֵ
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
	//��û�����ǵ����ϵ�����Ԫ�ض���ȥ����һСֵ
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
	std::cout << "����һ����δ�����ǵ�Ԫ�����ҵ���Сֵ��Ȼ�����ЩԪ�ض���ȥ����һСֵ�������ǽ�����Ԫ�ؼ�����һ��Сֵ" << std::endl;
	for (int i = 0; i < nOfRows; i++)
	{
		for (int j = 0; j < nOfColumns; j++)
			std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
		std::cout << std::endl;
	}
	step3(assignment, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}
