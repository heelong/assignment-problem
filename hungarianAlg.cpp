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
			std::cout << "�޹켣����ƥ��" << std::endl;
		}
		if (nOfColumns == 0)
		{
			std::cout << "�޶������ƥ��" << std::endl;
		}
		return 0;
	}
	distMatrix_t distMatrix_Square(nOfElements,0);
	std::vector<int> assignment_new(maxDim, -1);
	//���ɾ�����󣬲����Ԫ��
	track_t value;
	for (size_t row = 0; row < nOfRows; row++)
	{
		for (size_t col = 0; col < nOfColumns; col++)
		{
			value = distMatrixIn[row*nOfColumns + col];
			if (value < 0.0)
				std::cout << "��������������" << std::endl;
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
		assignment[row] = (assignment_new[row] >= nOfColumns) ? -1 : assignment_new[row];//�켣��Ӧ-���ϰ���
	}
	return cost;
}
// --------------------------------------------------------------------------
// ͨ���������㷨������С����ƥ��
// --------------------------------------------------------------------------
void AssignmentProblemSolver::assignmentoptimal(assignments_t& assignment, track_t& cost, const distMatrix_t& distMatrixIn, const size_t& nOfRows, const size_t& nOfColumns)
{

	//std::cout << "ԭʼ���󣬹켣���У���" << nOfRows << "    ����Ŀ�꣨�У���" << nOfColumns << std::endl;

	const size_t& nOfElements = nOfRows * nOfColumns;

	distMatrix_t distMatrix(nOfElements);
	// ָ�����һ��Ԫ��
	track_t* distMatrixEnd = distMatrix.data() + nOfElements;

	//���ɾ�����󣬲����Ԫ��
	track_t value;
	for (size_t row = 0; row < nOfElements; ++row)
	{
		value = distMatrixIn[row];
		if (value < 0.0)
			std::cout << "��������������" << std::endl;
		distMatrix[row] = value;
	}

	//for (int i = 0; i < nOfRows; i++)
	//{
	//	for (int j = 0; j < nOfColumns; j++)
	//		std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
	//	std::cout << std::endl;
	//}

	BoolVec coveredColumns(nOfColumns, 0);//��Ƕ��󱻷����ˣ�������0��
	BoolVec coveredRows(nOfRows, 0);//��ǹ켣�Ѿ��ж���ƥ����
	BoolVec starMatrix(nOfElements, 0);//��� ������켣�������
	BoolVec primeMatrix(nOfElements, 0);
	BoolVec newStarMatrix(nOfElements, 0); /* used in step4 */

	 markRows.resize(nOfRows, 0);//�����δ�������Ĺ켣
	 markColumns.resize(nOfColumns, 0);//������� **�ձ���ǹ�������** 0���ڵ���
	 CircleZero.resize(nOfElements, 0);//Ȧס��0Ԫ��
	 SlashZero.resize(nOfElements, 0); //������0Ԫ��

	//std::cout << " 1������С����ʱ���ҳ�ÿһ����ֵ��С��Ԫ�أ�Ȼ��Ѹ�������Ԫ�ض���ȥ��һ��Сֵ" << std::endl;
	//std::cout << " 2�����д�����ʱ���ҳ�ÿһ����ֵ��С��Ԫ�أ�Ȼ��Ѹ�������Ԫ�ض���ȥ��һ��Сֵ" << std::endl;
	if (nOfRows <= nOfColumns)//��С����
	{
		//std::cout << "����Ŀ����ڹ켣" << std::endl;
		track_t  minValue;
		track_t* distMatrixTemp;
		track_t* rowEnd;
		track_t* distMatrix_ptr = distMatrix.data();
		for (size_t row = 0; row < nOfRows; row++)
		{
			/* �ҳ�ÿһ����ֵ��С��Ԫ�� */
			distMatrixTemp = distMatrix_ptr + row*nOfColumns;
			rowEnd = distMatrixTemp + nOfColumns;//�õ��е�ĩβ
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
						//coveredRows[row] = true;//���������
						break;
					}
				}
			}
		}
	}
	else /* �н϶࣬����Ϊ��׼���б���*/
	{
		track_t* distMatrixTemp;
		track_t* columnEnd;
		track_t  minValue;
		track_t* distMatrix_ptr = distMatrix.data();
		int detal_size = (nOfRows - 1)*nOfColumns+1;//��ĩβ
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

		for (size_t col = 0; col < nOfColumns; col++)
		{
			for (size_t row = 0; row < nOfRows; row++)
			{
				if (distMatrix[row*nOfColumns + col] < 0.00001)
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
		for (size_t row = 0; row < nOfRows; row++)
		{
			coveredRows[row] = false;
		}
	}
	//std::cout << "��ȥÿ��/����Сֵ"<<std::endl;
	//for (int i = 0; i < nOfRows; i++)
	//{
	//	for (int j = 0; j < nOfColumns; j++)
	//		std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
	//	std::cout << std::endl;
	//}
	step2b(assignment,distMatrixIn ,distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, (nOfRows <= nOfColumns) ? nOfRows : nOfColumns);
	/* �������ɱ�*/
	computeassignmentcost(assignment, cost, distMatrixIn, nOfRows, nOfColumns);

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
// ������ɷ������2a
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step2a(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t&  nOfColumns, const size_t& minDim)
{
	uint *starMatrixTemp, *columnEnd;
	/* ���ÿһ���Ƿ��������0Ԫ��*/
	uint* starMatrix_ptr = starMatrix.data();
	int detal_size = (nOfRows - 1)*nOfColumns+1;//��ĩβ
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
	step2b(assignment, distMatrixIn, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}

// --------------------------------------------------------------------------
// �ж��Ƿ���ɷ��䣬û����ɷ�������벽��3���þ����ٵĺ��ߡ����߸��Ǿ����е�����0��
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step2b(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
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
		//std::cout << "��ɷ���" << std::endl;
	}
	else//���ڶ���û����ɷ���
	{
		step3(assignment,distMatrixIn, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
	}
}

// --------------------------------------------------------------------------
// �þ����ٵĺ��ߡ����߸��Ǿ����е�����0
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step3(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	//�þ����ٵĺ�������ȥ���Ǿ��󣬵������������ھ�����Сγ��ʱ�����ƥ��
	const size_t& nOfElements = nOfRows * nOfColumns;
	while (1)
	{
		for (size_t row = 0; row < nOfRows; row++)
		{
			markRows[row] = false;
			for (size_t col = 0; col < nOfColumns; col++)//����row_O�е�������,�ҵ�ԭʼ0
			{
				CircleZero[row*nOfColumns + col] = 0;
				SlashZero[row*nOfColumns + col] = 0;
				markColumns[col] = false;
			}
		}
		size_t CircleZeroNum = 0;          //Ȧס��0Ԫ�ص���Ŀ
		//�ҳ�����0Ԫ��
		step3b(assignment, distMatrixIn, distMatrix, starMatrix, CircleZero, SlashZero, CircleZeroNum, nOfRows, nOfColumns, minDim);

		if (CircleZeroNum == minDim)//�ж϶���0Ԫ�صĸ������������������Сά����ȣ���
		{
			step2a(assignment, distMatrixIn, distMatrix, CircleZero, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
			return;
		}

		//������δ����ǹ�����Ԫ�أ���Ȧ0�ĸ���<n
		//������ֱ�߸��ǵ�ǰ������Ԫ�أ������²����Ӷ�����Ԫ�صĸ�����

		//2����
		//2.1����Ȧ0�����ϴ򹳣�������һ����ɣ�
		for (size_t row = 0; row < nOfRows; row++)
		{
			markRows[row] = true;
			for (size_t col = 0; col < nOfColumns; col++)//����row_O�е�������,�ҵ�ԭʼ0
				if (CircleZero[row*nOfColumns + col])
				{
					markRows[row] = false;
					break;
				}
		}
		//2.2���򹳵����ϱ�������0�����д򹳣�������һ��
		//2.3���򹳵�������Ȧ0���д򹳣�������һ��
		bool zeorFrond = true;
		while (zeorFrond)
		{
			zeorFrond = false;
			//2.2���򹳵����ϱ�������0�����д�----------------
			for (size_t row = 0; row < nOfRows; row++)
			{
				if (markRows[row])
					for (size_t col = 0; col < nOfColumns; col++)//������
					{
						if (!markColumns[col] && (distMatrix[row*nOfColumns + col] < 0.00001) && SlashZero[row*nOfColumns + col])
						{
							markColumns[col] = true;
							zeorFrond = true;
						}
					}
			}
			//2.3���򹳵�������Ȧ0���д�-----------------------
			for (size_t col = 0; col < nOfColumns; col++)//������
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
		//4�������б�ǹ����к�δ��ǵ����ϻ��� 
		for (size_t col = 0; col < nOfColumns; col++)//������
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
		//��û�����ǵ�Ԫ��ͳһ��ȥ��Сֵ������ʮ�ֽ��渲�ǵ�Ԫ�ؼ�����Сֵ-----------------------

		//�������ǵ���������Ԫ�ؼ�����һ��Сֵ
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
		//��û�����ǵ����ϵ�����Ԫ�ض���ȥ����һСֵ
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
// �����Է��䣬�ж��Ƿ����n��������Ԫ��
// ���Զ�������Ԫ������ǣ�ȷ��������Ԫ�ء�
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step3b(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& CircleZero, BoolVec& SlashZero, size_t& CircleZeroNum, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
{
	size_t notMatchedZeroNum = 0;//û�б���ǵ�0����Ŀ
	for (size_t row_O = 0; row_O < nOfRows; row_O++)
	{
		for (size_t col = 0; col < nOfColumns; col++)//����row_O�е�������,�ҵ�ԭʼ0
			if ((distMatrix[row_O*nOfColumns + col] < 0.00001) && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
				notMatchedZeroNum++;
	}
	//�þ����ٵĺ�������ȥ���Ǿ��󣬵������������ھ�����Сά��ʱ�����ƥ��
	do
	{
		bool zeorFrond = true;
		//1.1�����м�顢�ҵ�ÿ�� >>>ֻ��һ��<<< û�б�Ȧס�򻮵���0Ԫ�أ�����Ȧס�����������������е�����δ��ǵ�0Ԫ�ػ���---------------------
		while (zeorFrond&&notMatchedZeroNum>0)
		{
			zeorFrond = false;
			for (size_t row_O = 0; row_O < nOfRows; row_O++)
			{
				size_t col_O = 0;
				size_t notMatchedZeroNumInRow_O = 0;
				for (size_t col = 0; col < nOfColumns; col++)//����row_O�е�������,�ҵ�ԭʼ0
					if ((distMatrix[row_O*nOfColumns + col] < 0.00001) && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
					{
						notMatchedZeroNumInRow_O++;
						col_O = col;
					}
				if (notMatchedZeroNumInRow_O == 0 || notMatchedZeroNumInRow_O > 1)
					continue;
				else
				{
					CircleZero[row_O*nOfColumns + col_O] = 1;//�ҵ�Ȧס��0
					notMatchedZeroNum--;
					starMatrix[row_O*nOfColumns + col_O] = 1;
					CircleZeroNum++;
					zeorFrond = true;
				}
				//����Ȧ�����Ԫ�����������е�����δ��ǵ���Ԫ���üǺ� / ��ȥ
				//for (size_t col = 0; col < nOfColumns; col++)////����row_O�е�������,��δ��ǵ�0����
				//{
				//	if (distMatrix[row_O*nOfColumns + col] < 0.00001 && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
				//		SlashZero[row_O*nOfColumns + col] = 1;
				//}
				for (size_t row = 0; row < nOfRows; row++)////����col_O�е�������,��δ��ǵ�0����
				{
					if (distMatrix[row*nOfColumns + col_O] < 0.00001 && (!CircleZero[row*nOfColumns + col_O] && !SlashZero[row*nOfColumns + col_O]))
					{
						SlashZero[row*nOfColumns + col_O] = 1;
						notMatchedZeroNum--;
					}
				}
			}
		}
		//1.2�����м�顢�ҵ�ÿ�� >>>ֻ��һ��<<< û�б�Ȧס�򻮵���0Ԫ�أ�����Ȧס�����������������е�����δ��ǵ�0Ԫ�ػ���---------------------
		zeorFrond = true;
		while (zeorFrond&&notMatchedZeroNum>0)
		{
			zeorFrond = false;
			for (size_t col_O = 0; col_O < nOfColumns; col_O++)
			{
				size_t row_O = 0;
				size_t notMatchedZeroNumInCol_O = 0;//û�б�ǵ�0
				for (size_t row = 0; row < nOfRows; row++)//����row_O�е�������,�ҵ�ԭʼ0
					if ((distMatrix[row*nOfColumns + col_O] < 0.00001) && (!CircleZero[row*nOfColumns + col_O] && !SlashZero[row*nOfColumns + col_O]))
					{
						notMatchedZeroNumInCol_O++;
						row_O = row;
					}
				if (notMatchedZeroNumInCol_O == 0 || notMatchedZeroNumInCol_O > 1)
					continue;
				else
				{
					CircleZero[row_O*nOfColumns + col_O] = 1;//�ҵ�Ȧס��0
					notMatchedZeroNum--;
					starMatrix[row_O*nOfColumns + col_O] = 1;
					CircleZeroNum++;
					zeorFrond = true;
				}
				//����Ȧ�����Ԫ�����������е�����δ��ǵ���Ԫ���üǺ� / ��ȥ
				for (size_t col = 0; col < nOfColumns; col++)////����row_O�е�������,��δ��ǵ�0����
				{
					if (distMatrix[row_O*nOfColumns + col] < 0.00001 && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
					{
						SlashZero[row_O*nOfColumns + col] = 1;
						notMatchedZeroNum--;
					}
				}
				//for (size_t row = 0; row < nOfRows; row++)////����col_O�е�������,��δ��ǵ�0����
				//{
				//	if (distMatrix[row*nOfColumns + col_O] < 0.00001 && (!CircleZero[row*nOfColumns + col_O] && !SlashZero[row*nOfColumns + col_O]))
				//		SlashZero[row*nOfColumns + col_O] = 1;
				//}
			}
		}
		//if ((CircleZeroNum < minDim))//ÿһ�о���Ȧ0���֣�Ȧ0�ĸ���ǡ�õ���n
		//	break;
		if (notMatchedZeroNum>0)
		{
			size_t row_O = -1;
			size_t col_O = -1;
			float minCost = (std::numeric_limits<float>::max)();//�������ҵ�����
			for (size_t row=0; row < nOfRows; row++)
			{
				for (size_t col=0; col < nOfColumns; col++)//����row_O�е�������,�ҵ�ԭʼ0
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
				//����Ȧ�����Ԫ�����������е�����δ��ǵ���Ԫ���üǺ� / ��ȥ
				for (size_t col = 0; col < nOfColumns; col++)////����row_O�е�������,��δ��ǵ�0����
				{
					if (distMatrix[row_O*nOfColumns + col] < 0.00001 && (!CircleZero[row_O*nOfColumns + col] && !SlashZero[row_O*nOfColumns + col]))
					{
						notMatchedZeroNum--;
						SlashZero[row_O*nOfColumns + col] = 1;
					}
				}
				for (size_t row = 0; row < nOfRows; row++)////����col_O�е�������,��δ��ǵ�0����
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
// ����һ����δ�����ǵ�Ԫ�����ҵ���Сֵ��Ȼ�����ЩԪ�ض���ȥ����һСֵ�������ǽ�����Ԫ�ؼ�����һ��Сֵ
// �ظ�����3��4
// ������Ԫ���е���Сֵʵ���������������������в��ɱ���Ŀ����� 
// ��һ�������������ӿ���������0�ĸ�����ʹ��������׷���
// --------------------------------------------------------------------------
void AssignmentProblemSolver::step4(assignments_t& assignment, const distMatrix_t& distMatrixIn, distMatrix_t& distMatrix, BoolVec& starMatrix, BoolVec& newStarMatrix, BoolVec& primeMatrix, BoolVec& coveredColumns, BoolVec& coveredRows, const size_t& nOfRows, const size_t& nOfColumns, const size_t& minDim)
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
	//��û�����ǵ�Ԫ��ͳһ��ȥ��Сֵ������ʮ�ֽ��渲�ǵ�Ԫ�ؼ�����Сֵ-----------------------

	//�������ǵ���������Ԫ�ؼ�����һ��Сֵ
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
	//��û�����ǵ����ϵ�����Ԫ�ض���ȥ����һСֵ
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
	//std::cout << "����һ����δ�����ǵ�Ԫ�����ҵ���Сֵ��Ȼ�����ЩԪ�ض���ȥ����һСֵ�������ǽ�����Ԫ�ؼ�����һ��Сֵ" << std::endl;
	//for (int i = 0; i < nOfRows; i++)
	//{
	//	for (int j = 0; j < nOfColumns; j++)
	//		std::cout << std::setw(10) << distMatrix[i*nOfColumns + j] << "   ";
	//	std::cout << std::endl;
	//}

	step3(assignment,distMatrixIn, distMatrix, starMatrix, newStarMatrix, primeMatrix, coveredColumns, coveredRows, nOfRows, nOfColumns, minDim);
}
