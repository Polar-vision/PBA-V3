
#include "IBA.h"
#include <stdio.h>

#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))

#elif defined(_WIN32)
	#include "stdafx.h"	
	
	Eigen::SparseMatrix<double> IBA::buildSparseMatrixCOO(int rows,int cols,
		const std::vector<int> & row_indices,const std::vector<int> & col_indices, const std::vector<double> & values){
		Eigen::SparseMatrix<double> mat(rows,cols);
		std::vector<Triplet<double>> triplets;
		for(size_t i=0;i<values.size();i++){
			triplets.push_back(Triplet<double>(row_indices[i],col_indices[i],values[i]));
		}
		mat.setFromTriplets(triplets.begin(),triplets.end());
		// mat.makeCompressed();
		return mat;
	}
	Eigen::SparseMatrix<double> IBA::buildSparseMatrixCOO_safe(
    int rows,
    int cols,
    const std::vector<int> &row_indices,
    const std::vector<int> &col_indices,
    const std::vector<double> &values)
	{
		// 创建稀疏矩阵
		Eigen::SparseMatrix<double> mat(rows, cols);

		// 提前估计非零元素数量，避免频繁扩容
		mat.reserve(values.size());

		// 插入元素
		for (size_t i = 0; i < values.size(); i++) {
			int r = row_indices[i];
			int c = col_indices[i];
			double v = values[i];

			// 安全检查
			if (r < 0 || r >= rows) continue;
			if (c < 0 || c >= cols) continue;
			if (!std::isfinite(v)) continue;  // 忽略 NaN 或 Inf

			mat.insert(r, c) = v;
		}

		// 压缩矩阵
		mat.makeCompressed();

		return mat;
	}
	//cite from SBA code in order to search sub matrix of S according to (i,j)
	void IBA::sba_crsm_alloc(struct sba_crsm* sm, int nr, int nc, int nnz)
	{
		int msz;
		sm->nr = nr;
		sm->nc = nc;
		sm->nnz = nnz;
		msz = 2 * nnz + nr + 1;
		sm->val = (int*)malloc(msz * sizeof(int));  /* required memory is allocated in a single step */
		if (!sm->val) {
			fprintf(stderr, "memory allocation request failed in sba_crsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
			exit(1);
		}
		sm->colidx = sm->val + nnz;
		sm->rowptr = sm->colidx + nnz;
	}

	void IBA::sba_crsm_free(struct sba_crsm* sm)
	{
		sm->nr = sm->nc = sm->nnz = -1;
		free(sm->val);
		sm->val = sm->colidx = sm->rowptr = NULL;
	}

	/* returns the index of the (i, j) element. No bounds checking! */
	int IBA::sba_crsm_elmidx(struct sba_crsm* sm, int i, int j)//����i��j��Ԫ�ص�������
	{
		int low, high, mid, diff;

		low = sm->rowptr[i];
		high = sm->rowptr[i + 1] - 1;

		/* binary search for finding the element at column j */
		while (low <= high)
		{
			mid = (low + high) >> 1; //(low+high)/2;
			diff = j - sm->colidx[mid];
			if (diff < 0)//��j�Ҳ�
				high = mid - 1;
			else if (diff > 0)//��j���
				low = mid + 1;
			else
				return mid;
		}

		return -1; /* not found */
	}
	double IBA::nrmL2xmy(double* const e, const double* const x, const double* const y, const int n)
	{
		const int blocksize = 8, bpwr = 3; /* 8=2^3 */
		int i;
		int j1, j2, j3, j4, j5, j6, j7;
		int blockn;
		double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

		/* n may not be divisible by blocksize,
		* go as near as we can first, then tidy up.
		*/
		blockn = (n >> bpwr) << bpwr; /* (n / blocksize) * blocksize; */

		/* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
		for (i = blockn - 1; i > 0; i -= blocksize) {
			e[i] = x[i] - y[i]; sum0 += e[i] * e[i];
			j1 = i - 1; e[j1] = x[j1] - y[j1]; sum1 += e[j1] * e[j1];
			j2 = i - 2; e[j2] = x[j2] - y[j2]; sum2 += e[j2] * e[j2];
			j3 = i - 3; e[j3] = x[j3] - y[j3]; sum3 += e[j3] * e[j3];
			j4 = i - 4; e[j4] = x[j4] - y[j4]; sum0 += e[j4] * e[j4];
			j5 = i - 5; e[j5] = x[j5] - y[j5]; sum1 += e[j5] * e[j5];
			j6 = i - 6; e[j6] = x[j6] - y[j6]; sum2 += e[j6] * e[j6];
			j7 = i - 7; e[j7] = x[j7] - y[j7]; sum3 += e[j7] * e[j7];

			//e[i] = x[i] - y[i]; sum0 += e[i];
			//j1 = i - 1; e[j1] = x[j1] - y[j1]; sum1 += abs(e[j1]);
			//j2 = i - 2; e[j2] = x[j2] - y[j2]; sum2 += abs(e[j2]);
			//j3 = i - 3; e[j3] = x[j3] - y[j3]; sum3 += abs(e[j3]);
			//j4 = i - 4; e[j4] = x[j4] - y[j4]; sum0 += abs(e[j4]);
			//j5 = i - 5; e[j5] = x[j5] - y[j5]; sum1 += abs(e[j5]);
			//j6 = i - 6; e[j6] = x[j6] - y[j6]; sum2 += abs(e[j6]);
			//j7 = i - 7; e[j7] = x[j7] - y[j7]; sum3 += abs(e[j7]);
		}

		i = blockn;
		if (i < n) {
			switch (n - i) {
			case 7: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
			case 6: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
			case 5: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
			case 4: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
			case 3: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
			case 2: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;
			case 1: e[i] = x[i] - y[i]; sum0 += e[i] * e[i]; ++i;

			//case 7: e[i] = x[i] - y[i]; sum0 += abs(e[i]); ++i;
			//case 6: e[i] = x[i] - y[i]; sum0 += abs(e[i]); ++i;
			//case 5: e[i] = x[i] - y[i]; sum0 += abs(e[i]); ++i;
			//case 4: e[i] = x[i] - y[i]; sum0 += abs(e[i]); ++i;
			//case 3: e[i] = x[i] - y[i]; sum0 += abs(e[i]); ++i;
			//case 2: e[i] = x[i] - y[i]; sum0 += abs(e[i]); ++i;
			//case 1: e[i] = x[i] - y[i]; sum0 += abs(e[i]); ++i;
			}
		}

		return sum0 + sum1 + sum2 + sum3;
	}
	void IBA::ba_saveTriangulatedxyz(const char* sz3Dpt, double* p)
	{
		static int i = 0;
		double x, y, z;
		FILE* fp = nullptr;
		//save features xyz
		if (sz3Dpt != NULL)
		{
			fopen_s(&fp, sz3Dpt, "w");
			for (i = 0; i < m_n3Dpts; i++)
			{
				x = *(p + m_ncams * 6 + i * 3);
				y = *(p + m_ncams * 6 + i * 3 + 1);
				z = *(p + m_ncams * 6 + i * 3 + 2);
				//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1, x, y, z);
				fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", x, y, z);
				//fprintf(fp, "%0.5lf\n", parallax * 180 / PI);
			}
			fclose(fp);
		}
		free(p);
	}

	int IBA::ba_ConstructSmask(sba_crsm& Sidxij, sba_crsm& Uidxij)//�ֱ���S�����U�����ϡ������
	{
		int i, j, k, ii, jj;
		int nuis, m = m_ncams;//nuis��U����ķ���Ԫ�ظ���
		//compute total smask
		for (i = 0; i < m; i++) for (j = 0; j < m; j++)
		{
			//��δ���δִ��
			if (m_umask[i * m + j] == 1 && m_smask[i * m + j] == 0)//m_umask[i*m+j] == 1 ��ʾ U �����ڸ�λ���Ƿ���ģ���� S ��û��ǣ��Ͱ�������ϡ�
			{
				m_smask[i * m + j] = 1;//S �����ϡ������ (Sidxij �Ľṹ),m_smask[i*m+j] ��һ�� m��m ά�Ķ�ά���飨�� 1D �洢������ʾ S ������ (i,j) λ���Ƿ��з���Ԫ�ء�
				m_nS += 1;
			}
		}
		//���� S ����Ĵ洢�ռ�,Sidxij��S ����Ľṹ��,m, mΪ�����С
		sba_crsm_alloc(&Sidxij, m, m, m_nS);//m_nS��S����ķ���Ԫ�ظ���
		for (i = k = 0; i < m; ++i)
		{
			Sidxij.rowptr[i] = k;// ��¼ S ����� i �е���ʼ����
			ii = i * m;
			for (j = 0; j < m; ++j)
				if (m_smask[ii + j])// �����λ���Ƿ���Ԫ��
				{
					Sidxij.val[k] = k;// ������ֵ������ k
					Sidxij.colidx[k++] = j; // ��¼�÷���Ԫ�ص�������
				}
		}
		Sidxij.rowptr[m] = m_nS;// ���һ��ָ��ĩβ

		for (i = nuis = 0, jj = m * m; i < jj; ++i)
			nuis += (m_umask[i] != 0);//U �����ϡ������ (Uidxij �Ľṹ)

		sba_crsm_alloc(&Uidxij, m, m, nuis);
		for (i = k = 0; i < m; ++i)
		{
			Uidxij.rowptr[i] = k;
			ii = i * m;
			for (j = 0; j < m; ++j)
				if (m_umask[ii + j])
				{
					Uidxij.val[k] = k;
					Uidxij.colidx[k++] = j;
				}
		}
		Uidxij.rowptr[m] = nuis;

		//����
		// FILE* fp;
		// const char* fn = "E:/zuo/projects/PBA_U.txt";
		// errno_t err = fopen_s(&fp, fn, "w");
		// if (err == 0 && fp != NULL) {  // ˫�ؼ��
		// 	for (i = 0; i < m; ++i)
		// 	{
		// 		for (j = 0; j < m; ++j)
		// 		{
		// 			fprintf(fp, "%d", m_umask[i * m + j]);
		// 		}
		// 		fprintf(fp, "\n");
		// 	}
		// 	fclose(fp);
		// }
		// else {
		// 	printf("��ʧ�ܣ�������: %d\n", err);
		// }

		//fn = "E:/zuo/projects/PBA-Colmap/data/colmap/odm_zoo/U-aff.txt";
		//err = fopen_s(&fp, fn, "w");
		//if (err == 0 && fp != NULL) {  // ˫�ؼ��
		//	for (i = 0; i <= m; ++i)
		//		fprintf(fp, "%d ", Uidxij.rowptr[i]);
		//	fprintf(fp, "\n\n\n");
		//	for (i = 0; i < nuis; ++i)
		//		fprintf(fp, "%d ", Uidxij.val[i]);
		//	fprintf(fp, "\n\n\n");
		//	for (i = 0; i < nuis; ++i)
		//		fprintf(fp, "%d ", Uidxij.colidx[i]);
		//	fprintf(fp, "\n\n\n");
		//	fclose(fp);
		//}
		//else {
		//	printf("��ʧ�ܣ�������: %d\n", err);
		//}


		////����
		//fn = "E:/zuo/projects/PBA-Colmap/data/colmap/odm_zoo/S.txt";
		//err = fopen_s(&fp, fn, "w");
		//if (err == 0 && fp != NULL) {  // ˫�ؼ��
		//	for (i = 0; i < m; ++i)
		//	{
		//		for (j = 0; j < m; ++j)
		//		{
		//			fprintf(fp, "%d", m_smask[i * m + j]);
		//		}
		//		fprintf(fp, "\n");
		//	}
		//	fclose(fp);
		//}
		//else {
		//	printf("��ʧ�ܣ�������: %d\n", err);
		//}

		//fn = "E:/zuo/projects/PBA-Colmap/data/colmap/odm_zoo/S-aff.txt";
		//err = fopen_s(&fp, fn, "w");
		//if (err == 0 && fp != NULL) {  // ˫�ؼ��
		//	for (i = 0; i <= m; ++i)
		//		fprintf(fp, "%d ", Sidxij.rowptr[i]);
		//	fprintf(fp, "\n\n\n");
		//	for (i = 0; i < m_nS; ++i)
		//		fprintf(fp, "%d ", Sidxij.val[i]);
		//	fprintf(fp, "\n\n\n");
		//	for (i = 0; i < m_nS; ++i)
		//		fprintf(fp, "%d ", Sidxij.colidx[i]);
		//	fprintf(fp, "\n\n\n");
		//	fclose(fp);
		//}
		//else {
		//	printf("��ʧ�ܣ�������: %d\n", err);
		//}


		return nuis;

	}
	//IV=(V+mu*I)^-1
	void IBA::ba_inverseVLM(double* V, double* IV, sba_crsm& Uidxij, double mu)
	{
		int i, j;
		int m = m_ncams, n = m_n3Dpts;
		int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
		double* ptr1, * ptr2;
		Matrix3d MatInv;

		//IV save inverse V matrix, V must unchange for the next step
		memcpy(IV, V, n * Vsz * sizeof(double));//IV是目标地址，V是源地址，n*Vsz是要拷贝的字节数
		for (i = 0; i < n; ++i)
		{
			ptr1 = V + i * Vsz; //original V
			ptr2 = IV + i * Vsz;//damped V

			for (j = 0; j < pnp; ++j){
				ptr2[j * pnp + j] += mu;

				// double diag_val = ptr1[j * pnp + j];
            	// ptr2[j * pnp + j] = diag_val * (1.0 + mu);
			}
				
			Eigen::Matrix3d matV(ptr2);
			MatInv = matV.inverse();
			ptr2[0] = MatInv(0, 0);
			ptr2[4] = MatInv(1, 1);
			ptr2[8] = MatInv(2, 2);
			ptr2[1] = ptr2[3] = MatInv(0, 1);
			ptr2[2] = ptr2[6] = MatInv(0, 2);
			ptr2[5] = ptr2[7] = MatInv(1, 2);
		}
	}
	void IBA::ba_inverseVGN(double* U, double* V, sba_crsm& Uidxij)
	{
		int i;
		int m = m_ncams, n = m_n3Dpts;
		int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
		double* ptr1;
		Matrix3d matV, MatInv;

		//compute V inverse matrix using Eigen that has better performance than Lapack
		for (i = 0; i < n; ++i)
		{
			ptr1 = V + i * Vsz; // set ptr1 to point to V_i
			matV << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];
			MatInv = matV.inverse();
			ptr1[0] = MatInv(0, 0);
			ptr1[4] = MatInv(1, 1);
			ptr1[8] = MatInv(2, 2);
			ptr1[1] = ptr1[3] = MatInv(0, 1);
			ptr1[2] = ptr1[6] = MatInv(0, 2);
			ptr1[5] = ptr1[7] = MatInv(1, 2);
		}
	}

	double IBA::ba_computeInitialmu(double* U, double* V, sba_crsm& Uidxij, double tau, int nvars)
	{
		int i, j;
		int pos, m = m_ncams, n = m_n3Dpts, cnp = 6, pnp = 3, Usz = 36, Vsz = 9;
		double tmp = 0;
		double* ptr1, * ptr2;
		double mu;

		double* diagUV = (double*)malloc(nvars * sizeof(double));

		double* diagU = diagUV;
		double* diagV = diagUV + m * cnp;

		for (j = 0; j < m; ++j)
		{
			pos = sba_crsm_elmidx(&Uidxij, j, j);
			ptr1 = U + pos * Usz;
			ptr2 = diagU + j * cnp;
			for (i = 0; i < cnp; ++i)
				ptr2[i] = ptr1[i * cnp + i];
		}
		for (i = 0; i < n; ++i)
		{
			ptr1 = V + i * Vsz; // set ptr1 to point to V_i
			ptr2 = diagV + i * pnp; // set ptr2 to point to diagV_i
			for (j = 0; j < pnp; ++j)
				ptr2[j] = ptr1[j * pnp + j];
		}

		/* find max diagonal element */
		for (i = 0, tmp = DBL_MIN; i < m * cnp; ++i)
			if (diagUV[i] > tmp)
				tmp = diagUV[i];
		for (i = m * cnp; i < nvars; ++i) /* tmp is not re-initialized! */
			if (diagUV[i] > tmp)
				tmp = diagUV[i];

		mu = m_Tau * tmp;

		free(diagUV);
		diagU = diagV = NULL;

		return mu;
	}


	void IBA::ba_solveFeatures(double* W, double* IV, double* ea, double* eb, double* dpa, double* dpb)
	{
		int i, j, ii, jj, pos, numfea;
		int nP1, cnp = 6, pnp = 3;
		double* ptr1, * ptr2, * ptr3, * ptr4, * ptr5;
		double sum, eb2[6];
		pos = 0;
		for (i = 0; i < m_n3Dpts; i++)
		{
			ptr1 = eb + i * 3;//ָ����ά���������
			ptr2 = IV + i * 3 * 3;//ָ��V����
			ptr5 = dpb + i * 3;//ָ����ά��δ֪����
			memset(eb2, 0, sizeof(double) * cnp);
			numfea = m_archor[i * 3];

			for (j = 0; j < numfea; j++)
			{
				nP1 = m_photo[pos];//��ͼ���
				ptr3 = W + pos * cnp * 3;//ָ��W����
				ptr4 = dpa + nP1 * cnp;//ָ�����δ֪����
				//Wta
				for (ii = 0; ii < pnp; ++ii)
				{
					for (jj = 0, sum = 0; jj < cnp; ++jj)
						sum += ptr3[jj * 3 + ii] * ptr4[jj];
					eb2[ii] += sum;
				}
				pos++;
			}

			//V*(eb-Wta)
			for (ii = 0; ii < pnp; ++ii)
			{
				for (jj = 0, sum = 0; jj < pnp; jj++)
					sum += ptr2[ii * 3 + jj] * (ptr1[jj] - eb2[jj]);
				ptr5[ii] = sum;
			}
		}
	}
	//��Ч���GN�г��ֵ�ϡ������ϵͳ
	//ʹ��Cholmod����о���Cholesky���ֽ�����
	//Ap����ָ�룬Aii������������CSC��ʽ
	bool IBA::ba_solveCholmodGN(int* Ap, int* Aii, bool init, bool ordering)
	{
		int i, j;
		int m = m_ncams;
		VectorXi scalarPermutation, blockPermutation;

		ordering = true;
		if (!init)//�״γ�ʼ��
		{
			if (!ordering)//����Ҫ�ֶ�����
			{
				m_cS.nmethods = 1;
				m_cS.method[0].ordering = CHOLMOD_AMD; //���򷽷�����ΪApproximately Minimum Degree��������С������
				m_cholFactorS = cholmod_analyze(m_cholSparseS, &m_cS); // symbolic factorization�����ɷֽ�������ڲ����ݽṹ
			}
			else
			{
				// get the ordering for the block matrix
				if (blockPermutation.size() == 0)
					blockPermutation.resize(m_ncams - 1);//���ڴ洢ÿ���������˳��

				// prepare AMD call via CHOLMOD
				cholmod_sparse auxCholmodSparse;//����cholmod_sparse���͵ĸ�������auxcholmodSparse����ϡ�����Ԫ���ݣ���Ap��Aii��ӳ�䵽Cholmod��ʽ
				auxCholmodSparse.nzmax = m_nS;
				auxCholmodSparse.nrow = auxCholmodSparse.ncol = m - 1;//������1֡
				auxCholmodSparse.p = Ap;
				auxCholmodSparse.i = Aii;
				auxCholmodSparse.nz = 0;
				auxCholmodSparse.x = 0;
				auxCholmodSparse.z = 0;
				auxCholmodSparse.stype = 1;
				auxCholmodSparse.xtype = CHOLMOD_PATTERN;
				auxCholmodSparse.itype = CHOLMOD_INT;
				auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
				auxCholmodSparse.sorted = 1;
				auxCholmodSparse.packed = 1;
				//AMD�������ڼ���ϡ�����ֽ��е����
				int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);//ִ��AMD���򣬲�������洢��blockPermutation
				if (!amdStatus) {
					return false;
				}

				// blow up the permutation to the scalar matrix
				//����������չΪ��������
				if (scalarPermutation.size() == 0)
					scalarPermutation.resize(m_cholSparseS->ncol);
				size_t scalarIdx = 0;

				int a = 0;
				for (i = 0; i < m_ncams - 1; ++i)
				{
					const int& pp = blockPermutation(i);
					int base = (pp == 0) ? 0 : pp * 6 - 1;
					int nCols = (pp == 0) ? 5 : 6;
					for (j = 0; j < nCols; ++j)
						scalarPermutation(scalarIdx++) = base++;
				}
				assert(scalarIdx == m_cholSparseS->ncol);//m_cholSparseS->ncol=90*6-7=533

				// apply the ordering
				m_cS.nmethods = 1;
				m_cS.method[0].ordering = CHOLMOD_GIVEN;//����CHOLMODʹ�ø���������˳��
				m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);//����cholmod_analyze_p���з��ŷ���
			}
			init = true;//��ɳ�ʼ��
		}

		//Cholmod package for solving sparse linear equation              
		cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS); //��ϡ����������ֵ�ֽ�
		m_cholSparseR = cholmod_solve(CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS);//������Է�����

		return true;
	}
	//learning this skill from G2O
	bool IBA::ba_solveCholmodLM(int* Ap, int* Aii, bool init, bool ordering)
	{
		int i, j;
		VectorXi scalarPermutation, blockPermutation;

		ordering = true;
		if (!init)
		{
			if (!ordering)
				m_cholFactorS = cholmod_analyze(m_cholSparseS, &m_cS); // symbolic factorization
			else
			{
				// get the ordering for the block matrix
				if (blockPermutation.size() == 0)
					blockPermutation.resize(m_ncams);
				if (blockPermutation.size() < m_ncams) // double space if resizing
					blockPermutation.resize(2 * m_ncams);

				// prepare AMD call via CHOLMOD
				cholmod_sparse auxCholmodSparse;
				auxCholmodSparse.nzmax = m_nS;
				auxCholmodSparse.nrow = auxCholmodSparse.ncol = m_ncams;
				auxCholmodSparse.p = Ap;
				auxCholmodSparse.i = Aii;
				auxCholmodSparse.nz = 0;
				auxCholmodSparse.x = 0;
				auxCholmodSparse.z = 0;
				auxCholmodSparse.stype = 1;
				auxCholmodSparse.xtype = CHOLMOD_PATTERN;
				auxCholmodSparse.itype = CHOLMOD_INT;
				auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
				auxCholmodSparse.sorted = 1;
				auxCholmodSparse.packed = 1;
				int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);
				if (!amdStatus)
					return false;


				// blow up the permutation to the scalar matrix
				if (scalarPermutation.size() == 0)
					scalarPermutation.resize(m_cholSparseS->ncol);
				if (scalarPermutation.size() < (int)m_cholSparseS->ncol)
					scalarPermutation.resize(2 * m_cholSparseS->ncol);
				size_t scalarIdx = 0;

				for (i = 0; i < m_ncams; ++i)
				{
					const int& pp = blockPermutation(i);
					int base = (pp == 0) ? 0 : pp * 6;
					int nCols = 6;

					for (j = 0; j < nCols; ++j)
						scalarPermutation(scalarIdx++) = base++;

				}
				assert(scalarIdx == m_cholSparseS->ncol);

				// apply the ordering
				m_cS.nmethods = 1;
				m_cS.method[0].ordering = CHOLMOD_GIVEN;
				m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);
			}
		}

		cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS);
		m_cholSparseR = cholmod_solve(CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS);

		return true;
	}
	void IBA::ba_constructCSSLM(int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init)
	{
		int ii, jj, jjj, k;
		int pos1, m = m_ncams;
		//Copy S matrix and E matrix to specific format structure for Cholmod 
		double* ptr5;
		int nZ = 0;

		Sx = (double*)m_cholSparseS->x;
		if (!init)
		{
			for (ii = 0; ii < m; ii++)  //colum
			{
				for (k = 0; k < 6; k++)
				{
					*Sp = nZ;
					for (jj = 0; jj <= ii; jj++)	//row
					{
						if (m_smask[jj * m + ii] == 1)
						{
							pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
							ptr5 = S + pos1 * 36;

							if (ii == jj)
							{
								for (jjj = 0; jjj <= k; jjj++)
								{
									*Si++ = jj * 6 + jjj;
									*Sx++ = ptr5[jjj * 6 + k];
									nZ++;
								}
							}
							else
							{
								for (jjj = 0; jjj < 6; jjj++)
								{
									*Si++ = jj * 6 + jjj;
									*Sx++ = ptr5[jjj * 6 + k];
									nZ++;
								}
							}
						}
					}
					Sp++;
				}
			}
			*Sp = nZ;
		}
		else
		{
			for (ii = 0; ii < m; ii++)  //colum
			{
				for (k = 0; k < 6; k++)
				{
					for (jj = 0; jj <= ii; jj++)	//row
					{
						if (m_smask[jj * m + ii] == 1)
						{
							pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
							ptr5 = S + pos1 * 36;

							if (ii == jj)
							{
								for (jjj = 0; jjj <= k; jjj++)
									*Sx++ = ptr5[jjj * 6 + k];
							}
							else
							{
								for (jjj = 0; jjj < 6; jjj++)
									*Sx++ = ptr5[jjj * 6 + k];
							}
						}
					}
				}
			}
		}
	}
	void IBA::ba_constructCSSGN(int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init, int nft)
	{
		int ii, jj, jjj, k;
		int pos1, m = m_ncams;
		//Copy S matrix and E matrix to specific format structure for Cholmod 
		double* ptr5;
		int nZ = 0;
		Sx = (double*)m_cholSparseS->x;//�洢ϡ����������з���Ԫ�ص�ʵ��ֵ
		//printf("\n\n\n\n");
		if (!init)
		{
			for (ii = 1; ii < m; ii++)  //column����0����ͼ�ǲο�֡��ͨ�����̶�������Ϊ�Ż�����
			{
				for (k = 0; k < 6; k++)
				{
					*Sp = nZ;//��ʾϡ�������ÿһ�е���ʼλ�ã�����Ԫ�ص�������,0,1,3,6,
					//printf("%d ", *Sp);
					if ((ii * 6 + k) == (9 + nft))//k=3����Xc��4����Yc��5����Zc��nft=0�����̶�Xc��1�����̶�Yc��2�����̶�Zc
						continue;//�̶���2֡��Xc��Ϊ�߶�Լ��

					for (jj = 1; jj <= ii; jj++)	//row
					{
						if ((m_smask[jj * m + ii] == 1))
						{
							pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
							//printf("%d ", pos1);
							ptr5 = S + pos1 * 36;//ָ��S�����Ӧ��jj,ii�����Ӿ�������

							if (ii == jj)//�Խ��߿�
							{
								for (jjj = 0; jjj <= k; jjj++)
								{
									if ((jj * 6 + jjj) != (9 + nft))//����jj=1��jjj=3��������2�������Xc�̶����⣬
									{
										if (jj * 6 + jjj < 9 + nft)//��2���������̬��
											*Si++ = jj * 6 + jjj - 6;
										else//��2�������Yc��Zc
											*Si++ = jj * 6 + jjj - 7;

										*Sx++ = ptr5[jjj * 6 + k];
										//printf("%d ", ptr5[jjj * 6 + k]);
										nZ++;
									}
								}
							}
							else
							{
								for (jjj = 0; jjj < 6; jjj++)
								{
									if ((jj * 6 + jjj) != (9 + nft))
									{
										if (jj * 6 + jjj < 9 + nft)
											*Si++ = jj * 6 + jjj - 6;
										else
											*Si++ = jj * 6 + jjj - 7;

										*Sx++ = ptr5[jjj * 6 + k];
										nZ++;
									}
								}
							}
						}
					}
					Sp++;
				}
			}
			*Sp = nZ;
		}
		else
		{
			for (ii = 1; ii < m; ii++)  //column
			{
				for (k = 0; k < 6; k++)
				{
					if ((ii * 6 + k) == (9 + nft))
						continue;

					for (jj = 1; jj <= ii; jj++)	//row
					{
						if ((m_smask[jj * m + ii] == 1))
						{
							pos1 = sba_crsm_elmidx(&Sidxij, jj, ii);
							ptr5 = S + pos1 * 36;

							if (ii == jj)
							{
								for (jjj = 0; jjj <= k; jjj++)
								{
									if ((jj * 6 + jjj) != (9 + nft))
										*Sx++ = ptr5[jjj * 6 + k];
								}
							}
							else
							{
								for (jjj = 0; jjj < 6; jjj++)
								{
									if ((jj * 6 + jjj) != (9 + nft))
										*Sx++ = ptr5[jjj * 6 + k];
								}
							}
						}
					}
				}
			}
		}
	}
	void IBA::ba_constructAuxCSSLM(int* Ap, int* Aii)
	{
		int* Cp = Ap;
		int* Ci = Aii;
		int ii, jj;
		int m = m_ncams, nZ = 0;
		for (ii = 0; ii < m; ii++)
		{
			*Cp = nZ;
			for (jj = 0; jj <= ii; jj++)
			{
				if (m_smask[jj * m + ii] == 1)
				{
					*Ci++ = jj;
					nZ++;
				}
			}
			Cp++;
		}
		*Cp = nZ;
	}
	//����ϡ�����������ṹ����1��1�п�ʼ���洢ÿ�е���ʼ����,�洢����Ԫ�ص�������
	void IBA::ba_constructAuxCSSGN(int* Ap, int* Aii)
	{
		//��0��0�п�ʼ����Ϊ�̶���һ֡
		int* Cp = Ap;//�洢ÿ�е���ʼ����
		int* Ci = Aii;//�洢����Ԫ�ص�������
		int ii, jj;
		int m = m_ncams, nZ = 0;
		//printf("\n\n\n");
		for (ii = 1; ii < m; ii++) //��
		{
			*Cp = nZ;
			for (jj = 1; jj <= ii; jj++)//��
			{
				if (m_smask[jj * m + ii] == 1)//���б���S���󣨶Գƾ��󣩵�������
				{
					*Ci++ = jj - 1;//���ܴ�0�п�ʼ
					//printf("%d ", jj-1);
					nZ++;
				}
			}
			Cp++;
		}
		*Cp = nZ;

		//����
		//Cp = Ap;//��ָ���Ƶ��ļ�ͷ
		//Ci = Aii;
		//printf("\n\n\n");
		//for (int i = 0; i < m; ++i)
		//	printf("%d ", *Cp++);
		//printf("\n\n\n");
		//for (int i = 0; i < m_nS; ++i)
		//	printf("%d ", *Ci++);
		//printf("\n\n\n");
	}
	std::pair<double, double> svd_max_min_singular_values_dense(MatrixXd denseS) {
		// 使用JacobiSVD计算奇异值
		JacobiSVD<MatrixXd> svd(denseS, ComputeThinU | ComputeThinV) ;
		
		if (svd.singularValues().size() == 0) {
		    std::cout << "SVD computation failed!" << std::endl;
		    return {-1.0, -1.0};
		}
		
		VectorXd singularValues = svd.singularValues();
		double sigma_max = singularValues[0];
		double sigma_min = singularValues[singularValues.size() - 1];
		
		return {sigma_max, sigma_min};
	}
	// VectorXd svd_singular_values_dense(const MatrixXd& denseS) {
	// 	JacobiSVD<MatrixXd> svd(denseS, Eigen::ComputeThinU | Eigen::ComputeThinV);

	// 	if (svd.singularValues().size() == 0) {
	// 		std::cerr << "SVD computation failed!" << std::endl;
	// 		return Eigen::VectorXd();
	// 	}

	// 	return svd.singularValues();
	// }
	VectorXd svd_singular_values_dense(const MatrixXd& denseS) {
		// BDCSVD 使用分治算法，对于大规模矩阵比 JacobiSVD 快很多
		// BDCSVD<MatrixXd> svd(denseS, Eigen::ComputeThinU | Eigen::ComputeThinV);
		BDCSVD<MatrixXd> svd(denseS);//只计算奇异值
		
		if (svd.singularValues().size() == 0) {
			std::cerr << "SVD computation failed!" << std::endl;
			return Eigen::VectorXd();
		}
		
		return svd.singularValues();
	}
	// VectorXd svd_singular_values_dense(const MatrixXd& denseS) {
	// 	int rows = denseS.rows();
	// 	int cols = denseS.cols();
	// 	int size = std::min(rows, cols);
		
	// 	if (size > 1000) {
	// 		// 对于大规模矩阵，只计算前 k 个最大奇异值
	// 		int k = std::min(50, size);  // 前50个最大奇异值
			
	// 		// 对于方阵，可以使用特征值分解 A^T A
	// 		MatrixXd ATA = denseS.transpose() * denseS;
			
	// 		// 使用 Spectra 计算最大的 k 个特征值
	// 		Spectra::DenseSymMatProd<double> op(ATA);
	// 		Spectra::SymEigsSolver<double> solver(op, k, std::min(2*k, size));
			
	// 		solver.init();
	// 		int nconv = solver.compute(Spectra::SortRule::LargestMagn);
			
	// 		if (solver.info() == Spectra::CompInfo::Successful) {
	// 			VectorXd eigenvalues = solver.eigenvalues();
	// 			VectorXd singular_values = eigenvalues.array().sqrt();
	// 			return singular_values;
	// 		} else {
	// 			std::cerr << "Spectra SVD computation failed!" << std::endl;
	// 			return Eigen::VectorXd();
	// 		}
	// 	} else {
	// 		// 小矩阵使用 BDCSVD
	// 		BDCSVD<MatrixXd> svd(denseS);
	// 		return svd.singularValues();
	// 	}
	// }
	//std::tuple<double,double,double>
	VectorXd IBA::convertStoDenseMatrix(double* S,sba_crsm& Sidxij,
		const char * name1,const char * name2){
		int blockSize = 6;
		int nCams = m_ncams;
		int fullSize = blockSize * nCams;
		Eigen::MatrixXd S_full = Eigen::MatrixXd::Zero(fullSize,fullSize);
		// Eigen::SparseMatrix<double> S_sparse(fullSize, fullSize);
		// std::vector<Eigen::Triplet<double>> triplets;
		for (int ci = 0; ci < nCams; ci++) {
			for (int cj = ci; cj < nCams; cj++) {
				int pos = sba_crsm_elmidx(&Sidxij, ci, cj);
				if (pos < 0) continue; 
				double* ppSpa = S + pos * blockSize * blockSize;

				if(cj==ci){
					for (int ii = 0; ii < blockSize; ii++) {
						for (int jj = ii; jj < blockSize; jj++) {
							double val = ppSpa[ii * blockSize + jj];
							int i=ci * blockSize + ii;
							int j=cj * blockSize + jj;
							S_full(i, j) = val;
							// triplets.emplace_back(i, j, val);
							if(i!=j){
								S_full(j,i)=val;
								// triplets.emplace_back(j, i, val);
							}

						}
					}
				}else{
					for (int ii = 0; ii < blockSize; ii++) {
						for (int jj = 0; jj < blockSize; jj++) {
							double val = ppSpa[ii * blockSize + jj];
							int i=ci * blockSize + ii;
							int j=cj * blockSize + jj;
							S_full(i, j) = val;
							// triplets.emplace_back(i, j, val);
							if(i!=j){
								S_full(j,i)=val;
								// triplets.emplace_back(j, i, val);
							}
						}
					}
				}

			}
		}
		
		// std::pair<double, double> sigma;
		// sigma = svd_max_min_singular_values_dense(S_full);
		VectorXd sigma_spectral = svd_singular_values_dense(S_full);
	
		// double sigma_max = sigma.first;
		// double sigma_min = sigma.second;

		// std::cout << "sigma_max: " << sigma_max << std::endl;
		// std::cout << "sigma_min: " << sigma_min << std::endl;
		// double eps = 1e-12;
		// if (sigma_min < eps)
		// {
		// 	std::cout << "Matrix is near singular!" << std::endl;
		// }
		// double cond = sigma_max / sigma_min;
		// std::cout << "Condition number: " << cond << std::endl;
		
		// for (int r = 0; r < fullSize; r++) {
		// 	for (int c = r + 1; c < fullSize; c++) {
		// 		if (S_full(r, c) != 0)
		// 			S_full(c, r) = S_full(r, c);
		// 	}
		// }
		// S_sparse.setFromTriplets(triplets.begin(), triplets.end());

		// Eigen::JacobiSVD<Eigen::MatrixXd> svd(S_full);
		// double lambda_max, lambda_min;
		// lambda_max = svd.singularValues()(0);
		// lambda_min = svd.singularValues().tail(1)(0);
		// double cond1 = log10(lambda_max / lambda_min);
		// printf("%s %f %f %f\n\n","lambda_max lambda_min cond = ",lambda_max,lambda_min,cond1);

		// double lambda_max, lambda_min;
		// // SparseSymMatProd<double> op(S_sparse);
		// DenseSymMatProd<double> op(S_full);
		// SymEigsSolver<decltype(op)> eigs_max(op, 1, 10);
		// eigs_max.init();
		// eigs_max.compute(SortRule::LargestAlge);
		// if (eigs_max.info() == Spectra::CompInfo::Successful) {
		// 	lambda_max = eigs_max.eigenvalues()[0];
		// 	printf("lambda_max = %f\n", lambda_max);
		// } else {
		// 	printf("Largest eigenvalue computation did not converge!\n");
		// }
		
		// SymEigsSolver<decltype(op)> eigs_min(op, 1, 10);
		// eigs_min.init();
		// eigs_min.compute(SortRule::SmallestAlge);
		// if (eigs_min.info() == Spectra::CompInfo::Successful) {
		// 	lambda_min = eigs_min.eigenvalues()[0];
		// 	printf("lambda_min = %f\n", lambda_min);
		// } else {
		// 	printf("Smallest eigenvalue computation did not converge!\n");
		// }
		// double cond = sqrt(lambda_max / lambda_min);
		// printf("%s %f\n","cond = ",cond);

		// FILE* fp1 = fopen(name1, "wb");
		// fprintf(fp1, "P5\n%d %d\n255\n", fullSize, fullSize);
		// for(int i=0;i<fullSize;i++){
		// 	for(int j=0;j<fullSize;j++){
		// 		unsigned char pixel = S_full(i, j) ? 255 : 0;
		// 		fwrite(&pixel, 1, 1, fp1);
		// 	}
		// }
		// fclose(fp1);

		// FILE* fp2;
		// errno_t err2 = fopen_s(&fp2, name2, "w");
		// for(int i=0;i<fullSize;i++){
		// 	for(int j=0;j<fullSize;j++){
		// 		fprintf(fp2,"%f ",S_full(i,j));
		// 	}
		// 	fprintf(fp2,"\n");
		// }
		// fclose(fp2);

		return sigma_spectral;
		// return {sigma_max, sigma_min, cond};
	}
	void IBA::convertUtoDenseMatrix(double* U,sba_crsm* Uidxij,const char * name1,const char * name2){
		int blockSize = 6;
		int nCams = m_ncams;
		int fullSize = blockSize * nCams;
		Eigen::MatrixXd U_full = Eigen::MatrixXd::Zero(fullSize,fullSize);
		for (int ci = 0; ci < nCams; ci++) {
			for (int cj = ci; cj < nCams; cj++) {
				int pos = sba_crsm_elmidx(Uidxij, ci, cj);
				if (pos < 0) continue; // 没有该块
				double* ppUpa = U + pos * blockSize * blockSize;

				if(cj==ci){
					// 每个ppUpa是6×6上三角块
					for (int ii = 0; ii < blockSize; ii++) {
						for (int jj = ii; jj < blockSize; jj++) {
							double val = ppUpa[ii * blockSize + jj];
							int i= ci * blockSize + ii;
							int j= cj * blockSize + jj;
							U_full(i, j) = val;
							if(i!=j){
								U_full(j, i)=val;
							}
						}
					}
				}else{
					for (int ii = 0; ii < blockSize; ii++) {
						for (int jj = 0; jj < blockSize; jj++) {
							double val = ppUpa[ii * blockSize + jj];
							int i= ci * blockSize + ii;
							int j= cj * blockSize + jj;
							U_full(i, j) = val;
							if(i!=j){
								U_full(j, i)=val;
							}
						}
					}
				}

			}
		}
		// for (int r = 0; r < fullSize; r++) {
		// 	for (int c = r + 1; c < fullSize; c++) {
		// 		if (U_full(r, c) != 0)
		// 			U_full(c, r) = U_full(r, c);
		// 	}
		// }
		FILE* fp1 = fopen(name1, "wb");
		fprintf(fp1, "P5\n%d %d\n255\n", fullSize, fullSize);
		for(int i=0;i<fullSize;i++){
			for(int j=0;j<fullSize;j++){
				unsigned char pixel = U_full(i, j) ? 255 : 0;
				fwrite(&pixel, 1, 1, fp1);
			}
			// fprintf(fp1, "\n");
		}
		fclose(fp1);
		FILE* fp2;
		// const char* fn2 = "E:/zuo/projects/PBA_U.txt";
		errno_t err2 = fopen_s(&fp2, name2, "w");
		for(int i=0;i<fullSize;i++){
			for(int j=0;j<fullSize;j++){
				fprintf(fp2,"%f ",U_full(i,j));
			}
			fprintf(fp2,"\n");
		}
		fclose(fp2);
	}
	double IBA::spectra_lanczos_max_lamuda(SpMat S){
		SparseSymMatProd<double> op(S);//Spectra的矩阵-向量乘算子 y=Sx
		int k = 1, ncv = 6;
		SymEigsSolver<SparseSymMatProd<double>> eigs_max(op, k, ncv);//求1个特征值，Krylov子空间维度是4，Lanczos方法
		eigs_max.init();
		eigs_max.compute(SortRule::LargestAlge);
		if (eigs_max.info() != CompInfo::Successful)
		{
			std::cout << "Max eigenvalue computation failed!" << std::endl;
			return -1.0;
		}
		double lambda_max = eigs_max.eigenvalues()[0];
		return lambda_max;
	}
	double IBA::spectra_lanczos_min_lamuda(SpMat S){
		SparseSymMatProd<double> op(S);//Spectra的矩阵-向量乘算子 y=Sx
		int k = 1, ncv = 6;
		SymEigsSolver<SparseSymMatProd<double>> eigs_min(op, k, ncv);
		eigs_min.init();
		eigs_min.compute(SortRule::SmallestAlge);
		if (eigs_min.info() != CompInfo::Successful)
		{
			std::cout << "Min eigenvalue computation failed!" << std::endl;
			return -1.0;
		}
		double lambda_min = eigs_min.eigenvalues()[0];
		return lambda_min;
	}
	double IBA::shift_invert_min_lamuda(SpMat S){
		// ===== Shift-Invert operator =====
		double sigma = 0.0;
		int k = 1, ncv = 6;
		SparseSymShiftSolve<double> op(S);
		SymEigsShiftSolver<SparseSymShiftSolve<double>>eigs_min(op, k, ncv, sigma);
		eigs_min.init();
		eigs_min.compute(SortRule::LargestMagn, 1000, 1e-10);
		if (eigs_min.info() != CompInfo::Successful)
		{
			std::cout << "Min eigenvalue computation failed!" << std::endl;
			return -1.0;
		}
		double lambda_min = eigs_min.eigenvalues()[0];
		return lambda_min;
	}
	void IBA::computeHMaxSingularValue(double* U, double* V, double* W, double& lambda_max, double& lambda_min, double& cond) {
		int m_cnp = 6, m_pnp = 3;
		int total_dim = m_ncams * m_cnp + m_n3Dpts * m_pnp;
		
		// 构建稀疏矩阵H
		Eigen::SparseMatrix<double> H(total_dim, total_dim);
		std::vector<Eigen::Triplet<double>> triplets;
		
		// 预估非零元素数量以优化内存
		int estimated_nnz = 0;
		
		// U块非零元素估算
		int pos_U = 0;
		for (int i = 0; i < m_ncams; i++) {
			for (int j = 0; j < m_ncams; j++) {
				if (m_umask[i * m_ncams + j] == 1) {
					estimated_nnz += 36;
					pos_U++;
				}
			}
		}
		
		// W块非零元素估算
		int total_obs = 0;
		for (int i = 0; i < m_n3Dpts; i++) {
			total_obs += m_archor[i * 3];
		}
		estimated_nnz += total_obs * 18;
		
		// V块非零元素估算
		estimated_nnz += m_n3Dpts * 9;
		
		triplets.reserve(estimated_nnz);
		
		// 添加U块（对称矩阵）
		pos_U = 0;
		for (int i = 0; i < m_ncams; i++) {
			for (int j = 0; j < m_ncams; j++) {
				if (m_umask[i * m_ncams + j] == 1) {
					double* ptrU = U + pos_U * 36;
					
					for (int ii = 0; ii < m_cnp; ii++) {
						for (int jj = 0; jj < m_cnp; jj++) {
							double val = ptrU[ii * m_cnp + jj];
							if (val != 0.0) {  // 只添加非零元素
								triplets.push_back(Eigen::Triplet<double>(
									i * m_cnp + ii, 
									j * m_cnp + jj, 
									val
								));
								
								// 对称位置（H是对称的）
								if (i != j || ii != jj) {
									triplets.push_back(Eigen::Triplet<double>(
										j * m_cnp + jj, 
										i * m_cnp + ii, 
										val
									));
								}
							}
						}
					}
					pos_U++;
				}
			}
		}
		
		// 添加W和W^T块
		int pos_W = 0;
		for (int i = 0; i < m_n3Dpts; i++) {
			int numfea = m_archor[i * 3];
			for (int j = 0; j < numfea; j++) {
				int nF1 = m_feature[pos_W];
				int nP1 = m_photo[pos_W];
				
				double* ptrW = W + pos_W * m_cnp * m_pnp;
				
				// W块
				for (int ii = 0; ii < m_cnp; ii++) {
					for (int jj = 0; jj < m_pnp; jj++) {
						double val = ptrW[ii * m_pnp + jj];
						if (val != 0.0) {
							triplets.push_back(Eigen::Triplet<double>(
								nP1 * m_cnp + ii,
								m_ncams * m_cnp + nF1 * m_pnp + jj,
								val
							));
						}
					}
				}
				pos_W++;
			}
		}
		
		// 添加V块（对称）
		for (int i = 0; i < m_n3Dpts; i++) {
			double* ptrV = V + i * 9;
			int base = m_ncams * m_cnp + i * m_pnp;
			
			for (int ii = 0; ii < m_pnp; ii++) {
				for (int jj = ii; jj < m_pnp; jj++) {  // 只添加上三角
					double val = ptrV[ii * m_pnp + jj];
					if (val != 0.0) {
						triplets.push_back(Eigen::Triplet<double>(
							base + ii, base + jj, val
						));
						if (ii != jj) {
							triplets.push_back(Eigen::Triplet<double>(
								base + jj, base + ii, val
							));
						}
					}
				}
			}
		}
		
		// 构建稀疏矩阵
		H.setFromTriplets(triplets.begin(), triplets.end());
		lambda_max = spectra_lanczos_max_lamuda(H);
		// lambda_min = spectra_lanczos_min_lamuda(H);
		lambda_min = shift_invert_min_lamuda(H);
		cond = lambda_max / lambda_min;
	}
	void IBA::convertHtoDenseMatrix(double* U,double *V,double *W,const char * name1,const char * name2){
		int m_cnp = 6, m_pnp = 3;
		int total_dim = m_ncams * m_cnp + m_n3Dpts * m_pnp;
		Eigen::SparseMatrix<double> H(total_dim, total_dim);
		std::vector<Eigen::Triplet<double>> triplets;

		// 添加U块（对称矩阵）
		int pos_U = 0;
		for (int i = 0; i < m_ncams; i++) {
			for (int j = 0; j < m_ncams; j++) {
				if (m_umask[i * m_ncams + j] == 1) {
					double* ptrU = U + pos_U * 36;
					
					// 对于对称矩阵，填充所有元素
					for (int ii = 0; ii < m_cnp; ii++) {
						for (int jj = 0; jj < m_cnp; jj++) {
							double val = ptrU[ii * m_cnp + jj];
							
							// 上三角元素
							triplets.push_back(Eigen::Triplet<double>(
								i * m_cnp + ii, 
								j * m_cnp + jj, 
								val
							));
							
							// 如果是非对角块，还需要添加下三角元素
							// 但要注意：对角块只需要添加一次（对称矩阵，上下三角相同）
							if (i != j || ii != jj) {
								// 下三角元素（由于对称性，与上三角相同）
								triplets.push_back(Eigen::Triplet<double>(
									j * m_cnp + jj, 
									i * m_cnp + ii, 
									val
								));
							}
						}
					}
					pos_U++;
				}
			}
		}

		// 添加W和W^T块
        int pos_W = 0;
        for (int i = 0; i < m_n3Dpts; i++) {
            int numfea = m_archor[i * 3];
            for (int j = 0; j < numfea; j++) {
                int nF1 = m_feature[pos_W];
                int nP1 = m_photo[pos_W];
                
                double* ptrW = W + pos_W * m_cnp * m_pnp;
                
                // W块
                for (int ii = 0; ii < m_cnp; ii++) {
                    for (int jj = 0; jj < m_pnp; jj++) {
                        triplets.push_back(Eigen::Triplet<double>(
                            nP1 * m_cnp + ii,
                            m_ncams * m_cnp + nF1 * m_pnp + jj,
                            ptrW[ii * m_pnp + jj]
                        ));
                    }
                }
                pos_W++;
            }
        }
        
        // 添加V块
        for (int i = 0; i < m_n3Dpts; i++) {
            double* ptrV = V + i * 9;
            for (int ii = 0; ii < m_pnp; ii++) {
                for (int jj = 0; jj < m_pnp; jj++) {
                    triplets.push_back(Eigen::Triplet<double>(
                        m_ncams * m_cnp + i * m_pnp + ii,
                        m_ncams * m_cnp + i * m_pnp + jj,
                        ptrV[ii * m_pnp + jj]
                    ));
                }
            }
        }
		H.setFromTriplets(triplets.begin(), triplets.end());
		Eigen::MatrixXd U_full = Eigen::MatrixXd(H);
		
		FILE* fp1 = fopen(name1, "wb");
		fprintf(fp1, "P5\n%d %d\n255\n", total_dim, total_dim);
		for(int i=0;i<total_dim;i++){
			for(int j=0;j<total_dim;j++){
				unsigned char pixel = U_full(i, j) ? 255 : 0;
				fwrite(&pixel, 1, 1, fp1);
			}
			// fprintf(fp1, "\n");
		}
		fclose(fp1);
		FILE* fp2;
		// const char* fn2 = "E:/zuo/projects/PBA_U.txt";
		errno_t err2 = fopen_s(&fp2, name2, "w");
		for(int i=0;i<total_dim;i++){
			for(int j=0;j<total_dim;j++){
				fprintf(fp2,"%f ",U_full(i,j));
			}
			fprintf(fp2,"\n");
		}
		fclose(fp2);
	}
	double IBA::compute_H_inf_norm(double* U, double* V, double* W, double mu)
	{
		int m   = m_ncams;
		int n   = m_n3Dpts;
		int cnp = 6;
		int pnp = 3;
		int Usz = 36;

		// 行累加
		double* row_sum_cam = new double[m * cnp];
		double* row_sum_pt  = new double[n * pnp];

		std::memset(row_sum_cam, 0, sizeof(double) * m * cnp);
		std::memset(row_sum_pt,  0, sizeof(double) * n * pnp);

		int i, j, ii, jj;

		// ==============================
		// 处理 U（相机块）
		// ==============================
		int pos = 0;
		for (i = 0; i < m; i++){//遍历稀疏block U
			for (j = 0; j < m; j++){
				if (m_umask[i * m + j] == 1)//非零元
				{
					double* ptrU = U + pos * Usz;//取到某个block的指针
					for (ii = 0; ii < cnp; ++ii){
						for (jj = 0; jj < cnp; ++jj){
							double val = ptrU[ii * cnp + jj];

							// 对角加阻尼
							// if (i == j && ii == jj)
							// 	val += mu;

							row_sum_cam[i * cnp + ii] += std::fabs(val);
						}
					}
					pos++;
				}
			}
		}

		// ==============================
		// 处理 V（点块）
		// ==============================
		for (i = 0; i < n; i++)
		{
			double* ptrV = V + i * 9;

			for (ii = 0; ii < pnp; ++ii)
			{
				for (jj = 0; jj < pnp; ++jj)
				{
					double val = ptrV[ii * pnp + jj];

					// 对角加阻尼
					// if (ii == jj)
					// 	val += mu;

					row_sum_pt[i * pnp + ii] += std::fabs(val);
				}
			}
		}

		// ==============================
		// 处理 W（最关键）
		// ==============================
		pos = 0;
		for (i = 0; i < m_n3Dpts; i++)
		{
			int numfea = m_archor[i * 3];

			for (j = 0; j < numfea; j++)
			{
				int nF1 = m_feature[pos];  // 点索引
				int nP1 = m_photo[pos];    // 相机索引

				double* ptrW = W + pos * cnp * pnp;

				for (ii = 0; ii < cnp; ++ii)
				{
					for (jj = 0; jj < pnp; ++jj)
					{
						double val = std::fabs(ptrW[ii * pnp + jj]);

						// 相机行贡献
						row_sum_cam[nP1 * cnp + ii] += val;

						// 点行贡献（W^T）
						row_sum_pt[nF1 * pnp + jj] += val;
					}
				}

				pos++;
			}
		}

		// ==============================
		// 求最大行和
		// ==============================
		double max_norm = 0.0;

		for (i = 0; i < m * cnp; i++)
			max_norm = (std::max)(max_norm, row_sum_cam[i]);

		for (i = 0; i < n * pnp; i++)
			max_norm = (std::max)(max_norm, row_sum_pt[i]);

		// 释放
		delete[] row_sum_cam;
		delete[] row_sum_pt;

		return mu / max_norm;
	}
	void IBA::ba_constructSLM(double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij, double mu)
	{
		int i, j, ii, jj, k, l;
		int m = m_ncams, cnp = 6, pnp = 3, Usz = 36, ebsz = 3;
		int pos, pos1, numfea;
		int nF1, nP1, nP2;
		double* ptr1, * ptr2, * ptr3, * ptr4, * ptr5, * ptrS, * ptrE;
		double WV[6 * 3], sum;

		//Copy U matrix to S matrix 
		pos = 0;
		for (i = 0; i < m; i++) for (j = 0; j < m; j++)
		{
			if (m_umask[i * m + j] == 1)// save upper triangle for diagonal element S
			{
				pos1 = sba_crsm_elmidx(&Sidxij, i, j);
				ptr2 = S + pos1 * 36;
				if (i == j)
				{
					ptr1 = U + pos * Usz;
					for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
					{
						ptr2[ii] = ptr1[ii * cnp + ii] + mu;
						// double diag_val = ptr1[ii * cnp + ii];
						// ptr2[ii] = diag_val * (1.0 + mu);
						for (jj = ii + 1; jj < cnp; ++jj)
							ptr2[jj] = ptr1[ii * cnp + jj];
					}
					pos++;
				}
				else
				{
					ptr1 = U + pos * Usz;
					for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
						for (jj = 0; jj < cnp; ++jj)
							ptr2[jj] = ptr1[ii * cnp + jj];
					pos++;
				}
			}
		}

		for (i = 0; i < m * cnp; i++)
			E[i] = ea[i];

		//Create integrated S matrix, S = U - W(V^-1)W^T
		pos = 0;
		for (i = 0; i < m_n3Dpts; i++)
		{
			numfea = m_archor[i * 3];
			for (j = 0; j < numfea; j++)
			{
				nF1 = m_feature[pos];
				nP1 = m_photo[pos];
				memset(WV, 0, sizeof(double) * cnp * 3);

				ptr1 = W + pos * cnp * 3;
				ptr2 = V + nF1 * 3 * 3;
				ptrE = E + nP1 * cnp;


				//WV
				for (ii = 0; ii < cnp; ++ii)
				{
					ptr3 = ptr1 + ii * pnp;
					for (jj = 0; jj < pnp; ++jj)
					{
						for (k = 0, sum = 0.0; k <= jj; ++k)
							sum += ptr3[k] * ptr2[jj * pnp + k];
						for (; k < pnp; ++k)
							sum += ptr3[k] * ptr2[k * pnp + jj];
						for (k = 0, sum = 0.0; k < pnp; k++)
							sum += ptr3[k] * ptr2[jj * pnp + k];
						WV[ii * pnp + jj] = sum;
					}
				}

				for (k = j; k < numfea; k++)
				{
					nP2 = m_photo[pos + (k - j)];

					//W(V^-1)W^T
					ptr3 = W + (pos + (k - j)) * cnp * 3;
					//ptrS = S + (nP1*m*36) + nP2*cnp;

					if (nP1 == nP2)
					{
						pos1 = sba_crsm_elmidx(&Sidxij, nP1, nP2);
						ptrS = S + pos1 * 36;
						for (ii = 0; ii < cnp; ++ii, ptrS += 6)
						{
							ptr5 = WV + ii * pnp;
							for (jj = ii; jj < cnp; ++jj)
							{
								ptr4 = ptr3 + jj * pnp;

								for (l = 0, sum = 0.0; l < pnp; ++l)
									sum += ptr5[l] * ptr4[l];

								ptrS[jj] -= sum;
							}
						}
					}
					else
					{
						pos1 = sba_crsm_elmidx(&Sidxij, nP1, nP2);
						ptrS = S + pos1 * 36;
						for (ii = 0; ii < cnp; ++ii, ptrS += 6)
						{
							ptr5 = WV + ii * pnp;
							for (jj = 0; jj < cnp; ++jj)
							{
								ptr4 = ptr3 + jj * pnp;

								for (l = 0, sum = 0.0; l < pnp; ++l)
									sum += ptr5[l] * ptr4[l];

								ptrS[jj] -= sum;
							}
						}
					}
				}
				//-W^tb
				ptr5 = eb + nF1 * ebsz;
				for (ii = 0; ii < cnp; ++ii)
				{
					ptr4 = WV + ii * pnp;
					for (jj = 0, sum = 0.0; jj < pnp; ++jj)
						sum += ptr4[jj] * ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
					ptrE[ii] -= sum;
				}
				pos++;
			}
		}
	}

	void IBA::ba_constructSGN(double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij)
	{
		int i, j, ii, jj, k, l;
		int m = m_ncams, cnp = 6, pnp = 3, Usz = 36, ebsz = 3;
		int pos, pos1, numfea;
		int nF1, nP1, nP2;
		double* ptr1, * ptr2, * ptr3, * ptr4, * ptr5, * ptrS, * ptrE;
		double WV[6 * 3], sum;

		//Copy U matrix to S matrix 
		pos = 0;
		for (i = 0; i < m; i++) for (j = 0; j < m; j++)
		{
			if (m_umask[i * m + j] == 1)// save upper triangle for diagonal element S
			{
				pos1 = sba_crsm_elmidx(&Sidxij, i, j);
				//printf("%d ", pos1);
				ptr2 = S + pos1 * 36;
				if (i == j)
				{
					ptr1 = U + pos * Usz;
					for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
					{
						ptr2[ii] = ptr1[ii * cnp + ii];
						for (jj = ii + 1; jj < cnp; ++jj)
							ptr2[jj] = ptr1[ii * cnp + jj];
					}
					pos++;
				}
				else
				{
					ptr1 = U + pos * Usz;
					for (ii = 0; ii < cnp; ++ii, ptr2 += 6)
						for (jj = 0; jj < cnp; ++jj)
							ptr2[jj] = ptr1[ii * cnp + jj];
					pos++;
				}
			}
		}

		for (i = 0; i < m * cnp; i++)
			E[i] = ea[i];//

		//Create integrated S matrix, S = U - W(V^-1)W^T
		pos = 0;
		for (i = 0; i < m_n3Dpts; i++)
		{
			numfea = m_archor[i * 3];
			for (j = 0; j < numfea; j++)
			{
				nF1 = m_feature[pos];
				nP1 = m_photo[pos];
				memset(WV, 0, sizeof(double) * cnp * 3);

				ptr1 = W + pos * cnp * 3;
				ptr2 = V + nF1 * 3 * 3;
				ptrE = E + nP1 * cnp;

				//
				for (ii = 0; ii < cnp; ++ii)
				{
					ptr3 = ptr1 + ii * pnp;//
					for (jj = 0; jj < pnp; ++jj)
					{
						for (k = 0, sum = 0.0; k <= jj; ++k)
							sum += ptr3[k] * ptr2[jj * pnp + k];
						for (; k < pnp; ++k)
							sum += ptr3[k] * ptr2[k * pnp + jj];
						for (k = 0, sum = 0.0; k < pnp; k++)//
							sum += ptr3[k] * ptr2[jj * pnp + k];//
						WV[ii * pnp + jj] = sum;
					}
				}

				for (k = j; k < numfea; k++)
				{
					nP2 = m_photo[pos + (k - j)];

					//W(V^-1)W^T
					ptr3 = W + (pos + (k - j)) * cnp * 3;
					pos1 = sba_crsm_elmidx(&Sidxij, nP1, nP2);
					ptrS = S + pos1 * 36;

					if (nP1 == nP2)
					{
						for (ii = 0; ii < cnp; ++ii, ptrS += 6)
						{
							ptr5 = WV + ii * pnp;							
							for (jj = ii; jj < cnp; ++jj)
							{
								ptr4 = ptr3 + jj * pnp;

								for (l = 0, sum = 0.0; l < pnp; ++l)
									sum += ptr5[l] * ptr4[l]; //-W(V^-1)W^T

								ptrS[jj] -= sum; //S-W(V^-1)W^T
							}
						}
					}
					else
					{
						for (ii = 0; ii < cnp; ++ii, ptrS += 6)
						{
							ptr5 = WV + ii * pnp;
							for (jj = 0; jj < cnp; ++jj)
							{
								ptr4 = ptr3 + jj * pnp;

								for (l = 0, sum = 0.0; l < pnp; ++l)
									sum += ptr5[l] * ptr4[l];

								ptrS[jj] -= sum;
							}
						}
					}
				}
				//-W^tb 
				ptr5 = eb + nF1 * ebsz;//
				for (ii = 0; ii < cnp; ++ii)
				{
					ptr4 = WV + ii * pnp;//W(V^-1)
					for (jj = 0, sum = 0.0; jj < pnp; ++jj)
						sum += ptr4[jj] * ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];   W(V^-1)eb
					ptrE[ii] -= sum; //ea-W(V^-1)eb
				}
				pos++;
			}
		}
	}
	void IBA::rotationMatrixToEulerAngles(double* R, double* eulerAngles)
	{
		//assert(isRotationMatrix(R));
		double sy = sqrt(R[0] * R[0] + R[1] * R[1]);

		bool singular = sy < 1e-6;

		double phi, omega, kappa;
		if (!singular)
		{
			phi = atan2(-R[2], sy);
			omega = atan2(R[5], R[8]);
			kappa = atan2(R[1], R[0]);
		}
		else
		{
			phi = 0;
			omega = atan2(R[5], R[8]);
			kappa = atan2(R[1], R[0]);
		}
		eulerAngles[0] = kappa;
		eulerAngles[1] = phi;
		eulerAngles[2] = omega;
	}
	void IBA::eulerAnglesToRotationMatrix(double* eulerAngles, double* R)
	{
		double kappa = eulerAngles[0], phi = eulerAngles[1], omega = eulerAngles[2];
		double R_x[9], R_y[9], R_z[9];
		R_x[0] = 1;   R_x[1] = 0;               R_x[2] = 0;
		R_x[3] = 0;   R_x[4] = cos(omega);      R_x[5] = sin(omega);
		R_x[6] = 0;   R_x[7] = -sin(omega);     R_x[8] = cos(omega);

		R_y[0] = cos(phi);   R_y[1] = 0;       R_y[2] = -sin(phi);
		R_y[3] = 0;          R_y[4] = 1;       R_y[5] = 0;
		R_y[6] = sin(phi);   R_y[7] = 0;       R_y[8] = cos(phi);

		R_z[0] = cos(kappa);  R_z[1] = sin(kappa); R_z[2] = 0;
		R_z[3] = -sin(kappa); R_z[4] = cos(kappa); R_z[5] = 0;
		R_z[6] = 0;           R_z[7] = 0;          R_z[8] = 1;

		double tmp[9];
		tmp[0] = R_y[0] * R_z[0] + R_y[1] * R_z[3] + R_y[2] * R_z[6];
		tmp[1] = R_y[0] * R_z[1] + R_y[1] * R_z[4] + R_y[2] * R_z[7];
		tmp[2] = R_y[0] * R_z[2] + R_y[1] * R_z[5] + R_y[2] * R_z[8];

		tmp[3] = R_y[3] * R_z[0] + R_y[4] * R_z[3] + R_y[5] * R_z[6];
		tmp[4] = R_y[3] * R_z[1] + R_y[4] * R_z[4] + R_y[5] * R_z[7];
		tmp[5] = R_y[3] * R_z[2] + R_y[4] * R_z[5] + R_y[5] * R_z[8];

		tmp[6] = R_y[6] * R_z[0] + R_y[7] * R_z[3] + R_y[8] * R_z[6];
		tmp[7] = R_y[6] * R_z[1] + R_y[7] * R_z[4] + R_y[8] * R_z[7];
		tmp[8] = R_y[6] * R_z[2] + R_y[7] * R_z[5] + R_y[8] * R_z[8];

		R[0] = R_x[0] * tmp[0] + R_x[1] * tmp[3] + R_x[2] * tmp[6];
		R[1] = R_x[0] * tmp[1] + R_x[1] * tmp[4] + R_x[2] * tmp[7];
		R[2] = R_x[0] * tmp[2] + R_x[1] * tmp[5] + R_x[2] * tmp[8];

		R[3] = R_x[3] * tmp[0] + R_x[4] * tmp[3] + R_x[5] * tmp[6];
		R[4] = R_x[3] * tmp[1] + R_x[4] * tmp[4] + R_x[5] * tmp[7];
		R[5] = R_x[3] * tmp[2] + R_x[4] * tmp[5] + R_x[5] * tmp[8];

		R[6] = R_x[6] * tmp[0] + R_x[7] * tmp[3] + R_x[8] * tmp[6];
		R[7] = R_x[6] * tmp[1] + R_x[7] * tmp[4] + R_x[8] * tmp[7];
		R[8] = R_x[6] * tmp[2] + R_x[7] * tmp[5] + R_x[8] * tmp[8];
	}
	//bool IBA::ba_initialize(char* szCamera, char* szFeature, char* szCalib, char* szXYZ)
	//{
	//	printf("SparseBA: Sparse Bundle Adjustment Version 1.0\n");
	//	FILE* fp;
	//
	//	//must input initial initial camera pose file and projection image points file
	//	fp = fopen(szCamera, "r");
	//	if (fp == NULL)
	//	{
	//		fprintf(stderr, "SparseBA: Missing initial camera poses file! \n");
	//		exit(1);
	//	}
	//	else
	//		fclose(fp);
	//
	//	fp = fopen(szFeature, "r");
	//	if (fp == NULL)
	//	{
	//		fprintf(stderr, "SparseBA: Missing feature projection points file! \n");
	//		exit(1);
	//	}
	//	else
	//		fclose(fp);
	//
	//	if (szCalib != NULL)
	//	{
	//		m_bFocal = false;
	//		m_K = (double*)malloc(9 * sizeof(double));
	//		ba_readCameraPoseration(szCalib, m_K);
	//	}
	//
	//	if (szXYZ != NULL)
	//		m_bProvideXYZ = true;
	//
	//	savePara = 0;
	//	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	//	sba_readAndInitialize(szCamera, szFeature, &m_ncams, &m_n3Dpts, &m_n2Dprojs, &m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
	//		&m_imgpts, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort);
	//
	//	std::string originalPath(szCamera);
	//	size_t pos = originalPath.find_last_of("/\\");
	//	std::string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	//	std::string newPath = parentPath + "/" + "Triangulatedxyz.txt";
	//
	//	if (savePara == 1)
	//		ba_saveTriangulatedxyz(newPath.c_str(), m_motstruct);
	//
	//	printf("Number of cameras: %d\n", m_ncams);
	//	printf("Number of points: %d\n", m_n3Dpts);
	//	printf("Number of projections: %d\n", m_n2Dprojs);
	//	return true;
	//}
#else
    #error "Unsupported platform!"
#endif

void IBA::ba_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams)
{
	int n;
	int nframes;
	int ptno = 0;
	int frameno;
	double* ptr1 = projs;//��¼ÿ��ͶӰ���xy��������ŵģ�2��2��������¼�ʹ洢��

	int i;
	//, j;
	//read all projection point, initialize three feature angle at the same time
	//��ȡ����ͶӰ�㣬ͬʱ��ʼ������������***�ص�***
	while (!feof(fp))//һ��һ�ж�Match��
	{
		n = readNInts(fp, &nframes, 1);  //��ȡMatch-FeaturePoint�ļ�ÿһ�еĵ�һ�У���ͬ����ĸ���Ϊnframes/3D���2DͶӰ����
		if (n != 1)
			break;//һ�㶼��1����ʾ�ɹ�

		Eigen::MatrixXd A(2 * nframes, 4);
		//if (nframes > 3)
		//	nframes = 3;


		for (i = 0; i < nframes; ++i)//һ��ͶӰ��һ��ͶӰ��Ķ�ȡ��������nframes
		{
			n = readNInts(fp, &frameno, 1); //�ڶ��Σ���ȡ��һ��ͶӰ������-image index

			if (frameno >= ncams)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "BA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles(fp, ptr1, 2); //�ڶ��Σ���ȡ��һ��ͶӰ���������� ptr1[0] ptr1[1]

			//Ҳ����1��id+2������
			if (n != 3)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "BA:reading image projections wrong!\n");
				return;
			}

			const Eigen::Vector2d pt(ptr1[0], ptr1[1]);//��������
			double* ptr = m_P + frameno * 12;
			const Eigen::Matrix<double, 3, 4> Pmat = (Eigen::Matrix<double, 3, 4>() <<
				ptr[0], ptr[1], ptr[2], ptr[3],
				ptr[4], ptr[5], ptr[6], ptr[7],
				ptr[8], ptr[9], ptr[10], ptr[11]).finished();
			A.row(2 * i) = pt(0) * Pmat.row(2) - Pmat.row(0);
			A.row(2 * i + 1) = pt(1) * Pmat.row(2) - Pmat.row(1);
		}

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
		Eigen::Vector4d X = svd.matrixV().col(3);
		X /= X(3);

		(m_motstruct + m_ncams * 6 + ptno * 3)[0] = X[0];
		(m_motstruct + m_ncams * 6 + ptno * 3)[1] = X[1];
		(m_motstruct + m_ncams * 6 + ptno * 3)[2] = X[2];
		ptno++;//3D���������ţ�point NO. OR photo NO.
	}
}

//�����ļ��з�ע���е�����
int IBA::findNcameras(FILE* fp)
{
	int lineno, ncams, ch;

	lineno = ncams = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#') { /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);

		SKIP_LINE(fp);
		++lineno;
		if (ferror(fp))
		{
			fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		++ncams;
	}
	return ncams;
}

IBA::IBA(void)
{
}

IBA::~IBA(void)
{
}

int IBA::countNDoubles(FILE* fp)
{
	int lineno, ch, np, i;
	char buf[MAXSTRLEN], * s;
	double dummy;

	lineno = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) return 0;

		ungetc(ch, fp);
		++lineno;
		if (!fgets(buf, MAXSTRLEN - 1, fp)) { /* read the line found... */
			fprintf(stderr, "countNDoubles(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		/* ...and count the number of doubles it has */
		for (np = i = 0, s = buf; 1; ++np, s += i) {
			ch = sscanf_s(s, "%lf%n", &dummy, &i);
			if (ch == 0 || ch == EOF) break;
		}

		rewind(fp);
		return np;
	}
	return 0; // should not reach this point
}

int IBA::skipNDoubles(FILE* fp, int nvals)
{
	int i;
	int j;

	for (i = 0; i < nvals; ++i)
	{
		j = fscanf_s(fp, "%*f");
		if (j == EOF) return EOF;

		if (ferror(fp)) return EOF - 1;
	}

	return nvals;
}

void IBA::readNpointsAndNprojections(FILE* fp, int* n3Dpts, int pnp, int* nprojs, int mnp)
{
	int nfirst, lineno, npts, nframes, ch, n;

	/* #parameters for the first line */
	nfirst = countNDoubles(fp);

	*n3Dpts = *nprojs = lineno = npts = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n = readNInts(fp, &nframes, 1);
		
		if (n != 1)
			exit(1);

		//printf("%d ", nframes);

		SKIP_LINE(fp);
		*nprojs += nframes;
		++npts;
	}

	*n3Dpts = npts;
}

void IBA::ba_readCablibration(FILE* fp, double* K)
{
	int n = fscanf_s(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &K[0], &K[3], &K[6], &K[1], &K[4], &K[7], &K[2], &K[5], &K[8]);

	if (n != 9)
	{
		fprintf(stderr, "BA error: Format of Calibaration is wrong\n");
		exit(1);
	}
}

void IBA::ba_readCameraPose(FILE* fp, double* params, int* m_v)
{
	int n, num, lineno = 0;
	double* tofilter;
	double* pPrams = params;
	int* pm_C = m_v;
	int jjg = 0;
	//the number of element per line is 8, it represents that focal length vary, or it is constant
	num = countNDoubles(fp);
	if (num == 8)
	{
		m_bFocal = true;
		m_K = (double*)malloc(m_ncams * 2 * sizeof(double));
		tofilter = (double*)malloc(8 * sizeof(double));
	}
	else if (num == 7) {
		tofilter = (double*)malloc(7 * sizeof(double));
	}
	else
		tofilter = (double*)malloc(6 * sizeof(double));

	while (!feof(fp))
	{
		if (num == 6) {
			n = readNDoubles(fp, tofilter, 6);
			if (n == -1) {
				//printf("%d %d\n", lineno, jjg);
				break;
			}
			#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
				Camera cam;
				cam.euler_angle[0] = tofilter[0];
				cam.euler_angle[1] = tofilter[1];
				cam.euler_angle[2] = tofilter[2];
				cam.camera_center[0] = tofilter[3];
				cam.camera_center[1] = tofilter[4];
				cam.camera_center[2] = tofilter[5];
				cam.camidx = 1;
				cams.push_back(cam);
			#elif defined(_WIN32)
				
			#else
				#error "Unsupported platform!"
			#endif

			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];
			*pm_C = 1;
			//m_v[lineno++] = 1;
		}
		if (num == 7){
			n = readNDoubles(fp, tofilter, 7);
			if (n == -1) {
				//printf("%d %d\n", lineno, jjg);
				break;
			}
			#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
				Camera cam;
				cam.euler_angle[0] = tofilter[0];
				cam.euler_angle[1] = tofilter[1];
				cam.euler_angle[2] = tofilter[2];
				cam.camera_center[0] = tofilter[3];
				cam.camera_center[1] = tofilter[4];
				cam.camera_center[2] = tofilter[5];
				cam.camidx = static_cast<int>(tofilter[6]);
				cams.push_back(cam);
			#elif defined(_WIN32)
				
			#else
				#error "Unsupported platform!"
			#endif

			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];
			*pm_C = static_cast<int>(tofilter[6]);
			++jjg;
			//m_v[lineno++] = tofilter[6];
			//printf("%d\n", *pm_C);
			/*------------------------------------������-------------------------------------------*/
			//for (int i = 0;i < 6;i++)
			//	printf("%f ", tofilter[i]);
			//printf("\n");
			/*-------------------------------------������------------------------------------------*/
		}
		if (num == 8){
			n = readNDoubles(fp, tofilter, 8);
			if (n == -1) {
				//printf("%d %d\n", lineno, jjg);
				break;
			}
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];

			m_K[lineno * 2] = tofilter[6];
			m_K[lineno * 2 + 1] = tofilter[7];
		}
			
		pPrams += 6;
		if (num == 7 || num == 6) {
			pm_C += 1;
		}
		++lineno;
	}
	if (tofilter != NULL) {
		free(tofilter);
		tofilter = NULL;
	}
}

int IBA::readNInts(FILE* fp, int* vals, int nvals)
{
	int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i) {//��һ�Σ�nvals����1
		j = fscanf_s(fp, "%d", vals + i);//��һ�Σ���ȡ��һ����Ҳ����1����ά�������2DͶӰ����
		if (j == EOF) return EOF;//��һ�Σ�һ�㲻���ǣ�����j=1����ʾ�ɹ�

		if (j != 1 || ferror(fp)) return EOF - 1;//��һ�Σ�һ�㶼�ǣ�j=1����ʾ�ɹ�

		n += j;//��һ����Ϊn=0�����Ծ���j��Ҳ���ǳɹ���־j=1
	}//֮��Ļ���fscanfÿ��һ����ô�ͻ��ƶ�ָ�룬Ȼ��%d��������ո񡢻��з��������հ��ַ���ֱ��������%d����������
	//����������EOF����j!=0��
	return n;
}

int IBA::readNDoubles(FILE* fp, double* vals, int nvals)
{
	int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i)//nvals��2����ȡxy����
	{
		j = fscanf_s(fp, "%lf", vals + i);//���￪ʼ������float��Ҳ��������ֵ
		if (j == EOF) return EOF;

		if (j != 1 || ferror(fp)) return EOF - 1;

		n += j;//����֮���2��Ϊ����2��
	}

	return n;//����2
}
//void IBA::ba_readCameraPose_(char* fname)
//{
//	FILE* fp;
//	int  ch = EOF;
//
//	if ((fp = fopen(fname, "r")) == NULL)
//	{
//		fprintf(stderr, "BA: Cannot open calbration file %s, exiting\n", fname);
//		return;
//	}
//	if (tpe == bal) {
//		//double* ptr;
//		double f;
//		double k0, k1, k2, k3, k4, k5, k6, k7, k8;
//		for (int i = 0; i < m_ncams; i++)
//		{
//			//ptr = ical + 9 * i;
//			int num = fscanf(fp, "%lf", &f);
//			if (num != 1)
//			{
//				fprintf(stderr, "BA error: Format of Calibration file is wrong");
//				return;
//			}
//			k4 = k0 = f;
//			k1 = k2 = k3 = k5 = k6 = k7 = 0;
//			k8 = 1;
//			m_K.push_back(k0); m_K.push_back(k3); m_K.push_back(k6);
//			m_K.push_back(k1); m_K.push_back(k4); m_K.push_back(k7);
//			m_K.push_back(k2); m_K.push_back(k5); m_K.push_back(k8);
//			//printf("%f %f %f\n", m_K[9 * i + 0], m_K[9 * i + 3], m_K[9 * i + 6]);
//			//printf("%f %f %f\n", m_K[9 * i + 1], m_K[9 * i + 4], m_K[9 * i + 7]);
//			//printf("%f %f %f\n", m_K[9 * i + 2], m_K[9 * i + 5], m_K[9 * i + 8]);
//		}
//	}
//	else if (tpe == colmap) {
//		double k0, k1, k2, k3, k4, k5, k6, k7, k8;
//		for (int i = 0; i < nc_; i++) {
//			//���ж�ȡ1-2-3 || 4-5-6 || 7-8-9
//			int num = fscanf_s(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", \
//				& k0, &k3, &k6, \
//				& k1, &k4, &k7, \
//				& k2, &k5, &k8);
//			m_K.push_back(k0); m_K.push_back(k3); m_K.push_back(k6);
//			m_K.push_back(k1); m_K.push_back(k4); m_K.push_back(k7);
//			m_K.push_back(k2); m_K.push_back(k5); m_K.push_back(k8);
//			printf("%f %f %f\n", m_K[9 * i + 0], m_K[9 * i + 3], m_K[9 * i + 6]);
//			printf("%f %f %f\n", m_K[9 * i + 1], m_K[9 * i + 4], m_K[9 * i + 7]);
//			printf("%f %f %f\n", m_K[9 * i + 2], m_K[9 * i + 5], m_K[9 * i + 8]);
//		}
//		if (m_K.size() != 9 * nc_)
//		{
//			fprintf(stderr, "BA error: Format of Calibration file is wrong");
//			return;
//		}
//	}
//
//	fclose(fp);
//	fp = NULL;
//}
void IBA::ba_readCameraPoseration(char* fname, double* ical)
{
	FILE* fp = nullptr;
	int  ch = EOF;
	fopen_s(&fp, fname, "r");
	if (fp == nullptr)
	{
		fprintf(stderr, "BA: Cannot open calbration file %s, exiting\n", fname);
		return;
	}
	if (tpe == bal) {
		double* ptr;
		for (int i = 0; i < m_ncams; i++)
		{
			ptr = ical + 9 * i;
			int num = fscanf_s(fp, "%lf", &ptr[0]);
			if (num != 1)
			{
				fprintf(stderr, "BA error: Format of Calibration file is wrong");
				return;
			}
			ptr[4] = ptr[0];
			ptr[1] = ptr[2] = ptr[3] = ptr[5] = ptr[6] = ptr[7] = 0;
			ptr[8] = 1;
			//printf("%f %f %f\n", ptr[0], ptr[1], ptr[2]);
			//printf("%f %f %f\n", ptr[3], ptr[4], ptr[5]);
			//printf("%f %f %f\n", ptr[6], ptr[7], ptr[8]);
		}
	}
	else if (tpe == colmap) {
		int s = 0;
		for (int i = 0; i < nc_; i++) {
			//���ж�ȡ1-2-3 || 4-5-6 || 7-8-9
			int num = fscanf_s(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", \
				& ical[9 * i + 0], &ical[9 * i + 3], &ical[9 * i + 6], \
				& ical[9 * i + 1], &ical[9 * i + 4], &ical[9 * i + 7], \
				& ical[9 * i + 2], &ical[9 * i + 5], &ical[9 * i + 8]);
			#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
				Intrinsic intr;
				intr.fx = ical[9 * i + 0];
				intr.fy = ical[9 * i + 4];
				intr.cx = ical[9 * i + 6];
				intr.cy = ical[9 * i + 7];
				intrs.push_back(intr);
			#elif defined(_WIN32)
				
			#else
				#error "Unsupported platform!"
			#endif

			// intrs.push_back({ical[9 * i + 0], ical[9 * i + 4], ical[9 * i + 6], ical[9 * i + 7]});
			s += num;
			//printf("%f %f %f\n", ical[9 * i + 0], ical[9 * i + 3], ical[9 * i + 6]);
			//printf("%f %f %f\n", ical[9 * i + 1], ical[9 * i + 4], ical[9 * i + 7]);
			//printf("%f %f %f\n", ical[9 * i + 2], ical[9 * i + 5], ical[9 * i + 8]);
		}
		if (s != 9 * nc_)
		{
			fprintf(stderr, "BA error: Format of Calibration file is wrong");
			return;
		}
	}

	fclose(fp);
	fp = NULL;
}

void IBA::ba_updateKR(double* KR, double* KdA, double* KdB, double* KdG, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pKR, * pKdA, * pKdB, * pKdG, * pK;
		double matR[9];//�����ת����
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];
		//for (i = 0; i < m_ncams; i++) {
		//	printf("%d\n", m_V[i]);
		//}
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;//ָ�������ξ���
			/*kappa phi omegaϵͳ*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			//omega��ת����
			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);

			//phi��ת����
			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);

			//kappa��ת����
			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			//matRG�������omega��һ�׵�
			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[2]);	matDRG[5] = cos(ptAngle[2]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[2]);	matDRG[8] = -sin(ptAngle[2]);

			//matRB�������phi��һ�׵�
			matDRB[0] = -sin(ptAngle[1]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[1]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[1]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[1]);

			//matRA�������kappa��һ�׵�
			matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//pKR=KR*matR
			pKR = KR + i * 9;
			if (tpe == bal) {
				pK = K + i * 9;//����BAL����
			}
			else if (tpe == colmap) {
				pK = K + (m_V[i] - 1) * 9;//���ڷ�BAL����
			}
			//printf("%d\n", m_V[i]);
			//printf("%f %f %f\n", pK[0], pK[3], pK[6]);
			//printf("%f %f %f\n", pK[1], pK[4], pK[7]);
			//printf("%f %f %f\n", pK[2], pK[5], pK[8]);
			pKR[0] = pK[0] * matR[0] + pK[3] * matR[3] + pK[6] * matR[6];
			pKR[1] = pK[0] * matR[1] + pK[3] * matR[4] + pK[6] * matR[7];
			pKR[2] = pK[0] * matR[2] + pK[3] * matR[5] + pK[6] * matR[8];
			pKR[3] = pK[1] * matR[0] + pK[4] * matR[3] + pK[7] * matR[6];
			pKR[4] = pK[1] * matR[1] + pK[4] * matR[4] + pK[7] * matR[7];
			pKR[5] = pK[1] * matR[2] + pK[4] * matR[5] + pK[7] * matR[8];
			pKR[6] = pK[2] * matR[0] + pK[5] * matR[3] + pK[8] * matR[6];
			pKR[7] = pK[2] * matR[1] + pK[5] * matR[4] + pK[8] * matR[7];
			pKR[8] = pK[2] * matR[2] + pK[5] * matR[5] + pK[8] * matR[8];

			//pKR�������omega��һ�׵�
			pKdG = KdG + i * 9;
			tmp1[0] = pK[0] * matDRG[0] + pK[3] * matDRG[3] + pK[6] * matDRG[6];
			tmp1[1] = pK[1] * matDRG[0] + pK[4] * matDRG[3] + pK[7] * matDRG[6];
			tmp1[2] = pK[2] * matDRG[0] + pK[5] * matDRG[3] + pK[8] * matDRG[6];
			tmp1[3] = pK[0] * matDRG[1] + pK[3] * matDRG[4] + pK[6] * matDRG[7];
			tmp1[4] = pK[1] * matDRG[1] + pK[4] * matDRG[4] + pK[7] * matDRG[7];
			tmp1[5] = pK[2] * matDRG[1] + pK[5] * matDRG[4] + pK[8] * matDRG[7];
			tmp1[6] = pK[0] * matDRG[2] + pK[3] * matDRG[5] + pK[6] * matDRG[8];
			tmp1[7] = pK[1] * matDRG[2] + pK[4] * matDRG[5] + pK[7] * matDRG[8];
			tmp1[8] = pK[2] * matDRG[2] + pK[5] * matDRG[5] + pK[8] * matDRG[8];

			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdG[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdG[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdG[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdG[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdG[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdG[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdG[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdG[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdG[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//pKR�������phi��һ�׵�
			pKdB = KdB + i * 9;
			tmp1[0] = pK[0] * matRG[0] + pK[3] * matRG[3] + pK[6] * matRG[6];
			tmp1[1] = pK[1] * matRG[0] + pK[4] * matRG[3] + pK[7] * matRG[6];
			tmp1[2] = pK[2] * matRG[0] + pK[5] * matRG[3] + pK[8] * matRG[6];
			tmp1[3] = pK[0] * matRG[1] + pK[3] * matRG[4] + pK[6] * matRG[7];
			tmp1[4] = pK[1] * matRG[1] + pK[4] * matRG[4] + pK[7] * matRG[7];
			tmp1[5] = pK[2] * matRG[1] + pK[5] * matRG[4] + pK[8] * matRG[7];
			tmp1[6] = pK[0] * matRG[2] + pK[3] * matRG[5] + pK[6] * matRG[8];
			tmp1[7] = pK[1] * matRG[2] + pK[4] * matRG[5] + pK[7] * matRG[8];
			tmp1[8] = pK[2] * matRG[2] + pK[5] * matRG[5] + pK[8] * matRG[8];

			tmp2[0] = tmp1[0] * matDRB[0] + tmp1[3] * matDRB[3] + tmp1[6] * matDRB[6];
			tmp2[1] = tmp1[1] * matDRB[0] + tmp1[4] * matDRB[3] + tmp1[7] * matDRB[6];
			tmp2[2] = tmp1[2] * matDRB[0] + tmp1[5] * matDRB[3] + tmp1[8] * matDRB[6];
			tmp2[3] = tmp1[0] * matDRB[1] + tmp1[3] * matDRB[4] + tmp1[6] * matDRB[7];
			tmp2[4] = tmp1[1] * matDRB[1] + tmp1[4] * matDRB[4] + tmp1[7] * matDRB[7];
			tmp2[5] = tmp1[2] * matDRB[1] + tmp1[5] * matDRB[4] + tmp1[8] * matDRB[7];
			tmp2[6] = tmp1[0] * matDRB[2] + tmp1[3] * matDRB[5] + tmp1[6] * matDRB[8];
			tmp2[7] = tmp1[1] * matDRB[2] + tmp1[4] * matDRB[5] + tmp1[7] * matDRB[8];
			tmp2[8] = tmp1[2] * matDRB[2] + tmp1[5] * matDRB[5] + tmp1[8] * matDRB[8];

			pKdB[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdB[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdB[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdB[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdB[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdB[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdB[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdB[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdB[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//pKR�������kappa��һ�׵�
			pKdA = KdA + i * 9;
			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdA[0] = tmp2[0] * matDRA[0] + tmp2[3] * matDRA[3] + tmp2[6] * matDRA[6];
			pKdA[3] = tmp2[1] * matDRA[0] + tmp2[4] * matDRA[3] + tmp2[7] * matDRA[6];
			pKdA[6] = tmp2[2] * matDRA[0] + tmp2[5] * matDRA[3] + tmp2[8] * matDRA[6];
			pKdA[1] = tmp2[0] * matDRA[1] + tmp2[3] * matDRA[4] + tmp2[6] * matDRA[7];
			pKdA[4] = tmp2[1] * matDRA[1] + tmp2[4] * matDRA[4] + tmp2[7] * matDRA[7];
			pKdA[7] = tmp2[2] * matDRA[1] + tmp2[5] * matDRA[4] + tmp2[8] * matDRA[7];
			pKdA[2] = tmp2[0] * matDRA[2] + tmp2[3] * matDRA[5] + tmp2[6] * matDRA[8];
			pKdA[5] = tmp2[1] * matDRA[2] + tmp2[4] * matDRA[5] + tmp2[7] * matDRA[8];
			pKdA[8] = tmp2[2] * matDRA[2] + tmp2[5] * matDRA[5] + tmp2[8] * matDRA[8];

		}
	}
	else
	{
		int i = 0;
		double* ptAngle;
		double* pKR, * pKdA, * pKdB, * pKdG;
		double matR[9];
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];
		double K[9];
		memset(K, 0, 9 * sizeof(double));
		K[8] = 1;
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;

			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			//����omega����ת����
			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);
			//����phi����ת����
			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);
			//����kappa����ת����
			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;
			//��omega��һ�׵�
			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[2]);	matDRG[5] = cos(ptAngle[2]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[2]);	matDRG[8] = -sin(ptAngle[2]);
			//��phi��һ�׵�
			matDRB[0] = -sin(ptAngle[1]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[1]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[1]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[1]);
			//��kappa��һ�׵�
			matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//KR

			K[0] = m_K[i * 3];
			K[4] = m_K[i * 3];

			pKR = KR + i * 9;
			pKR[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pKR[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pKR[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pKR[3] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pKR[4] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pKR[5] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pKR[6] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pKR[7] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pKR[8] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];

			//KdG
			pKdG = KdG + i * 9;
			tmp1[0] = K[0] * matDRG[0] + K[3] * matDRG[3] + K[6] * matDRG[6];
			tmp1[1] = K[1] * matDRG[0] + K[4] * matDRG[3] + K[7] * matDRG[6];
			tmp1[2] = K[2] * matDRG[0] + K[5] * matDRG[3] + K[8] * matDRG[6];
			tmp1[3] = K[0] * matDRG[1] + K[3] * matDRG[4] + K[6] * matDRG[7];
			tmp1[4] = K[1] * matDRG[1] + K[4] * matDRG[4] + K[7] * matDRG[7];
			tmp1[5] = K[2] * matDRG[1] + K[5] * matDRG[4] + K[8] * matDRG[7];
			tmp1[6] = K[0] * matDRG[2] + K[3] * matDRG[5] + K[6] * matDRG[8];
			tmp1[7] = K[1] * matDRG[2] + K[4] * matDRG[5] + K[7] * matDRG[8];
			tmp1[8] = K[2] * matDRG[2] + K[5] * matDRG[5] + K[8] * matDRG[8];

			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdG[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdG[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdG[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdG[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdG[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdG[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdG[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdG[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdG[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//KdB
			pKdB = KdB + i * 9;
			tmp1[0] = K[0] * matRG[0] + K[3] * matRG[3] + K[6] * matRG[6];
			tmp1[1] = K[1] * matRG[0] + K[4] * matRG[3] + K[7] * matRG[6];
			tmp1[2] = K[2] * matRG[0] + K[5] * matRG[3] + K[8] * matRG[6];
			tmp1[3] = K[0] * matRG[1] + K[3] * matRG[4] + K[6] * matRG[7];
			tmp1[4] = K[1] * matRG[1] + K[4] * matRG[4] + K[7] * matRG[7];
			tmp1[5] = K[2] * matRG[1] + K[5] * matRG[4] + K[8] * matRG[7];
			tmp1[6] = K[0] * matRG[2] + K[3] * matRG[5] + K[6] * matRG[8];
			tmp1[7] = K[1] * matRG[2] + K[4] * matRG[5] + K[7] * matRG[8];
			tmp1[8] = K[2] * matRG[2] + K[5] * matRG[5] + K[8] * matRG[8];

			tmp2[0] = tmp1[0] * matDRB[0] + tmp1[3] * matDRB[3] + tmp1[6] * matDRB[6];
			tmp2[1] = tmp1[1] * matDRB[0] + tmp1[4] * matDRB[3] + tmp1[7] * matDRB[6];
			tmp2[2] = tmp1[2] * matDRB[0] + tmp1[5] * matDRB[3] + tmp1[8] * matDRB[6];
			tmp2[3] = tmp1[0] * matDRB[1] + tmp1[3] * matDRB[4] + tmp1[6] * matDRB[7];
			tmp2[4] = tmp1[1] * matDRB[1] + tmp1[4] * matDRB[4] + tmp1[7] * matDRB[7];
			tmp2[5] = tmp1[2] * matDRB[1] + tmp1[5] * matDRB[4] + tmp1[8] * matDRB[7];
			tmp2[6] = tmp1[0] * matDRB[2] + tmp1[3] * matDRB[5] + tmp1[6] * matDRB[8];
			tmp2[7] = tmp1[1] * matDRB[2] + tmp1[4] * matDRB[5] + tmp1[7] * matDRB[8];
			tmp2[8] = tmp1[2] * matDRB[2] + tmp1[5] * matDRB[5] + tmp1[8] * matDRB[8];

			pKdB[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdB[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdB[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdB[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdB[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdB[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdB[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdB[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdB[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//KdA
			pKdA = KdA + i * 9;
			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdA[0] = tmp2[0] * matDRA[0] + tmp2[3] * matDRA[3] + tmp2[6] * matDRA[6];
			pKdA[3] = tmp2[1] * matDRA[0] + tmp2[4] * matDRA[3] + tmp2[7] * matDRA[6];
			pKdA[6] = tmp2[2] * matDRA[0] + tmp2[5] * matDRA[3] + tmp2[8] * matDRA[6];
			pKdA[1] = tmp2[0] * matDRA[1] + tmp2[3] * matDRA[4] + tmp2[6] * matDRA[7];
			pKdA[4] = tmp2[1] * matDRA[1] + tmp2[4] * matDRA[4] + tmp2[7] * matDRA[7];
			pKdA[7] = tmp2[2] * matDRA[1] + tmp2[5] * matDRA[4] + tmp2[8] * matDRA[7];
			pKdA[2] = tmp2[0] * matDRA[2] + tmp2[3] * matDRA[5] + tmp2[6] * matDRA[8];
			pKdA[5] = tmp2[1] * matDRA[2] + tmp2[4] * matDRA[5] + tmp2[7] * matDRA[8];
			pKdA[8] = tmp2[2] * matDRA[2] + tmp2[5] * matDRA[5] + tmp2[8] * matDRA[8];

		}
	}
}
void IBA::ba_constructP(double* P, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9];
		double matT[3];

		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;//ָ�������ξ���
			/*kappa phi omegaϵͳ*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];

			//test
			//double t1 = pP[0] * (2 - ptAngle[3]) + pP[1] * (2 - ptAngle[4]) + pP[2] * (2 - ptAngle[5]);
			//double t2 = pP[4] * (2 - ptAngle[3]) + pP[5] * (2 - ptAngle[4]) + pP[6] * (2 - ptAngle[5]);
			//double t3 = pP[8] * (2 - ptAngle[3]) + pP[9] * (2 - ptAngle[4]) + pP[10] * (2 - ptAngle[5]);

			//double u1 = t1 / t3;
			//double v1 = t2 / t3;

			//t1 = pP[0] * 2 + pP[1] * 2 + pP[2] * 2 + pP[3];
			//t2 = pP[4] * 2 + pP[5] * 2 + pP[6] * 2 + pP[7];
			//t3 = pP[8] * 2 + pP[9] * 2 + pP[10] * 2 + pP[11];
			//double u2 = t1 / t3;
			//double v2 = t2 / t3;
		}
	}
	else
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9], matT[3];
		double K[9];
		memset(K, 0, 9 * sizeof(double));
		K[8] = 1;
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;

			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//KR

			K[0] = m_K[i * 2];
			K[4] = m_K[i * 2 + 1];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];
		}
	}
}



void IBA::readNpointsAndNprojectionsFromProj(FILE* fp, int& n3Dpts, int& nprojs)
{
	int nfirst, lineno, npts, nframes, ch, n;
	nprojs = 0;
	n3Dpts = 0;
	npts = 0;

	/* #parameters for the first line */
	nfirst = countNDoubles(fp);

	//*n3Dpts=*nprojs=lineno=npts=0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n = readNInts(fp, &nframes, 1);
		if (n != 1)
		{
			fprintf(stderr, "readNpointsAndNprojections(): error reading input file, line %d: "
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		SKIP_LINE(fp);
		nprojs += nframes;
		++npts;
	}

	n3Dpts = npts;
}

void IBA::readPointProjections(FILE* fp, double* imgpts, int* photo, int* imgptsSum, int n3Dpts, int n2Dprojs)
{
	int nframes, ch, lineno, ptno, frameno, n;
	int i;
	int nproj2D = 0;

	lineno = ptno = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			lineno++;

			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);

		n = readNInts(fp, &nframes, 1);  /* read in number of image projections */
		if (n != 1)
		{
			fprintf(stderr, "sba_readProjectionAndInitilizeFeature(): error reading input file, line %d:\n"
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		imgptsSum[ptno] = nframes;

		for (i = 0; i < nframes; ++i)
		{
			n = readNInts(fp, &frameno, 1); /* read in frame number... */

			photo[nproj2D] = frameno;

			n += readNDoubles(fp, imgpts + nproj2D * 2, 2); /* ...and image projection */

			nproj2D++;
		}
		fscanf_s(fp, "\n"); // consume trailing newline

		lineno++;
		ptno++;
	}
}
void IBA::readImagePts(const char* szProj, double** imgpts, int** photo, int** imgptsSum, int& n3Dpts, int& n2Dprojs)
{
	FILE* fpp = nullptr;
	fopen_s(&fpp, szProj, "r");
	if (fpp == nullptr) {
		fprintf(stderr, "cannot open file %s, exiting\n", szProj);
		exit(1);
	}
	readNpointsAndNprojectionsFromProj(fpp, n3Dpts, n2Dprojs);

	*imgpts = (double*)malloc(n2Dprojs * 2 * sizeof(double));
	if (*imgpts == NULL) {
		fprintf(stderr, "memory allocation for 'imgpts' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*photo = (int*)malloc(n2Dprojs * sizeof(int));
	if (*photo == NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*imgptsSum = (int*)malloc(n3Dpts * sizeof(int));
	if (*imgptsSum == NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	rewind(fpp);
	readPointProjections(fpp, *imgpts, *photo, *imgptsSum, n3Dpts, n2Dprojs);

	fclose(fpp);
}
bool IBA::ba_parseArgs(int argc, char* argv[])
{
	int i;
	string param;
	bool bSuccess, bRKF;

	for (i = 1; i < argc; i++)
	{
		bSuccess = false;
		string name = argv[i];

		if (name[0] != '-') { // each param has to start with at least one dash
			return false;
		}

		string::size_type dashPos = name.find_first_not_of('-');
		if (dashPos != string::npos)
			name = name.substr(dashPos);

		if (strcmp(name.c_str(), "help") == 0)
		{
			ba_printHelp(pba);
			return false;
		}

		if (strcmp(name.c_str(), "cam") == 0)
		{
			i++;
			param = argv[i];
			m_szCameraInit = (char*)malloc(param.length());
			strcpy_s(m_szCameraInit,sizeof(m_szCameraInit), param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "fea") == 0)
		{
			i++;
			param = argv[i];
			m_szFeatures = (char*)malloc(param.length());
			strcpy_s(m_szFeatures, sizeof(m_szFeatures), param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "calib") == 0)
		{
			i++;
			param = argv[i];
			m_szCalibration = (char*)malloc(param.length());
			strcpy_s(m_szCalibration, sizeof(m_szCalibration), param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "pose") == 0)
		{
			i++;
			param = argv[i];
			m_szCamePose = (char*)malloc(param.length());
			strcpy_s(m_szCamePose, sizeof(m_szCamePose), param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "3D") == 0)
		{
			i++;
			param = argv[i];
			m_sz3Dpts = (char*)malloc(param.length());
			strcpy_s(m_sz3Dpts, sizeof(m_sz3Dpts), param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "report") == 0)
		{
			i++;
			param = argv[i];
			m_szReport = (char*)malloc(param.length());
			strcpy_s(m_szReport, sizeof(m_szReport), param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "xyz") == 0)
		{
			i++;
			param = argv[i];
			m_szXYZ = (char*)malloc(param.length());
			strcpy_s(m_szXYZ, sizeof(m_szXYZ), param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "i") == 0)
		{
			i++;
			param = argv[i];
			m_nMaxIter = atoi(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "robustKernel") == 0)
		{
			m_bRobustKernel = true;
			bRKF = false;
			i++;
			param = argv[i];
			if (strcmp(param.c_str(), "Huber") == 0)
			{
				m_nRobustType = 2;
				bRKF = true;
			}

			if (strcmp(param.c_str(), "Cauchy") == 0)
			{
				m_nRobustType = 1;
				bRKF = true;
			}

			bSuccess = true;

			if (!bRKF)
			{
				printf("BA: Must input right robust kernel function!\n");
				return false;
			}
		}

		if (strcmp(name.c_str(), "solve") == 0)
		{
			i++;
			param = argv[i];
			if (strcmp(param.c_str(), "LM") == 0)
				m_bsolverLM = true;

			if (strcmp(param.c_str(), "GN") == 0)
				m_bsolverGN = true;

			bSuccess = true;
		}

		if (strcmp(name.c_str(), "t") == 0)
		{
			i++;
			param = argv[i];

			m_Tau = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e1") == 0)
		{
			i++;
			param = argv[i];

			m_e1 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e2") == 0)
		{
			i++;
			param = argv[i];

			m_e2 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e3") == 0)
		{
			i++;
			param = argv[i];

			m_e3 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "e4") == 0)
		{
			i++;
			param = argv[i];

			m_e4 = atof(param.c_str());
			bSuccess = true;
		}

		if (strcmp(name.c_str(), "robustKernelWidth") == 0)
		{
			i++;
			param = argv[i];

			m_delt = atof(param.c_str());
			bSuccess = true;
		}

		if (!bSuccess)
		{
			printf("BA error: %s command is wrong!\n", name.c_str());
			return false;
		}

	}

	return true;

}

void IBA::ba_printHelp(BAType ba)
{
	string name;
	if (ba == pba)
		name = "Parallax";
	else if (ba == sba)
		name = "Sparse";
	printf(name.append("Sparse Bundle Adjustment General Options\n").c_str());
	printf("\n");

	printf("-cam			Provide initial camera pose.\n");
	printf("-fea			Provide features.\n");
	printf("-calib			Provide calibration.\n");
	printf("-xyz			Provide initial XYZ.\n");
	printf("-pose			Output optimal camera pose.\n");
	printf("-3D			Output optimal 3D point cloud.\n");
	printf("-report			Output report.\n");
	printf("-solve			Solve method including LevenbergMarquart(LM) and Gauss-Newton(GN).\n");
	printf("-i			Max Iteration.\n");
	printf("-t			LevenbergMarquart parameters.\n");
	printf("-robustKernel		use Cauchy Robust Kernel Function.\n");
	printf("-robustKernelWidth		width for the robust Kernel.\n");
}

