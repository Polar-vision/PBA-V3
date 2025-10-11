#pragma once
#include "BAExporter.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/StdVector>
#include <Eigen/Cholesky>

#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
	#include "ceres/ceres.h"
	#include "ceres/rotation.h"
#elif defined(_WIN32)
	#include "cholmod/x86/cholmod.h"
#else
    #error "Unsupported platform!"
#endif


#define PI  3.1415926535898 
#define MAXARCHOR 0.5
using namespace std;
using namespace Eigen;
#define MAXSTRLEN  2048 /* 2K */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
	while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}

////use sba_crsm structure from SBA (http://www.ics.forth.gr/~lourakis/sba/) to store sparse matrix
//struct sba_crsm
//{
//    int nr, nc;   //稀疏矩阵的行列
//    int nnz;      //非零元素个数
//    int* val;     //存储非零元素
//    int* colidx;  //非零元素的列号
//    int* rowptr;  //指向每行第一个非零元素的索引数组 (size: nr+1)
//};

class IBA
{
public:
    IBA(void);
    ~IBA(void);

	virtual bool ba_run(bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL,
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6)=0;

	virtual	bool ba_run(int argc, char** argv)=0;

	virtual bool ba_initialize(char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL)=0;

	bool ba_parseArgs(int argc, char* argv[]);
	int findNcameras(FILE* fp);
	void ba_readCameraPoseration(char* fname, double* ical);
	void readNpointsAndNprojections(FILE* fp, int* n3Dpts, int pnp, int* nprojs, int mnp);
	void ba_readCameraPose(FILE* fp, double* params, int* m_v);
	void ba_updateKR(double* KR, double* KdA, double* KdB, double* KdG, double* K, double* p);
	void ba_constructP(double* P, double* K, double* p);
	void ba_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams);
	int readNInts(FILE* fp, int* vals, int nvals);
	int readNDoubles(FILE* fp, double* vals, int nvals);
	int countNDoubles(FILE* fp);
	int skipNDoubles(FILE* fp, int nvals);
	void ba_readCablibration(FILE* fp, double* K);
	void readNpointsAndNprojectionsFromProj(FILE* fp, int& n3Dpts, int& nprojs);
	void readPointProjections(FILE* fp, double* imgpts, int* photo, int* imgptsSum, int n3Dpts, int n2Dprojs);
	void readImagePts(const char* szProj, double** imgpts, int** photo, int** imgptsSum, int& n3Dpts, int& n2Dprojs);
	void ba_printHelp(BAType ba);
	int		m_ncams, m_n3Dpts, m_n2Dprojs, m_nS, nc_;  //number of camera, 3D points, 2D projection points, non-zero element of S matrix
	int* m_archor;
	int* m_photo, * m_feature;
	double* m_motstruct, * m_imgpts;			  //6 camera pose and 3 feature parameters/PBA parameter,
	double* m_XYZ;								  //initial XYZ provided 	
	double* m_K;								  //calibration parameters
	//std::vector<double> m_K;
	int* m_V;//camera id
	char* m_vmask, * m_umask, * m_smask;
	int* m_imgptsSum, * m_struct, * m_pnt2main, * m_archorSort;
	double* m_KR, * m_KdA, * m_KdB, * m_KdG, * m_P;

	bool    m_bProvideXYZ, m_bFocal;
	char* m_szCameraInit;
	char* m_szFeatures;
	char* m_szCalibration;
	char* m_szXYZ;
	char* m_szCamePose;
	char* m_sz3Dpts;
	char* m_szReport;
	int		m_nMaxIter;
	double  m_Tau, m_e1, m_e2, m_e3, m_e4;
	bool	m_bRobustKernel;
	bool	m_bsolverLM;
	bool	m_bsolverGN;
	int		m_nRobustType;
	double  m_delt;
	int savePara = 0;
	//Solve Sparse Matrix using CHOLMOD (http://www.cise.ufl.edu/research/sparse/SuiteSparse/) 
#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
	struct Intrinsic{
		double fx, fy, cx, cy;
	};
	struct Camera{
		double euler_angle[3];
		double camera_center[3];
		int camidx;
	};
	struct Point3D{
		double xyz[3];
		double aep[3];//视差角参数化的三维点
		int nM;//主锚点
		int nA;//副锚点
	};
	struct Observation{
		int view_idx;
		double u, v;
	};
	struct Track{
		int nview;
		std::vector<Observation> obss;
	};
	std::vector<Intrinsic> intrs;
	std::vector<Camera> cams;
	std::vector<Point3D> points;
	std::vector<Track> tracks;
#elif defined(_WIN32)
	virtual bool ba_motstr_levmar()=0;
	virtual bool ba_motstr_gn(FixType ft = BA_FixDefault)=0;
	void sba_crsm_alloc(struct sba_crsm* sm, int nr, int nc, int nnz);
	void sba_crsm_free(struct sba_crsm* sm);
	int sba_crsm_elmidx(struct sba_crsm* sm, int i, int j);

	//void ba_readCameraPose_(char* fname);
	double nrmL2xmy(double* const e, const double* const x, const double* const y, const int n);
	void ba_saveTriangulatedxyz(const char* sz3Dpt, double* p);
	int ba_ConstructSmask(sba_crsm& Sidxij, sba_crsm& Uidxij);//分别是S矩阵和U矩阵的稀疏索引
	void ba_inverseVLM(double* V, double* IV, sba_crsm& Uidxij, double mu);
	void ba_inverseVGN(double* U, double* V, sba_crsm& Uidxij);
	double ba_computeInitialmu(double* U, double* V, sba_crsm& Uidxij, double tau, int nvars);
	
	void ba_solveFeatures(double* W, double* IV, double* ea, double* eb, double* dpa, double* dpb);
	bool ba_solveCholmodGN(int* Ap, int* Aii, bool init, bool ordering);
	bool ba_solveCholmodLM(int* Ap, int* Aii, bool init, bool ordering);
	//set CSS format
	void ba_constructCSSLM(int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init);
	void ba_constructCSSGN(int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init, int nft);
	//prepare for ordering of CHOLMOD
	void ba_constructAuxCSSLM(int* Ap, int* Aii);
	void ba_constructAuxCSSGN(int* Ap, int* Aii);
	//construct S matrix, S = U - W*V^-1*W^T
	void ba_constructSLM(double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij, double mu);
	void ba_constructSGN(double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij);
	void eulerAnglesToRotationMatrix(double* eulerAngles, double* R);
	void rotationMatrixToEulerAngles(double* R, double* eulerAngles);
	cholmod_sparse* m_cholSparseS;
	cholmod_factor* m_cholFactorS;
	cholmod_common m_cS;
	cholmod_dense* m_cholSparseR, * m_cholSparseE;
#else
	#error "Unsupported platform!"
#endif

	//
	DataType tpe = colmap;
};
