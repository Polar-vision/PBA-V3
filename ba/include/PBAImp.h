#pragma once
#include "BAExporter.h"
#include "IBA.h"
#include <map>
#include <vector>
#include <set>
using namespace std;
//using namespace Eigen;

class PBA : public IBA
{
public:
	PBA(void);
	~PBA(void);
	virtual bool ba_run(bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL,
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6);

	virtual	bool ba_run( int argc, char** argv );

	virtual bool ba_initialize( char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL );

	void	pba_readAndInitialize( char *camsfname, char *ptsfname, char *calibfname, int *ncams, int *n3Dpts, int *n2Dprojs,double **motstruct, 
				double **imgpts, int **archor, char **vmask, char **umask,int **nphoto, int** nfeature, int** archorSort );
	void	pba_readProjectionAndInitilizeFeature(	FILE *fp, double *params, double *projs, char *vmask, int ncams, 
				int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort );
	//initialize feature. You can provide xyz or system also provide them by itself;
	bool    pba_initializeMainArchor( double* imgpts, double* camera,double* K,double* feature, int nP, int FID, double* KR );
	bool    pba_initializeAssoArchor( double* imgpts, int* photo, double* camera,double* K,double* feature,int nMI, int nAI, int FID, bool bLast );
	bool	pba_initializeOtheArchors( double* imgpts, int* photo, double* camera,double* K,double* feature,int* archorSort,int nfeacout, int nOI, int FID );
	int zu = 0;

	#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
		//SBA残差方程
		struct SBA_ReprojectionError {
		SBA_ReprojectionError(double observed_u, double observed_v,
						double fx, double fy, double cx, double cy)
			: u(observed_u), v(observed_v),
			fx(fx), fy(fy), cx(cx), cy(cy) {}

		mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
		template <typename T>
		bool operator()(const T* const euler_angles,
						const T* const camera_center,
						const T* const point3D,
						T* residuals) const {
			// T matR[9];
			// T eu[3];
			// eu[0] = euler_angles[0] * T(180) / T(PI);
			// eu[1] = euler_angles[1] * T(180) / T(PI);
			// eu[2] = euler_angles[2] * T(180) / T(PI);
			// ceres::EulerAnglesToRotationMatrix(eu, 3, matR);//以度为单位
			T R[9];
			T c0, c1, c2, s0, s1, s2;
			c0 = cos(euler_angles[0]);
			c1 = cos(euler_angles[1]);
			c2 = cos(euler_angles[2]);
			s0 = sin(euler_angles[0]);
			s1 = sin(euler_angles[1]);
			s2 = sin(euler_angles[2]);
			R[0] = c1 * c0;
			R[1] = c1 * s0;
			R[2] = -s1;
			R[3] = s2 * s1 * c0 - c2 * s0;
			R[4] = s2 * s1 * s0 + c2 * c0;
			R[5] = s2 * c1;
			R[6] = c2 * s1 * c0 + s2 * s0;
			R[7] = c2 * s1 * s0 - s2 * c0;
			R[8] = c2 * c1;
			// R[0] = cos(euler_angles[1]) * cos(euler_angles[0]);
			// R[1] = cos(euler_angles[1]) * sin(euler_angles[0]);
			// R[2] = -sin(euler_angles[1]);
			// R[3] = sin(euler_angles[2]) * sin(euler_angles[1]) * cos(euler_angles[0]) - cos(euler_angles[2]) * sin(euler_angles[0]);
			// R[4] = sin(euler_angles[2]) * sin(euler_angles[1]) * sin(euler_angles[0]) + cos(euler_angles[2]) * cos(euler_angles[0]);
			// R[5] = sin(euler_angles[2]) * cos(euler_angles[1]);
			// R[6] = cos(euler_angles[2]) * sin(euler_angles[1]) * cos(euler_angles[0]) + sin(euler_angles[2]) * sin(euler_angles[0]);
			// R[7] = cos(euler_angles[2]) * sin(euler_angles[1]) * sin(euler_angles[0]) - sin(euler_angles[2]) * cos(euler_angles[0]);
			// R[8] = cos(euler_angles[2]) * cos(euler_angles[1]);
			// X - C
			T Xc[3];
			Xc[0] = point3D[0] - camera_center[0];
			Xc[1] = point3D[1] - camera_center[1];
			Xc[2] = point3D[2] - camera_center[2];

			// P = R * (X - C)
			T p[3];
			p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
			p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
			p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

			// Normalize
			T xp = p[0] / p[2];
			T yp = p[1] / p[2];

			// Project
			T predicted_u = fx * xp + cx;
			T predicted_v = fy * yp + cy;

			// Residual
			residuals[0] = predicted_u - T(u);
			residuals[1] = predicted_v - T(v);

			// if (!printed_once_) {
				// std::cout << "Euler: " << euler_angles[0] << ", "
				// 					<< euler_angles[1] << ", "
				// 					<< euler_angles[2] << std::endl;
				// std::cout << "Camera center: " << camera_center[0] << ", "
				// 							<< camera_center[1] << ", "
				// 							<< camera_center[2] << std::endl;
				// std::cout << "Point3D: " << point3D[0] << ", "
				// 						<< point3D[1] << ", "
				// 						<< point3D[2] << std::endl;
				// printed_once_ = true;
				// std::cout << predicted_u << ", " << u << std::endl;
				// std::cout << predicted_v << ", " << v << std::endl;
				// std::cout << residuals[0] << ", " << residuals[1] << std::endl;
			// }
			return true;
		}

		static ceres::CostFunction* Create(double u, double v,
										double fx, double fy, double cx, double cy) {	
			// printf("%s\n","ReprojectionError entered");	
			return new ceres::AutoDiffCostFunction<SBA_ReprojectionError, 2, 3, 3, 3>(
				new SBA_ReprojectionError(u, v, fx, fy, cx, cy));
		}

		double u, v;
		double fx, fy, cx, cy;
		};


		//PBA_nM残差方程
		struct PBA_nM_ReprojectionError {
		PBA_nM_ReprojectionError(double observed_u, double observed_v,
						double fx, double fy, double cx, double cy)
			: u(observed_u), v(observed_v),
			fx(fx), fy(fy), cx(cx), cy(cy) {}

		mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
		template <typename T>
		bool operator()(const T* const euler_angles,
						const T* const point3D,
						T* residuals) const {
			T R[9];
			T c0, c1, c2, s0, s1, s2;
			c0 = cos(euler_angles[0]);
			c1 = cos(euler_angles[1]);
			c2 = cos(euler_angles[2]);
			s0 = sin(euler_angles[0]);
			s1 = sin(euler_angles[1]);
			s2 = sin(euler_angles[2]);
			R[0] = c1 * c0;
			R[1] = c1 * s0;
			R[2] = -s1;
			R[3] = s2 * s1 * c0 - c2 * s0;
			R[4] = s2 * s1 * s0 + c2 * c0;
			R[5] = s2 * c1;
			R[6] = c2 * s1 * c0 + s2 * s0;
			R[7] = c2 * s1 * s0 - s2 * c0;
			R[8] = c2 * c1;

			T p[3];
			T ptXj[3];
			ptXj[0] = sin(point3D[0]) * cos(point3D[1]);
			ptXj[1] = sin(point3D[1]);
			ptXj[2] = cos(point3D[0]) * cos(point3D[1]);
			p[0] = R[0] * ptXj[0] + R[1] * ptXj[1] + R[2] * ptXj[2];
			p[1] = R[3] * ptXj[0] + R[4] * ptXj[1] + R[5] * ptXj[2];
			p[2] = R[6] * ptXj[0] + R[7] * ptXj[1] + R[8] * ptXj[2];

			// Normalize
			T xp = p[0] / p[2];
			T yp = p[1] / p[2];

			// Project
			T predicted_u = fx * xp + cx;
			T predicted_v = fy * yp + cy;

			// Residual
			residuals[0] = predicted_u - T(u);
			residuals[1] = predicted_v - T(v);

				// std::cout << "Euler: " << euler_angles[0] << ", "
				// 					<< euler_angles[1] << ", "
				// 					<< euler_angles[2] << std::endl;
				// std::cout << "Camera center: " << camera_center[0] << ", "
				// 							<< camera_center[1] << ", "
				// 							<< camera_center[2] << std::endl;
				// std::cout << "Point3D: " << point3D[0] << ", "
				// 						<< point3D[1] << ", "
				// 						<< point3D[2] << std::endl;
				// printed_once_ = true;
				// std::cout << predicted_u << ", " << u << std::endl;
				// std::cout << predicted_v << ", " << v << std::endl;
				// std::cout << residuals[0] << ", " << residuals[1] << std::endl;
			// }
			return true;
		}

		static ceres::CostFunction* Create(double u, double v,
										double fx, double fy, double cx, double cy) {	
			// printf("%s\n","ReprojectionError entered");	
			return new ceres::AutoDiffCostFunction<PBA_nM_ReprojectionError, 2, 3, 3>(
				new PBA_nM_ReprojectionError(u, v, fx, fy, cx, cy));
		}

		double u, v;
		double fx, fy, cx, cy;
		};



		//副锚点残差方程
		struct PBA_nA_ReprojectionError {
		PBA_nA_ReprojectionError(double observed_u, double observed_v,
						double fx, double fy, double cx, double cy)
			: u(observed_u), v(observed_v),
			fx(fx), fy(fy), cx(cx), cy(cy) {}

		mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
		template <typename T>
		bool operator()(const T* const euler_angles,
						const T* const camera_center_nP,
						const T* const camera_center_nM,
						const T* const point3D,
						T* residuals) const {
			T R[9];
			T c0, c1, c2, s0, s1, s2;
			c0 = cos(euler_angles[0]);
			c1 = cos(euler_angles[1]);
			c2 = cos(euler_angles[2]);
			s0 = sin(euler_angles[0]);
			s1 = sin(euler_angles[1]);
			s2 = sin(euler_angles[2]);
			R[0] = c1 * c0;
			R[1] = c1 * s0;
			R[2] = -s1;
			R[3] = s2 * s1 * c0 - c2 * s0;
			R[4] = s2 * s1 * s0 + c2 * c0;
			R[5] = s2 * c1;
			R[6] = c2 * s1 * c0 + s2 * s0;
			R[7] = c2 * s1 * s0 - s2 * c0;
			R[8] = c2 * c1;

			T p[3], pti2k[3], ptXUnit[3], ptXk[3];
			//主锚点到副锚点的平移向量
			pti2k[0] = camera_center_nP[0] - camera_center_nM[0];	
			pti2k[1] = camera_center_nP[1] - camera_center_nM[1];	
			pti2k[2] = camera_center_nP[2] - camera_center_nM[2];	

			//主锚点到特征点的单位向量
			ptXUnit[0] = sin(point3D[0]) * cos(point3D[1]);
			ptXUnit[1] = sin(point3D[1]);
			ptXUnit[2] = cos(point3D[0]) * cos(point3D[1]);

			//compute angle w2
			T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
			T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
			T dW2;
			if (dDot/dDisi2k > T(1))
				dW2 = T(0);
			if (dDot/dDisi2k < T(-1))
				dW2 = T(PI);
			else
				dW2  = acos(dDot/dDisi2k);

			//compute Xk vector according sin theory
			ptXk[0] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[0] - sin(point3D[2]) * pti2k[0];
			ptXk[1] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[1] - sin(point3D[2]) * pti2k[1];
			ptXk[2] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[2] - sin(point3D[2]) * pti2k[2];

			p[0] = R[0] * ptXk[0] + R[1] * ptXk[1] + R[2] * ptXk[2];
			p[1] = R[3] * ptXk[0] + R[4] * ptXk[1] + R[5] * ptXk[2];
			p[2] = R[6] * ptXk[0] + R[7] * ptXk[1] + R[8] * ptXk[2];

			// Normalize
			T xp = p[0] / p[2];
			T yp = p[1] / p[2];

			// Project
			T predicted_u = fx * xp + cx;
			T predicted_v = fy * yp + cy;

			// Residual
			residuals[0] = predicted_u - T(u);
			residuals[1] = predicted_v - T(v);
			return true;
		}

		static ceres::CostFunction* Create(double u, double v,
										double fx, double fy, double cx, double cy) {	
			// printf("%s\n","ReprojectionError entered");	
			return new ceres::AutoDiffCostFunction<PBA_nA_ReprojectionError, 2, 3, 3, 3, 3>(
				new PBA_nA_ReprojectionError(u, v, fx, fy, cx, cy));
		}

		double u, v;
		double fx, fy, cx, cy;
		};


		//其它锚点残差方程
		struct PBA_nP_ReprojectionError {
		PBA_nP_ReprojectionError(double observed_u, double observed_v,
						double fx, double fy, double cx, double cy)
			: u(observed_u), v(observed_v),
			fx(fx), fy(fy), cx(cx), cy(cy) {}

		mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
		template <typename T>
		bool operator()(const T* const euler_angles,
						const T* const camera_center_nP,
						const T* const camera_center_nM,
						const T* const camera_center_nA,
						const T* const point3D,
						T* residuals) const {
			T R[9];
			T c0, c1, c2, s0, s1, s2;
			c0 = cos(euler_angles[0]);
			c1 = cos(euler_angles[1]);
			c2 = cos(euler_angles[2]);
			s0 = sin(euler_angles[0]);
			s1 = sin(euler_angles[1]);
			s2 = sin(euler_angles[2]);
			R[0] = c1 * c0;
			R[1] = c1 * s0;
			R[2] = -s1;
			R[3] = s2 * s1 * c0 - c2 * s0;
			R[4] = s2 * s1 * s0 + c2 * c0;
			R[5] = s2 * c1;
			R[6] = c2 * s1 * c0 + s2 * s0;
			R[7] = c2 * s1 * s0 - s2 * c0;
			R[8] = c2 * c1;

			T p[3], pti2k[3], pti2l[3], ptXUnit[3], ptXk[3];

			//主锚点到副锚点的平移向量
			pti2k[0] = camera_center_nA[0] - camera_center_nM[0];		
			pti2k[1] = camera_center_nA[1] - camera_center_nM[1];		
			pti2k[2] = camera_center_nA[2] - camera_center_nM[2];
			//主锚点到当且锚点的平移向量
			pti2l[0] = camera_center_nP[0] - camera_center_nM[0];	
			pti2l[1] = camera_center_nP[1] - camera_center_nM[1];		
			pti2l[2] = camera_center_nP[2] - camera_center_nM[2];
			
			//XUnit 主锚点到特征点的单位向量
			ptXUnit[0] = sin(point3D[0]) * cos(point3D[1]);
			ptXUnit[1] = sin(point3D[1]);
			ptXUnit[2] = cos(point3D[0]) * cos(point3D[1]);

			//compute angle w2
			T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
			T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
			T dW2;
			//dW2  = acos( dDot/dDisi2k );
			if (dDot/dDisi2k > T(1))
				dW2 = T(0);
			if ( dDot/dDisi2k < T(-1))
				dW2 = T(PI);
			else
				dW2  = acos(dDot/dDisi2k);

			//compute Xl vector according sin theory
			ptXk[0] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[0] - sin(point3D[2]) * pti2l[0];
			ptXk[1] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[1] - sin(point3D[2]) * pti2l[1];
			ptXk[2] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[2] - sin(point3D[2]) * pti2l[2];
			
			p[0] = R[0] * ptXk[0] + R[1] * ptXk[1] + R[2] * ptXk[2];
			p[1] = R[3] * ptXk[0] + R[4] * ptXk[1] + R[5] * ptXk[2];
			p[2] = R[6] * ptXk[0] + R[7] * ptXk[1] + R[8] * ptXk[2];

			// Normalize
			T xp = p[0] / p[2];
			T yp = p[1] / p[2];

			// Project
			T predicted_u = fx * xp + cx;
			T predicted_v = fy * yp + cy;

			// Residual
			residuals[0] = predicted_u - T(u);
			residuals[1] = predicted_v - T(v);
			
			return true;
		}

		static ceres::CostFunction* Create(double u, double v,
										double fx, double fy, double cx, double cy) {	
			// printf("%s\n","ReprojectionError entered");	
			return new ceres::AutoDiffCostFunction<PBA_nP_ReprojectionError, 2, 3, 3, 3, 3, 3>(
				new PBA_nP_ReprojectionError(u, v, fx, fy, cx, cy));
		}

		double u, v;
		double fx, fy, cx, cy;
		};
	#elif defined(_WIN32)
		public:
			virtual bool ba_motstr_levmar( );
			virtual bool ba_motstr_gn(  FixType ft = BA_FixDefault  );
		private:
			//compute reprojection error
			void	pba_cost(double *p, double *hx, int* archor );
			void pba_reprojectEachPts(double* KR, double* pa, double* pb, int nM, int nN, int nP, double n[2]);
			//compute Jacobian, not save Jp, Jc etc
			void	pba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );
			void	pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );
			void pba_jacobianEachPts(double* KR, double* KdA, double* KdB, double* KdG, double* pa, double* ppt, int nM, int nN, int nP, double* pAM, double* pAA, double* pPA, double* pPB);

			//Add By Zuo
			bool pba_initializeOtheArchors_Mindw(double* imgpts, int* photo, double* camera, double* K, double* feature, int* archorSort, int nfeacout, int nOI, int FID);
			void pba_saveInitialXYZ(const char* sz3Dpt, double* p);
			void pba_saveInitialParallax(const char* sz3Dpt, double* p);
			double* pba_angle2xyz(double* p);

			//transform angle into XYZ
			int		pba_angle2xytGN(double* p);
			int		pba_angle2xytLM(double* p);
			void	pba_saveXYZ(const char* camera, const char* sz3Dpt, double* p, bool gn = true);
			
	#else
		#error "Unsupported platform!"
	#endif

};

