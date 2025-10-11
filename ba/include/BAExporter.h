#include <cstddef>
#include <iostream>
#ifndef _BA_H
#define _BA_H

#ifdef BA_EXPORTS
    #define BAapi __declspec(dllexport)
#elif defined(BA_STATIC)
    #define BAapi  // 静态库无修饰符
#else
    #define BAapi __declspec(dllimport)
#endif


//use sba_crsm structure from SBA (http://www.ics.forth.gr/~lourakis/sba/) to store sparse matrix
struct sba_crsm
{
	int nr, nc;   //稀疏矩阵的行列
	int nnz;      //非零元素个数
	int* val;     //存储非零元素
	int* colidx;  //非零元素的列号
	int* rowptr;  //指向每行第一个非零元素的索引数组 (size: nr+1)
};

typedef enum BA_FixType
{	
	BA_FixDefault = -1, //fix the longest axis  
	BA_FixX = 0 ,	    // fix X axis	
	BA_FixY = 1 ,	    // fix Y axis
	BA_FixZ = 2 ,	    // fix Z axis
} FixType ; 

typedef enum BA_Type
{
	pba = 1,
	sba = 2,
	balpba = 3,
    balsba = 4,
}BAType;

typedef enum Data_Type
{
	bal = 1,
	colmap = 2, 
}DataType;

class IBA;//PBA和SBA的基类
// class IBALBA;//BALPBA和BALSBA的基类
class BAapi BAExporter
{
public:
	bool ba_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 );

	bool ba_run( int argc, char** argv );

	bool ba_initialize( char* szCamera, char* szFeature, char* szCalib =  NULL, char* szXYZ = NULL );

#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))

#elif defined(_WIN32)
	bool ba_motstr_levmar( );
	bool ba_motstr_gn( FixType ft = BA_FixDefault );
#else
    #error "Unsupported platform!"
#endif

	BAExporter(BAType);
	~BAExporter(void);

private:
	IBA* ptr;
	// IBALBA* balptr;
};

#endif