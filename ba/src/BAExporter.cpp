
#include "BAExporter.h"
#include "PBAImp.h"
#include "SBAImp.h"
// #include "BALPBAImp.h"
// #include "BALSBAImp.h"
BAType batype;
#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))

#elif defined(_WIN32)
	#include "stdafx.h"
	bool BAExporter::ba_motstr_levmar()
	{
		if (batype == pba || batype == sba)
			return ptr->ba_motstr_levmar();
		// else if (batype == balpba || batype == balsba)
		// 	return balptr->ba_motstr_levmar();
		else {
			std::cerr << "Unknown BA type in BAExporter::ba_run" << std::endl;
			return false;
		}
	}
	bool BAExporter::ba_motstr_gn(FixType ft)
	{
		if (batype == pba || batype == sba)
			return ptr->ba_motstr_gn(ft);
		// else if (batype == balpba || batype == balsba)
		// 	return balptr->ba_motstr_gn(ft);
		else {
			std::cerr << "Unknown BA type in BAExporter::ba_run" << std::endl;
			return false;
		}
	}
#else
    #error "Unsupported platform!"
#endif


BAExporter::BAExporter(BAType ba)
{
	batype = ba;
	if (ba == pba)
		ptr = new PBA;
	else if (ba == sba)
		ptr = new SBA;
	// else if (ba == balpba)
	// 	balptr = new BALPBA;
	// else if (ba == balsba)
	// 	balptr = new BALSBA;
}

BAExporter::~BAExporter(void)
{
	if (batype == pba || batype == sba)
		delete ptr;
	// else if (batype == balpba || batype == balsba)
	// 	delete balptr;
}

bool BAExporter::ba_run(bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ, char* szCalib, char* szReport,
	char* szPose, char* sz3D, double Tau)
{
	if (batype == pba || batype == sba)
		return ptr->ba_run(bRobust, bLM, nMaxIter, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D, Tau);
	// else if (batype == balpba || batype == balsba)
	// 	return balptr->ba_run(bRobust, bLM, nMaxIter, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D, Tau);
	else {
		std::cerr << "Unknown BA type in BAExporter::ba_run" << std::endl;
		return false;
	}
}
bool BAExporter::ba_run(int argc, char** argv)
{
	if (batype == pba || batype == sba)
		return ptr->ba_run(argc, argv);
	// else if (batype == balpba || batype == balsba)
	// 	return balptr->ba_run(argc, argv);
	else {
		std::cerr << "Unknown BA type in BAExporter::ba_run" << std::endl;
		return false;
	}
}
bool BAExporter::ba_initialize(char* szCamera, char* szFeature, char* szCalib, char* szXYZ)
{
	if (batype == pba || batype == sba)
		return	ptr->ba_initialize(szCamera, szFeature, szCalib, szXYZ);
	// else if (batype == balpba || batype == balsba)
	// 	return	balptr->ba_initialize(szCamera, szFeature, szCalib, szXYZ);
	else {
		std::cerr << "Unknown BA type in BAExporter::ba_run" << std::endl;
		return false;
	}
}
