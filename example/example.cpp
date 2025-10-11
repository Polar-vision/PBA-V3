#include "stdafx.h"
#include "BAExporter.h"
#include "dataPath.h"
#include <vector>
#include <map>
#include <string>
#include <filesystem>
using namespace std;

string get_current_dir(){
	return filesystem::path(__FILE__).parent_path().string();
}
string get_dataset_path(){
	filesystem::path current_dir=get_current_dir();
	filesystem::path dataset_path =current_dir.parent_path() / "datasets/";
	return dataset_path.string();
}

int main(int argc, char* argv[] )
{
	argv[1] = const_cast<char*>(dt);
	argv[2] = "pba_lm";
	bool bLM = true;   //true is Levenberg-Marquardt
	BAType ba = BAType::pba;
	string dt = argv[1];
	string ba_ = argv[2];
	string pCheck;
	if (ba_ == "pba_lm")
	{ 
		ba = BAType::pba;
		bLM = true;
		pCheck = "PBA-LM";
	}
	if (ba_ == "sba_lm")
	{
		ba = BAType::sba;
		bLM = true;
		pCheck = "SBA-LM";
	}
	if (ba_ == "pba_gn")
	{
		ba = BAType::pba;
		bLM = false;
		pCheck = "PBA-GN";
	}
	if (ba_ == "sba_gn")
	{
		ba = BAType::sba;
		bLM = false;
		pCheck = "SBA-GN";
	}

	string zsw(dt);
	string data_file = get_dataset_path()+zsw;
	size_t pos = data_file.find_last_of("/\\");
	string parentPath = (pos != string::npos) ? data_file.substr(0, pos) : "";
	string pP = parentPath + "/";
	string p1 = pP + "Cam.txt";
	string p2 = pP + "Feature.txt";
	string p3 = pP + "XYZ.txt";
	string p4 = pP + "cal.txt";
	//string p1 = pP + "Cam_new.txt";
	//string p2 = pP + "Feature_new.txt";
	//string p3 = pP + "XYZ_new.txt";
	//string p4 = pP + "cal_new.txt";
	//string p1 = pP + "PBA-LM-FinalPose.txt";
	//string p3 = pP + "PBA-LM-Final3D.txt";// "XYZ1.txt";
	string pReport = "-report.txt";
	string pPose = "-FinalPose.txt";
	string p3D = "-Final3D.ply";
	string p5, p6, p7, pInit3D;
	pInit3D = pP + "XYZ.ply";

/*������BAL��98�����ݼ����ԣ�https://grail.cs.washington.edu/projects/bal/ */

	//if (ba == BAType::sba && bLM)
	//	pCheck = "SBA-LM";
	//if (ba == BAType::sba && !bLM)
	//	pCheck = "SBA-GN";
	//if (ba == BAType::pba && bLM)
	//	pCheck = "PBA-LM";
	//if (ba == BAType::pba && !bLM)
	//	pCheck = "PBA-GN";
	p5 = pP + pCheck + pReport;
	p6 = pP + pCheck + pPose;
	p7 = pP + pCheck + p3D;

	char* szCam = const_cast<char*>(p1.c_str());
	char* szFea = const_cast<char*>(p2.c_str());
	//char* szXYZ = NULL;
	char* szXYZ = const_cast<char*>(p3.c_str());
	char* szCalib = const_cast<char*>(p4.c_str());
	char* szReport = const_cast<char*>(p5.c_str()); 
	char* szPose = const_cast<char*>(p6.c_str()); 
	char* sz3D = const_cast<char*>(p7.c_str()); 


	printf("%s\n",szXYZ);
	BAExporter BA(ba);
	BA.ba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*������BAL��98�����ݼ�����*/


	//BeiJing dataset
	//BAExporter BA(pba);
	//char* szCam = "../../data/BeiJing/Cam40.txt";
	//char* szFea = "../../data/BeiJing/Feature40.txt";
	//char* szCalib = "../../data/BeiJing/cal40.txt";
	//char* szXYZ = "../../data/BeiJing/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/BeiJing/report.txt";
	//char* szPose = "../../data/BeiJing/FinalPose.txt";
	//char* sz3D = "../../data/BeiJing/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)       PBA(GN)          SBA(LM)       PBA(LM)
#Iteration          4             4                11             8
Initial MSE      232.408        232.408          232.408       232.408
Final MSE        0.210473       0.210473         0.210473      0.210473
Runtime          0.111          0.135            0.244         0.208
Reason              2             2                 2             2
*/

	//Toronto dataset
	//char* szCam = "../../data/Toronto/Cam13.txt";
	//char* szFea = "../../data/Toronto/Feature13.txt";
	//char* szCalib = "../../data/Toronto/cal13.txt";
	//char* szXYZ = "../../data/Toronto/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/Toronto/report.txt";
	//char* szPose = "../../data/Toronto/FinalPose.txt";
	//char* sz3D = "../../data/Toronto/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 500, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
				   SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            9                      9                    487                   18
Initial MSE    3326138.80221379       3326138.80221379      3326138.80221379      3326138.80221379
Final MSE         0.041856               0.041856              50.083779             0.041856
Runtime           0.717                  0.652                 37.665                1.178
Reason                2                      2                     2                     2
*/

	//WuLong dataset
	//BAExporter BA(sba);
	//char* szCam = "../../data/WuLong/Cam42.txt";
	//char* szFea = "../../data/WuLong/Feature42.txt";
	//char* szCalib = "../../data/WuLong/cal42.txt";
	//char* szXYZ = "../../data/WuLong/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/WuLong/report.txt";
	//char* szPose = "../../data/WuLong/FinalPose.txt";
	//char* sz3D = "../../data/WuLong/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       5                    321                   13
Initial MSE       16.3765                  16.3765               16.3765              16.3765
Final MSE             -                    0.277196              0.277196             0.277196
Runtime               -                    0.139                 6.796                0.295
Reason                -                       2                     2                     2
*/

	//Taian dataset
	//char* szCam = "../../data/Taian/Cam737.txt";
	//char* szFea = "../../data/Taian/Feature737.txt";
	//char* szCalib = "../../data/Taian/cal737.txt";
	//char* szXYZ = "../../data/Taian/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/Taian/report.txt";
	//char* szPose = "../../data/Taian/FinalPose.txt";
	//char* sz3D = "../../data/Taian/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 500, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       11                    427                  61
Initial MSE    747896.41362525        747896.41362525        747896.41362525       747896.41362525
Final MSE             -                    0.038748              5.895712             0.038748
Runtime               -                    20.394                 960.195              123.876
Reason                -                       2                     2                     2
*/

	//DunHuang dataset
	//char* szCam = "../../data/DunHuang/Cam63.txt";
	//char* szFea = "../../data/DunHuang/Feature63.txt";
	//char* szCalib = "../../data/DunHuang/cal63.txt";
	//char* szXYZ = "../../data/DunHuang/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/DunHuang/report.txt";
	//char* szPose = "../../data/DunHuang/FinalPose.txt";
	//char* sz3D = "../../data/DunHuang/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 2000, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            8                       7                    651                  11
Initial MSE    261763.35960011        261763.35960011        261763.35960011       261763.35960011
Final MSE          0.168251                0.168251              0.793898             0.168251
Runtime            1.631                   1.411                 128.398              1.939
Reason                2                       2                     2                     2
*/

	//Vaihingen dataset
	//char* szCam = "../../data/Vaihingen/Cam20.txt";
	//char* szFea = "../../data/Vaihingen/Feature20.txt";
	//char* szCalib = "../../data/Vaihingen/cal20.txt";
	//char* szXYZ = "../../data/Vaihingen/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szPose = "../../data/Vaihingen/FinalPose.txt";
	//char* sz3D = "../../data/Vaihingen/Final3D.txt";
	//char* szReport = "../../data/Vaihingen/report.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 500, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            6                       6                    36                    13
Initial MSE    142089.50626664        142089.50626664        142089.50626664       142089.50626664
Final MSE          0.119365                0.119365              0.119365             0.119365
Runtime            2.560                   2.383                 14.319               4.475
Reason                2                       2                     2                     2
*/

	//College dataset
	//char* szCam = "../../data/College/Cam468.txt";
	//char* szFea = "../../data/College/Feature468.txt";
	//char* szCalib = "../../data/College/cal468.txt";
	//char* szXYZ = "../../data/College/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/College/report.txt";
	//char* szPose = "../../data/College/FinalPose.txt";
	//char* sz3D = "../../data/College/Final3D.txt";
	//bool bLM = false;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       12                   1001                  17
Initial MSE    202329.44856826         202329.44856826        202329.44856826       202329.44856826
Final MSE             -                   0.734738               16.598062            0.734738
Runtime               -                   12.240                 1765.791             16.395
Reason                -                       2                     3                     2
*/

	//Village dataset
	//BAExporter BA(pba);
	//char* szCam = "../../data/Village/Cam90.txt";
	//char* szFea = "../../data/Village/Feature90.txt";
	//char* szCalib = "../../data/Village/cal90.txt";
 //   char* szXYZ = "../../data/Village/XYZ.txt";
 //   //char* szXYZ = NULL;
	//char* szReport = "../../data/Village/report.txt";	
 //   char* szPose = "../../data/Village/FinalPose.txt";
 //   char* sz3D = "../../data/Village/Final3D.txt";
	//bool bLM     = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            6                       6                    15                    11
Initial MSE    28170.98109754           28170.98109754          28170.98109754       28170.98109754
Final MSE          0.083716                 0.083716              0.083716            0.083716
Runtime            1.653                    1.653                 3.751               2.606
Reason                2                       2                     2                     2
*/

	//NewCollege dataset
	//char* szCam = "../../data/NewCollege/Cam3500.txt";
	//char* szFea = "../../data/NewCollege/Feature3500.txt";
	//char* szCalib = "../../data/NewCollege/cal3500.txt";
	//char* szXYZ = "../../data/NewCollege/XYZ.txt";
	////char* szXYZ = NULL;
	//char* szReport = "../../data/NewCollege/report.txt";
	//char* szPose = "../../data/NewCollege/FinalPose.txt";
	//char* sz3D = "../../data/NewCollege/Final3D.txt";
	//bool bLM = true;   //true is Levenberg-Marquardt
	//BA.ba_run(false, bLM, 1000, szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D);

/*
					SBA(GN)                 PBA(GN)               SBA(LM)              PBA(LM)
#Iteration            -                       -                    1001                  33
Initial MSE      295.79746736           295.79746736           295.79746736          295.79746736
Final MSE             -                       -                  0.991214              1.090911
Runtime               -                       -                  3692.488              92.964
Reason                -                       -                     3                     1
*/


	return 0;
}

