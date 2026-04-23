
#include "PBAImp.h"

#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))

#elif defined(_WIN32)
	#include "stdafx.h"
	//compute reprojection image coordinates of each image points
	void PBA::pba_reprojectEachPts(double* KR, double* pa, double* pb, int nM, int nN, int nP, double n[2])
	{
		double *posei, *posek;
		double pti2k[3];
		double ptXUnit[3];
		double dDot, dDisi2k, dW2;
		double ptXk[3];
		double *posel;
		double pti2l[3];
		double ptXj[3];
		double *pKR;
		int bs = 1;
		if (tpe == colmap) {
			bs = 1;
		}
		else if (tpe == bal) {
			bs = -1;
		}
		if ( nP == nM )	//Main archor
		{
			ptXj[0] = sin( pb[0] ) * cos( pb[1] );
			ptXj[1] = sin( pb[1] );
			ptXj[2] = cos( pb[0] ) * cos( pb[1] );

			pKR = KR + nP*9;

			n[0] = bs*(pKR[0] * ptXj[0] + pKR[1] * ptXj[1] + pKR[2] * ptXj[2]) /
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			n[1] = bs*(pKR[3] * ptXj[0] + pKR[4] * ptXj[1] + pKR[5] * ptXj[2]) /
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);
		}else 
			if ( nP == nN )	//Associate archor
			{
				posei = pa+6*nM+3;
				posek = pa+6*nN+3;

				pti2k[0] = posek[0]-posei[0];	pti2k[1] = posek[1]-posei[1];	pti2k[2] = posek[2]-posei[2];	

				ptXUnit[0] = sin( pb[0] ) * cos( pb[1] );
				ptXUnit[1] = sin( pb[1] );
				ptXUnit[2] = cos( pb[0] ) * cos( pb[1] );

				//compute angle w2
				dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
				dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
				if (dDot/dDisi2k > 1)
					dW2 = 0;
				if ( dDot/dDisi2k < -1)
					dW2 = PI;
				else
					dW2  = acos( dDot/dDisi2k );

				//compute Xk vector according sin theory
				ptXk[0] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[0] - sin( pb[2] ) * pti2k[0];
				ptXk[1] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[1] - sin( pb[2] ) * pti2k[1];
				ptXk[2] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[2] - sin( pb[2] ) * pti2k[2];

				pKR = KR + nN*9;	

				n[0] = bs*(pKR[0] * ptXk[0] + pKR[1] * ptXk[1] + pKR[2] * ptXk[2]) /
					(pKR[6] * ptXk[0] + pKR[7] * ptXk[1] + pKR[8] * ptXk[2]);

				n[1] = bs*(pKR[3] * ptXk[0] + pKR[4] * ptXk[1] + pKR[5] * ptXk[2]) /
					(pKR[6] * ptXk[0] + pKR[7] * ptXk[1] + pKR[8] * ptXk[2]);
			}
			else
			{
				posel = pa + nP*6 + 3;
				posek = pa + nN *6 + 3;
				posei = pa + nM *6 + 3;

				pti2k[0] = posek[0] - posei[0];		pti2k[1] = posek[1] - posei[1];		pti2k[2] = posek[2] - posei[2];
				pti2l[0] = posel[0] - posei[0];		pti2l[1] = posel[1] - posei[1];		pti2l[2] = posel[2] - posei[2];
				
				ptXUnit[0] = sin( pb[0] ) * cos( pb[1] );
				ptXUnit[1] = sin( pb[1] );
				ptXUnit[2] = cos( pb[0] ) * cos( pb[1] );

				//compute angle w2
				dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
				dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
				//dW2  = acos( dDot/dDisi2k );
				if (dDot/dDisi2k > 1)
					dW2 = 0;
				if ( dDot/dDisi2k < -1)
					dW2 = PI;
				else
					dW2  = acos( dDot/dDisi2k );

				//compute Xl vector according sin theory
				ptXk[0] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[0] - sin( pb[2] ) * pti2l[0];
				ptXk[1] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[1] - sin( pb[2] ) * pti2l[1];
				ptXk[2] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[2] - sin( pb[2] ) * pti2l[2];

				pKR = KR + nP*9;	

				n[0] = bs*(pKR[0] * ptXk[0] + pKR[1] * ptXk[1] + pKR[2] * ptXk[2])/
					(pKR[6] * ptXk[0] + pKR[7] * ptXk[1] + pKR[8] * ptXk[2]);

				n[1] = bs*(pKR[3] * ptXk[0] + pKR[4] * ptXk[1] + pKR[5] * ptXk[2]) /
					(pKR[6] * ptXk[0] + pKR[7] * ptXk[1] + pKR[8] * ptXk[2]);
			}
	}

	//compute Jacobian for each image point
	void PBA::pba_jacobianEachPts(double* KR, double* KdA, double *KdB, double* KdG, double* pa, double* ppt, int nM, int nN, int nP,double* pAM, double* pAA, double* pPA, double* pPB)
	{
		double matXj[3];
		double matxyt[3];
		double matDuvDxyt[6];
		double matDxytDRA[3], matDxytDRB[3], matDxytDRG[3];
		double matDuvDRA[2], matDuvDRB[2], matDuvDRG[2];
		double matDXjDDH[6];
		double matDxytDDH[6];
		double matDuvDDH[4];
		double *ptAngle = NULL;

		double *posei, *posek, *posel;
		double matPosei2k[3], matPosei2l[3];
		double	dDisTi2k, dDot, dArcCosIn, dW;

		double matXk[3],matXl[3];
		double matDDotDDH[2];
		double matDArcCosInDDH[2];
		double dDWDArcCosIn;
		double matDsinwWDDH[2];
		double matDXkDpa[3];
		double matDuvDpa[2];
		double tmp1[6];
		double tmp2[6];
		
		double matDdisDTi[3];
		double matDDotDTi[3];
		double matDArcCosInDTi[3];
		double matDsinwWDTi[3];
		double matDTi2kDTi[9], matDTi2kDTk[9];
		double matDXkDTi[9], matDXkDTk[9], matDXlDTl[9];
		double matDuvDTi[6], matDuvDTk[6], matDuvDTl[6];
		double matDdisDTk[3];
		double matDDotDTk[3], matDArcCosInDTk[3], matDsinwWDTk[3];
		double *pKR, *pKdA, *pKdB, *pKdG;
		int bs = 1;
		if (tpe == colmap) {
			bs = 1;
		}
		else if (tpe == bal) {
			bs = -1;
		}
		if ( nP == nM )	//main archor
		{
			matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
			matXj[1] = sin( ppt[1] );
			matXj[2] = cos( ppt[0] ) * cos( ppt[1] );

			pKR = KR + nP*9;
			matxyt[0] = pKR[0]*matXj[0] +pKR[1]*matXj[1] +pKR[2]*matXj[2];
			matxyt[1] = pKR[3]*matXj[0] +pKR[4]*matXj[1] +pKR[5]*matXj[2];
			matxyt[2] = pKR[6]*matXj[0] +pKR[7]*matXj[1] +pKR[8]*matXj[2];

			matDuvDxyt[0] = bs * 1 / matxyt[2]; 
			matDuvDxyt[1] = 0;
			matDuvDxyt[2] = -bs * matxyt[0] / (matxyt[2] * matxyt[2]);
			matDuvDxyt[3] = 0;	
			matDuvDxyt[4] = bs * 1 / matxyt[2]; 
			matDuvDxyt[5] = -bs * matxyt[1] / (matxyt[2] * matxyt[2]);
			
			//camera angles
			pKdG = KdG + nP*9;
			matDxytDRG[0] =  pKdG[0]*matXj[0] + pKdG[1]*matXj[1] + pKdG[2]*matXj[2];//omega
			matDxytDRG[1] =  pKdG[3]*matXj[0] + pKdG[4]*matXj[1] + pKdG[5]*matXj[2];
			matDxytDRG[2] =  pKdG[6]*matXj[0] + pKdG[7]*matXj[1] + pKdG[8]*matXj[2];

			pKdB = KdB + nP*9;
			matDxytDRB[0] =  pKdB[0]*matXj[0] + pKdB[1]*matXj[1] + pKdB[2]*matXj[2];//phi
			matDxytDRB[1] =  pKdB[3]*matXj[0] + pKdB[4]*matXj[1] + pKdB[5]*matXj[2];
			matDxytDRB[2] =  pKdB[6]*matXj[0] + pKdB[7]*matXj[1] + pKdB[8]*matXj[2];

			pKdA = KdA + nP*9;
			matDxytDRA[0] =  pKdA[0]*matXj[0] + pKdA[1]*matXj[1] + pKdA[2]*matXj[2];//kappa
			matDxytDRA[1] =  pKdA[3]*matXj[0] + pKdA[4]*matXj[1] + pKdA[5]*matXj[2];
			matDxytDRA[2] =  pKdA[6]*matXj[0] + pKdA[7]*matXj[1] + pKdA[8]*matXj[2];

			matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv????kappa??????
			matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];

			matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv????phi??????
			matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];

			matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv????omega??????
			matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];

			pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv????��???????????????
			pPA[3] = 0;						pPA[4] = 0;						pPA[5] = 0;
			pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
			pPA[9] = 0;						pPA[10] = 0;					pPA[11] = 0;
			
			//matXj[0] = sin(ppt[0]) * cos(ppt[1]);
			//matXj[1] = sin(ppt[1]);
			//matXj[2] = cos(ppt[0]) * cos(ppt[1]);
			//matxyt[0] = pKR[0] * matXj[0] + pKR[1] * matXj[1] + pKR[2] * matXj[2];
			//matxyt[1] = pKR[3] * matXj[0] + pKR[4] * matXj[1] + pKR[5] * matXj[2];
			//matxyt[2] = pKR[6] * matXj[0] + pKR[7] * matXj[1] + pKR[8] * matXj[2];
			//azimuth and elevation angle
			matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
			matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
			matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);

			matDxytDDH[0] = pKR[0]*matDXjDDH[0]	+ pKR[1]*matDXjDDH[2] + pKR[2]*matDXjDDH[4];//x????azimuth?????
			matDxytDDH[1] = pKR[0]*matDXjDDH[1] + pKR[1]*matDXjDDH[3] + pKR[2]*matDXjDDH[5];//x????elevation?????

			matDxytDDH[2] = pKR[3]*matDXjDDH[0]	+ pKR[4]*matDXjDDH[2] + pKR[5]*matDXjDDH[4];//y????azimuth?????
			matDxytDDH[3] = pKR[3]*matDXjDDH[1] + pKR[4]*matDXjDDH[3] + pKR[5]*matDXjDDH[5];//y????elevation?????

			matDxytDDH[4] = pKR[6]*matDXjDDH[0]	+ pKR[7]*matDXjDDH[2] + pKR[8]*matDXjDDH[4];//t????azimuth?????
			matDxytDDH[5] = pKR[6]*matDXjDDH[1] + pKR[7]*matDXjDDH[3] + pKR[8]*matDXjDDH[5];//t????elevation?????

			matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//uv????azimuth?????
			matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];
			matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//uv????elevation?????
			matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];

			pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = 0;//uv??azimuth??elevation??????
			pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = 0;
			
		}else 
			if ( nP == nN )	//associate archor
			{			
				ptAngle = pa + nN*6;
				matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
				matXj[1] = sin( ppt[1] );
				matXj[2] = cos( ppt[0] ) * cos( ppt[1] );
				
				pKR = KR + nN*9;
				pKdA= KdA+ nN*9;
				pKdB= KdB+ nN*9;
				pKdG= KdG+ nN*9;
				
				posei = pa + nM*6 + 3;//??��??Xc
				posek = pa + nN*6 + 3;//??��??Xc
				matPosei2k[0] = posek[0]-posei[0];		matPosei2k[1] = posek[1]-posei[1];		matPosei2k[2] = posek[2]-posei[2];//??��????��???????

				dDisTi2k = sqrt( matPosei2k[0]*matPosei2k[0] + matPosei2k[1]*matPosei2k[1] + matPosei2k[2]*matPosei2k[2] );//??��????��????????????
				dDot     = matXj[0]*matPosei2k[0] + matXj[1]*matPosei2k[1] + matXj[2]*matPosei2k[2];
				dArcCosIn= dDot / dDisTi2k;
				dW       = acos( dArcCosIn );

				matXk[0] = dDisTi2k * sin( dW + ppt[2] ) * matXj[0] - sin(ppt[2])*matPosei2k[0];
				matXk[1] = dDisTi2k * sin( dW + ppt[2] ) * matXj[1] - sin(ppt[2])*matPosei2k[1];
				matXk[2] = dDisTi2k * sin( dW + ppt[2] ) * matXj[2] - sin(ppt[2])*matPosei2k[2];

				matxyt[0] = pKR[0]*matXk[0] + pKR[1]*matXk[1] + pKR[2]*matXk[2];
				matxyt[1] = pKR[3]*matXk[0] + pKR[4]*matXk[1] + pKR[5]*matXk[2];
				matxyt[2] = pKR[6]*matXk[0] + pKR[7]*matXk[1] + pKR[8]*matXk[2];

				matDuvDxyt[0] = bs * 1 / matxyt[2]; 
				matDuvDxyt[1] = 0;
				matDuvDxyt[2] = -bs * matxyt[0] / (matxyt[2] * matxyt[2]);
				matDuvDxyt[3] = 0;	
				matDuvDxyt[4] = bs * 1 / matxyt[2]; 
				matDuvDxyt[5] = -bs * matxyt[1] / (matxyt[2] * matxyt[2]);
				
	//pKdG:pKR???????omega??????
	//pKdB:pKR???????phi??????
	//pKdA:pKR???????kappa??????
				//camera angles
				matDxytDRG[0] =  pKdG[0]*matXk[0] + pKdG[1]*matXk[1] + pKdG[2]*matXk[2];
				matDxytDRG[1] =  pKdG[3]*matXk[0] + pKdG[4]*matXk[1] + pKdG[5]*matXk[2];
				matDxytDRG[2] =  pKdG[6]*matXk[0] + pKdG[7]*matXk[1] + pKdG[8]*matXk[2];

				
				matDxytDRB[0] =  pKdB[0]*matXk[0] + pKdB[1]*matXk[1] + pKdB[2]*matXk[2];
				matDxytDRB[1] =  pKdB[3]*matXk[0] + pKdB[4]*matXk[1] + pKdB[5]*matXk[2];
				matDxytDRB[2] =  pKdB[6]*matXk[0] + pKdB[7]*matXk[1] + pKdB[8]*matXk[2];

				matDxytDRA[0] =  pKdA[0]*matXk[0] + pKdA[1]*matXk[1] + pKdA[2]*matXk[2];
				matDxytDRA[1] =  pKdA[3]*matXk[0] + pKdA[4]*matXk[1] + pKdA[5]*matXk[2];
				matDxytDRA[2] =  pKdA[6]*matXk[0] + pKdA[7]*matXk[1] + pKdA[8]*matXk[2];

				matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv????kappa??????
				matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];

				matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv????phi??????
				matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];

				matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv????omega??????
				matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];
				
				//Xj????azimuth??elevation??????
				//azimuth and elevation angle
				matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
				matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
				matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);

	//matPosei2k[0] = posek[0] - posei[0];		matPosei2k[1] = posek[1] - posei[1];		matPosei2k[2] = posek[2] - posei[2];//??��????��???????

	//dDisTi2k = sqrt(matPosei2k[0] * matPosei2k[0] + matPosei2k[1] * matPosei2k[1] + matPosei2k[2] * matPosei2k[2]);//??��????��????????????
	//dDot = matXj[0] * matPosei2k[0] + matXj[1] * matPosei2k[1] + matXj[2] * matPosei2k[2];
	//dArcCosIn = dDot / dDisTi2k;
	//dW = acos(dArcCosIn);

				matDDotDDH[0] = cos(ppt[0])*cos(ppt[1])*matPosei2k[0] - sin(ppt[0])*cos(ppt[1])*matPosei2k[2];//dDot????azimuth??????
				matDDotDDH[1] = -sin(ppt[0])*sin(ppt[1])*matPosei2k[0] + cos(ppt[1])*matPosei2k[1]-cos(ppt[0])*sin(ppt[1])*matPosei2k[2];//dDot????elevation??????

				matDArcCosInDDH[0] = matDDotDDH[0]/dDisTi2k;//ArcCosIn????azimuth??????
				matDArcCosInDDH[1] = matDDotDDH[1]/dDisTi2k;//ArcCosIn????elevation??????

				dDWDArcCosIn = -1/sqrt(1-dArcCosIn*dArcCosIn);//dW????ArccosIn??????

	//matXk[0] = dDisTi2k * sin(dW + ppt[2]) * matXj[0] - sin(ppt[2]) * matPosei2k[0];
	//matXk[1] = dDisTi2k * sin(dW + ppt[2]) * matXj[1] - sin(ppt[2]) * matPosei2k[1];
	//matXk[2] = dDisTi2k * sin(dW + ppt[2]) * matXj[2] - sin(ppt[2]) * matPosei2k[2];

				matDsinwWDDH[0] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[0];//sin(dW + ppt[2])????azimuth??????
				matDsinwWDDH[1] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[1];//sin(dW + ppt[2])????elevation??????

				tmp2[0] =  dDisTi2k*matXj[0]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[0];//matXk[0]????azimuth??????
				tmp2[1] =  dDisTi2k*matXj[0]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[1];//matXk[0]????elevation??????
				tmp2[2] =  dDisTi2k*matXj[1]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[2];//matXk[1]????azimuth??????
				tmp2[3] =  dDisTi2k*matXj[1]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[3];//matXk[1]????elevation??????
				tmp2[4] =  dDisTi2k*matXj[2]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[4];//matXk[2]????azimuth??????
				tmp2[5] =  dDisTi2k*matXj[2]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[5];//matXk[2]????elevation??????

				matDxytDDH[0] = pKR[0]*tmp2[0] + pKR[1]*tmp2[2] + pKR[2]*tmp2[4];//x????azimuth??????
				matDxytDDH[1] = pKR[0]*tmp2[1] + pKR[1]*tmp2[3] + pKR[2]*tmp2[5];//x????elevation??????
				matDxytDDH[2] = pKR[3]*tmp2[0] + pKR[4]*tmp2[2] + pKR[5]*tmp2[4];//y????azimuth??????
				matDxytDDH[3] = pKR[3]*tmp2[1] + pKR[4]*tmp2[3] + pKR[5]*tmp2[5];//y????elevation??????
				matDxytDDH[4] = pKR[6]*tmp2[0] + pKR[7]*tmp2[2] + pKR[8]*tmp2[4];//t????azimuth??????
				matDxytDDH[5] = pKR[6]*tmp2[1] + pKR[7]*tmp2[3] + pKR[8]*tmp2[5];//t????elevation??????

				matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//u????azimuth??????
				matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];//u????elevation??????
				matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//v????azimuth??????
				matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];//v????elevation??????
				
				//parallax angle
				matDXkDpa[0] = dDisTi2k*cos(dW+ppt[2])*matXj[0] - cos(ppt[2])*matPosei2k[0];//matXk[0]????parallax??????
				matDXkDpa[1] = dDisTi2k*cos(dW+ppt[2])*matXj[1] - cos(ppt[2])*matPosei2k[1];//matXk[1]????parallax??????
				matDXkDpa[2] = dDisTi2k*cos(dW+ppt[2])*matXj[2] - cos(ppt[2])*matPosei2k[2];//matXk[2]????parallax??????

				tmp1[0] = pKR[0]*matDuvDxyt[0] + pKR[3]*matDuvDxyt[1] + pKR[6]*matDuvDxyt[2];
				tmp1[1] = pKR[1]*matDuvDxyt[0] + pKR[4]*matDuvDxyt[1] + pKR[7]*matDuvDxyt[2];
				tmp1[2] = pKR[2]*matDuvDxyt[0] + pKR[5]*matDuvDxyt[1] + pKR[8]*matDuvDxyt[2];
				tmp1[3] = pKR[0]*matDuvDxyt[3] + pKR[3]*matDuvDxyt[4] + pKR[6]*matDuvDxyt[5];
				tmp1[4] = pKR[1]*matDuvDxyt[3] + pKR[4]*matDuvDxyt[4] + pKR[7]*matDuvDxyt[5];
				tmp1[5] = pKR[2]*matDuvDxyt[3] + pKR[5]*matDuvDxyt[4] + pKR[8]*matDuvDxyt[5];

				matDuvDpa[0] = tmp1[0]*matDXkDpa[0] + tmp1[1]*matDXkDpa[1] + tmp1[2]*matDXkDpa[2];//u????parallax??????
				matDuvDpa[1] = tmp1[3]*matDXkDpa[0] + tmp1[4]*matDXkDpa[1] + tmp1[5]*matDXkDpa[2];//v????parallax??????

				//uv????azimuth??elevation??parallax??????
				pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = matDuvDpa[0];
				pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = matDuvDpa[1];

	//posei = pa + nM * 6 + 3;//??��??Xc
	//posek = pa + nN * 6 + 3;//??��??Xc
	//matPosei2k[0] = posek[0] - posei[0];		matPosei2k[1] = posek[1] - posei[1];		matPosei2k[2] = posek[2] - posei[2];//??��????��???????
	//dDisTi2k = sqrt(matPosei2k[0] * matPosei2k[0] + matPosei2k[1] * matPosei2k[1] + matPosei2k[2] * matPosei2k[2]);//??��????��????????????
	//dDot = matXj[0] * matPosei2k[0] + matXj[1] * matPosei2k[1] + matXj[2] * matPosei2k[2];
	//dArcCosIn = dDot / dDisTi2k;
	//dW = acos(dArcCosIn);
	//matXk[0] = dDisTi2k * sin(dW + ppt[2]) * matXj[0] - sin(ppt[2]) * matPosei2k[0];
	//matXk[1] = dDisTi2k * sin(dW + ppt[2]) * matXj[1] - sin(ppt[2]) * matPosei2k[1];
	//matXk[2] = dDisTi2k * sin(dW + ppt[2]) * matXj[2] - sin(ppt[2]) * matPosei2k[2];
	// 
				//Ti
				matDdisDTi[0] = -matPosei2k[0]/dDisTi2k;//dDisTi2k????��??????Ti??????
				matDdisDTi[1] = -matPosei2k[1]/dDisTi2k;
				matDdisDTi[2] = -matPosei2k[2]/dDisTi2k;


				matDDotDTi[0] = -sin(ppt[0])*cos(ppt[1]);//dDot????��??????Ti??????
				matDDotDTi[1] = -sin(ppt[1]);
				matDDotDTi[2] = -cos(ppt[0])*cos(ppt[1]);

				
				matDArcCosInDTi[0] = (dDisTi2k*matDDotDTi[0]-dDot*matDdisDTi[0])/(dDisTi2k*dDisTi2k);//dArcCosIn????��??????Ti??????
				matDArcCosInDTi[1] = (dDisTi2k*matDDotDTi[1]-dDot*matDdisDTi[1])/(dDisTi2k*dDisTi2k);
				matDArcCosInDTi[2] = (dDisTi2k*matDDotDTi[2]-dDot*matDdisDTi[2])/(dDisTi2k*dDisTi2k);

				matDsinwWDTi[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[0];//sin(dW + ppt[2])????��??????Ti??????
				matDsinwWDTi[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[1];
				matDsinwWDTi[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[2];

				matDTi2kDTi[0] = matDTi2kDTi[4] = matDTi2kDTi[8] = -1;//matPosei2k????��??????Ti??????
				matDTi2kDTi[1] = matDTi2kDTi[2] = matDTi2kDTi[3] = 0;
				matDTi2kDTi[5] = matDTi2kDTi[6] = matDTi2kDTi[7] = 0;

				matDXkDTi[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[0] + dDisTi2k*matXj[0]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[0];//matXk????��??????Ti??????
				matDXkDTi[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[1] + dDisTi2k*matXj[0]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[1];
				matDXkDTi[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[2] + dDisTi2k*matXj[0]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[2];
				matDXkDTi[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[0] + dDisTi2k*matXj[1]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[3];
				matDXkDTi[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[1] + dDisTi2k*matXj[1]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[4];
				matDXkDTi[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[2] + dDisTi2k*matXj[1]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[5];
				matDXkDTi[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[0] + dDisTi2k*matXj[2]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[6];
				matDXkDTi[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[1] + dDisTi2k*matXj[2]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[7];
				matDXkDTi[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[2] + dDisTi2k*matXj[2]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[8];

				matDuvDTi[0] = tmp1[0]*matDXkDTi[0] + tmp1[1]*matDXkDTi[3] + tmp1[2]*matDXkDTi[6];//uv????��??????Ti??????
				matDuvDTi[1] = tmp1[0]*matDXkDTi[1] + tmp1[1]*matDXkDTi[4] + tmp1[2]*matDXkDTi[7];
				matDuvDTi[2] = tmp1[0]*matDXkDTi[2] + tmp1[1]*matDXkDTi[5] + tmp1[2]*matDXkDTi[8];
				matDuvDTi[3] = tmp1[3]*matDXkDTi[0] + tmp1[4]*matDXkDTi[3] + tmp1[5]*matDXkDTi[6];
				matDuvDTi[4] = tmp1[3]*matDXkDTi[1] + tmp1[4]*matDXkDTi[4] + tmp1[5]*matDXkDTi[7];
				matDuvDTi[5] = tmp1[3]*matDXkDTi[2] + tmp1[4]*matDXkDTi[5] + tmp1[5]*matDXkDTi[8];

				pAM[0] = matDuvDTi[0];		pAM[1] = matDuvDTi[1];		pAM[2] = matDuvDTi[2];
				pAM[3] = matDuvDTi[3];		pAM[4] = matDuvDTi[4];		pAM[5] = matDuvDTi[5];				 
				
				//Tk
				matDdisDTk[0] = matPosei2k[0]/dDisTi2k;
				matDdisDTk[1] = matPosei2k[1]/dDisTi2k;
				matDdisDTk[2] = matPosei2k[2]/dDisTi2k;

				matDDotDTk[0] = sin(ppt[0])*cos(ppt[1]);
				matDDotDTk[1] = sin(ppt[1]);
				matDDotDTk[2] = cos(ppt[0])*cos(ppt[1]);

				matDArcCosInDTk[0] = (dDisTi2k*matDDotDTk[0] - dDot*matDdisDTk[0])/(dDisTi2k*dDisTi2k);
				matDArcCosInDTk[1] = (dDisTi2k*matDDotDTk[1] - dDot*matDdisDTk[1])/(dDisTi2k*dDisTi2k);
				matDArcCosInDTk[2] = (dDisTi2k*matDDotDTk[2] - dDot*matDdisDTk[2])/(dDisTi2k*dDisTi2k);

				matDsinwWDTk[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[0];
				matDsinwWDTk[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[1];
				matDsinwWDTk[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[2];

				matDTi2kDTk[0] = matDTi2kDTk[4] = matDTi2kDTk[8] = 1;
				matDTi2kDTk[1] = matDTi2kDTk[2] = matDTi2kDTk[3] = 0;
				matDTi2kDTk[5] = matDTi2kDTk[6] = matDTi2kDTk[7] = 0;


				matDXkDTk[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[0] + dDisTi2k*matXj[0]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[0];
				matDXkDTk[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[1] + dDisTi2k*matXj[0]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[1];
				matDXkDTk[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[2] + dDisTi2k*matXj[0]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[2];
				matDXkDTk[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[0] + dDisTi2k*matXj[1]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[3];
				matDXkDTk[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[1] + dDisTi2k*matXj[1]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[4];
				matDXkDTk[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[2] + dDisTi2k*matXj[1]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[5];
				matDXkDTk[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[0] + dDisTi2k*matXj[2]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[6];
				matDXkDTk[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[1] + dDisTi2k*matXj[2]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[7];
				matDXkDTk[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[2] + dDisTi2k*matXj[2]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[8];

				matDuvDTk[0] = tmp1[0]*matDXkDTk[0] + tmp1[1]*matDXkDTk[3] + tmp1[2]*matDXkDTk[6];//uv???��??????Tk??????
				matDuvDTk[1] = tmp1[0]*matDXkDTk[1] + tmp1[1]*matDXkDTk[4] + tmp1[2]*matDXkDTk[7];
				matDuvDTk[2] = tmp1[0]*matDXkDTk[2] + tmp1[1]*matDXkDTk[5] + tmp1[2]*matDXkDTk[8];
				matDuvDTk[3] = tmp1[3]*matDXkDTk[0] + tmp1[4]*matDXkDTk[3] + tmp1[5]*matDXkDTk[6];
				matDuvDTk[4] = tmp1[3]*matDXkDTk[1] + tmp1[4]*matDXkDTk[4] + tmp1[5]*matDXkDTk[7];
				matDuvDTk[5] = tmp1[3]*matDXkDTk[2] + tmp1[4]*matDXkDTk[5] + tmp1[5]*matDXkDTk[8];
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv???��???????????
				pPA[3] = matDuvDTk[0];			pPA[4] = matDuvDTk[1];			pPA[5] = matDuvDTk[2];
				pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
				pPA[9] = matDuvDTk[3];			pPA[10] = matDuvDTk[4];			pPA[11] = matDuvDTk[5];
			}
			else
			{
				ptAngle = pa + nP*6;
				matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
				matXj[1] = sin( ppt[1] );
				matXj[2] = cos( ppt[0] ) * cos( ppt[1] );

				pKR = KR + nP*9;
				pKdA= KdA+ nP*9;
				pKdB= KdB+ nP*9;
				pKdG= KdG+ nP*9;
				
				posei = pa + nM*6 + 3;//??��??
				posek = pa + nN*6 + 3;//??��??
				posel = pa + nP*6  + 3;//????��??

				matPosei2k[0] = posek[0]-posei[0];		matPosei2k[1] = posek[1]-posei[1];		matPosei2k[2] = posek[2]-posei[2];//??��????��??
				matPosei2l[0] = posel[0]-posei[0];		matPosei2l[1] = posel[1]-posei[1];		matPosei2l[2] = posel[2]-posei[2];//??��??????��??

				dDisTi2k = sqrt( matPosei2k[0]*matPosei2k[0] + matPosei2k[1]*matPosei2k[1] + matPosei2k[2]*matPosei2k[2] );
				dDot     = matXj[0]*matPosei2k[0] + matXj[1]*matPosei2k[1] + matXj[2]*matPosei2k[2];
				dArcCosIn= dDot / dDisTi2k;
				dW       = acos( dArcCosIn );

				matXl[0] = dDisTi2k * sin( dW + ppt[2] ) * matXj[0] - sin(ppt[2])*matPosei2l[0];
				matXl[1] = dDisTi2k * sin( dW + ppt[2] ) * matXj[1] - sin(ppt[2])*matPosei2l[1];
				matXl[2] = dDisTi2k * sin( dW + ppt[2] ) * matXj[2] - sin(ppt[2])*matPosei2l[2];

				matxyt[0] = pKR[0]*matXl[0] + pKR[1]*matXl[1] + pKR[2]*matXl[2];
				matxyt[1] = pKR[3]*matXl[0] + pKR[4]*matXl[1] + pKR[5]*matXl[2];
				matxyt[2] = pKR[6]*matXl[0] + pKR[7]*matXl[1] + pKR[8]*matXl[2];

				matDuvDxyt[0] = bs*1/matxyt[2];	
				matDuvDxyt[1] = 0;
				matDuvDxyt[2] = -bs*matxyt[0]/(matxyt[2]*matxyt[2]);
				matDuvDxyt[3] = 0;	
				matDuvDxyt[4] = bs*1/matxyt[2];		
				matDuvDxyt[5] = -bs*matxyt[1]/(matxyt[2]*matxyt[2]);	
			
				//camera angle
				matDxytDRG[0] =  pKdG[0]*matXl[0] + pKdG[1]*matXl[1] + pKdG[2]*matXl[2];
				matDxytDRG[1] =  pKdG[3]*matXl[0] + pKdG[4]*matXl[1] + pKdG[5]*matXl[2];
				matDxytDRG[2] =  pKdG[6]*matXl[0] + pKdG[7]*matXl[1] + pKdG[8]*matXl[2];

				matDxytDRB[0] =  pKdB[0]*matXl[0] + pKdB[1]*matXl[1] + pKdB[2]*matXl[2];
				matDxytDRB[1] =  pKdB[3]*matXl[0] + pKdB[4]*matXl[1] + pKdB[5]*matXl[2];
				matDxytDRB[2] =  pKdB[6]*matXl[0] + pKdB[7]*matXl[1] + pKdB[8]*matXl[2];

				matDxytDRA[0] =  pKdA[0]*matXl[0] + pKdA[1]*matXl[1] + pKdA[2]*matXl[2];
				matDxytDRA[1] =  pKdA[3]*matXl[0] + pKdA[4]*matXl[1] + pKdA[5]*matXl[2];
				matDxytDRA[2] =  pKdA[6]*matXl[0] + pKdA[7]*matXl[1] + pKdA[8]*matXl[2];			

				matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv??????��??kappa??????????
				matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];

				matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv??????��??phi??????????
				matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];

				matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv??????��??omega??????????
				matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];
			
				//azimuth and elevation angle
				matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
				matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
				matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);

				matDDotDDH[0] = cos(ppt[0])*cos(ppt[1])*matPosei2k[0] - sin(ppt[0])*cos(ppt[1])*matPosei2k[2];
				matDDotDDH[1] = -sin(ppt[0])*sin(ppt[1])*matPosei2k[0] + cos(ppt[1])*matPosei2k[1]-cos(ppt[0])*sin(ppt[1])*matPosei2k[2];

				matDArcCosInDDH[0] = matDDotDDH[0]/dDisTi2k;
				matDArcCosInDDH[1] = matDDotDDH[1]/dDisTi2k;

				dDWDArcCosIn = -1/sqrt(1-dArcCosIn*dArcCosIn);

				matDsinwWDDH[0] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[0];
				matDsinwWDDH[1] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[1];

				tmp2[0] =  dDisTi2k*matXj[0]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[0];
				tmp2[1] =  dDisTi2k*matXj[0]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[1];
				tmp2[2] =  dDisTi2k*matXj[1]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[2];
				tmp2[3] =  dDisTi2k*matXj[1]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[3];
				tmp2[4] =  dDisTi2k*matXj[2]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[4];
				tmp2[5] =  dDisTi2k*matXj[2]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[5];

				matDxytDDH[0] = pKR[0]*tmp2[0] + pKR[1]*tmp2[2] + pKR[2]*tmp2[4];
				matDxytDDH[1] = pKR[0]*tmp2[1] + pKR[1]*tmp2[3] + pKR[2]*tmp2[5];
				matDxytDDH[2] = pKR[3]*tmp2[0] + pKR[4]*tmp2[2] + pKR[5]*tmp2[4];
				matDxytDDH[3] = pKR[3]*tmp2[1] + pKR[4]*tmp2[3] + pKR[5]*tmp2[5];
				matDxytDDH[4] = pKR[6]*tmp2[0] + pKR[7]*tmp2[2] + pKR[8]*tmp2[4];
				matDxytDDH[5] = pKR[6]*tmp2[1] + pKR[7]*tmp2[3] + pKR[8]*tmp2[5];

				matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//u??azimuth??????
				matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];//u??elevation??????
				matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//v??azimuth??????
				matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];//v??elevation??????
			
				//parallax angle
				matDXkDpa[0] = dDisTi2k*cos(dW+ppt[2])*matXj[0] - cos(ppt[2])*matPosei2l[0];
				matDXkDpa[1] = dDisTi2k*cos(dW+ppt[2])*matXj[1] - cos(ppt[2])*matPosei2l[1];
				matDXkDpa[2] = dDisTi2k*cos(dW+ppt[2])*matXj[2] - cos(ppt[2])*matPosei2l[2];	

				tmp1[0] = pKR[0]*matDuvDxyt[0] + pKR[3]*matDuvDxyt[1] + pKR[6]*matDuvDxyt[2];
				tmp1[1] = pKR[1]*matDuvDxyt[0] + pKR[4]*matDuvDxyt[1] + pKR[7]*matDuvDxyt[2];
				tmp1[2] = pKR[2]*matDuvDxyt[0] + pKR[5]*matDuvDxyt[1] + pKR[8]*matDuvDxyt[2];
				tmp1[3] = pKR[0]*matDuvDxyt[3] + pKR[3]*matDuvDxyt[4] + pKR[6]*matDuvDxyt[5];
				tmp1[4] = pKR[1]*matDuvDxyt[3] + pKR[4]*matDuvDxyt[4] + pKR[7]*matDuvDxyt[5];
				tmp1[5] = pKR[2]*matDuvDxyt[3] + pKR[5]*matDuvDxyt[4] + pKR[8]*matDuvDxyt[5];

				matDuvDpa[0] = tmp1[0]*matDXkDpa[0] + tmp1[1]*matDXkDpa[1] + tmp1[2]*matDXkDpa[2];//u??parallax??????
				matDuvDpa[1] = tmp1[3]*matDXkDpa[0] + tmp1[4]*matDXkDpa[1] + tmp1[5]*matDXkDpa[2];//v??parallax??????

				pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = matDuvDpa[0];
				pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = matDuvDpa[1];
			
				//Ti
				matDdisDTi[0] = -matPosei2k[0]/dDisTi2k;
				matDdisDTi[1] = -matPosei2k[1]/dDisTi2k;
				matDdisDTi[2] = -matPosei2k[2]/dDisTi2k;

				matDDotDTi[0] = -sin(ppt[0])*cos(ppt[1]);
				matDDotDTi[1] = -sin(ppt[1]);
				matDDotDTi[2] = -cos(ppt[0])*cos(ppt[1]);

				matDArcCosInDTi[0] = (dDisTi2k*matDDotDTi[0]-dDot*matDdisDTi[0])/(dDisTi2k*dDisTi2k);
				matDArcCosInDTi[1] = (dDisTi2k*matDDotDTi[1]-dDot*matDdisDTi[1])/(dDisTi2k*dDisTi2k);
				matDArcCosInDTi[2] = (dDisTi2k*matDDotDTi[2]-dDot*matDdisDTi[2])/(dDisTi2k*dDisTi2k);

				matDsinwWDTi[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[0];
				matDsinwWDTi[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[1];
				matDsinwWDTi[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[2];

				matDTi2kDTi[0] = matDTi2kDTi[4] = matDTi2kDTi[8] = -1;
				matDTi2kDTi[1] = matDTi2kDTi[2] = matDTi2kDTi[3] = 0;
				matDTi2kDTi[5] = matDTi2kDTi[6] = matDTi2kDTi[7] = 0;

				matDXkDTi[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[0] + dDisTi2k*matXj[0]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[0];
				matDXkDTi[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[1] + dDisTi2k*matXj[0]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[1];
				matDXkDTi[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[2] + dDisTi2k*matXj[0]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[2];
				matDXkDTi[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[0] + dDisTi2k*matXj[1]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[3];
				matDXkDTi[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[1] + dDisTi2k*matXj[1]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[4];
				matDXkDTi[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[2] + dDisTi2k*matXj[1]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[5];
				matDXkDTi[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[0] + dDisTi2k*matXj[2]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[6];
				matDXkDTi[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[1] + dDisTi2k*matXj[2]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[7];
				matDXkDTi[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[2] + dDisTi2k*matXj[2]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[8];

				matDuvDTi[0] = tmp1[0]*matDXkDTi[0] + tmp1[1]*matDXkDTi[3] + tmp1[2]*matDXkDTi[6];
				matDuvDTi[1] = tmp1[0]*matDXkDTi[1] + tmp1[1]*matDXkDTi[4] + tmp1[2]*matDXkDTi[7];
				matDuvDTi[2] = tmp1[0]*matDXkDTi[2] + tmp1[1]*matDXkDTi[5] + tmp1[2]*matDXkDTi[8];
				matDuvDTi[3] = tmp1[3]*matDXkDTi[0] + tmp1[4]*matDXkDTi[3] + tmp1[5]*matDXkDTi[6];
				matDuvDTi[4] = tmp1[3]*matDXkDTi[1] + tmp1[4]*matDXkDTi[4] + tmp1[5]*matDXkDTi[7];
				matDuvDTi[5] = tmp1[3]*matDXkDTi[2] + tmp1[4]*matDXkDTi[5] + tmp1[5]*matDXkDTi[8];

				pAM[0] = matDuvDTi[0];		pAM[1] = matDuvDTi[1];		pAM[2] = matDuvDTi[2];//uv????��???????????????
				pAM[3] = matDuvDTi[3];		pAM[4] = matDuvDTi[4];		pAM[5] = matDuvDTi[5];	
			
				//Tk
				matDdisDTk[0] = matPosei2k[0]/dDisTi2k;
				matDdisDTk[1] = matPosei2k[1]/dDisTi2k;
				matDdisDTk[2] = matPosei2k[2]/dDisTi2k;

				matDDotDTk[0] = sin(ppt[0])*cos(ppt[1]);
				matDDotDTk[1] = sin(ppt[1]);
				matDDotDTk[2] = cos(ppt[0])*cos(ppt[1]);

				matDArcCosInDTk[0] = (dDisTi2k*matDDotDTk[0] - dDot*matDdisDTk[0])/(dDisTi2k*dDisTi2k);
				matDArcCosInDTk[1] = (dDisTi2k*matDDotDTk[1] - dDot*matDdisDTk[1])/(dDisTi2k*dDisTi2k);
				matDArcCosInDTk[2] = (dDisTi2k*matDDotDTk[2] - dDot*matDdisDTk[2])/(dDisTi2k*dDisTi2k);

				matDsinwWDTk[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[0];
				matDsinwWDTk[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[1];
				matDsinwWDTk[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[2];

				matDTi2kDTk[0] = matDTi2kDTk[4] = matDTi2kDTk[8] = 1;
				matDTi2kDTk[1] = matDTi2kDTk[2] = matDTi2kDTk[3] = 0;
				matDTi2kDTk[5] = matDTi2kDTk[6] = matDTi2kDTk[7] = 0;

				matDXkDTk[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[0] + dDisTi2k*matXj[0]*matDsinwWDTk[0];
				matDXkDTk[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[1] + dDisTi2k*matXj[0]*matDsinwWDTk[1];
				matDXkDTk[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[2] + dDisTi2k*matXj[0]*matDsinwWDTk[2];
				matDXkDTk[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[0] + dDisTi2k*matXj[1]*matDsinwWDTk[0];
				matDXkDTk[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[1] + dDisTi2k*matXj[1]*matDsinwWDTk[1];
				matDXkDTk[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[2] + dDisTi2k*matXj[1]*matDsinwWDTk[2];
				matDXkDTk[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[0] + dDisTi2k*matXj[2]*matDsinwWDTk[0];
				matDXkDTk[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[1] + dDisTi2k*matXj[2]*matDsinwWDTk[1];
				matDXkDTk[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[2] + dDisTi2k*matXj[2]*matDsinwWDTk[2];

				matDuvDTk[0] = tmp1[0]*matDXkDTk[0] + tmp1[1]*matDXkDTk[3] + tmp1[2]*matDXkDTk[6];
				matDuvDTk[1] = tmp1[0]*matDXkDTk[1] + tmp1[1]*matDXkDTk[4] + tmp1[2]*matDXkDTk[7];
				matDuvDTk[2] = tmp1[0]*matDXkDTk[2] + tmp1[1]*matDXkDTk[5] + tmp1[2]*matDXkDTk[8];
				matDuvDTk[3] = tmp1[3]*matDXkDTk[0] + tmp1[4]*matDXkDTk[3] + tmp1[5]*matDXkDTk[6];
				matDuvDTk[4] = tmp1[3]*matDXkDTk[1] + tmp1[4]*matDXkDTk[4] + tmp1[5]*matDXkDTk[7];
				matDuvDTk[5] = tmp1[3]*matDXkDTk[2] + tmp1[4]*matDXkDTk[5] + tmp1[5]*matDXkDTk[8];
				
				pAA[0] = matDuvDTk[0];		pAA[1] = matDuvDTk[1];		pAA[2] = matDuvDTk[2];//uv???��???????????????
				pAA[3] = matDuvDTk[3];		pAA[4] = matDuvDTk[4];		pAA[5] = matDuvDTk[5];
			
				//Tl
				matDXlDTl[0] = matDXlDTl[4] = matDXlDTl[8] = -sin(ppt[2]);
				matDXlDTl[1] = matDXlDTl[2] = matDXlDTl[3] = 0;
				matDXlDTl[5] = matDXlDTl[6] = matDXlDTl[7] = 0;	

				matDuvDTl[0] = tmp1[0]*matDXlDTl[0] + tmp1[1]*matDXlDTl[3] + tmp1[2]*matDXlDTl[6];
				matDuvDTl[1] = tmp1[0]*matDXlDTl[1] + tmp1[1]*matDXlDTl[4] + tmp1[2]*matDXlDTl[7];
				matDuvDTl[2] = tmp1[0]*matDXlDTl[2] + tmp1[1]*matDXlDTl[5] + tmp1[2]*matDXlDTl[8];
				matDuvDTl[3] = tmp1[3]*matDXlDTl[0] + tmp1[4]*matDXlDTl[3] + tmp1[5]*matDXlDTl[6];
				matDuvDTl[4] = tmp1[3]*matDXlDTl[1] + tmp1[4]*matDXlDTl[4] + tmp1[5]*matDXlDTl[7];
				matDuvDTl[5] = tmp1[3]*matDXlDTl[2] + tmp1[4]*matDXlDTl[5] + tmp1[5]*matDXlDTl[8];
			
				pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv??????��???????????
				pPA[3] = matDuvDTl[0];			pPA[4] = matDuvDTl[1];			pPA[5] = matDuvDTl[2];
				pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
				pPA[9] = matDuvDTl[3];			pPA[10] = matDuvDTl[4];			pPA[11] = matDuvDTl[5];
			}
	}
	double *PBA::pba_angle2xyz(double* p)
	{
		double *n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));
		static int i, j;
		double* pAngle;
		double xj[3], xk[3];
		double Tik[3];
		int nM, nA;
		double Dik;
		double w, w2;
		int cnp = 6, pnp = 3;
		for (i = 0; i < m_n3Dpts; i++)
		{
			pAngle = p + m_ncams * 6 + i * 3;//azimuth angle
			w = *(p + m_ncams * 6 + i * 3 + 2);//parallax angle

			xj[0] = sin(*(pAngle)) * cos(*(pAngle + 1));//x
			xj[1] = sin(*(pAngle + 1));//y
			xj[2] = cos(*(pAngle)) * cos(*(pAngle + 1));//z

			nM = *(m_archor + i * 3 + 1);
			nA = *(m_archor + i * 3 + 2);

			Tik[0] = -*(p + nM * 6 + 3) + *(p + nA * 6 + 3);
			Tik[1] = -*(p + nM * 6 + 4) + *(p + nA * 6 + 4);
			Tik[2] = -*(p + nM * 6 + 5) + *(p + nA * 6 + 5);

			Dik = sqrt(Tik[0] * Tik[0] + Tik[1] * Tik[1] + Tik[2] * Tik[2]);

			w2 = acos((xj[0] * Tik[0] + xj[1] * Tik[1] + xj[2] * Tik[2]) / Dik);

			//?????????????��???xyz????
			xk[0] = (Dik * sin(w2 + w) * xj[0]) / sin(w);
			xk[1] = (Dik * sin(w2 + w) * xj[1]) / sin(w);
			xk[2] = (Dik * sin(w2 + w) * xj[2]) / sin(w);

			//??????????xyz????
			double* ptr = n3DptsT + i * 3;
			double* ptr1 = p + nM * 6 + 3;
			// ????? ptr ????��??
			if (ptr != nullptr) {
				*ptr = *ptr1 + xk[0];
				*(ptr + 1) = *(ptr1 + 1) + xk[1];
				*(ptr + 2) = *(ptr1 + 2) + xk[2];
				//????
				//printf("%f %f %f\n", *(n3DptsT + i * 3), *(n3DptsT + i * 3 + 1), *(n3DptsT + i * 3 + 2));
			}
			else {
				throw std::runtime_error("????????: ptr ????? NULL");
			}
		}

		return n3DptsT;
	}
	void PBA::pba_saveInitialXYZ(const char* sz3Dpt, double* p)
	{
		static int i = 0;
		double dx, dy, dz;
		FILE* fp = nullptr;
		// int nM, nA;
		//save features xyz
		if (sz3Dpt != NULL)
		{
			fopen_s(&fp, sz3Dpt, "w");
			for (i = 0; i < m_n3Dpts; i++)
			{
				dx = *(p + i * 3);
				dy = *(p + i * 3 + 1);
				dz = *(p + i * 3 + 2);
				//nM = *(m_archor + i * 3 + 1);
				//nA = *(m_archor + i * 3 + 2);
				//fprintf(fp, "%0.5lf %0.5lf %0.5lf %d %d\n", dx, dy, dz,nM, nA);
				//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1,dx, dy, dz);
				fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", dx, dy, dz);
				//if (i == 2)
				//{
				//	printf("%0.5lf %0.5lf %0.5lf\n", dx, dy, dz);
				//	printf("%d\n", m_n3Dpts);
				//}
					
			}
			//printf("%d\n", i);
			fclose(fp);
		}
		free(p);
	}
	void PBA::pba_saveInitialParallax(const char* sz3Dpt, double* p)
	{
		//double* n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));

		static int i = 0;
		double azimuth, elevation, parallax;
		FILE* fp = nullptr, * fpc = nullptr;
		double* pAngle;
		//save features xyz
		if (sz3Dpt != NULL)
		{
			fopen_s(&fp, sz3Dpt, "w");
			for (i = 0; i < m_n3Dpts; i++)
			{
				pAngle = p + m_ncams * 6 + i * 3;

				azimuth = *pAngle;
				elevation = *(pAngle + 1);
				parallax = *(pAngle + 2);
				fprintf(fp, "%0.5lf     %0.5lf     %0.5lf\n", azimuth, elevation, parallax);
				//fprintf(fp, "%0.5lf\n", parallax*180/PI);
			}
			fclose(fp);
		}
		free(p);
	}

	bool PBA::ba_motstr_levmar( )
	{
		/*printf("%f\n", m_e1);*/
		//if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
		//{	
		//	fprintf( stderr, "ParallaxBA: Missing related source file, please input them first\n");	
		//	return false;
		//}
		if (m_motstruct == NULL || m_imgpts == NULL || m_K == NULL)
		{
			fprintf(stderr, "ParallaxBA: Missing related source file, please input them first\n");
			return false;
		}
		
		string originalPath(m_szReport);
		size_t pos = originalPath.find_last_of("/\\");
		string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
		string cond_log = parentPath + "/" + "pba_cond_log.txt";
		string spectral_log = parentPath + "/" + "pba_sigma_spectral_log.txt";
		string Upgm = parentPath + "/" + "PBA_U.pgm";
		string Utxt = parentPath + "/" + "PBA_U.txt";
		string Spgm = parentPath + "/" + "PBA_S.pgm";
		string Stxt = parentPath + "/" + "PBA_S.txt";
		string Hpgm = parentPath + "/" + "PBA_H.pgm";
		string Htxt = parentPath + "/" + "PBA_H.txt";

		FILE *fpRe = nullptr;
		if (m_szReport!= NULL)
		{
			fopen_s(&fpRe, m_szReport, "w ");
			fprintf(fpRe, "%d  poses, %d 3D features, %d projection\n", m_ncams, m_n3Dpts, m_n2Dprojs );
			fprintf(fpRe, "Levenberg-Marquardt is used\n" );
		}

		int i, ii;
		int n = m_n3Dpts, m = m_ncams, cnp = 6, pnp = 3, mnp = 2;
		int nvis, nuis, itno, issolved, nobs, nvars, nMaxS, nu=2, stop=0;

		double *p = m_motstruct, *x = m_imgpts;	//p pointer refers to unknown parameters, x pointer refers to image coordinate	
		double *U, *V, *W, *e, *eab, *eab_pre, *E, *S, *dp, *IV; // pointers into U V W E S IV
		double *pa = NULL, *pb = NULL, *ea, *eb, *dpa, *dpb, *hx, *Ex, *rx, *pdp; // pointers into p, jac, eab and dp respectively 	double *Ex, *rx; 
		double initialerror = 0, error = 0;
		int    nIter = 0, nLinear = 0;
		int Usz, Vsz, Wsz, easz, esz, ebsz, Sblsz, Sdim; 
		static double mu, tmp, p_eL2, eab_inf, pdp_eL2, x2, delt2, init_p_eL2 = 0; 
		double p_L2, dp_L2=DBL_MAX, dF, dL;
		double tau=fabs(m_Tau), eps1=fabs(m_e1),eps3_sq=m_e3*m_e3;
		bool init = false, ordering = true;
		double tStart, tEnd, tTimeUse, t0, t1, t11, t2, t3, t4, t5;
		struct sba_crsm Sidxij, Uidxij;		// S mask and U mask
		double tCons = 0, tSolve = 0, tCost = 0, tTimeIter = 0, tTimeCal = 0;
		int nLmIterations = 0;

		Usz=cnp*cnp;												//6*6
		Vsz=pnp * pnp;												//3*3
		Wsz=cnp * pnp;												//6*3
		esz=mnp; easz=cnp; ebsz=pnp;								//2 ?? 6 ?? 3
		Sblsz=cnp*cnp;												//6*6
		Sdim=m * cnp;												//m_ncams??????????*6
		nvis = m_n2Dprojs;											//?????????2D????????
		mu=eab_inf=0.0;												//0.0

		nuis = ba_ConstructSmask( Sidxij, Uidxij );
		nobs=nvis*mnp;
		nvars=m*cnp + n*pnp;
				
		S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));
		W	=	(double *)malloc(nvis*Wsz*sizeof(double));
		U	=	(double *)malloc(nuis*Usz*sizeof(double));
		V	=	(double *)malloc(n*Vsz*sizeof(double));
		IV	=	(double *)malloc(n*Vsz*sizeof(double));
		e	=	(double *)malloc(nobs*sizeof(double));
		eab	=	(double *)malloc(nvars*sizeof(double));
		eab_pre = (double *)malloc(nvars*sizeof(double));
		E	=	(double *)malloc(m*cnp*sizeof(double));		
		dp	=	(double *)malloc(nvars*sizeof(double));
		hx	=	(double *)malloc(nobs*sizeof(double));	
		pdp	=	(double *)malloc(nvars*sizeof(double));	
		pa=p; pb=p+m*cnp; ea=eab; eb=eab+m*cnp;	dpa=dp; dpb=dp+m*cnp;	

		cholmod_start (&m_cS) ;  	
		//m_cS.print_function = NULL;	
		int *Ap  = (int*)malloc((m+1)*sizeof(int));
		int * Aii = (int*)malloc(m_nS*sizeof(int));
		ba_constructAuxCSSLM( Ap, Aii );

		m_cholSparseE = cholmod_zeros( cnp*m, 1, CHOLMOD_REAL, &m_cS);
		Ex = (double*)m_cholSparseE->x;
		nMaxS = (m_nS-m_ncams)*36+m_ncams*21;	//maximum non-zero element in S matrix 

		m_cholSparseS = cholmod_allocate_sparse(m_ncams*6,m_ncams*6,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
		int *Sp, *Si;
		double* Sx = NULL;
		Sp = (int*)m_cholSparseS->p;		//column pointer
		Si = (int*)m_cholSparseS->i;		//row pointer
		
		//Compute initial error and initial reprojection error 
		tStart = clock();
		pba_cost(p, hx, m_archor ); 

		p_eL2 = nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */

		initialerror = p_eL2/nvis;
		printf("Initial Error %0.1lf [%0.8lf]\n", p_eL2, initialerror);

		if( m_szReport != NULL )
			fprintf( fpRe, "Initial Error  %lf\n", initialerror );

		init_p_eL2 = p_eL2;

		//Iteration 
		FILE* fp3 = nullptr, *fp4 = nullptr;
		errno_t err3 = fopen_s(&fp3, cond_log.c_str(), "w");
		errno_t err4 = fopen_s(&fp4, spectral_log.c_str(), "w");

		VectorXd sigma_spectral;
		double S_sigma_max,S_sigma_min,S_cond,rho_ = 0, relative_state_change = 0, relative_mse_change = 0;
		double Lips = 0, descent_dir_valid = 0;

		// m_nMaxIter = 2;
		for(itno=0; itno<m_nMaxIter && !stop; ++itno)
		{
			//Setup S matrix, include two step
			memset( U, 0, nuis*Usz*sizeof(double) );
			memset( ea, 0,m*easz*sizeof(double) );
			memset( V, 0, n*Vsz*sizeof(double));
			memset( IV, 0, n*Vsz*sizeof(double));
			memset( eb, 0, n*ebsz*sizeof(double));
			memset( W, 0, nvis*Wsz*sizeof(double));	
							
			//Step one; compute W V U directly, don't save each projection image Jacobian 
			t0 = clock();
			if ( m_bRobustKernel)
				pba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
			else
				pba_jacobian(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
				
			
			// convertUtoDenseMatrix(U,&Uidxij,Upgm.c_str(),Utxt.c_str());
			// convertHtoDenseMatrix(U,V,W,Hpgm.c_str(),Htxt.c_str());

			t1 = clock();
			if( itno == 0)
				mu = ba_computeInitialmu( U, V, Uidxij, tau, nvars );//damping factor

			nLmIterations = 0;
			while(1) //determine increment using adaptive damping
			{
				nLmIterations++;
				t11 = clock();

				ba_inverseVLM( V, IV, Uidxij, mu ); //compute inverse matrix of V
				
				//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
				memset( E, 0, m*easz*sizeof(double));
				memset(S, 0, m_nS * Sblsz * sizeof(double));
				// Eigen::MatrixXd S_full = Eigen::MatrixXd::Zero(fullSize,fullSize);
	
				ba_constructSLM(S, E, U, IV, W, ea, eb, Sidxij, mu);
				t2 = clock();	
				
				//Solve linear equation
				//set CSS format using S matrix
				ba_constructCSSLM(Si, Sp, Sx, S, m_cholSparseS, Sidxij, init); 
				for ( ii = 0; ii < cnp*m; ii++  )
					Ex[ii] = E[ii];	
				
				ba_solveCholmodLM( Ap, Aii, init, ordering);
				nLinear++;
				
				init = true;
				rx = (double*)m_cholSparseR->x;

				if (m_cS.status != CHOLMOD_NOT_POSDEF )
				{
					for ( ii = 0; ii < cnp*m; ii++ )
						dpa[ii] = rx[ii];
					issolved = 1;
				}
				else
					issolved = 1;
				
				t3 = clock();
							
				if(issolved)
				{
					//Solve features
					ba_solveFeatures( W, IV, ea, eb, dpa, dpb );
				
					// Compute ||J^T e||_inf and ||p||^2 
					for(i=0, p_L2=eab_inf=0.0; i<nvars; ++i)
					{
						if(eab_inf < (tmp=fabs(eab[i])))//eab_inf:最大梯度分量
							eab_inf=tmp;
						
						p_L2+=p[i]*p[i];
					}								
		
					//update
					for (i = 0, dp_L2 = 0.0; i < nvars; ++i)// compute p's new estimate and ||dp||^2 
					{
						tmp = dp[i];
						pdp[i]=p[i] + dp[i];
						dp_L2 += tmp*tmp;
						// dp_L2 += abs(tmp)*tmp;
					}
					//sqrt(dp_L2)：状态向量的实际变化幅度
					//sqrt(p_L2)：当前状态向量的模长
					relative_state_change = sqrt(dp_L2) / (sqrt(p_L2) + m_e1);
					if (sqrt(dp_L2) <= m_e1 * (sqrt(p_L2) + m_e1))
					{
						stop = 1;	break;//relative change of state vector is small
					}

					ba_updateKR(m_KR, m_KdA, m_KdB, m_KdG, m_K, pdp);
					t4 = clock();
					
					pba_cost(pdp, hx, m_archor );
					pdp_eL2=nrmL2xmy(hx, x, hx, nobs); 	
					error = pdp_eL2/nvis;
					t5 = clock();
									
					if ( m_bRobustKernel )
					{
						pdp_eL2 = 0;
						delt2 = m_delt*m_delt;
						if ( m_nRobustType==1)						//Cauchy Kernel Function
						{
							for ( i = 0; i < m_n2Dprojs; i++ )
							{					
								x2 = hx[i*2]*hx[i*2]+hx[i*2+1]*hx[i*2+1];
								x2 = delt2*log( x2/delt2 + 1 );
								pdp_eL2 += x2;
							}
						}
						else										//Huber Kernel Function
						{
							for ( i = 0; i < m_n2Dprojs; i++ )
							{					
								x2 = hx[i*2]*hx[i*2]+hx[i*2+1]*hx[i*2+1];
								
								if (x2 <= delt2)  // inlier
									x2 = x2;
								else  // outliers
									x2 = 2*sqrt(x2)*m_delt - delt2;
									
								pdp_eL2 += x2;
							}
						}
						error = pdp_eL2/nvis;
					}				
					for(i=0, dL=0.0; i<nvars; ++i){
						dL+=dp[i]*(mu*dp[i]+eab[i]);  //low  predicted reduction
						//attention: eab = -g
					}
					dF=p_eL2-pdp_eL2;//actual reduction

					if(itno == 0 && nLmIterations == 1){
						sigma_spectral = convertStoDenseMatrix(S,Sidxij,Spgm.c_str(),Stxt.c_str());
						S_sigma_max = sigma_spectral[0];
						S_sigma_min = sigma_spectral[sigma_spectral.size() - 1];
						S_cond = S_sigma_max / S_sigma_min;
						rho_ = dF/dL;
						fprintf(fp3,"%f %f %f %f %f %f %f %f %f %f ", initialerror, rho_, mu, S_sigma_max, S_sigma_min, S_cond, 
							eab_inf, relative_state_change, relative_mse_change, Lips);
						for(int k=0;k<sigma_spectral.size()-1;++k){
							fprintf(fp4,"%f ", sigma_spectral[k]);
						}
						fprintf(fp4,"%f\n", sigma_spectral[sigma_spectral.size()-1]);
					}

					if((dF/dL)>0.0)
					{ 
						//sqrt(p_eL2):当前误差的MSE
						//sqrt(pdp_eL2):上一次误差的MSE
						relative_mse_change = (sqrt(p_eL2)-sqrt(pdp_eL2)) / sqrt(p_eL2);
						if((sqrt(p_eL2)-sqrt(pdp_eL2))<m_e3*sqrt(p_eL2)) 
						{	stop=2;		break;	}//relative change of projection error is small
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
						for(i=0; i<nvars; ++i) 
							p[i]=pdp[i];
						for(i=0; i<nobs; ++i) 
							e[i]=hx[i];		
				

						p_eL2=pdp_eL2;
						if((eab_inf <= eps1))
						{	dp_L2=0.0; 		stop=4;		break;	}//maximum value of b is small

						tmp=(2.0*dF/dL-1.0);
						tmp=1.0-tmp*tmp*tmp;
						mu=mu*( (tmp>=1.0/3.0)? tmp : 1.0/3.0 );
						nu=2;

						tTimeIter = (t5 - t0)*0.001;
						tTimeCal += tTimeIter;

						printf( "Iteration=%d  MSE=%0.8lf   LmIters=%d  Pertime=%0.2lf TotalTime=%0.2lf\n", itno, pdp_eL2/nvis, nLmIterations, tTimeIter, tTimeCal ); 					
						if( m_szReport!= NULL )
							fprintf( fpRe, "Iteration %d  Error  %0.8lf\n", itno, pdp_eL2/nvis );
						nIter++;

						//计算梯度的Lipschitz连续性与梯度方向有效性
						double dot_grad_dp = 0.0, grad_sq = 0.0, dp_sq = 0.0, grad_norm, dp_norm;
						for(i=0; i<nvars; ++i){
							dot_grad_dp += eab[i] * dp[i];
							grad_sq += eab[i] * eab[i];
							dp_sq += dp[i] * dp[i];
						}
						grad_norm = sqrt(grad_sq);
						dp_norm = sqrt(dp_sq);
						descent_dir_valid = dot_grad_dp / (grad_norm * dp_norm);

						//已经算出下一个点处的梯度eab_pre
						if(itno > 0){
							double grad_diff_sq = 0.0, grad_diff_norm;
							for(i=0; i<nvars; ++i){
								double diff = eab[i] - eab_pre[i];//相邻两次迭代的负梯度差
								grad_diff_sq += diff * diff;
								//attention: eab = -g
							}
							grad_diff_norm = sqrt(grad_diff_sq);
							Lips = grad_diff_norm / dp_norm;
						}
						if(itno == 0)
							fprintf(fp3,"%f %d\n", descent_dir_valid, nLmIterations);
						break;
					}
					else
					{
						mu*=nu;
						nu *= 2;
					}
				} 		
			}
				
			if(p_eL2<=eps3_sq) stop=5; //total reprojection error is small

			if(!stop && itno > 0){
				sigma_spectral = convertStoDenseMatrix(S,Sidxij,Spgm.c_str(),Stxt.c_str());
				S_sigma_max = sigma_spectral[0];
				S_sigma_min = sigma_spectral[sigma_spectral.size() - 1];
				S_cond = S_sigma_max / S_sigma_min;
				rho_ = dF/dL;
				// double H_inf = compute_H_inf_norm(U, V, W, mu);
				// double H_lambda_max, H_lambda_min, H_cond;
				// computeHMaxSingularValue(U, V, W, H_lambda_max, H_lambda_min, H_cond);
				fprintf(fp3,"%f %f %f %f %f %f %f %f %f %f %f %d\n", pdp_eL2/nvis, rho_, mu, S_sigma_max, S_sigma_min, S_cond,
					eab_inf, relative_state_change, relative_mse_change, Lips, descent_dir_valid, nLmIterations);
				for(int k=0;k<sigma_spectral.size()-1;++k){
					fprintf(fp4,"%f ", sigma_spectral[k]);
				}
				fprintf(fp4,"%f\n", sigma_spectral[sigma_spectral.size()-1]);
			}
			for(i=0; i<nvars; ++i){
				eab_pre[i] = eab[i];
			}
		}
		if(itno>=m_nMaxIter)
			stop=3;	//maximum iteration
		
		//clear memory and print
	// iterstop:
		
		sba_crsm_free(&Uidxij);
		sba_crsm_free(&Sidxij);

		cholmod_free_factor(&m_cholFactorS, &m_cS) ;              
		cholmod_l_free_dense(&m_cholSparseE, &m_cS);
		cholmod_l_free_dense(&m_cholSparseR, &m_cS);
		cholmod_finish (&m_cS) ;  
		free(Ap);
		free(Aii);
		cholmod_free_sparse(&m_cholSparseS, &m_cS) ;
		
		tEnd  = clock();
		tTimeUse = tEnd - tStart;

		//save optimal camera pose and feature
		pba_saveXYZ( m_szCamePose, m_sz3Dpts, p, false );

		free(S);	free(W);	free(U);	free(V);	free(IV);
		free(e);	free(eab);	free(E);   	free(dp);	free(hx);	free(pdp);	
		free(m_KR); free(m_KdA); free(m_KdB); free(m_KdG);	
		//free(m_V);

		printf( "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, p_eL2/nvis, initialerror, nIter+1, nLinear, tTimeUse / CLOCKS_PER_SEC );
		printf( "ParallaxBA reasons is listed as following:\n" );
		printf( "reason 1: relative change of state vector is small\n" );
		printf( "reason 2: relative change of projection error is small\n" );
		printf( "reason 3: maximum iteration\n");
		printf( "reason 4: maximum value of b is small\n");
		printf( "reason 5: total reprojection error is small\n");

		if( m_szReport!=NULL )
		{
			fprintf( fpRe, "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
				nvars, m_n2Dprojs*2, stop, error, initialerror, nIter+1, nLinear, tTimeUse*0.001 );

			fprintf( fpRe, "ParallaxBA reasons is listed as following:\n" );
			fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
			fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
			fprintf( fpRe, "reason 3: maximum iteration\n");
			fprintf( fpRe, "reason 4: maximum value of b is small\n");
			fprintf( fpRe, "reason 5: total reprojection error is small\n");

			fclose(fpRe);
		}	
		fclose(fp3);
		fclose(fp4);
		return true;
	}

	bool PBA::ba_motstr_gn( FixType ft )
	{	
		if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
		{	
			printf( "ParallaxBA: Missing related file, please input them first\n");	
			return false;
		}

		FILE *fpRe = nullptr;
		if( m_szReport != NULL )
		{
			fopen_s(&fpRe, m_szReport, "w ");
			fprintf(fpRe, "%d  poses, %d 3D features, %d projection\n", m_ncams, m_n3Dpts, m_n2Dprojs );
			fprintf(fpRe, "Gauss-Newton is used\n" );
		}

		int i, ii, nft;
		int n = m_n3Dpts, m = m_ncams, cnp = 6, pnp = 3, mnp = 2;
		int nvis, nuis, itno, nobs, nvars, nMaxS;

		double *p = m_motstruct, *x = m_imgpts;	//p pointer refers to unknown parameters, x pointer refers to image coordinate	
		double *U, *V, *W, *e, *eab, *E, *S, *dp; // pointers into U V W E S IV
		double *pa, *pb, *ea, *eb, *dpa, *dpb, *hx, *Ex, *rx; // pointers into p, jac, eab and dp respectively 	double *Ex, *rx; 
		double initialerror = 0, error = 0, lasterror = 0;
		int Usz, Vsz, Wsz, esz, easz, ebsz, Sblsz, Sdim; 

		double tmp, tmp1, init_p_eL2, p_eL2, eab_inf, pdp_eL2, x2, delt2; 
		double dp_L2=DBL_MAX;
		int nu=2, stop=0;
		eab_inf=0.0;		

		bool init = false;
		bool ordering = true;
		double tStart, tEnd, tTimeUse, tPerTime, t0, t1, t2, t3, t4, t5, tTotalTime = 0;
		struct sba_crsm Sidxij, Uidxij;		// S mask and U mask          
		nuis = ba_ConstructSmask( Sidxij, Uidxij );   

		Usz=cnp*cnp; Vsz=pnp * pnp; Wsz=cnp * pnp; 	
		esz=mnp; easz=cnp; ebsz=pnp; Sblsz=cnp*cnp;	Sdim=m * cnp;	
		nvis = m_n2Dprojs;                                            
		nobs=nvis*mnp;                                               
		nvars=m*cnp + n*pnp;                                         

		S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));         
		W	=	(double *)malloc(nvis*Wsz*sizeof(double));            
		U	=	(double *)malloc(nuis*Usz*sizeof(double));            
		V	=	(double *)malloc(n*Vsz*sizeof(double));               
		e	=	(double *)malloc(nobs*sizeof(double));                 
		eab	=	(double *)malloc(nvars*sizeof(double));              
		E	=	(double *)malloc(m*cnp*sizeof(double));	            
		dp	=	(double *)malloc(nvars*sizeof(double));              
		hx	=	(double *)malloc(nobs*sizeof(double));	              
		pa=p; 
		pb=p+m*cnp; 

		ea=eab;
		eb=eab+m*cnp;

		dpa=dp; 
		dpb=dp+m*cnp;
		
		//Select fix axis for Gauss-Newton
		if ( ft == BA_FixDefault )
		{
			tmp = abs(*(m_motstruct+9)-*(m_motstruct+3));
			ft = BA_FixX;									
			tmp1 = abs(*(m_motstruct+10)-*(m_motstruct+4));
			if ( tmp1 > tmp )
			{ft= BA_FixY; tmp = tmp1;	}
			tmp1 = abs(*(m_motstruct+11)-*(m_motstruct+5));
			if ( tmp1 > tmp )
			{ft= BA_FixZ; }
		}
		nft = ft;
		
		cholmod_start (&m_cS) ;  
		int *Ap  = (int*)malloc((m)*sizeof(int));
		int * Aii = (int*)malloc(m_nS*sizeof(int));
		ba_constructAuxCSSGN( Ap, Aii );

		m_cholSparseE = cholmod_zeros( cnp*m-7, 1, CHOLMOD_REAL, &m_cS);
		Ex = (double*)m_cholSparseE->x;
		nMaxS = (m_nS-m_ncams)*36+m_ncams*21;	//maximum non-zero element in S matrix 
		m_cholSparseS = cholmod_allocate_sparse(m_ncams*6-7,m_ncams*6-7,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
		int *Sp, *Si;
		double* Sx = NULL;
		Sp = (int*)m_cholSparseS->p;		//column pointer
		Si = (int*)m_cholSparseS->i;		//row pointer
		
		//Compute initial error and initial reprojection error 
		tStart = clock();
		//dpa[0] = dpa[1] = dpa[2] = dpa[3] = dpa[4] = dpa[5] = 0;//��???
		//dpa[9+nft] = 0;
		pba_cost(p, hx, m_archor );
		p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */

		initialerror = p_eL2/nvis;
		lasterror = initialerror;//Add by zuo
		printf("Initial Error %0.1lf [%0.8lf]\n", p_eL2, initialerror);

		if( m_szReport!= NULL )
			fprintf( fpRe, "Initial Error  %0.8lf\n", initialerror );
		init_p_eL2=p_eL2;	
		
		//Iteration 
		for(itno=0; itno<m_nMaxIter && !stop; ++itno)
		{
			//Setup S matrix, include two step
			memset( U, 0, nuis*Usz*sizeof(double) );
			memset( ea, 0,m*easz*sizeof(double) );
			memset( V, 0, n*Vsz*sizeof(double));
			memset( eb, 0, n*ebsz*sizeof(double));
			memset( W, 0, nvis*Wsz*sizeof(double));

			//Step one; compute W V U directly, don't save each projection image Jacobian 
			t0 = clock();
			if ( m_bRobustKernel)
				pba_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
			else
				pba_jacobian(p, m_archor, &Uidxij, e, U, ea, V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature);

			t1 = clock();

			ba_inverseVGN( U, V, Uidxij ); //compute inverse matrix of V

			//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
			memset( E, 0, m*easz*sizeof(double));
			memset( S, 0, m_nS*Sblsz*sizeof(double) );
			ba_constructSGN( S, E, U, V, W, ea, eb, Sidxij );//S????????E????????
			t2 = clock();
		
			//Solve equation
			ba_constructCSSGN( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init, nft ); //set CSS format using S matrix//break
			for ( ii = 0; ii < 3+nft; ii++ )
				Ex[ii] = E[6+ii];//??2????????

			for ( ii = 0; ii < 6*m-(10+nft); ii++  )
				Ex[3+nft+ii] = E[10+nft+ii];//??????2???Xc????��Yc??Zc

			//Sx????Ex????

			ba_solveCholmodGN( Ap, Aii, init, ordering);//break
			init = true;
			rx = (double*)m_cholSparseR->x;//????????????

			if (m_cS.status == CHOLMOD_NOT_POSDEF)
			{
				printf ( "Cholesky failure, writing debug.txt (Hessian loadable by Octave)" );
				goto iterstop;
			}
			else
			{
				for ( ii = 0; ii < 6*m-(10+nft); ii++ )
					dpa[6+ii] = rx[ii];
				for ( ii = 0; ii < 6*m-(10+nft); ii++ )
					dpa[ii+10+nft] = rx[3+nft+ii];//?��????????????
			}					
			t3 = clock();
			
			//Solve features
			ba_solveFeatures( W, V, ea, eb, dpa, dpb );//break
			
			//update
			for(i=0, dp_L2=0.0; i<nvars; ++i)//nvars??��?????????
			{
				p[i]=p[i] + (tmp=dp[i]);//dpa???dp?��???��???????dpb???dp?��?????????????p???��?????
				dp_L2+=tmp*tmp;//?��?????????????
			}

			//dp_L2=sqrt(dp_L2);
			if (dp_L2<1E-8*1E-8)//break
			{	stop = 1;	goto iterstop;	}	//??????��?????��						

			ba_updateKR( m_KR, m_KdA, m_KdB, m_KdG, m_K, p );
			t4 = clock();
			
			pba_cost(p, hx, m_archor ); 
			pdp_eL2=nrmL2xmy(e, x, hx, nobs); 				
			error = pdp_eL2/nvis;//break ????mse
			t5 = clock();
					
			if ( m_bRobustKernel )
			{
				pdp_eL2 = 0;
				delt2 = m_delt*m_delt;
				if ( m_nRobustType==1)						//Cauchy Kernel Function
				{
					for ( i = 0; i < m_n2Dprojs; i++ )
					{					
						x2 = e[i*2]*e[i*2]+e[i*2+1]*e[i*2+1];
						x2 = delt2*log( x2/delt2 + 1 );
						pdp_eL2 += x2;
					}
				}
				else										//Huber Kernel Function
				{
					for ( i = 0; i < m_n2Dprojs; i++ )
					{					
						x2 = e[i*2]*e[i*2]+e[i*2+1]*e[i*2+1];

						if (x2 <= delt2)  // inlier
							x2 = x2;
						else  // outlier
							x2 = 2*sqrt(x2)*m_delt - delt2;

						pdp_eL2 += x2;
					}
				}
				error = pdp_eL2/nvis;
			}
			
			if (abs(error-lasterror) < 1E-9)//change by zuo
			{	stop = 2;	goto iterstop;		}//?????mse???????????????
			
			tPerTime = (t5-t0)*0.001;//0.001?????????????????
			tTotalTime += tPerTime;//???��???????

			printf( "Iteration=%d  MSE=%0.8lf Pertime=%0.2lf TotalTime=%0.2lf\n", itno, pdp_eL2/nvis, tPerTime, tTotalTime );
			//printf("setup: %0.2lf  solve: %0.2lf  update: %0.2lf cost: %0.2lf totoal: %0.2lf\n",
			//(t2-t0)/1000.0,(t3-t2)/1000.0, (t4-t3)/1000.0, (t5-t4)/1000.0,(t5-t0)/1000.0);
			//printf("\n");

			if( m_szReport != NULL )
				fprintf( fpRe, "Iteration %d  Error  %0.8lf\n", itno, pdp_eL2/nvis );

			lasterror = error;		//change by zuo
		}
		if(itno>=m_nMaxIter) stop=3;
		
		//clear memory and print
	iterstop:
		cholmod_finish (&m_cS) ;  
		sba_crsm_free(&Uidxij);
		sba_crsm_free(&Sidxij);

		cholmod_free_factor(&m_cholFactorS, &m_cS) ;              
		cholmod_l_free_dense(&m_cholSparseE, &m_cS);
		cholmod_l_free_dense(&m_cholSparseR, &m_cS);
		free(Ap);
		free(Aii);
		cholmod_free_sparse(&m_cholSparseS, &m_cS) ;

		tEnd  = clock();
		tTimeUse = tEnd - tStart;
		

		//save optimal camera pose and feature
		pba_saveXYZ( m_szCamePose, m_sz3Dpts, p );

		free(S);	free(W);	free(U);	free(V);	
		free(e);	free(eab);	free(E);   	free(dp);	free(hx);	
		free(m_KR); free(m_KdA);free(m_KdB);free(m_KdG);

		printf( "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
		printf( "ParallaxBA reasons is listed as following:\n " );
		printf( "reason 1: relative change of state vector is small\n" );
		printf( "reason 2: relative change of projection error is small\n " );
		printf( "reason 3: maximum iteration\n");

		if( m_szReport != NULL )
		{
			fprintf( fpRe, "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
				nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
			fprintf( fpRe, "ParallaxBA reasons is listed as following:\n" );
			fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
			fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
			fprintf( fpRe, "reason 3: maximum iteration\n");
			fclose(fpRe);
		}
		
		return true;
	}

	void Jocobian_pre(std::vector<int> & row_indices1,std::vector<int> & col_indices1,std::vector<double> & values1,
		std::vector<int> & row_indices2,std::vector<int> & col_indices2,std::vector<double> & values2,
		Eigen::MatrixXd & J_dense,
		int & J_row,int nP,int nM,int nN, int nF,
		double* pAM, double* pAA, double* pPA, double* pPB){
		J_dense(J_row*2+0,nP*6+0)=pPA[0];
		J_dense(J_row*2+0,nP*6+1)=pPA[1];
		J_dense(J_row*2+0,nP*6+2)=pPA[2];
		J_dense(J_row*2+0,nP*6+3)=pPA[3];
		J_dense(J_row*2+0,nP*6+4)=pPA[4];
		J_dense(J_row*2+0,nP*6+5)=pPA[5];

		J_dense(J_row*2+1,nP*6+0)=pPA[6];
		J_dense(J_row*2+1,nP*6+1)=pPA[7];
		J_dense(J_row*2+1,nP*6+2)=pPA[8];
		J_dense(J_row*2+1,nP*6+3)=pPA[9];
		J_dense(J_row*2+1,nP*6+4)=pPA[10];
		J_dense(J_row*2+1,nP*6+5)=pPA[11];
		if(nP==nM){
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+0);values1.push_back(pPA[0]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+1);values1.push_back(pPA[1]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+2);values1.push_back(pPA[2]);

			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+0);values1.push_back(pPA[6]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+1);values1.push_back(pPA[7]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+2);values1.push_back(pPA[8]);

			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+0);values2.push_back(pPB[0]);
			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+1);values2.push_back(pPB[1]);

			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+0);values2.push_back(pPB[3]);
			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+1);values2.push_back(pPB[4]);
		}else if(nP==nN){

			J_dense(J_row*2+0,nM*6+3)=pAM[0];
			J_dense(J_row*2+0,nM*6+4)=pAM[1];
			J_dense(J_row*2+0,nM*6+5)=pAM[2];
			J_dense(J_row*2+1,nM*6+3)=pAM[3];
			J_dense(J_row*2+1,nM*6+4)=pAM[4];
			J_dense(J_row*2+1,nM*6+5)=pAM[5];
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+0);values1.push_back(pPA[0]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+1);values1.push_back(pPA[1]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+2);values1.push_back(pPA[2]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+3);values1.push_back(pPA[3]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+4);values1.push_back(pPA[4]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+5);values1.push_back(pPA[5]);

			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nM*6+3);values1.push_back(pAM[0]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nM*6+4);values1.push_back(pAM[1]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nM*6+5);values1.push_back(pAM[2]);

			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+0);values1.push_back(pPA[6]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+1);values1.push_back(pPA[7]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+2);values1.push_back(pPA[8]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+3);values1.push_back(pPA[9]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+4);values1.push_back(pPA[10]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+5);values1.push_back(pPA[11]);

			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nM*6+3);values1.push_back(pAM[3]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nM*6+4);values1.push_back(pAM[4]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nM*6+5);values1.push_back(pAM[5]);

			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+0);values2.push_back(pPB[0]);
			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+1);values2.push_back(pPB[1]);
			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+2);values2.push_back(pPB[2]);

			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+0);values2.push_back(pPB[3]);
			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+1);values2.push_back(pPB[4]);
			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+2);values2.push_back(pPB[5]);
		}else{

			J_dense(J_row*2+0,nM*6+3)=pAM[0];
			J_dense(J_row*2+0,nM*6+4)=pAM[1];
			J_dense(J_row*2+0,nM*6+5)=pAM[2];
			J_dense(J_row*2+1,nM*6+3)=pAM[3];
			J_dense(J_row*2+1,nM*6+4)=pAM[4];
			J_dense(J_row*2+1,nM*6+5)=pAM[5];

			J_dense(J_row*2+0,nN*6+3)=pAA[0];
			J_dense(J_row*2+0,nN*6+4)=pAA[1];
			J_dense(J_row*2+0,nN*6+5)=pAA[2];
			J_dense(J_row*2+1,nN*6+3)=pAA[3];
			J_dense(J_row*2+1,nN*6+4)=pAA[4];
			J_dense(J_row*2+1,nN*6+5)=pAA[5];
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+0);values1.push_back(pPA[0]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+1);values1.push_back(pPA[1]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+2);values1.push_back(pPA[2]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+3);values1.push_back(pPA[3]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+4);values1.push_back(pPA[4]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nP*6+5);values1.push_back(pPA[5]);

			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nM*6+3);values1.push_back(pAM[0]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nM*6+4);values1.push_back(pAM[1]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nM*6+5);values1.push_back(pAM[2]);

			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nN*6+3);values1.push_back(pAA[0]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nN*6+4);values1.push_back(pAA[1]);
			row_indices1.push_back(J_row*2+0);col_indices1.push_back(nN*6+5);values1.push_back(pAA[2]);

			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+0);values1.push_back(pPA[6]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+1);values1.push_back(pPA[7]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+2);values1.push_back(pPA[8]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+3);values1.push_back(pPA[9]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+4);values1.push_back(pPA[10]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nP*6+5);values1.push_back(pPA[11]);

			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nM*6+3);values1.push_back(pAM[3]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nM*6+4);values1.push_back(pAM[4]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nM*6+5);values1.push_back(pAM[5]);

			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nN*6+3);values1.push_back(pAA[3]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nN*6+4);values1.push_back(pAA[4]);
			row_indices1.push_back(J_row*2+1);col_indices1.push_back(nN*6+5);values1.push_back(pAA[5]);

			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+0);values2.push_back(pPB[0]);
			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+1);values2.push_back(pPB[1]);
			row_indices2.push_back(J_row*2+0);col_indices2.push_back(3*nF+2);values2.push_back(pPB[2]);

			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+0);values2.push_back(pPB[3]);
			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+1);values2.push_back(pPB[4]);
			row_indices2.push_back(J_row*2+1);col_indices2.push_back(3*nF+2);values2.push_back(pPB[5]);
		}
		++J_row;
	}

	void PBA::pba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
	{
		int nM, nN;	  
		int i, j, ii, jj, k;
		int cnp, pnp, mnp;
		double *pa, *pb, *ppt;	
		double *ppUpa, *pea, *peb, *pe, *pV, *pW;
		int pos, pos2, nP, nF, numfea, cur;
		double sum;
		double tmp = 0;

		double pAM[6], pAA[6], pPA[12], pPB[6];

		cnp = 6; pnp = 3; mnp = 2;
		pa=p; pb=p+m*cnp;

		pos2 = 0;

		// std::vector<int> row_indices1,col_indices1, row_indices2,col_indices2;
		// std::vector<double> values1, values2;
		// int J_row=0;
		// Eigen::MatrixXd J_dense = Eigen::MatrixXd::Zero(IBA::m_n2Dprojs*2, IBA::m_ncams*6);


		// Eigen::MatrixXd U_full = Eigen::MatrixXd::Zero(fullSize,fullSize);
		

		for ( i = 0; i < m_n3Dpts; i++ )
		{
			numfea = m_archor[i*3];//??????????????
			cur = 0;
			nF = i;
			//printf("%d %d\n", m_archorSort[i * 2], m_archorSort[i * 2 + 1]);
			
			for ( j = 0; j < numfea; j++ )
			{
				nP = nphoto[pos2];//????????
				ppt=pb + nF*pnp;//?????????  
				pe= e + pos2*2;		//error

				nM = archor[nF*3+1];
				nN = archor[nF*3+2];

				// fprintf(fp, "\n");

				
				pba_jacobianEachPts(m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB);
				// Jocobian_pre(row_indices1, col_indices1, values1, row_indices2, col_indices2, values2, J_dense, J_row, nP, nM, nN, nF, pAM, pAA, pPA, pPB);

				//U
				pos = sba_crsm_elmidx(Uidxij, nP, nP);//????U??????nP??nP?��???????????��???????
				//printf("%d  ", pos);
				ppUpa = U + pos*(6*6);//?????????????????��????????????????????????????????????????????

				//printf("\n\n\n\n");
				for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//????????????????????
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pPA[k*6+ii]*pPA[k*6+jj];//U^T*U
					ppUpa[ii*6+jj] += sum;
					// U_full(nP*6+ii,nP*6+jj)+=sum;
					//printf("%d ", ppUpa[ii * 6 + jj]);
				}
				
				//ea
				pea = ea + nP*6;
				for ( ii = 0; ii < 6; ii++ )
				{
					for( jj = 0, sum = 0; jj < 2; jj++ )
						sum += pPA[jj*6+ii] * pe[jj];
					pea[ii] += sum;
				}
				//V
				pV = V + nF*3*3;
				for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pPB[k*3+ii]*pPB[k*3+jj];
					pV[ii*3+jj] += sum;

				}
				//eb
				peb = eb + nF*3;
				for ( ii = 0; ii < 3; ii++ )
				{
					for( k = 0, sum = 0; k < 2; k++ )
						sum += pPB[k*3+ii] * pe[k];
					peb[ii] += sum;
				}
				//W
				pW = W + pos2 * 3*6;
				for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pPA[k*6+ii] * pPB[k*3+jj];
					pW[ii*3+jj] += sum;
				}
				if ( nP == nM )
				{
					cur++;
					pos2++;
					continue;
				}
				else
				if( nP == nN )
				{
					//U		  
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						ppUpa[(ii+3)*6+3+jj] += sum;
						// U_full(nM*6+3+ii,nM*6+3+jj)+=sum;
					}
					if( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum;
							// U_full(nM*6+3+ii,nP*6+jj)+=sum;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum;
							// U_full(nP*6+ii,nM*6+3+jj)+=sum;
						}
					}

					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];//?????????????��?????????��??��?????????? * ???????????????��?
						pea[ii+3] += sum;
					}
					//W
					pos = pos2 + (m_archorSort[i*2]-j);//???��??pos2++????????????-j
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum;
					}
				}
				else
				{
					//U
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						}
						
						ppUpa[(ii+3)*6+3+jj] += sum;
						// U_full(nM*6+3+ii,nM*6+3+jj)+=sum;
					}

					pos = sba_crsm_elmidx(Uidxij, nN, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAA[k*3+ii] * pAA[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum;
						// U_full(nN*6+3+ii,nN*6+3+jj)+=sum;
					}
					if ( nM < nN )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAM[k*3+ii] * pAA[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum;
							// U_full(nM*6+3+ii,nN*6+3+jj)+=sum;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAA[k*3+ii] * pAM[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum;
							// U_full(nN*6+3+ii,nM*6+3+jj)+=sum;
						}
					}
					if ( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum;
							// U_full(nM*6+3+ii,nP*6+jj)+=sum;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum;
							// U_full(nP*6+ii,nM*6+3+jj)+=sum;
						}
					}
					
					if ( nN < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum;
							// U_full(nN*6+3+ii,nP*6+jj)+=sum;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAA[k*3+jj];
							ppUpa[ii*6+jj+3] += sum;
							// U_full(nP*6+ii,nN*6+3+jj)+=sum;
						}
					}
					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						pea[ii+3] += sum;
					}
					pea = ea + nN*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pe[k];
						pea[ii+3] += sum;
					}  
					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum;
					}
					pos = pos2 + (m_archorSort[i*2+1]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum;
					}
				}
				pos2++;
				cur++;
			}
		}
		
		// int blockSize = 6;
		// int nCams = m_ncams;
		// int fullSize = blockSize * nCams;
		// Eigen::MatrixXd U_full1 = Eigen::MatrixXd::Zero(fullSize,fullSize);
		// for (int ci = 0; ci < nCams; ci++) {
		// 	for (int cj = ci; cj < nCams; cj++) {
		// 		int pos = sba_crsm_elmidx(Uidxij, ci, cj);
		// 		if (pos < 0) continue; // 没有该块
		// 		double* ppUpa = U + pos * blockSize * blockSize;

		// 		if(cj==ci){
		// 			// 每个ppUpa�?6×6上三角块
		// 			for (int ii = 0; ii < blockSize; ii++) {
		// 				for (int jj = ii; jj < blockSize; jj++) {
		// 					double val = ppUpa[ii * blockSize + jj];
		// 					U_full1(ci * blockSize + ii, cj * blockSize + jj) = val;
		// 				}
		// 			}
		// 		}else{
		// 			for (int ii = 0; ii < blockSize; ii++) {
		// 				for (int jj = 0; jj < blockSize; jj++) {
		// 					double val = ppUpa[ii * blockSize + jj];
		// 					U_full1(ci * blockSize + ii, cj * blockSize + jj) = val;
		// 				}
		// 			}
		// 		}

		// 	}
		// }
		// for (int r = 0; r < fullSize; r++) {
		// 	for (int c = r + 1; c < fullSize; c++) {
		// 		if (U_full1(r, c) != 0)
		// 			U_full1(c, r) = U_full1(r, c);
		// 	}
		// }
		// FILE* fp1 = fopen("E:/zuo/projects/PBA_U.pgm", "wb");
		// fprintf(fp1, "P5\n%d %d\n255\n", fullSize, fullSize);
		// for(int i=0;i<fullSize;i++){
		// 	for(int j=0;j<fullSize;j++){
		// 		unsigned char pixel = U_full1(i, j) ? 0 : 255;
		// 		fwrite(&pixel, 1, 1, fp1);
		// 	}
		// 	// fprintf(fp1, "\n");
		// }
		// fclose(fp1);
		// FILE* fp2;
		// const char* fn2 = "E:/zuo/projects/PBA_U.txt";
		// errno_t err2 = fopen_s(&fp2, fn2, "w");
		// for(int i=0;i<fullSize;i++){
		// 	for(int j=0;j<fullSize;j++){
		// 		fprintf(fp2,"%f ",U_full1(i,j));
		// 	}
		// 	fprintf(fp2,"\n");
		// }
		// fclose(fp2);
		// int rows = J_row*2;
		// int cols1 = IBA::m_ncams * 6;
		// int cols2 = IBA::m_n3Dpts * 3;
        // int hhh = IBA::m_n2Dprojs*2;



		// for (int r = 0; r < fullSize; r++) {
		// 	for (int c = r + 1; c < fullSize; c++) {
		// 		if (U_full(r, c) != 0)
		// 			U_full(c, r) = U_full(r, c);
		// 	}
		// }



		// Eigen::SparseMatrix<double> Jc=buildSparseMatrixCOO(rows, cols1, row_indices1, col_indices1, values1);
		// Eigen::SparseMatrix<double> U_ = Jc.transpose() * Jc;
		// Eigen::MatrixXd U_dense = Eigen::MatrixXd(U_);
		// // Eigen::MatrixXd U_dense = J_dense.transpose()*J_dense;
		// FILE* fp;
		// const char* fn = "E:/zuo/projects/U1.txt";
		// errno_t err = fopen_s(&fp, fn, "w");
		// for(int i=0;i<cols1;i++){
		// 	for(int j=0;j<cols1;j++){
		// 		fprintf(fp,"%f ",U_dense(i,j));
		// 	}
		// 	fprintf(fp,"\n");
		// }
		// fclose(fp);

		// Eigen::SparseMatrix<double> Jf=buildSparseMatrixCOO(rows, cols2, row_indices2, col_indices2, values2);
		// Eigen::SparseMatrix<double> V_ = Jf.transpose() * Jf;
		// Eigen::SparseMatrix<double> W_ = Jc.transpose() * Jf;
		// Eigen::SparseMatrix<double> Wt_ = W_.transpose();
		// Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
		// cg.compute(V_);
		// Eigen::SparseMatrix<double> V_inv = cg.solve(Wt_);
		// Eigen::SparseMatrix<double> SS_ = W_ * V_inv;
		// Eigen::SparseMatrix<double> S_ = U_ - SS_;
		// Eigen::MatrixXd S_dense = Eigen::MatrixXd(S_);
		// FILE* fp;
		// const char* fn = "E:/zuo/projects/S2.txt";
		// errno_t err = fopen_s(&fp, fn, "w");
		// for(int i=0;i<cols1;i++){
		// 	for(int j=0;j<cols1;j++){
		// 		fprintf(fp,"%d",S_dense(i,j)?1:0);
		// 	}
		// 	fprintf(fp,"\n");
		// }
		// fclose(fp);


		


	}

	void PBA::pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
	{
		int nM, nN;	  
		int i, j, ii, jj, k;
		int cnp, pnp, mnp;
		double *pa, *pb, *ppt;	
		double *ppUpa, *pea, *peb, *pe, *pV, *pW;
		int pos, pos2, nP, nF, numfea, cur;
		double sum, x2;
		double tmp = 0;

		double pAM[6], pAA[6], pPA[12], pPB[6];

		cnp = 6; pnp = 3; mnp = 2;
		pa=p; pb=p+m*cnp;

		pos2 = 0;
		for ( i = 0; i < m_n3Dpts; i++ )
		{
			numfea = m_archor[i*3];
			cur = 0;
			nF = i;
			for ( j = 0; j < numfea; j++ )
			{
				nP = nphoto[pos2];
				ppt=pb + nF*pnp;      
				pe= e + pos2*2;		//error

				nM = archor[nF*3+1];
				nN = archor[nF*3+2];

				double delt2 = m_delt*m_delt;
				if( m_nRobustType == 1)
				{
					x2 = pe[0]*pe[0]+pe[1]*pe[1];
					x2 = 1.0/( x2/delt2+1 );
				}
				else
				{
					x2 = pe[0]*pe[0]+pe[1]*pe[1];

					if (sqrt(x2) < m_delt) //inlier
						x2 = 1;
					else // outliers 
						x2 = m_delt/sqrt(x2);
				}			

				pba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );

				//U
				pos = sba_crsm_elmidx(Uidxij, nP, nP);		
				ppUpa = U + pos*(6*6);

				for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//diag
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pPA[k*6+ii]*pPA[k*6+jj];
					//ppUpa[ii*6+jj] += sum;
					ppUpa[ii*6+jj] += sum*x2;
				}

				//ea
				pea = ea + nP*6;
				for ( ii = 0; ii < 6; ii++ )
				{
					for( jj = 0, sum = 0; jj < 2; jj++ )
						sum += pPA[jj*6+ii] * pe[jj];
					//pea[ii] += sum;
					pea[ii] += sum*x2;
				}
				//V
				pV = V + nF*3*3;
				for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pPB[k*3+ii]*pPB[k*3+jj];
					//pV[ii*3+jj] += sum;
					pV[ii*3+jj] += sum*x2;

				}
				//eb
				peb = eb + nF*3;
				for ( ii = 0; ii < 3; ii++ )
				{
					for( k = 0, sum = 0; k < 2; k++ )
						sum += pPB[k*3+ii] * pe[k];
					//peb[ii] += sum;
					peb[ii] += sum*x2;
				}
				//W
				pW = W + pos2 * 3*6;
				for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
				{
					for ( k = 0, sum = 0; k < 2; k++ )
						sum += pPA[k*6+ii] * pPB[k*3+jj];
					//pW[ii*3+jj] += sum;
					pW[ii*3+jj] += sum*x2;
				}

				if ( nP == nM )
				{
					cur++;
					pos2++;
					continue;
				}
				else
					if( nP == nN )
					{
						//U		  
						pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pAM[k*3+jj];
							//ppUpa[(ii+3)*6+3+jj] += sum;
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}

						if( nM < nP )
						{
							pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
									sum += pAM[k*3+ii] * pPA[k*6+jj];
								//ppUpa[(ii+3)*6+jj] += sum;
								ppUpa[(ii+3)*6+jj] += sum*x2;
							}
						}
						else
						{
							pos = sba_crsm_elmidx(Uidxij, nP, nM );		
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
									sum +=  pPA[k*6+ii]*pAM[k*3+jj];
								//ppUpa[ii*6+jj+3] += sum;
								ppUpa[ii*6+jj+3] += sum*x2;
							}
						}

						//ea
						pea = ea + nM*6;
						for ( ii = 0; ii < 3; ii++ )
						{
							for( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pe[k];
							//pea[ii+3] += sum;
							pea[ii+3] += sum*x2;
						}
						//W
						pos = pos2 + (m_archorSort[i*2]-j);
						pW = W + pos*6*3;
						for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPB[k*3+jj];
							//pW[(ii+3)*3+jj] += sum;
							pW[(ii+3)*3+jj] += sum*x2;
						}
					}
					else
					{
						//U
						pos = sba_crsm_elmidx(Uidxij, nM, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAM[k*3+ii] * pAM[k*3+jj];
							}
							//ppUpa[(ii+3)*6+3+jj] += sum;
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}

						pos = sba_crsm_elmidx(Uidxij, nN, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAA[k*3+ii] * pAA[k*3+jj];
							}
							//ppUpa[(ii+3)*6+3+jj] += sum;
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}

						if ( nM < nN )
						{
							pos = sba_crsm_elmidx(Uidxij, nM, nN );		
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
								{
									sum += pAM[k*3+ii] * pAA[k*3+jj];
								}
								//ppUpa[(ii+3)*6+3+jj] += sum;
								ppUpa[(ii+3)*6+3+jj] += sum*x2;
							}
						}	
						else
						{
							pos = sba_crsm_elmidx(Uidxij, nN, nM );		
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
								{
									sum += pAA[k*3+ii] * pAM[k*3+jj];
								}
								//ppUpa[(ii+3)*6+3+jj] += sum;
								ppUpa[(ii+3)*6+3+jj] += sum*x2;
							}
						}

						if ( nM < nP )
						{
							pos = sba_crsm_elmidx(Uidxij, nM, nP );		
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
									sum += pAM[k*3+ii] * pPA[k*6+jj];
								//ppUpa[(ii+3)*6+jj] += sum;
								ppUpa[(ii+3)*6+jj] += sum*x2;
							}
						}
						else
						{
							pos = sba_crsm_elmidx(Uidxij, nP, nM );		
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
									sum +=  pPA[k*6+ii]*pAM[k*3+jj];
								//ppUpa[ii*6+jj+3] += sum;
								ppUpa[ii*6+jj+3] += sum*x2;
							}
						}

						if ( nN < nP )
						{
							pos = sba_crsm_elmidx(Uidxij, nN, nP );		
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
									sum += pAA[k*3+ii] * pPA[k*6+jj];
								//ppUpa[(ii+3)*6+jj] += sum;
								ppUpa[(ii+3)*6+jj] += sum*x2;
							}
						}	
						else
						{
							pos = sba_crsm_elmidx(Uidxij, nP, nN );		
							ppUpa = U + pos*(6*6);
							for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
							{
								for ( k = 0, sum = 0; k < 2; k++ )
									sum +=  pPA[k*6+ii]*pAA[k*3+jj];
								//ppUpa[ii*6+jj+3] += sum;
								ppUpa[ii*6+jj+3] += sum*x2;
							}
						}


						//ea
						pea = ea + nM*6;
						for ( ii = 0; ii < 3; ii++ )
						{
							for( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pe[k];
							//pea[ii+3] += sum;
							pea[ii+3] += sum*x2;
						}

						pea = ea + nN*6;
						for ( ii = 0; ii < 3; ii++ )
						{
							for( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pe[k];
							//pea[ii+3] += sum;
							pea[ii+3] += sum*x2;
						}  

						//W
						pos = pos2 + (m_archorSort[i*2]-j);
						pW = W + pos*6*3;
						for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPB[k*3+jj];
							//pW[(ii+3)*3+jj] += sum;
							pW[(ii+3)*3+jj] += sum*x2;
						}

						pos = pos2 + (m_archorSort[i*2+1]-j);
						pW = W + pos*6*3;
						for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pPB[k*3+jj];
							//pW[(ii+3)*3+jj] += sum;
							pW[(ii+3)*3+jj] += sum*x2;
						}
					}

					pos2++;
					cur++;

			}
		}
	}
	void PBA::pba_cost(double* p, double* hx, int* archor)
	{
		int i;
		int cnp, pnp, mnp;
		double* pa, * pb, * ppt, * pmeas;
		int nF, nP, nM, nN;
		int m;

		cnp = 6, pnp = 3, mnp = 2;
		m = m_ncams;
		pa = p; pb = p + m * cnp;

		for (i = 0; i < m_n2Dprojs; i++)
		{
			nF = m_feature[i];
			nP = m_photo[i];

			ppt = pb + nF * pnp;
			pmeas = hx + i * mnp; // set pmeas to point to hx_ij ?��?????????

			nM = archor[nF * 3 + 1];
			nN = archor[nF * 3 + 2];

			pba_reprojectEachPts(m_KR, pa, ppt, nM, nN, nP, pmeas);
			//printf("%f %f\n", pmeas[0], pmeas[1]);
		}
	}
	bool PBA::pba_initializeOtheArchors_Mindw(
		double* imgpts,
		int* photo,
		double* camera,
		double* K,
		double* feature,
		int* archorSort,
		int nfeacout,
		int nOI,
		int FID)
	{
		static int i = 0;
		double dw = feature[2];                   //????
		double dwNew;
		double dminw = dw;
		int   nNewI = 0;
		bool bAdjust = false;
		double dDot, dDisM, dDisA;

		if (dw < MAXARCHOR)
		{
			//current archor vector 
			int nO = photo[nOI];

			Vector3d  xO;
			if (m_bProvideXYZ)
			{
				xO(0) = m_XYZ[FID * 3] - *(camera + nO * 6 + 3);
				xO(1) = m_XYZ[FID * 3 + 1] - *(camera + nO * 6 + 4);
				xO(2) = m_XYZ[FID * 3 + 2] - *(camera + nO * 6 + 5);
			}
			else
			{
				double* ptr1 = m_KR + nO * 9;
				Matrix3d  AO;
				AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

				double matxO[3];
				matxO[0] = *(imgpts + nOI * 2);
				matxO[1] = *(imgpts + nOI * 2 + 1);
				matxO[2] = 1;

				Vector3d  bO(matxO);
				xO = AO.colPivHouseholderQr().solve(bO);
			}

			double dDAngle = atan2(xO(0), xO(2));
			double dHAngle = atan2(xO(1), sqrt(xO(0) * xO(0) + xO(2) * xO(2)));

			for (i = 0; i < nfeacout; i++)
			{
				//Main Archor Vector
				int nM = photo[i];                             //?????
				Vector3d  xM;

				if (m_bProvideXYZ)
				{
					xM(0) = m_XYZ[FID * 3] - *(camera + nM * 6 + 3);
					xM(1) = m_XYZ[FID * 3 + 1] - *(camera + nM * 6 + 4);
					xM(2) = m_XYZ[FID * 3 + 2] - *(camera + nM * 6 + 5);
				}
				else
				{
					double* ptr2 = m_KR + nM * 9;
					Matrix3d  AM;
					AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

					double matxM[3];
					matxM[0] = *(imgpts + i * 2);
					matxM[1] = *(imgpts + i * 2 + 1);
					matxM[2] = 1;

					Vector3d  bM(matxM);
					xM = AM.colPivHouseholderQr().solve(bM);
				}

				//Parallax angle between current archor and main archor
				dDot = xM(0) * xO(0) + xM(1) * xO(1) + xM(2) * xO(2);
				dDisM = sqrt(xM(0) * xM(0) + xM(1) * xM(1) + xM(2) * xM(2));
				dDisA = sqrt(xO(0) * xO(0) + xO(1) * xO(1) + xO(2) * xO(2));

				if (dDot / (dDisM * dDisA) > 1)
					dwNew = 0;
				else if (dDot / (dDisM * dDisA) < -1)
					dwNew = PI;
				else
					dwNew = acos(dDot / (dDisM * dDisA));

				if (dwNew < dminw)
				{
					dminw = dwNew;
					archorSort[0] = nOI;   //?????��???
					archorSort[1] = i;     //????��???
					feature[0] = dDAngle;
					feature[1] = dHAngle;
					feature[2] = dminw;
					bAdjust = true;
				}
			}
		}
		return bAdjust;
	}

	//--------------------------------------------------------
	int PBA::pba_angle2xytGN( double *p )
	{
		static int i, j;
		double* pAngle;
		double xj[3], xk[3];
		double Tik[3];
		int nM, nA;
		double Dik;
		double w, w2;
		static double aa;
		for ( i = 0; i < m_n3Dpts; i++ )
		{
			pAngle = p + m_ncams*6 + i*3;
			w = *(p + m_ncams*6 + i*3 + 2);

			if ( i == 408741 )
			{
				aa = w;
			}

			xj[0] = sin(*(pAngle)) * cos(*(pAngle+1));
			xj[1] = sin(*(pAngle+1));
			xj[2] = cos(*(pAngle)) * cos(*(pAngle+1));

			nM = *(m_archor+i*3+1);
			nA = *(m_archor+i*3+2);

			Tik[0] = -*(p+nM*6+3)+*(p+nA*6+3);
			Tik[1] = -*(p+nM*6+4)+*(p+nA*6+4);
			Tik[2] = -*(p+nM*6+5)+*(p+nA*6+5);
			
			Dik = sqrt( Tik[0]*Tik[0] + Tik[1]*Tik[1] + Tik[2]*Tik[2] );

			w2 = acos( (xj[0]*Tik[0]+xj[1]*Tik[1]+xj[2]*Tik[2])/Dik );

			xk[0] = ( Dik * sin(w2+w) * xj[0] )/sin(w);
			xk[1] = ( Dik * sin(w2+w) * xj[1] )/sin(w);
			xk[2] = ( Dik * sin(w2+w) * xj[2] )/sin(w);

			*(p+m_ncams*6+i*3)   = *(p+nM*6+3) + xk[0];
			*(p+m_ncams*6+i*3+1) = *(p+nM*6+4) + xk[1];
			*(p+m_ncams*6+i*3+2) = *(p+nM*6+5) + xk[2];
		}

		return 1;
	}

	void PBA::pba_saveXYZ(const char* camera, const char* sz3Dpt, double* p, bool gn)
	{
		static int i = 0;
		double dx, dy, dz;
		FILE *fp = nullptr, *fpc = nullptr;
		
		//transform from angle to xyz
		if (gn)
			pba_angle2xytGN(p);
		else
			pba_angle2xytLM(p);
		
		//save camera poss
		if (camera != NULL)
		{
			fopen_s(&fpc, camera, "w");

			for (i = 0; i < m_ncams; i++)
			{
				fprintf( fpc, "%0.5lf     %0.5lf      %0.5lf     %0.5lf       %0.5lf     %0.5lf\n",
					*(p + i * 6), *(p + i * 6 + 1), *(p + i * 6 + 2), *(p + i * 6 + 3), *(p + i * 6 + 4), *(p + i * 6 + 5));;
			}
			fclose(fpc);
		}	
			
		//save features xyz
		if (sz3Dpt != NULL)
		{
			fopen_s(&fp, sz3Dpt, "w");
			fprintf(fp, "%s\n", "ply");
			fprintf(fp, "%s\n", "format ascii 1.0");
			fprintf(fp, "%s %d\n", "element vertex", m_n3Dpts);
			fprintf(fp, "%s\n", "property float x");
			fprintf(fp, "%s\n", "property float y");
			fprintf(fp, "%s\n", "property float z");
			fprintf(fp, "%s\n", "end_header");
			for (i = 0; i < m_n3Dpts; i++)
			{
				dx = *(p+m_ncams*6+i*3);
				dy = *(p+m_ncams*6+i*3+1);
				dz = *(p+m_ncams*6+i*3+2);
				fprintf(fp, "%0.5lf     %0.5lf     %0.5lf\n", dx, dy, dz);
			}
			fclose(fp);
		}

	}

	int	PBA::pba_angle2xytLM(double* p)
	{
		static int i, j;
		double* pAngle;
		double xj[3], xk[3];
		double Tik[3];
		int nM, nA;
		double Dik;
		double w, w2;
		int cnp = 6, pnp = 3;
		for ( i = 0; i < m_n3Dpts; i++ )
		{
			pAngle = p + m_ncams*6 + i*3;
			w = *(p + m_ncams*6 + i*3 + 2);

			xj[0] = sin(*(pAngle)) * cos(*(pAngle+1));//x
			xj[1] = sin(*(pAngle+1));//y
			xj[2] = cos(*(pAngle)) * cos(*(pAngle+1));//z

			nM = *(m_archor+i*3+1);
			nA = *(m_archor+i*3+2);

			Tik[0] = -*(p+nM*6+3)+*(p+nA*6+3);
			Tik[1] = -*(p+nM*6+4)+*(p+nA*6+4);
			Tik[2] = -*(p+nM*6+5)+*(p+nA*6+5);

			Dik = sqrt( Tik[0]*Tik[0] + Tik[1]*Tik[1] + Tik[2]*Tik[2] );

			w2 = acos( (xj[0]*Tik[0]+xj[1]*Tik[1]+xj[2]*Tik[2])/Dik );

			//?????????????��???xyz????
			xk[0] = ( Dik * sin(w2+w) * xj[0] )/sin(w);
			xk[1] = ( Dik * sin(w2+w) * xj[1] )/sin(w);
			xk[2] = ( Dik * sin(w2+w) * xj[2] )/sin(w);

			//??????????xyz????
			*(p+m_ncams*6+i*3)   = *(p+nM*6+3) + xk[0];
			*(p+m_ncams*6+i*3+1) = *(p+nM*6+4) + xk[1];
			*(p+m_ncams*6+i*3+2) = *(p+nM*6+5) + xk[2];
		}

		return 1;
	}
#else
	#error "Unsupported platform!"
#endif


PBA::PBA(void)
{
	m_e1 = m_e2 = m_e3 = 1E-6;
	m_e4 = 0;	
	m_bProvideXYZ = false;
	m_bFocal = false;
	m_szCameraInit = m_szFeatures = m_szCalibration = m_szXYZ = m_sz3Dpts = m_szCamePose = m_szReport = NULL;

	m_nMaxIter = 100;
	m_Tau      = 1E-6;	
	m_bRobustKernel = false;
	m_nRobustType = 1;
	m_bsolverLM = true;
	m_bsolverGN = false;
	m_delt      = 1;
}


PBA::~PBA(void)
{
	if ( m_szCameraInit!= NULL)	free(m_szCameraInit);
	if ( m_szFeatures!= NULL)	free(m_szFeatures);
	if ( m_szCalibration!= NULL)	free(m_szCalibration);
	if ( m_szXYZ!= NULL)		free(m_szXYZ);
	if ( m_szCamePose!= NULL)	free(m_szCamePose);
	if ( m_sz3Dpts!= NULL)		free(m_sz3Dpts);
	if ( m_szReport!= NULL)		free(m_szReport);
}




bool PBA::ba_run(bool bRobust,
	bool bLM,
	int nMaxIter,
	char* szCam,
	char* szFea,
	char* szXYZ,
	char* szCalib,
	char* szReport,
	char* szPose,
	char* sz3D,
	double Tau)
{
	m_szCameraInit = szCam;                                                       
	m_szFeatures   = szFea;                                                       
	m_szCalibration= szCalib;                                                     
	m_szXYZ        = szXYZ;                                                       
	m_bRobustKernel= bRobust;                                                    
	m_bsolverGN    = !bLM;                                                        
	m_nMaxIter     = nMaxIter;                                                   
	m_szCamePose   = szPose;                                                      
	m_sz3Dpts      = sz3D;                                                        
	m_szReport     = szReport;                                                   
	m_Tau          = Tau;                                                         
	//printf("%f\n", m_e1);
	BAType ba = pba;

	ba_initialize(m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ); 
	//printf("%f\n", m_e1);
	#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
		
		//???????????
		// for (int i = 0; i < cams.size();i++){
		// 	double euler_angles[3];
		// 	euler_angles[0] = cams[i].euler_angle[0] * 180 / PI;//pitch kappa
		// 	euler_angles[1] = cams[i].euler_angle[1] * 180 / PI;//roll phi
		// 	euler_angles[2] = cams[i].euler_angle[2] * 180 / PI;//yaw omega
		// 	double R[9];
		// 	double matR[9];
		// 	double p1[3], p2[3], Xc[3];
		// 	ceres::EulerAnglesToRotationMatrix(euler_angles, 3, R);
		// 	matR[0] = cos(cams[i].euler_angle[1]) * cos(cams[i].euler_angle[0]);
		// 	matR[1] = cos(cams[i].euler_angle[1]) * sin(cams[i].euler_angle[0]);
		// 	matR[2] = -sin(cams[i].euler_angle[1]);
		// 	matR[3] = sin(cams[i].euler_angle[2]) * sin(cams[i].euler_angle[1]) * cos(cams[i].euler_angle[0]) - cos(cams[i].euler_angle[2]) * sin(cams[i].euler_angle[0]);
		// 	matR[4] = sin(cams[i].euler_angle[2]) * sin(cams[i].euler_angle[1]) * sin(cams[i].euler_angle[0]) + cos(cams[i].euler_angle[2]) * cos(cams[i].euler_angle[0]);
		// 	matR[5] = sin(cams[i].euler_angle[2]) * cos(cams[i].euler_angle[1]);
		// 	matR[6] = cos(cams[i].euler_angle[2]) * sin(cams[i].euler_angle[1]) * cos(cams[i].euler_angle[0]) + sin(cams[i].euler_angle[2]) * sin(cams[i].euler_angle[0]);
		// 	matR[7] = cos(cams[i].euler_angle[2]) * sin(cams[i].euler_angle[1]) * sin(cams[i].euler_angle[0]) - sin(cams[i].euler_angle[2]) * cos(cams[i].euler_angle[0]);
		// 	matR[8] = cos(cams[i].euler_angle[2]) * cos(cams[i].euler_angle[1]);
		// 	// printf("%f %f %f %f %f %f %f %f %f\n", R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
		// 	// printf("%f %f %f %f %f %f %f %f %f\n", matR[0], matR[1], matR[2], matR[3], matR[4], matR[5], matR[6], matR[7], matR[8]);
		// 	Xc[0] = 1;  
		// 	Xc[1] = 1;
		// 	Xc[2] = 1;
		// 	p1[2] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];
		// 	p1[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];
		// 	p1[0] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];
		// 	p2[0] = matR[0] * Xc[0] + matR[1] * Xc[1] + matR[2] * Xc[2];//Z kappa
		// 	p2[1] = matR[3] * Xc[0] + matR[4] * Xc[1] + matR[5] * Xc[2];//Y phi
		// 	p2[2] = matR[6] * Xc[0] + matR[7] * Xc[1] + matR[8] * Xc[2];//X omega
		// 	printf("%f %f %f\n", p1[0], p1[1], p1[2]);//ZYX
		// 	printf("%f %f %f\n", p2[0], p2[1], p2[2]);//ZYX
		// }

		//????intrinsics
		// for (int i = 0; i < intrs.size();i++){
		// 	printf("%f %f %f %f\n", intrs[i].fx, intrs[i].fy, intrs[i].cx, intrs[i].cy);
		// }

		//????cameras
		// for(int i=0;i<cams.size();i++){
		// 	printf("%f %f %f %f %f %f %d\n", cams[i].euler_angle[0], cams[i].euler_angle[1], cams[i].euler_angle[2],
		// 		   cams[i].camera_center[0], cams[i].camera_center[1], cams[i].camera_center[2],
		// 		   cams[i].camidx);
		// }

		//???????????????????
		// printf("%d\n", tracks.size());
		// for (int i = 0; i < tracks.size();i++){
		// 	// printf("%f %f %f\n", points[i].xyz[0], points[i].xyz[1], points[i].xyz[2]);
		// 	// printf("%d %d %f %f %f\n", points[i].nM, points[i].nA, points[i].aep[0], points[i].aep[1], points[i].aep[2]);
		// 	double w = points[i].aep[2];//parallax angle
		// 	double xj[3], Tik[3], xk[3], p[3];
		// 	xj[0] = sin(points[i].aep[0]) * cos(points[i].aep[1]);//x
		// 	xj[1] = sin(points[i].aep[1]);//y
		// 	xj[2] = cos(points[i].aep[0]) * cos(points[i].aep[1]);//z
		// 	int nM = points[i].nM;
		// 	int nA = points[i].nA;
		// 	Tik[0] = -cams[nM].camera_center[0] + cams[nA].camera_center[0];
		// 	Tik[1] = -cams[nM].camera_center[1] + cams[nA].camera_center[1];
		// 	Tik[2] = -cams[nM].camera_center[2] + cams[nA].camera_center[2];
		// 	double Dik = sqrt(Tik[0] * Tik[0] + Tik[1] * Tik[1] + Tik[2] * Tik[2]);
		// 	double w2 = acos((xj[0] * Tik[0] + xj[1] * Tik[1] + xj[2] * Tik[2]) / Dik);
		// 	//?????????????��???xyz????
		// 	xk[0] = (Dik * sin(w2 + w) * xj[0]) / sin(w);
		// 	xk[1] = (Dik * sin(w2 + w) * xj[1]) / sin(w);
		// 	xk[2] = (Dik * sin(w2 + w) * xj[2]) / sin(w);
		// 	//??????????xyz????
		// 	p[0] = cams[nM].camera_center[0] + xk[0];
		// 	p[1] = cams[nM].camera_center[1] + xk[1];
		// 	p[2] = cams[nM].camera_center[2] + xk[2];
		// 	if(abs(points[i].xyz[0]-p[0])>0.01 || abs(points[i].xyz[1]-p[1])>0.01 || abs(points[i].xyz[2]-p[2])>0.01){
		// 		printf("%f %f %f\n", cams[nM].camera_center[0], cams[nM].camera_center[1], cams[nM].camera_center[2]);
		// 		printf("%f %f %f\n", xk[0], xk[1], xk[2]);
		// 		printf("%d %d %d %f %f %f\n", i+1, nM, nA, points[i].xyz[0], points[i].xyz[1], points[i].xyz[2]);
		// 		printf("%d %d %d %f %f %f\n", i+1, nM, nA, p[0], p[1], p[2]);
		// 	}
		// }
	    
		//????tracks
		// for (int i = 0; i < tracks.size();i++){
		// 	printf("%d ", tracks[i].nview);
		// 	for (int j = 0; j < tracks[i].nview;j++){
		// 		printf("%d %f %f  ", tracks[i].obss[j].view_idx, tracks[i].obss[j].u, tracks[i].obss[j].v);
		// 	}
		// 	printf("\n");
		// }
		ceres::Problem problem;
		ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
		int nobs = 0;
		// double sba_cost_v = 0, pba_cost_v = 0;
		// double sba_er = 0, pba_er = 0;
		for (int i = 0; i < tracks.size(); i++)
		{ // point
			// printf("%d ", tracks[i].nview);
			nobs += tracks[i].nview;
			int nM = points[i].nM;
			int nA = points[i].nA;
			// double camera_center_nM[3], camera_center_nA[3];
			// camera_center_nM[0] = cams[nM].camera_center[0];
			// camera_center_nM[1] = cams[nM].camera_center[1];
			// camera_center_nM[2] = cams[nM].camera_center[2];
			// camera_center_nA[0] = cams[nA].camera_center[0];
			// camera_center_nA[1] = cams[nA].camera_center[1];
			// camera_center_nA[2] = cams[nA].camera_center[2];

			for (int j = 0; j < tracks[i].nview; j++)
			{ // view
				// sba_er = 0;
				// pba_er = 0;
				int view_idx = tracks[i].obss[j].view_idx;
				double u = tracks[i].obss[j].u;
				double v = tracks[i].obss[j].v;
				double fx = 0, fy = 0, cx = 0, cy = 0;
				int cam_idx = cams[view_idx].camidx;
				// printf("%d %d\n", cam_idx, static_cast<int>(intrs.size()));
				// if (cam_idx <= static_cast<int>(intrs.size())){
				fx = intrs[cam_idx-1].fx;
				fy = intrs[cam_idx-1].fy;
				cx = intrs[cam_idx-1].cx;
				cy = intrs[cam_idx-1].cy;
				int nP = view_idx;
				//??????????
				// double R[9],KR[9],euler_angles[3], Xc[3],point3D[3],camera_center_nP[3],p[3],residuals[2];
				// double c0, c1, c2, s0, s1, s2;
				// euler_angles[0] = cams[view_idx].euler_angle[0];
				// euler_angles[1] = cams[view_idx].euler_angle[1];
				// euler_angles[2] = cams[view_idx].euler_angle[2];
				// camera_center_nP[0] = cams[view_idx].camera_center[0];
				// camera_center_nP[1] = cams[view_idx].camera_center[1];
				// camera_center_nP[2] = cams[view_idx].camera_center[2];
				// point3D[0] = points[i].xyz[0];
				// point3D[1] = points[i].xyz[1];
				// point3D[2] = points[i].xyz[2];
				// c0 = cos(euler_angles[0]);
				// c1 = cos(euler_angles[1]);
				// c2 = cos(euler_angles[2]);
				// s0 = sin(euler_angles[0]);
				// s1 = sin(euler_angles[1]);
				// s2 = sin(euler_angles[2]);
				// R[0] = c1 * c0;
				// R[1] = c1 * s0;
				// R[2] = -s1;
				// R[3] = s2 * s1 * c0 - c2 * s0;
				// R[4] = s2 * s1 * s0 + c2 * c0;
				// R[5] = s2 * c1;
				// R[6] = c2 * s1 * c0 + s2 * s0;
				// R[7] = c2 * s1 * s0 - s2 * c0;
				// R[8] = c2 * c1;
				// Xc[0] = point3D[0] - camera_center_nP[0];
				// Xc[1] = point3D[1] - camera_center_nP[1];
				// Xc[2] = point3D[2] - camera_center_nP[2];
				// p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
				// p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
				// p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z
				// // Normalize
				// double xp = p[0] / p[2];
				// double yp = p[1] / p[2];
				// // Project
				// double predicted_u = fx * xp + cx;
				// double predicted_v = fy * yp + cy;
				// // Residual
				// residuals[0] = predicted_u - u;
				// residuals[1] = predicted_v - v;
				// double err = residuals[0] * residuals[0] + residuals[1] * residuals[1];
				// sba_cost_v += err;
				// sba_er += err;
				// point3D[0] = points[i].aep[0];
				// point3D[1] = points[i].aep[1];
				// point3D[2] = points[i].aep[2];
				// if(nP==nM){
				// 	double ptXj[3];
				// 	ptXj[0] = sin(point3D[0]) * cos(point3D[1]);
				// 	ptXj[1] = sin(point3D[1]);
				// 	ptXj[2] = cos(point3D[0]) * cos(point3D[1]);
				// 	p[0] = R[0] * ptXj[0] + R[1] * ptXj[1] + R[2] * ptXj[2];
				// 	p[1] = R[3] * ptXj[0] + R[4] * ptXj[1] + R[5] * ptXj[2];
				// 	p[2] = R[6] * ptXj[0] + R[7] * ptXj[1] + R[8] * ptXj[2];
				// 	// Normalize
				// 	double xp = p[0] / p[2];
				// 	double yp = p[1] / p[2];
				// 	// Project
				// 	double predicted_u = fx * xp + cx;
				// 	double predicted_v = fy * yp + cy;
				// 	// Residual
				// 	residuals[0] = predicted_u - u;
				// 	residuals[1] = predicted_v - v;
				// }else if(nP==nA){
				// 	double pti2k[3], ptXUnit[3], ptXk[3];
				// 	//??��????��??????????
				// 	pti2k[0] = camera_center_nP[0] - camera_center_nM[0];	
				// 	pti2k[1] = camera_center_nP[1] - camera_center_nM[1];	
				// 	pti2k[2] = camera_center_nP[2] - camera_center_nM[2];	
				// 	//??��??????????��????
				// 	ptXUnit[0] = sin(point3D[0]) * cos(point3D[1]);
				// 	ptXUnit[1] = sin(point3D[1]);
				// 	ptXUnit[2] = cos(point3D[0]) * cos(point3D[1]);
				// 	//compute angle w2
				// 	double dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
				// 	double dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
				// 	double dW2;
				// 	if (dDot/dDisi2k > 1)
				// 		dW2 = 0;
				// 	if (dDot/dDisi2k < -1)
				// 		dW2 = PI;
				// 	else
				// 		dW2  = acos(dDot/dDisi2k);
				// 	//compute Xk vector according sin theory
				// 	ptXk[0] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[0] - sin(point3D[2]) * pti2k[0];
				// 	ptXk[1] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[1] - sin(point3D[2]) * pti2k[1];
				// 	ptXk[2] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[2] - sin(point3D[2]) * pti2k[2];
				// 	p[0] = R[0] * ptXk[0] + R[1] * ptXk[1] + R[2] * ptXk[2];
				// 	p[1] = R[3] * ptXk[0] + R[4] * ptXk[1] + R[5] * ptXk[2];
				// 	p[2] = R[6] * ptXk[0] + R[7] * ptXk[1] + R[8] * ptXk[2];
				// 	// Normalize
				// 	double xp = p[0] / p[2];
				// 	double yp = p[1] / p[2];
				// 	// Project
				// 	double predicted_u = fx * xp + cx;
				// 	double predicted_v = fy * yp + cy;
				// 	// Residual
				// 	residuals[0] = predicted_u - u;
				// 	residuals[1] = predicted_v - v;
				// }else{
				// 	double pti2k[3], pti2l[3], ptXUnit[3], ptXk[3];
				// 	//??��????��??????????
				// 	pti2k[0] = camera_center_nA[0] - camera_center_nM[0];		
				// 	pti2k[1] = camera_center_nA[1] - camera_center_nM[1];		
				// 	pti2k[2] = camera_center_nA[2] - camera_center_nM[2];
				// 	//??��??????��??????????
				// 	pti2l[0] = camera_center_nP[0] - camera_center_nM[0];	
				// 	pti2l[1] = camera_center_nP[1] - camera_center_nM[1];		
				// 	pti2l[2] = camera_center_nP[2] - camera_center_nM[2];
				// 	//XUnit ??��??????????��????
				// 	ptXUnit[0] = sin(point3D[0]) * cos(point3D[1]);
				// 	ptXUnit[1] = sin(point3D[1]);
				// 	ptXUnit[2] = cos(point3D[0]) * cos(point3D[1]);
				// 	//compute angle w2
				// 	double dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
				// 	double dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
				// 	double dW2;
				// 	//dW2  = acos( dDot/dDisi2k );
				// 	if (dDot/dDisi2k > 1)
				// 		dW2 = 0;
				// 	if ( dDot/dDisi2k < -1)
				// 		dW2 = PI;
				// 	else
				// 		dW2  = acos(dDot/dDisi2k);
				// 	//compute Xl vector according sin theory
				// 	ptXk[0] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[0] - sin(point3D[2]) * pti2l[0];
				// 	ptXk[1] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[1] - sin(point3D[2]) * pti2l[1];
				// 	ptXk[2] = dDisi2k * sin(dW2+point3D[2]) * ptXUnit[2] - sin(point3D[2]) * pti2l[2];	
				// 	p[0] = R[0] * ptXk[0] + R[1] * ptXk[1] + R[2] * ptXk[2];
				// 	p[1] = R[3] * ptXk[0] + R[4] * ptXk[1] + R[5] * ptXk[2];
				// 	p[2] = R[6] * ptXk[0] + R[7] * ptXk[1] + R[8] * ptXk[2];
				// 	// Normalize
				// 	double xp = p[0] / p[2];
				// 	double yp = p[1] / p[2];
				// 	// Project
				// 	double predicted_u = fx * xp + cx;
				// 	double predicted_v = fy * yp + cy;
				// 	// Residual
				// 	residuals[0] = predicted_u - u;
				// 	residuals[1] = predicted_v - v;
				// }
				// err = residuals[0] * residuals[0] + residuals[1] * residuals[1];
				// pba_cost_v += err;
				// pba_er += err;

				// ceres::CostFunction *cost_function = SBA_ReprojectionError::Create(u, v, fx, fy, cx, cy);
				// problem.AddResidualBlock(
				// 	cost_function,
				// 	nullptr, // ?????? Huber loss ?????? HuberLoss(1.0) nullptr,new ceres::HuberLoss(1.0) loss_function
				// 	&cams[view_idx].euler_angle[0],
				// 	&cams[view_idx].camera_center[0],
				// 	&points[i].xyz[0]);

				if(nP==nM){
					ceres::CostFunction *cost_function = PBA_nM_ReprojectionError::Create(u, v, fx, fy, cx, cy);
					problem.AddResidualBlock(
						cost_function,
						nullptr, //loss_function,
						&cams[nP].euler_angle[0],
						&points[i].aep[0]);
				}else if(nP==nA){
					ceres::CostFunction *cost_function = PBA_nA_ReprojectionError::Create(u, v, fx, fy, cx, cy);
					problem.AddResidualBlock(
						cost_function,
						nullptr, //loss_function, 
						&cams[nP].euler_angle[0],
						&cams[nP].camera_center[0],
						&cams[nM].camera_center[0],
						&points[i].aep[0]);
				}else{
					ceres::CostFunction *cost_function = PBA_nP_ReprojectionError::Create(u, v, fx, fy, cx, cy);
					problem.AddResidualBlock(
						cost_function,
						nullptr, //loss_function, 
						&cams[nP].euler_angle[0],
						&cams[nP].camera_center[0],
						&cams[nM].camera_center[0],
						&cams[nA].camera_center[0],
						&points[i].aep[0]);
				}

			}
			// if(abs(pba_er-sba_er)>0.01){
			// 	printf("%d %f %f\n", i, sba_er, pba_er);
			// } 
			// printf("%s%f\n", "sba_cost_v = ", sba_cost_v);
			// printf("%s%f\n", "pba_cost_v = ", pba_cost_v);
		}
		
		printf("Ceres Options\n");
		ceres::Solver::Options options;
		options.logging_type = ceres::SILENT;
		options.minimizer_progress_to_stdout = true;
		options.num_threads = 8;
		// options.trust_region_strategy_type = ceres::DOGLEG;         
		options.linear_solver_type = ceres::SPARSE_SCHUR;            
		options.sparse_linear_algebra_library_type = ceres::CUDA_SPARSE;
		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem, &summary);
		std::cout << summary.FullReport() << "\n";

		// problem.Evaluate(ceres::Problem::EvaluateOptions(), &initial_cost, nullptr, nullptr, nullptr);
		// printf("%s%f\n", "sba_cost_v = ", sba_cost_v/nobs);
		// printf("%s%f\n", "pba_cost_v = ", pba_cost_v/nobs);
		double initial_cost = summary.initial_cost;
		double final_cost = summary.final_cost;
		printf("%f %s %f\n", initial_cost*2/nobs, "??", final_cost*2/nobs);

#elif defined(_WIN32)
		if (m_bsolverGN) {
			ba_motstr_gn();
		}	
		else {
			ba_motstr_levmar();
		}	
	#else
		#error "Unsupported platform!"
	#endif

	return true;
}

bool PBA::ba_run( int argc, char** argv )
{
	bool bTrue = ba_parseArgs( argc, argv );
	if (!bTrue)
	{
		fprintf( stderr, "ParallaxBA: Input wrong commands, please check them!\n");
		return false;
	}
	
	ba_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );

	#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))

	#elif defined(_WIN32)
		if (m_bsolverGN) {
			ba_motstr_gn();
		}	
		else {
			ba_motstr_levmar();
		}	
	#else
		#error "Unsupported platform!"
	#endif

	return true;
}


bool PBA::ba_initialize( char* szCamera, char* szFeature,  char* szCalib, char* szXYZ )
{
	printf("ParallaxBA: Parallax Angle Bundle Adjustment Version 1.0\n");
	FILE* fp = nullptr;

	//must input initial initial camera pose file and projection image points file
	fopen_s(&fp, szCamera, "r" );
	if ( fp == NULL )
	{
		fprintf( stderr, "ParallaxBA: Missing initial camera poses file! \n");
		exit(1);
	}
	else
		fclose(fp);

	fopen_s(&fp, szFeature, "r" );
	if ( fp == NULL )
	{	
		fprintf( stderr, "ParallaxBA: Missing feature projection points file! \n");
		exit(1); 
	}
	else
		fclose(fp);

	if (tpe == bal) {

	}
	else if (tpe == colmap) {
		if (szCalib != NULL)
		{
			FILE* fpc = nullptr;
			fopen_s(&fpc, szCalib, "r");//????????
			nc_ = findNcameras(fpc) / 3;//????????
			fclose(fpc);
			fpc = nullptr;
			// printf("%d\n", nc_);
			m_bFocal = false;
			m_K = (double*)malloc(9 * nc_ * sizeof(double));
			ba_readCameraPoseration(szCalib, m_K);
			//ba_readCameraPose_(szCalib);
			// for (int i = 0; i < nc_;i++){
			// 	printf("%f %f %f %f\n", intrs[i].fx, intrs[i].fy, intrs[i].cx, intrs[i].cy);
			// }
			// int zsw = 0;
			printf("????????\n");
		}
	}

	if ( szXYZ != NULL )
		m_bProvideXYZ = true;	
	zu = 0;
	savePara = 0; 
	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	pba_readAndInitialize( szCamera, szFeature,szCalib, &m_ncams, &m_n3Dpts, &m_n2Dprojs,&m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
		&m_imgpts, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort);

	//string originalPath(szCamera); 
	//size_t pos = originalPath.find_last_of("/\\");
	//string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	//string newPath1 = parentPath + "/" + "InitialParallax.txt";
	//string newPath2 = parentPath + "/" + "Triangulatedxyz.txt";
	//string newPath3 = parentPath + "/" + "XYZ1.txt";

	//if(savePara == 1)
	//	pba_saveInitialParallax(newPath1.c_str(), m_motstruct);
	//if (zu == 1)
	//	ba_saveTriangulatedxyz(newPath2.c_str(), m_motstruct);
	//else if (zu == 2)
	//	pba_saveInitialXYZ(newPath3.c_str(), pba_angle2xyz(m_motstruct));//????xyz
	//else
	//	;//Other codes

	//printf("Number of cameras: %d\n", m_ncams);
	//printf("Number of points: %d\n", m_n3Dpts);
	//printf("Number of projections: %d\n", m_n2Dprojs);
	return true;
}



void PBA::pba_readProjectionAndInitilizeFeature(FILE *fp,
	double *params, double *projs, char *vmask, int ncams, 
	int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort )
{
	int n;
	int nframes;
	int ptno = 0, cur;

	int nproj2D = 0;
	
	int count = 0;

	int frameno;
	int feastart = 0;

	int nP, nP2;
	
	double* ptr1 = projs;

	int i, j; 
	int  sum, cnp = 6;

	int nFlag;
	
	int *ptr2;
	bool bAdjust;	

	m_smask = (char*)malloc(m_ncams*m_ncams*sizeof(char));
	memset( m_smask, 0, m_ncams*m_ncams*sizeof(char) );

	int* archorEx = new int[m_n3Dpts*2];

	bool bM, bN;
	int max_nframes = -1;
	//read all projection point, initialize three feature angle at the same time
	while(!feof(fp))//????track
	{
		nFlag = 0;
		n = readNInts( fp, &nframes, 1 );  
		if( n!=1 )
			break;

		#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
			Track track;
			track.nview = nframes;
		#elif defined(_WIN32)

		#else
			#error "Unsupported platform!"
		#endif

		archor[ptno*3] = nframes;
		cur = 0;
		bM = bN = false;
		for( i=0, sum = 0; i<nframes; ++i )
		{
			n = readNInts( fp, &frameno, 1 );
			nphoto[nproj2D] = frameno;
			nfeature[nproj2D] = ptno;
			nproj2D++;

			if(frameno>=ncams)
			{
				fprintf(stderr, "ParallaxBA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles( fp, ptr1, 2 ); 

			#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
				Observation obs;
				obs.view_idx = frameno;
				obs.u = ptr1[0];
				obs.v = ptr1[1];
				track.obss.push_back(obs);
			#elif defined(_WIN32)

			#else
				#error "Unsupported platform!"
			#endif


			ptr1+=2;
			if(n!=3)
			{
				fprintf(stderr, "ParallaxBA:reading image projections wrong!\n");
				return;
			}

			//??????
			if ( bM && bN )
			{
				ptr2 = archorSort+ptno*2;
				//bAdjust = pba_initializeOtheArchors_Mindw( //_Mindw
				bAdjust = pba_initializeOtheArchors( //Maxdw
					projs+feastart*2,
					nphoto+feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams*cnp + ptno * 3, 
					ptr2,
					sum,
					i,
					ptno );
				if ( bAdjust )
				{
					archor[ptno*3+1] = *(nphoto+feastart+ptr2[0]);
					archor[ptno*3+2] = *(nphoto+feastart+ptr2[1]);

					archorEx[ptno*2] = ptr2[0];
					archorEx[ptno*2+1] = ptr2[1];
				}
				sum++;
			}

			//?????
			if ( bM && !bN )
			{	
				bool bLast = (i == nframes-1);
				bool bT = pba_initializeAssoArchor( 
					projs+feastart*2,	
					nphoto+feastart,	
					m_motstruct,
					m_K,
					m_motstruct+m_ncams*cnp+ptno*3,
					0,
					1,
					ptno,
					bLast );

				if (bT)
				{
					archorSort[ptno*2+1] = i;
					archor[ptno*3+2] = nphoto[count];
					sum++;

					archorEx[ptno*2+1] = i;

					bN = true;
				}
			}


			//?????
			if ( !bM )
			{
				bool bLast = (i == nframes-2);
				bool bT = pba_initializeMainArchor( 
					projs+feastart*2,	
					m_motstruct,		
					m_K,				
					m_motstruct+m_ncams*cnp+ptno*3,
					nphoto[count],		
					ptno,				
					m_KR );				

				archorSort[ptno*2] = i;
				archor[ptno*3+1] = nphoto[count];
				sum++;

				archorEx[ptno*2] = i;
				bM = true;
			}	
			count++;	
		}
		#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
			tracks.push_back(track);
			// points[ptno].nM=archor[ptno * 3 + 1];
			// points[ptno].nA=archor[ptno * 3 + 2];
		#elif defined(_WIN32)

		#else
			#error "Unsupported platform!"
		#endif
		


		//set masks for U and S matrix             
		for( i = 0; i < nframes; i++ )
		{
			nP = nphoto[feastart+i];                         
			int nM_ = archor[ptno * 3 + 1];
			int nA_ = archor[ptno * 3 + 2];
			int tmp3 = archor[ptno * 3 + 3];
			
			if (nM_<nP)                         
				umask[nM_*(ncams)+nP] = 1;
			else                                             
				umask[nP*(ncams)+nM_] = 1;

			
			if (nA_<nP)                         
				umask[nA_*(ncams)+nP] = 1;
			else                                             
				umask[nP*(ncams)+nA_] = 1;

			umask[nP*ncams+nP] = 1;

			for ( j = i; j < nframes; j++  )
			{
				nP2 = nphoto[feastart+j];                

				if ( nP == nP2 )                              
					m_smask[nP*m_ncams+nP2] = 1;
				else if ( nP < nP2 )
					m_smask[nP*m_ncams+nP2] = 1;
				else
				{
					m_smask[nP2*m_ncams + nP] = 1;
				}
					
			}
		}					
		feastart += nframes;
		ptno++;
	}


	// count number of non-zero element in S matrix
	m_nS = 0;
	for ( i = 0; i < m_ncams; i++ ) 
	{
		for (j = 0; j < m_ncams; j++)
		{
			if (m_smask[i*m_ncams + j] == 1)
			{
				m_nS++;
			}
		}
	}

}

void PBA::pba_readAndInitialize(char *camsfname, char *ptsfname,char *calibfname, int *ncams,
	int *n3Dpts, int *n2Dprojs,
	double **motstruct, double **imgpts,
	int **archor, char **vmask,
	char **umask, int **nphoto,
	int** nfeature, int** archorSort)
{
	FILE *fpc = nullptr, *fpp = nullptr, *fpXYZ = nullptr;
	int i, tmp1, tmp2;
	double ptMain[3], ptA[3];
	double dW1, dW2;	

	//calculate number of cameras, 3D points and projection points
	fopen_s(&fpc, camsfname, "r" );//????��??????
	*ncams	=	findNcameras( fpc );//????????????
	m_ncams =	*ncams;
	m_V = (int*)malloc(sizeof(int) * m_ncams);
	
	if (tpe == bal) {
		if (calibfname != NULL)
		{
			m_bFocal = false;
			m_K = (double*)malloc(*ncams * 9 * sizeof(double));
			ba_readCameraPoseration(calibfname, m_K);
			//for (int i = 0; i < m_ncams; i++) {
			//	printf("%f\n", m_K[i * 9]);
			//}
		}
	}
	else if (tpe == colmap) {

	}


	fopen_s(&fpp, ptsfname, "r" );//??????????????
	readNpointsAndNprojections( fpp, n3Dpts, 3, n2Dprojs, 2 );//???3D???????????????????

	*motstruct = (double*)malloc((*ncams * 6 + *n3Dpts * 3) * sizeof(double));//????��???��???��??? ?? ???????????3D????
	if(	*motstruct==NULL )
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'motstruct' failed \n");
		exit(1);
	}

	*imgpts = (double*)malloc(*n2Dprojs * 2 * sizeof(double));//????��?????????2D????
	if(	*imgpts==NULL )
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'imgpts' failed\n");
		exit(1);
	}

	rewind(fpc);
	rewind(fpp);

	//allocate indicator of U
	*umask = (char*)malloc(*ncams * *ncams );
	memset(*umask, 0, *ncams * *ncams * sizeof(char));

	//allocate main and associate anchors
	*archor = (int*)malloc(*n3Dpts*3*sizeof(int));//
	memset(*archor, -1, *n3Dpts * 3 * sizeof(int)); //????????��???��?��?????-1

	*nphoto		= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*nfeature	= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*archorSort = (int*)malloc(*n3Dpts*3*sizeof(int));

	ba_readCameraPose(fpc, *motstruct, m_V);//m_V?????????
	printf("????????\n");
	// for (int i = 0;i < m_ncams; i++){
	// 	printf("%f %f %f %f %f %f %d\n", cams[i].euler_angle[0], cams[i].euler_angle[1], cams[i].euler_angle[2], cams[i].camera_center[0], cams[i].camera_center[1], cams[i].camera_center[2], cams[i].camidx);
	// }

	fclose(fpc);
	fpc = NULL;

	//for (int i = 0; i < m_ncams; i++) {
	//	printf("%d\n", m_V[i]);
	//}
	
	//if (m_K != NULL) {
	//	printf("m_K address: %p\n", (void*)m_K);
	//	free(m_K);
	//	m_K = NULL;
	//}
	//if (m_V != NULL) {
	//	printf("m_V address: %p\n", (void*)m_V);
	//	free(m_V);
	//	m_V = NULL;
	//}
	//if (m_motstruct != NULL) {
	//	printf("m_motstruct address: %p\n", (void*)m_motstruct);
	//	free(m_motstruct);
	//	m_motstruct = NULL;
	//}
	//if (m_imgpts != NULL) {
	//	printf("m_imgpts address: %p\n", (void*)m_imgpts);
	//	free(m_imgpts);
	//	m_imgpts = NULL;
	//}
	//if (m_umask != NULL) {
	//	printf("m_umask address: %p\n", (void*)m_umask);
	//	free(m_umask);
	//	m_umask = NULL;
	//}
	//if (m_archor != NULL) {
	//	printf("m_archor address: %p\n", (void*)m_archor);
	//	free(m_archor);
	//	m_archor = NULL;
	//}
	//if (m_photo != NULL) {
	//	printf("m_photo address: %p\n", (void*)m_photo);
	//	free(m_photo);
	//	m_photo = NULL;
	//}
	//if (m_feature != NULL) {
	//	printf("m_feature address: %p\n", (void*)m_feature);
	//	free(m_feature);
	//	m_feature = NULL;
	//}
	//if (m_archorSort != NULL) {
	//	printf("m_archorSort address: %p\n", (void*)m_archorSort);
	//	free(m_archorSort);
	//	m_archorSort = NULL;
	//}
	
	//Update KR
	m_KR  = (double*)malloc(m_ncams*9*sizeof(double));//?????��??? ?? ??????????
	m_KdA = (double*)malloc(m_ncams*9*sizeof(double));//???????M_x (??) M_y (??) M_z (??)??????
	m_KdB = (double*)malloc(m_ncams*9*sizeof(double));//???????M_x (??) M_y (??) M_z (??)??????
	m_KdG = (double*)malloc(m_ncams*9*sizeof(double));//???????M_x (??) M_y (??) M_z (??)??????
	
	if (zu == 1)
	{
		//????P????
		m_P = (double*)malloc(m_ncams * 12 * sizeof(double));
		ba_constructP(m_P, m_K, *motstruct);
	}
	else
		ba_updateKR(m_KR, m_KdA, m_KdB, m_KdG, m_K, *motstruct);

	//if (m_KR == NULL) {
	//	// ??????��?????
	//	fprintf(stderr, "Memory allocation failed!\n");
	//	exit(EXIT_FAILURE);
	//}
	//if (m_KR != NULL) {
	//	printf("m_KR address: %p\n", (void*)m_KR);
	//	free(m_KR);
	//	m_KR = NULL;
	//}
	//test
	//for (int i = 0; i < m_ncams; i++)
	//{
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9], m_KR[i * 9 + 1], m_KR[i * 9 + 2], m_P[i * 12], m_P[i * 12 + 1], m_P[i * 12 + 2]);
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9 + 3], m_KR[i * 9 + 4], m_KR[i * 9 + 5], m_P[i * 12 + 4], m_P[i * 12 + 5], m_P[i * 12 + 6]);
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9 + 6], m_KR[i * 9 + 7], m_KR[i * 9 + 8], m_P[i * 12 + 8], m_P[i * 12 + 9], m_P[i * 12 + 10]);
	//}


	//if XYZ are provided, we can use them as feature initialization.
	if (m_bProvideXYZ)
	{
		fopen_s(&fpXYZ, m_szXYZ, "r");
		m_XYZ = (double*)malloc(m_n3Dpts*3*sizeof(double));

		for( i = 0; i < m_n3Dpts; i++){
			fscanf_s(fpXYZ, "%lf  %lf  %lf", m_XYZ + i * 3, m_XYZ + i * 3 + 1, m_XYZ + i * 3 + 2);

			#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
				Point3D p3d;
				p3d.xyz[0] = m_XYZ[i * 3];
				p3d.xyz[1] = m_XYZ[i * 3 + 1];
				p3d.xyz[2] = m_XYZ[i * 3 + 2];
				points.push_back(p3d);
			#elif defined(_WIN32)
			
			#else
				#error "Unsupported platform!"
			#endif

		}
		printf("??????????\n");
		// for (int i = 0; i < m_n3Dpts;i++){
		// 	printf("%f %f %f\n", points[i].xyz[0], points[i].xyz[1], points[i].xyz[2]);
		// }

		fclose(fpXYZ);

		string originalPath(m_szXYZ);
		size_t pos = originalPath.find_last_of("/\\");
		string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
		string init3D = parentPath + "/" + "XYZ.ply";
		
		//save features xyz
		if (init3D.c_str() != NULL)
		{
			FILE *fp = nullptr;
			fopen_s(&fp, init3D.c_str(), "w");
			fprintf(fp, "%s\n", "ply");
			fprintf(fp, "%s\n", "format ascii 1.0");
			fprintf(fp, "%s %d\n", "element vertex", m_n3Dpts);
			fprintf(fp, "%s\n", "property float x");
			fprintf(fp, "%s\n", "property float y");
			fprintf(fp, "%s\n", "property float z");
			fprintf(fp, "%s\n", "end_header");
			for (i = 0; i < m_n3Dpts; i++)
				fprintf(fp, "%lf %lf %lf\n", m_XYZ[i * 3], m_XYZ[i * 3 + 1], m_XYZ[i * 3 + 2]);
			fclose(fp);
		}
	}

	if(zu==1)
		ba_readProjectionAndTriangulateFeature(fpp, *imgpts, *ncams);
	else
	{
		pba_readProjectionAndInitilizeFeature(fpp,
			*motstruct + *ncams * 6,
			*imgpts,
			*vmask,
			*ncams,
			*archor,
			*umask,
			*nphoto,
			*nfeature,
			*archorSort);
	}
	printf("?????????\n");
	fclose(fpp);

	
	int nCount = 0;
	double pti2k[3];
	int cur = 0;
	int bs = 1;//???????????????
	if (tpe == colmap) {
		bs = 1;
	}
	else if (tpe == bal) {
		bs = -1;
	}

	if (m_bProvideXYZ)
	{
		for (i = 0; i < m_n3Dpts; i++)
		{
			int nM = m_archor[i*3+1];
			int nN = m_archor[i*3+2];
			//printf("%d %d\n", nM, nN);

			//??��??��??
			ptMain[0] = *(*motstruct + nM*6 + 3);
			ptMain[1] = *(*motstruct + nM*6 + 4);
			ptMain[2] = *(*motstruct + nM*6 + 5);
			//??��??��??
			ptA[0] = *(*motstruct + nN*6 + 3);
			ptA[1] = *(*motstruct + nN*6 + 4);
			ptA[2] = *(*motstruct + nN*6 + 5);
			//??��????��???????
			pti2k[0] = ptA[0] - ptMain[0];
			pti2k[1] = ptA[1] - ptMain[1];
			pti2k[2] = ptA[2] - ptMain[2];

			double dispti2k;
			dispti2k = sqrt(pti2k[0] * pti2k[0] + pti2k[1] * pti2k[1] + pti2k[2] * pti2k[2]);//???????
			//??��??3D???????1
			ptMain[0] = m_XYZ[i * 3 + 0] - ptMain[0];
			ptMain[1] = m_XYZ[i * 3 + 1] - ptMain[1];
			ptMain[2] = m_XYZ[i * 3 + 2] - ptMain[2]; 
			//??��??3D???????2
			ptA[0] = m_XYZ[i * 3 + 0] - ptA[0];
			ptA[1] = m_XYZ[i * 3 + 1] - ptA[1];
			ptA[2] = m_XYZ[i * 3 + 2] - ptA[2];
			//????1????????
			dW1 = ptMain[0] * ptMain[0] + ptMain[1] * ptMain[1] + ptMain[2] * ptMain[2];
			//????2????????
			dW2 = ptA[0] * ptA[0] + ptA[1] * ptA[1] + ptA[2] * ptA[2];

			double disDot2;
			disDot2 = ptMain[0] * pti2k[0] + ptMain[1] * pti2k[1] + ptMain[2] * pti2k[2]; 
			double dww = disDot2 / (dispti2k * sqrt(dW1));//??-??��??????????-3D????????��?


			double* pKR = m_KR + nM * 9;
			double n[2], n2[2], ptXj[3];

			ptXj[0] = ptMain[0];	ptXj[1] = ptMain[1];	ptXj[2] = ptMain[2];
			//3D????????��???uv????
			n[0] = bs*(pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			n[1] = bs*(pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			pKR = m_KR + nN*9;
			//3D????????��???uv????
			ptXj[0] = ptA[0];	ptXj[1] = ptA[1];	ptXj[2] = ptA[2];
			n2[0] = bs*(pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			n2[1] = bs*(pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			//printf("%d %d %d\n", m_archorSort[i * 2], m_archorSort[i * 2 + 1], m_archor[i * 3]);
			//???m_archor??m_archorSort?????m_archor?????????��?????????m_archorSort?????????��??????track??????
			int id1 = cur + m_archorSort[i * 2];//3D??i????��???uv????
			int id2 = cur + m_archorSort[i * 2 + 1];//3D??i???��???uv????
			//printf("%f %f\n", m_imgpts[id1 * 2], m_imgpts[id1 * 2 + 1]);
			//printf("%f %f\n", m_imgpts[id2 * 2], m_imgpts[id2 * 2 + 1]);
			double err1 = (m_imgpts[id1 * 2] - n[0]) * (m_imgpts[id1 * 2] - n[0]) + (m_imgpts[id1 * 2 + 1] - n[1]) * (m_imgpts[id1 * 2 + 1] - n[1]);;
			double err2 = (m_imgpts[id2 * 2] - n2[0]) * (m_imgpts[id2 * 2] - n2[0]) + (m_imgpts[id2 * 2 + 1] - n2[1]) * (m_imgpts[id2 * 2 + 1] - n2[1]);

			//printf("%f %f\n", err1, err2);
			cur += m_archor[i * 3];
			/*
			* 1.??????????????????????
			* 2.???????????????��?????????>??��?? 
			*/
			if ((sqrt(dW1) / sqrt(dW2) > 30) || (err1 > err2) && (m_archor[3 * i] == 2))
			{
				nCount++;
				//????��????
				m_archor[i*3+1] = nN;
				m_archor[i*3+2] = nM;

				tmp1 = m_archorSort[i*2] ;
				tmp2 = m_archorSort[i*2+1] ;

				m_archorSort[i*2] = tmp2;
				m_archorSort[i*2+1] = tmp1;

				//??????????
				double dDAngle = atan2(ptA[0], ptA[2]);
				double dHAngle = atan2(ptA[1], sqrt(ptA[0] * ptA[0] + ptA[2] * ptA[2]));

				(*motstruct)[m_ncams * 6 + i * 3] = dDAngle;
				(*motstruct)[m_ncams * 6 + i * 3 + 1] = dHAngle;

				//double dwwDot = ptMain[0]*ptA[0] + ptMain[1]*ptA[1] + ptMain[2]*ptA[2];				
			}	
		}
	}
#if defined(_WIN64) || (defined(_MSC_VER) && defined(_M_X64))
	for (int i = 0; i < points.size();i++){
		points[i].nM = m_archor[i * 3 + 1];
		points[i].nA = m_archor[i * 3 + 2];
		points[i].aep[0] = (*motstruct)[m_ncams * 6 + i * 3];
		points[i].aep[1] = (*motstruct)[m_ncams * 6 + i * 3 + 1];
		points[i].aep[2] = (*motstruct)[m_ncams * 6 + i * 3 + 2];
	}
#elif defined(_WIN32)
	
#else
	#error "Unsupported platform!"
#endif


	if (m_bProvideXYZ)
		free(m_XYZ);

	fpXYZ = NULL;
}


/*void PBA::pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
{
	int nM, nN;	  
	register int i, j, ii, jj, k;
	int cnp, pnp, mnp;
	double *pa, *pb, *ppt;	
	double *ppUpa, *pea, *peb, *pe, *pV, *pW;
	int pos, pos2, nP, nF, numfea, cur;
	static double sum, x2, delt2;
	register double tmp = 0;

	double pAM[6], pAA[6], pPA[12], pPB[6];

	cnp = 6; pnp = 3; mnp = 2;
	pa=p; pb=p+m*cnp;
	delt2 = m_delt*m_delt;

	pos2 = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		cur = 0;
		nF = i;
		for ( j = 0; j < numfea; j++ )
		{
			nP = nphoto[pos2];
			ppt=pb + nF*pnp;      
			pe= e + pos2*2;		//error

			nM = archor[nF*3+1];
			nN = archor[nF*3+2];	

			if( m_nRobustType == 1)
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];
				x2 = 1.0/( x2/delt2+1 );
			}
			else
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];

				if (sqrt(x2) < m_delt) //inlier
					x2 = 1;
				else // outliers 
					x2 = m_delt/sqrt(x2);
			}			

			pba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );

			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);		
			ppUpa = U + pos*(6*6);

			for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//diag
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii]*pPA[k*6+jj];
				ppUpa[ii*6+jj] += sum*x2;
			}

			//ea
			pea = ea + nP*6;
			for ( ii = 0; ii < 6; ii++ )
			{
				for( jj = 0, sum = 0; jj < 2; jj++ )
					sum += pPA[jj*6+ii] * pe[jj];
				pea[ii] += sum*x2;
			}
			//V
			pV = V + nF*3*3;
			for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii]*pPB[k*3+jj];
				pV[ii*3+jj] += sum*x2;

			}
			//eb
			peb = eb + nF*3;
			for ( ii = 0; ii < 3; ii++ )
			{
				for( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii] * pe[k];
				peb[ii] += sum*x2;
			}
			//W
			pW = W + pos2 * 3*6;
			for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii] * pPB[k*3+jj];
				pW[ii*3+jj] += sum*x2;
			}

			if ( nP == nM )
			{
				cur++;
				pos2++;
				continue;
			}
			else
				if( nP == nN )
				{
					//U		  
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						pea[ii+3] += sum*x2;
					}
					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}
				else
				{
					//U
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					pos = sba_crsm_elmidx(Uidxij, nN, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAA[k*3+ii] * pAA[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if ( nM < nN )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAM[k*3+ii] * pAA[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAA[k*3+ii] * pAM[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}

					if ( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					if ( nN < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAA[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}


					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						pea[ii+3] += sum*x2;
					}

					pea = ea + nN*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
						pea[ii+3] += sum*x2;
					}  

					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}

					pos = pos2 + (m_archorSort[i*2+1]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}

				pos2++;
				cur++;
		}
	}
}
*/


bool PBA::pba_initializeMainArchor( 
	double* imgpts,	
	double* camera,	
	double* K,		
	double* feature,
	int nP,			
	int FID,		
	double* KR )	
{
	//solve  KRX = x
	Vector3d x;//???????????/??????��???XYZ????
	if (m_bProvideXYZ)//??????????3D???????
	{
		x(0) = m_XYZ[FID*3] - *(camera + nP*6+3);
		x(1) = m_XYZ[FID*3+1] - *(camera + nP*6+4);
		x(2) = m_XYZ[FID*3+2] - *(camera + nP*6+5);
	}
	
	else
	{
		double *ptr = m_KR + nP*9;	
		Matrix3d  A;
		A << ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7], ptr[8];//??????????KR????
		//?? KR ????? 3x3 ????A??
		double matx[3];				
		matx[0] = imgpts[0];
		matx[1] = imgpts[1];
		matx[2] = 1;
		//Matrix A??KR????????????????�?????2D??????
		Vector3d  b(matx);
		//Matrix A??KR?????Matrix b?????????
		x = A.colPivHouseholderQr().solve(b);//???Ax=b  ??????????????????��???XYZ????
	}

	double* pKR = KR + nP*9;
	double t = pKR[6]*x(0) + pKR[7]*x(1) + pKR[8]*x(2);	

	//compute azimuth and elevation angle
	double dDAngle = atan2( x(0), x(2) );				//\Psi???????Azimuth angle-??��??
	double dHAngle = atan2( x(1), sqrt(x(0)*x(0)+ x(2)*x(2)) );//\theta???????Elevation angle
	//feature?????????????
	feature[0] = dDAngle;
	feature[1] = dHAngle;
	feature[2] = 0;

	if ( t < 0 )//??
		return true;
	else
		return false;
}


bool PBA::pba_initializeAssoArchor( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int nMI,
	int nAI,
	int FID,
	bool bLast )
{
	int nM = photo[nMI];                              //??��???
	int nA = photo[nAI];                              //??��???

	Vector3d  xM, xA;

	if (m_bProvideXYZ)
	{
		xM[0] = m_XYZ[FID*3]   - *(camera + nM*6+3);
		xM[1] = m_XYZ[FID*3+1] - *(camera + nM*6+4);
		xM[2] = m_XYZ[FID*3+2] - *(camera + nM*6+5);

		xA[0] = m_XYZ[FID*3]   - *(camera + nA*6+3);
		xA[1] = m_XYZ[FID*3+1] - *(camera + nA*6+4);
		xA[2] = m_XYZ[FID*3+2] - *(camera + nA*6+5);
	}
	else
	{
		//Main anchor ray
		double *ptr1 = m_KR + nM*9;
		Matrix3d  AM;	
		AM << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

		double matxM[3];
		matxM[0] = *(imgpts+2*nMI);
		matxM[1] = *(imgpts+2*nMI+1);
		matxM[2] = 1;

		Vector3d  bM(matxM);
		xM = AM.colPivHouseholderQr().solve(bM);			//?????????????��???XYZ????

		//Associate archor ray
		double *ptr2 = m_KR + nA*9;
		Matrix3d  AA;	
		AA << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

		double matxA[3];
		matxA[0] = *(imgpts+2*nAI);
		matxA[1] = *(imgpts+2*nAI+1);
		matxA[2] = 1;

		Vector3d  bA(matxA);
		xA = AA.colPivHouseholderQr().solve(bA);			//????????????��???XYZ????
	}
		
	//Parallax Angle
	double dDot = xM(0)*xA(0) + xM(1)*xA(1) + xM(2)*xA(2);			//IJRR??????3.13???????

	double dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );		//IJRR??????3.13??????
	double dDisA = sqrt( xA(0)*xA(0)+xA(1)*xA(1)+xA(2)*xA(2) );		//IJRR??????3.13??????

	if (dDot/(dDisM*dDisA)>1)			
		feature[2] = 0;
	else if (dDot/(dDisM*dDisA)<-1)		
		feature[2] = PI;
	else
	{
		double dw = acos( dDot/(dDisM*dDisA) );
		feature[2] = dw;
	}

	
	double pti2k[3];
	//+3??X????????��??-??��???????t^m_a??????????��??m????��?????????X????XA-XM
	pti2k[0] = *(camera + nA*6+3) - *(camera + nM*6+3);    
	//YA-YM??Y????
	pti2k[1] = *(camera + nA*6+4) - *(camera + nM*6+4);    
	//ZA-ZM??Z????
	pti2k[2] = *(camera + nA*6+5) - *(camera + nM*6+5);    
	//????????pti2k????????t_m->a

	//dDot1??????��????????F_j??x????xM??t_m->a????
	double dDot1 = xM[0]*pti2k[0] + xM[1]*pti2k[1] + xM[2]*pti2k[2];
	double dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
	double tmp = dDot1/(dDisM*dDisi2k);//???????IJRR?2??/varphi??cos
	double dW2;
	if (tmp > 1)						
		dW2 = 0;
	if ( tmp < -1)						
		dW2 = PI;
	else
		dW2  = acos( tmp );//???????IJRR?2??/varphi

	return true;//dW2??/varphi_j;dw=feature[2]??/omega_j;????????????dW2?????????tmp?��????IJRR??BA????��?
}

bool PBA::pba_initializeOtheArchors( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int* archorSort,
	int nfeacout,
	int nOI,
	int FID )
{
	static int i = 0;
	double dw = feature[2];                   //????
	double dwNew;
	double dmaxw = dw;
	int   nNewI = 0;
	bool bAdjust = false;
	double dDot,dDisM,dDisA;

	if ( dw < MAXARCHOR   )
	{
		//current archor vector 
		int nO = photo[nOI];

		Vector3d  xO;
		if ( m_bProvideXYZ)
		{
			xO(0) = m_XYZ[FID*3] - *(camera + nO*6 + 3);
			xO(1) = m_XYZ[FID*3+1]-*(camera + nO*6 + 4);
			xO(2) = m_XYZ[FID*3+2]-*(camera + nO*6 + 5);
		}
		else
		{
			double *ptr1 = m_KR + nO*9;
			Matrix3d  AO;	
			AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

			double matxO[3];
			matxO[0] = *(imgpts+nOI*2);
			matxO[1] = *(imgpts+nOI*2+1);
			matxO[2] = 1;

			Vector3d  bO(matxO);
			xO = AO.colPivHouseholderQr().solve(bO);
		}

		double dDAngle = atan2( xO(0), xO(2) );
		double dHAngle = atan2( xO(1), sqrt(xO(0)*xO(0)+ xO(2)*xO(2)) );

		for ( i = 0; i < nfeacout; i++ )
		{
			//Main Archor Vector
			int nM = photo[i];                            
			Vector3d  xM;

			if (m_bProvideXYZ)
			{
				xM(0) = m_XYZ[FID*3]  -*(camera + nM*6+3);
				xM(1) = m_XYZ[FID*3+1]-*(camera + nM*6+4);
				xM(2) = m_XYZ[FID*3+2]-*(camera + nM*6+5);
			}
			else
			{
				double *ptr2 = m_KR + nM*9;
				Matrix3d  AM;	
				AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

				double matxM[3];
				matxM[0] = *(imgpts+i*2);
				matxM[1] = *(imgpts+i*2+1);
				matxM[2] = 1;

				Vector3d  bM(matxM);
				xM = AM.colPivHouseholderQr().solve(bM);
			}

			//Parallax angle between current archor and main archor
			dDot = xM(0)*xO(0) + xM(1)*xO(1) + xM(2)*xO(2);
			dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );
			dDisA = sqrt( xO(0)*xO(0)+xO(1)*xO(1)+xO(2)*xO(2) );

			if( dDot/(dDisM*dDisA) > 1 )
				dwNew = 0;
			else if(dDot/(dDisM*dDisA)<-1)
				dwNew = PI;
			else
				dwNew = acos( dDot/(dDisM*dDisA) );

			if ( dwNew > dmaxw )
			{
				dmaxw = dwNew;
				archorSort[0] = nOI;   
				archorSort[1] = i;     
				feature[0] = dDAngle;
				feature[1] = dHAngle;
				feature[2] = dmaxw;
				bAdjust = true;
			}
		}
	}
	return bAdjust;
}



