#include <vector>
#include <cmath>

using namespace std;

#include "ELossHeader.hpp"
#include "Arsenal.hpp"
#include "LinearInterpolation.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FdA INTEGRALS:

static int FdAMaxPoints2, FdAMaxPoints3, FdAMaxPoints4, FdAMaxPoints5; //number of points for dA integration
static vector<double> FdAHS2, FdAHS3, FdAHS4, FdAHS5;				   //vectors that store Halton sequences for dA integrals

//function that initializes Halton sequences for dA integrals calculations:
void FdAHaltonSeqInit(int FdAMaxPts)
{
	FdAMaxPoints2 = FdAMaxPts; 	  //setting values of dAMaxPoints
	FdAMaxPoints3 = FdAMaxPts-25;
	FdAMaxPoints4 = FdAMaxPts-50;
	FdAMaxPoints5 = FdAMaxPts-75;
	
	for (int i=0; i<FdAMaxPts; i++) //generating Halton sequences
	{
		FdAHS2.push_back(HaltonSequence((i+1)*409, 2));
		FdAHS3.push_back(HaltonSequence((i+1)*409, 3));
		FdAHS4.push_back(HaltonSequence((i+1)*409, 5));
		FdAHS5.push_back(HaltonSequence((i+1)*409, 7));
	}
}

//dAp410 integral
double dAp410(double ph, interpFun &normint)
{
	double res = 1.0 / exp(normint.interp(ph));

	return res;
}

//FdA411 integral
double FdA411(double ph, double dp, interpFun &normint, interpFun &dndxint)
{
	double res = 1.0 / exp(normint.interp(ph + dp))*dndxint.interp(ph + dp, 1.0 - ph/(ph + dp));

	return res;
}

//FdA412 integral
double FdA412(double ph, double dp, interpFun &normint, interpFun &dndxint)
{
	if (dp < 2.0*mgC / 2.0) return 0.0;

	double p = ph + dp;

	double yl, yh, yq, y;

	double sum = 0.0;

	for (int i=0; i<FdAMaxPoints2; i++)
	{
		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;

		sum += 1.0 / exp(normint.interp(p))*(1.0 / 2.0)*dndxint.interp(p, 1.0 - ph/p - y)*
			dndxint.interp(p, y)*(yh - yl);
	}

	return (sum/FdAMaxPoints2);
}

//FdA413 integral
double FdA413(double ph, double dp, interpFun &normint, interpFun &dndxint)
{
	if (dp < 3.0*mgC / 2.0) return 0.0;

	double p = ph + dp;

	double yl, yh, yq, y;

	double zl, zh, zq, z;

	double sum = 0.0;

	for (int i=0; i<FdAMaxPoints3; i++)
	{
		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + FdAHS3[i]*zq;

		sum += 1.0 / exp(normint.interp(p))*(1.0 / 2.0 / 3.0)*dndxint.interp(p, 1 - ph/p - y - z)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*(yh - yl)*(zh - zl);
	}

	return (sum/FdAMaxPoints3);
}

//FdA414 integral
double FdA414(double ph, double dp, interpFun &normint, interpFun &dndxint)
{
	if (dp < 4.0*mgC / 2.0) return 0.0;

	double p = ph + dp;

	double yl, yh, yq, y;

	double zl, zh, zq, z;

	double zzl, zzh, zzq, zz;

	double sum = 0.0;

	for (int i=0; i<FdAMaxPoints4; i++)
	{
		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 3.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + FdAHS3[i]*zq;

		zzl = mgC/(p + sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - mgC/(p + sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + FdAHS4[i]*zzq;

		sum += 1.0 / exp(normint.interp(p))*(1.0 / 2.0 / 3.0 / 4.0)*dndxint.interp(p, 1.0 - ph/p - y - z - zz)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*dndxint.interp(p, zz)*(yh - yl)*(zh - zl)*(zzh - zzl);
	}

	return (sum/FdAMaxPoints4);
}

//FdA415 integral
double FdA415(double ph, double dp, interpFun &normint, interpFun &dndxint)
{
	if (dp < 5.0*mgC / 2.0) return 0.0;

	double p = ph + dp;

	double yl, yh, yq, y;

	double zl, zh, zq, z;

	double zzl, zzh, zzq, zz;

	double zzzl, zzzh, zzzq, zzz;

	double sum = 0.0;

	for (int i=0; i<FdAMaxPoints5; i++)
	{
		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 4.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 3.0*mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + FdAHS3[i]*zq;

		zzl = mgC/(p + sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + FdAHS4[i]*zzq;

		zzzl = mgC/(p + sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - mgC/(p + sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + FdAHS5[i]*zzzq;

		sum += 1.0 / exp(normint.interp(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0)*dndxint.interp(p, 1.0 - ph/p - y - z - zz - zzz)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*dndxint.interp(p, zz)*dndxint.interp(p, zzz)*(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl);

	}

	return (sum/FdAMaxPoints5);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//dA INTEGRALS:

static int dAMaxPoints1, dAMaxPoints2, dAMaxPoints3, dAMaxPoints4, dAMaxPoints5, dAMaxPoints6, dAMaxPoints7; //number of points for dA integration
static vector<double> dAHS1, dAHS2, dAHS3, dAHS4, dAHS5, dAHS6, dAHS7; 										 //vectors that store Halton sequences for dA integrals

//function that initializes Halton sequences for dA integrals calculations:
void dAHaltonSeqInit(int dAMaxPts)
{
	dAMaxPoints1 = dAMaxPts;	 //setting values of dAMaxPoints
	dAMaxPoints2 = dAMaxPts-100;
	dAMaxPoints3 = dAMaxPts-200;
	dAMaxPoints4 = dAMaxPts-300;
	dAMaxPoints5 = dAMaxPts-400;
	dAMaxPoints6 = dAMaxPts-500;
	dAMaxPoints7 = dAMaxPts-600;

	for (int i=0; i<dAMaxPts; i++) //generating Halton sequences
	{
		dAHS1.push_back(HaltonSequence((i+1)*409, 2));
		dAHS2.push_back(HaltonSequence((i+1)*409, 3));
		dAHS3.push_back(HaltonSequence((i+1)*409, 5));
		dAHS4.push_back(HaltonSequence((i+1)*409, 7));
		dAHS5.push_back(HaltonSequence((i+1)*409, 11));
		dAHS6.push_back(HaltonSequence((i+1)*409, 13));
		dAHS7.push_back(HaltonSequence((i+1)*409, 17));
	}
}

//dA410 integral
double dA410(double ph, interpFun &normint)
{
	double res = dsdpti2.interp(ph)/exp(normint.interp(ph));

	return res;
}

//dA411 integral
double dA411(double ph, interpFun &normint, interpFun &dndxint)
{
	double p1 = ph + mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;

	double sum = 0.0;

	for (int i=0; i<dAMaxPoints1; i++)
	{
		p = p1 + dAHS1[i]*pq;

		sum += dsdpti2.interp(p) / p / exp(normint.interp(p))*dndxint.interp(p, 1.0 - ph/p);
	}

	return (sum*pq/dAMaxPoints1);
}

//dA412 integral
double dA412(double ph, interpFun &normint, interpFun &dndxint)
{
	double p1 = ph + 2.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;

	double sum = 0.0;

	for (int i=0; i<dAMaxPoints2; i++)
	{
		p = p1 + dAHS1[i]*pq;

		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;

		sum += dsdpti2.interp(p) / p / exp(normint.interp(p))*(1.0 / 2.0)*dndxint.interp(p, 1.0 - ph / p - y)*
			dndxint.interp(p, y)*(yh - yl);
	}

	return (sum*pq / dAMaxPoints2);
}

//dA413 integral
double dA413(double ph, interpFun &normint, interpFun &dndxint)
{
	double p1 = ph + 3.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;

	double zl, zh, zq, z;

	double sum = 0.0;

	for (int i=0; i<dAMaxPoints3; i++)
	{
		p = p1 + dAHS1[i]*pq;

		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;

		sum += dsdpti2.interp(p) / p / exp(normint.interp(p))*(1.0 / 2.0 / 3.0)*dndxint.interp(p, 1 - ph/p - y - z)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*(yh - yl)*(zh - zl);
	}

	return (sum*pq/dAMaxPoints3);
}

//dA414 integral
double dA414(double ph, interpFun &normint, interpFun &dndxint)
{
	double p1 = ph + 4.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;

	double zl, zh, zq, z;

	double zzl, zzh, zzq, zz;

	double sum = 0.0;

	for (int i=0; i<dAMaxPoints4; i++)
	{
		p = p1 + dAHS1[i]*pq;

		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 3.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;

		zzl = mgC/(p + sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - mgC/(p + sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;

		sum += dsdpti2.interp(p) / p / exp(normint.interp(p))*(1.0 / 2.0 / 3.0 / 4.0)*dndxint.interp(p, 1.0 - ph/p - y - z - zz)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*dndxint.interp(p, zz)*(yh - yl)*(zh - zl)*(zzh - zzl);
	}

	return (sum*pq/dAMaxPoints4);
}

//dA415 integral
double dA415(double ph, interpFun &normint, interpFun &dndxint)
{
	double p1 = ph + 5.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;

	double zl, zh, zq, z;

	double zzl, zzh, zzq, zz;

	double zzzl, zzzh, zzzq, zzz;

	double sum = 0.0;

	for (int i=0; i<dAMaxPoints5; i++)
	{
		p = p1 + dAHS1[i]*pq;

		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 4.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 3.0*mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;

		zzl = mgC/(p + sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;

		zzzl = mgC/(p + sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - mgC/(p + sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + dAHS5[i]*zzzq;

		sum += dsdpti2.interp(p) / p / exp(normint.interp(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0)*dndxint.interp(p, 1.0 - ph/p - y - z - zz - zzz)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*dndxint.interp(p, zz)*dndxint.interp(p, zzz)*(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl);

	}

	return (sum*pq/dAMaxPoints5);
}

//dA416 integral
double dA416(double ph, interpFun &normint, interpFun &dndxint)
{
	double p1 = ph + 6.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;

	double zl, zh, zq, z;

	double zzl, zzh, zzq, zz;

	double zzzl, zzzh, zzzq, zzz;

	double zzzzl, zzzzh, zzzzq, zzzz;

	double sum = 0.0;

	for (int i=0; i<dAMaxPoints6; i++)
	{
		p = p1 + dAHS1[i]*pq;

		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 5.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 4.0*mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;

		zzl = mgC/(p + sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 3.0*mgC/(p + sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;

		zzzl = mgC/(p + sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + dAHS5[i]*zzzq;

		zzzzl = mgC/(p + sqrt(MC*MC + p*p));
		zzzzh = 1.0 - ph/p - y - z - zz - zzz - mgC/(p + sqrt(MC*MC + p*p));
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + dAHS6[i]*zzzzq;

		sum += dsdpti2.interp(p) / p / exp(normint.interp(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0 / 6.0)*dndxint.interp(p, 1.0 - ph/p - y - z - zz - zzz - zzzz)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*dndxint.interp(p, zz)*dndxint.interp(p, zzz)*dndxint.interp(p, zzzz)*(yh - yl)*(zh - zl)*
			(zzh - zzl)*(zzzh - zzzl)*(zzzzh - zzzzl);

	}

	return (sum*pq/dAMaxPoints6);
}

//dA417 integral
double dA417(double ph, interpFun &normint, interpFun &dndxint)
{
	double p1 = ph + 7.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;

	double yl, yh, yq, y; //dA412

	double zl, zh, zq, z; //dA413

	double zzl, zzh, zzq, zz; //dA414

	double zzzl, zzzh, zzzq, zzz; //dA415

	double zzzzl, zzzzh, zzzzq, zzzz; //dA416

	double zzzzzl, zzzzzh, zzzzzq, zzzzz; //dA417

	double sum = 0.0;

	for (int i=0; i<dAMaxPoints7; i++)
	{
		p = p1 + dAHS1[i]*pq;

		yl = mgC/(p + sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 6.0*mgC/(p + sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;

		zl = mgC/(p + sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 5.0*mgC/(p + sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;

		zzl = mgC/(p + sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 4.0*mgC/(p + sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;

		zzzl = mgC/(p + sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - 3.0*mgC/(p + sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + dAHS5[i]*zzzq;

		zzzzl = mgC/(p + sqrt(MC*MC + p*p));
		zzzzh = 1.0 - ph/p - y - z - zz - zzz - 2.0*mgC/(p + sqrt(MC*MC + p*p));
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + dAHS6[i]*zzzzq;

		zzzzzl = mgC/(p + sqrt(MC*MC + p*p));
		zzzzzh = 1.0 - ph/p - y - z - zz - zzz - zzzz - mgC/(p + sqrt(MC*MC + p*p));
		zzzzzq = zzzzzh - zzzzzl;
		zzzzz = zzzzzl + dAHS7[i]*zzzzzq;

		sum += dsdpti2.interp(p) / p / exp(normint.interp(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0 / 6.0 / 7.0)*dndxint.interp(p, 1.0 - ph/p - y - z - zz - zzz - zzzz - zzzzz)*
			dndxint.interp(p, y)*dndxint.interp(p, z)*dndxint.interp(p, zz)*dndxint.interp(p, zzz)*dndxint.interp(p, zzzz)*dndxint.interp(p, zzzzz)*
			(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl)*(zzzzh - zzzzl)*(zzzzzh - zzzzzl);

	}

	return (sum*pq/dAMaxPoints7);
}