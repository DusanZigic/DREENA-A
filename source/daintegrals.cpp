#include "energyloss.hpp"
#include "linearinterpolation.hpp"

#include <vector>
#include <cmath>

double energyLoss::haltonSequence(int index, int base) const
{
	double f = 1.0;
	double res = 0.0;

	while (index > 0) {
		f = f / static_cast<double>(base);
		res += f * static_cast<double>(index % base);
		index = index / base; // integer division
	}

	return res;
}


void energyLoss::FdAHaltonSeqInit(size_t FdAMaxPts)
{
	m_FdAMaxPoints2 = FdAMaxPts; 	  //setting values of dAMaxPoints
	m_FdAMaxPoints3 = FdAMaxPts-25;
	m_FdAMaxPoints4 = FdAMaxPts-50;
	m_FdAMaxPoints5 = FdAMaxPts-75;
	
	for (size_t i=0; i<FdAMaxPts; i++) //generating Halton sequences
	{
		m_FdAHS2.push_back(haltonSequence((i+1)*409, 2));
		m_FdAHS3.push_back(haltonSequence((i+1)*409, 3));
		m_FdAHS4.push_back(haltonSequence((i+1)*409, 5));
		m_FdAHS5.push_back(haltonSequence((i+1)*409, 7));
	}
}

double energyLoss::dAp410(double ph, const interpolationF<double> &norm) const {
	return (1.0 / std::exp(norm.interpolation(ph)));
}

double energyLoss::FdA411(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	return (1.0 / std::exp(norm.interpolation(ph + dp))*dndx.interpolation(ph + dp, 1.0 - ph/(ph + dp)));
}

double energyLoss::FdA412(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const{
	if (dp < 2.0*m_mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double sum = 0.0;
	for (size_t i=0; i<m_FdAMaxPoints2; i++) {
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_FdAHS2[i]*yq;
		
		sum += 1.0 / std::exp(norm.interpolation(p))*(1.0 / 2.0)*dndx.interpolation(p, 1.0 - ph/p - y)*
			dndx.interpolation(p, y)*(yh - yl);
	}

	return (sum/static_cast<double>(m_FdAMaxPoints2));
}

double energyLoss::FdA413(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	if (dp < 3.0*m_mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double sum = 0.0;
	for (size_t i=0; i<m_FdAMaxPoints3; i++) {
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_FdAHS2[i]*yq;
		
		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_FdAHS3[i]*zq;
		
		sum += 1.0 / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0)*dndx.interpolation(p, 1.0 - ph/p - y - z)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*(yh - yl)*(zh - zl);
	}

	return (sum/static_cast<double>(m_FdAMaxPoints3));
}

double energyLoss::FdA414(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	if (dp < 4.0*m_mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double sum = 0.0;
	for (size_t i=0; i<m_FdAMaxPoints4; i++) {
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 3.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_FdAHS2[i]*yq;
		
		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_FdAHS3[i]*zq;
		
		zzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzh = 1.0 - ph/p - y - z - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + m_FdAHS4[i]*zzq;
		
		sum += 1.0 / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0)*dndx.interpolation(p, 1.0 - ph/p - y - z - zz)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*dndx.interpolation(p, zz)*(yh - yl)*(zh - zl)*(zzh - zzl);
	}

	return (sum/static_cast<double>(m_FdAMaxPoints4));
}

double energyLoss::FdA415(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	if (dp < 5.0*m_mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double sum = 0.0;
	for (size_t i=0; i<m_FdAMaxPoints5; i++) {
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 4.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;		
		y = yl + m_FdAHS2[i]*yq;
		
		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - 3.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_FdAHS3[i]*zq;

		zzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzh = 1.0 - ph/p - y - z - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + m_FdAHS4[i]*zzq;

		zzzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_FdAHS5[i]*zzzq;

		sum += 1.0 / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0)*dndx.interpolation(p, 1.0 - ph/p - y - z - zz - zzz)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*dndx.interpolation(p, zz)*dndx.interpolation(p, zzz)*(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl);

	}

	return (sum/static_cast<double>(m_FdAMaxPoints5));
}

double energyLoss::FdA(double ph, double dp, const interpolationF<double> &currnorm, const interpolationF<double> &currdndx) const {
	return (FdA411(ph, dp, currnorm, currdndx) + FdA412(ph, dp, currnorm, currdndx) + FdA413(ph, dp, currnorm, currdndx) +
				FdA414(ph, dp, currnorm, currdndx) + FdA415(ph, dp, currnorm, currdndx));
}


void energyLoss::dAHaltonSeqInit(size_t dAMaxPts)
{
	m_dAMaxPoints1 = dAMaxPts;	 //setting values of dAMaxPoints
	m_dAMaxPoints2 = dAMaxPts-100;
	m_dAMaxPoints3 = dAMaxPts-200;
	m_dAMaxPoints4 = dAMaxPts-300;
	m_dAMaxPoints5 = dAMaxPts-400;
	m_dAMaxPoints6 = dAMaxPts-500;
	m_dAMaxPoints7 = dAMaxPts-600;

	for (size_t i=0; i<dAMaxPts; i++) //generating Halton sequences
	{
		m_dAHS1.push_back(haltonSequence((i+1)*409, 2));
		m_dAHS2.push_back(haltonSequence((i+1)*409, 3));
		m_dAHS3.push_back(haltonSequence((i+1)*409, 5));
		m_dAHS4.push_back(haltonSequence((i+1)*409, 7));
		m_dAHS5.push_back(haltonSequence((i+1)*409, 11));
		m_dAHS6.push_back(haltonSequence((i+1)*409, 13));
		m_dAHS7.push_back(haltonSequence((i+1)*409, 17));
	}
}

double energyLoss::dA410(double ph, const interpolationF<double> &norm) const {
	return (m_dsdpti2.interpolation(ph)/exp(norm.interpolation(ph)));
}

double energyLoss::dA411(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	double p1 = ph + m_mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double sum = 0.0;
	for (size_t i=0; i<m_dAMaxPoints1; i++) {
		p = p1 + m_dAHS1[i]*pq;
		
		sum += m_dsdpti2.interpolation(p) / p / std::exp(norm.interpolation(p))*dndx.interpolation(p, 1.0 - ph/p);
	}

	return (sum*pq/static_cast<double>(m_dAMaxPoints1));
}

double energyLoss::dA412(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	double p1 = ph + 2.0*m_mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double sum = 0.0;
	for (size_t i=0; i<m_dAMaxPoints2; i++) {
		p = p1 + m_dAHS1[i]*pq;
		
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_dAHS2[i]*yq;
		
		sum += m_dsdpti2.interpolation(p) / p / std::exp(norm.interpolation(p))*(1.0 / 2.0)*dndx.interpolation(p, 1.0 - ph / p - y)*
			dndx.interpolation(p, y)*(yh - yl);
	}

	return (sum*pq/static_cast<double>(m_dAMaxPoints2));
}

double energyLoss::dA413(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	double p1 = ph + 3.0*m_mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double sum = 0.0;
	for (size_t i=0; i<m_dAMaxPoints3; i++) {
		p = p1 + m_dAHS1[i]*pq;
		
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_dAHS2[i]*yq;

		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_dAHS3[i]*zq;

		sum += m_dsdpti2.interpolation(p) / p / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0)*dndx.interpolation(p, 1.0 - ph/p - y - z)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*(yh - yl)*(zh - zl);
	}

	return (sum*pq/static_cast<double>(m_dAMaxPoints3));
}

double energyLoss::dA414(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	double p1 = ph + 4.0*m_mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double sum = 0.0;
	for (size_t i=0; i<m_dAMaxPoints4; i++) {
		p = p1 + m_dAHS1[i]*pq;
		
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 3.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_dAHS2[i]*yq;

		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_dAHS3[i]*zq;

		zzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzh = 1.0 - ph/p - y - z - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i]*zzq;

		sum += m_dsdpti2.interpolation(p) / p / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0)*dndx.interpolation(p, 1.0 - ph/p - y - z - zz)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*dndx.interpolation(p, zz)*(yh - yl)*(zh - zl)*(zzh - zzl);
	}

	return (sum*pq/static_cast<double>(m_dAMaxPoints4));
}

double energyLoss::dA415(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	double p1 = ph + 5.0*m_mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double sum = 0.0;
	for (size_t i=0; i<m_dAMaxPoints5; i++) {
		p = p1 + m_dAHS1[i]*pq;
		
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 4.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_dAHS2[i]*yq;

		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - 3.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_dAHS3[i]*zq;

		zzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzh = 1.0 - ph/p - y - z - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i]*zzq;

		zzzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[i]*zzzq;

		sum += m_dsdpti2.interpolation(p) / p / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0)*dndx.interpolation(p, 1.0 - ph/p - y - z - zz - zzz)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*dndx.interpolation(p, zz)*dndx.interpolation(p, zzz)*(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl);
	}

	return (sum*pq/static_cast<double>(m_dAMaxPoints5));
}

double energyLoss::dA416(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	double p1 = ph + 6.0*m_mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double zzzzl, zzzzh, zzzzq, zzzz;
	double sum = 0.0;
	for (size_t i=0; i<m_dAMaxPoints6; i++) {
		p = p1 + m_dAHS1[i]*pq;
		
		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 5.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_dAHS2[i]*yq;

		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - 4.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_dAHS3[i]*zq;

		zzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzh = 1.0 - ph/p - y - z - 3.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i]*zzq;

		zzzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[i]*zzzq;

		zzzzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzzh = 1.0 - ph/p - y - z - zz - zzz - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + m_dAHS6[i]*zzzzq;

		sum += m_dsdpti2.interpolation(p) / p / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0 / 6.0)*dndx.interpolation(p, 1.0 - ph/p - y - z - zz - zzz - zzzz)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*dndx.interpolation(p, zz)*dndx.interpolation(p, zzz)*dndx.interpolation(p, zzzz)*(yh - yl)*(zh - zl)*
			(zzh - zzl)*(zzzh - zzzl)*(zzzzh - zzzzl);
	}

	return (sum*pq/static_cast<double>(m_dAMaxPoints6));
}

double energyLoss::dA417(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const {
	double p1 = ph + 7.0*m_mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double zzzzl, zzzzh, zzzzq, zzzz;
	double zzzzzl, zzzzzh, zzzzzq, zzzzz;
	double sum = 0.0;
	for (size_t i=0; i<m_dAMaxPoints7; i++) {
		p = p1 + m_dAHS1[i]*pq;

		yl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yh = 1.0 - ph/p - 6.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		yq = yh - yl;
		y = yl + m_dAHS2[i]*yq;

		zl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zh = 1.0 - ph/p - y - 5.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zq = zh - zl;
		z = zl + m_dAHS3[i]*zq;

		zzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzh = 1.0 - ph/p - y - z - 4.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i]*zzq;

		zzzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - 3.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[i]*zzzq;

		zzzzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzzh = 1.0 - ph/p - y - z - zz - zzz - 2.0*m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + m_dAHS6[i]*zzzzq;

		zzzzzl = m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzzzh = 1.0 - ph/p - y - z - zz - zzz - zzzz - m_mgC/(p + std::sqrt(m_MC*m_MC + p*p));
		zzzzzq = zzzzzh - zzzzzl;
		zzzzz = zzzzzl + m_dAHS7[i]*zzzzzq;

		sum += m_dsdpti2.interpolation(p) / p / std::exp(norm.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0 / 6.0 / 7.0)*dndx.interpolation(p, 1.0 - ph/p - y - z - zz - zzz - zzzz - zzzzz)*
			dndx.interpolation(p, y)*dndx.interpolation(p, z)*dndx.interpolation(p, zz)*dndx.interpolation(p, zzz)*dndx.interpolation(p, zzzz)*dndx.interpolation(p, zzzzz)*
			(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl)*(zzzzh - zzzzl)*(zzzzzh - zzzzzl);
	}

	return (sum*pq/static_cast<double>(m_dAMaxPoints7));
}

double energyLoss::dA41(double ph, interpolationF<double> &currnorm, interpolationF<double> &currdndx) const {
	if (m_pName == "Gluon") { //gluon needs 7 dA integrals
		return (dA410(ph, currnorm) + dA411(ph, currnorm, currdndx) + dA412(ph, currnorm, currdndx) +dA413(ph, currnorm, currdndx) +
					dA414(ph, currnorm, currdndx) + dA415(ph, currnorm, currdndx) + dA416(ph, currnorm, currdndx) +
					dA417(ph, currnorm, currdndx));
	}
	else { //light quarks need 5 dA integrals
		return (dA410(ph, currnorm) + dA411(ph, currnorm, currdndx) + dA412(ph, currnorm, currdndx) + dA413(ph, currnorm, currdndx) +
					dA414(ph, currnorm, currdndx) + dA415(ph, currnorm, currdndx));
	}
}