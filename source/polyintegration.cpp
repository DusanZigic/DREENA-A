#include "polyintegration.hpp"

#include <vector>
#include <cmath>

template <typename T>
static size_t locatePoint(const std::vector<T> &data, T x, int interpolationOrder)
{
	int ju, jm, jl;
	int mm = interpolationOrder + 1;
	int n = data.size();
	bool ascnd = (data.back() >= data.front());
	jl = 0;
	ju = n - 1;
	while (ju - jl > 1)
	{
		jm = (ju + jl) >> 1;
		if ((x >= data[jm]) == ascnd) {
			jl = jm;
		}
		else {
			ju = jm;
		}
	}
	int pointLocation = std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
	return static_cast<size_t>(pointLocation);
}

template <typename T>
static void polynomialCoeff(const std::vector<T> &dataX, const std::vector<T> &dataF, std::vector<double> &coeff)
{
	size_t n = dataX.size();
	coeff.resize(n, 0.0);
	double phi, ff, b;
	std::vector<double> s(n, 0.0);	
	s[n-1] = -static_cast<double>(dataX[0]);

	for (size_t i=1; i<n; i++) {
		for (size_t j=n-1-i; j<n-1; j++)
			s[j] -= static_cast<double>(dataX[i]) * s[j+1];
		s[n-1] -= static_cast<double>(dataX[i]);
	}

	for (size_t j=0; j<n; j++) {
		phi = static_cast<double>(n);
		
		for (size_t k=n-1; k>0; k--)
			phi = static_cast<double>(k)*s[k] + static_cast<double>(dataX[j])*phi;
		
		ff = static_cast<double>(dataF[j])/phi;		
		b = 1.0;
		
		for (int k=n-1; k>=0; k--) {
			coeff[k] += b * ff;
			b = s[k] + static_cast<double>(dataX[j]) * b;
		}
	}
}

template <typename T>
T poly::linearIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata)
{
	if (xdata.size() < 2) return 0.0;

	std::vector<double> k, c;
	for (size_t i=0; i<(xdata.size()-1); i++)
	{
		k.push_back(static_cast<double>(fdata[i+1]-fdata[i])/static_cast<double>(xdata[i+1]-xdata[i]));
		c.push_back(static_cast<double>(fdata[i])-static_cast<double>(k.back())*static_cast<double>(xdata[i]));
	}

	double res = 0.0;

	for (size_t i=0; i<(xdata.size()-1); i++)
		res += 0.5*k[i]*(static_cast<double>(xdata[i+1]*xdata[i+1]) - static_cast<double>(xdata[i]*xdata[i])) + c[i]*static_cast<double>(xdata[i+1] - xdata[i]);

	return res;
}
template float poly::linearIntegrate<float>(const std::vector<float> &xdata, const std::vector<float> &fdata);
template double poly::linearIntegrate<double>(const std::vector<double> &xdata, const std::vector<double> &fdata);
template long double poly::linearIntegrate<long double>(const std::vector<long double> &xdata, const std::vector<long double> &fdata);

template <typename T>
T poly::linearIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata, T lowLimit, T highLimit)
{
	if (xdata.size() < 2) return 0.0;

	std::vector<double> k, c;
	for (size_t i=0; i<(xdata.size()-1); i++)
	{
		k.push_back(static_cast<double>(fdata[i+1]-fdata[i])/static_cast<double>(xdata[i+1]-xdata[i]));
		c.push_back(static_cast<double>(fdata[i]-k.back()*xdata[i]));
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of full integral (in it's whole range):
	double sum = 0.0L;

	for (size_t i=0; i<(xdata.size()-1); i++)
		sum += 0.5*k[i]*(static_cast<double>(xdata[i+1]*xdata[i+1]) - static_cast<double>(xdata[i]*xdata[i])) + 
							c[i]*static_cast<double>(xdata[i+1] - xdata[i]);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from lower range to lower limit:
	size_t lowLimitPos = locatePoint(xdata, lowLimit, 1);

	double lowSum = 0.0L;

	for (size_t i=0; i<lowLimitPos; i++)
		lowSum += 0.5*k[i]*(static_cast<double>(xdata[i+1]*xdata[i+1]) - static_cast<double>(xdata[i]*xdata[i])) + 
								c[i]*static_cast<double>(xdata[i+1] - xdata[i]);

	lowSum += 0.5*k[lowLimitPos]*(static_cast<double>(lowLimit*lowLimit) - static_cast<double>(xdata[lowLimitPos]*xdata[lowLimitPos])) + 
									c[lowLimitPos]*static_cast<double>(lowLimit - xdata[lowLimitPos]);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from higher limit to higer range:
	size_t highLimitPos = locatePoint(xdata, highLimit, 1);

	double highSum = 0.0L;

	highSum += 0.5*k[highLimitPos]*(static_cast<double>(xdata[highLimitPos+1]*xdata[highLimitPos+1]) - static_cast<double>(highLimit*highLimit)) + 
						c[highLimitPos]*static_cast<double>(xdata[highLimitPos+1] - highLimit);

	for (size_t i=highLimitPos+1; i<xdata.size()-1; i++)
		highSum += 0.5*k[i]*(static_cast<double>(xdata[i+1]*xdata[i+1]) - static_cast<double>(xdata[i]*xdata[i])) + 
								c[i]*static_cast<double>(xdata[i+1] - xdata[i]);	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//integral value is full-low-high
	return (sum - highSum - lowSum);
}
template float poly::linearIntegrate<float>(const std::vector<float> &xdata, const std::vector<float> &fdata, float lowLimit, float highLimit);
template double poly::linearIntegrate<double>(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit);
template long double poly::linearIntegrate<long double>(const std::vector<long double> &xdata, const std::vector<long double> &fdata, long double lowLimit, long double highLimit);

template <typename T>
T poly::cubicIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata)
{
	if (xdata.size() < 2) return 0;

	//calculating polynomial coefficients for each segment:
	std::vector<std::vector<double>> coefficents; coefficents.resize(xdata.size() - 1);
	for (size_t i=0; i<coefficents.size(); i++)
	{
		size_t pointLocation = locatePoint(xdata, xdata[i], 3);
		std::vector<T> xdatatemp(xdata.begin()+pointLocation, xdata.begin()+pointLocation+4);
		std::vector<T> fdatatemp(fdata.begin()+pointLocation, fdata.begin()+pointLocation+4);
		polynomialCoeff(xdatatemp, fdatatemp, coefficents[i]);
	}

	//calculating value of integral:
	double sum = 0.0L;
	for (size_t i=0; i<xdata.size()-1; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			sum += 1.0/static_cast<double>(j+1)*coefficents[i][j]*(std::pow(static_cast<double>(xdata[i+1]), static_cast<double>(j+1)) - std::pow(static_cast<double>(xdata[i]), static_cast<double>(j+1)));

	return sum;
}
template float poly::cubicIntegrate<float>(const std::vector<float> &xdata, const std::vector<float> &fdata);
template double poly::cubicIntegrate<double>(const std::vector<double> &xdata, const std::vector<double> &fdata);
template long double poly::cubicIntegrate<long double>(const std::vector<long double> &xdata, const std::vector<long double> &fdata);

template <typename T>
T poly::cubicIntegrate(const std::vector<T> &xdata, const std::vector<T> &fdata, T lowLimit, T highLimit)
{
	if (xdata.size() < 2) return 0;

	//calculating polynomial coefficients for each segment:
	std::vector<std::vector<double>> coefficents; coefficents.resize(xdata.size() - 1);
	for (size_t i=0; i<coefficents.size(); i++)
	{
		size_t pointLocation = locatePoint(xdata, xdata[i], 3);
		std::vector<T> xdatatemp(xdata.begin()+pointLocation, xdata.begin()+pointLocation+4);
		std::vector<T> fdatatemp(fdata.begin()+pointLocation, fdata.begin()+pointLocation+4);
		polynomialCoeff(xdatatemp, fdatatemp, coefficents[i]);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of full integral (in it's whole range):
	double sum = 0.0L;
	for (size_t i=0; i<xdata.size()-1; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			sum += 1.0L/static_cast<double>(j+1)*coefficents[i][j]*
                        (std::pow(static_cast<double>(xdata[i+1]), static_cast<double>(j+1)) - std::pow(static_cast<double>(xdata[i]), static_cast<double>(j+1)));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from lower range to lower limit:
	size_t lowLimitPos = locatePoint(xdata, lowLimit, 1);

	double lowSum = 0.0L;

	for (size_t i=0; i<lowLimitPos; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			lowSum += 1.0L/static_cast<double>(j+1)*coefficents[i][j]*
                            (std::pow(static_cast<double>(xdata[i+1]), static_cast<double>(j+1)) - std::pow(static_cast<double>(xdata[i]), static_cast<double>(j+1)));

	for (size_t j=0; j<coefficents[lowLimitPos].size(); j++)
		lowSum += 1.0L/static_cast<double>(j+1)*coefficents[lowLimitPos][j]*
                        (std::pow(static_cast<double>(lowLimit), static_cast<double>(j+1)) - std::pow(static_cast<double>(xdata[lowLimitPos]), static_cast<double>(j+1)));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from higher limit to higer range:
	size_t highLimitPos = locatePoint(xdata, highLimit, 1);

	double highSum = 0.0L;

	for (size_t j=0; j<coefficents[highLimitPos].size(); j++)
		highSum += 1.0L/static_cast<double>(j+1)*coefficents[highLimitPos][j]*
                        (std::pow(static_cast<double>(xdata[highLimitPos+1]), static_cast<double>(j+1)) - std::pow(static_cast<double>(highLimit), static_cast<double>(j+1)));

	for (size_t i=highLimitPos+1; i<xdata.size()-1; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			highSum += 1.0L/static_cast<double>(j+1)*coefficents[i][j]*
                            (std::pow(static_cast<double>(xdata[i+1]), static_cast<double>(j+1)) - std::pow(static_cast<double>(xdata[i]), static_cast<double>(j+1)));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//integral value is full-low-high
	return (sum - highSum - lowSum);
}
template float poly::cubicIntegrate<float>(const std::vector<float> &xdata, const std::vector<float> &fdata, float lowLimit, float highLimit);
template double poly::cubicIntegrate<double>(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit);
template long double poly::cubicIntegrate<long double>(const std::vector<long double> &xdata, const std::vector<long double> &fdata, long double lowLimit, long double highLimit);