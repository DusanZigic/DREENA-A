#include "linearinterpolation.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

//CONSTRUCTORS:
template <typename T>
interpolationF<T>::interpolationF() {}

//input is 2 1D arrays:
template <typename T>
interpolationF<T>::interpolationF(const T *xData, const T *fData, size_t NofElements)
{
	setData(xData, fData, NofElements);
}

template <typename T>
void interpolationF<T>::setData(const T *xData, const T *fData, size_t NofElements)
{
	m_variableN = 1;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(xData, xData + m_dataLength);
    m_data[1] = std::vector<T>(fData, fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 2 1D vectors:
template <typename T>
interpolationF<T>::interpolationF(const std::vector<T> &xData, const std::vector<T> &fData)
{
	setData(xData, fData);
}

template <typename T>
void interpolationF<T>::setData(const std::vector<T> &xData, const std::vector<T> &fData)
{
	m_variableN = 1;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(xData.begin(), xData.begin() + m_dataLength);
    m_data[1] = std::vector<T>(fData.begin(), fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 3 1D arrays:
template <typename T>
interpolationF<T>::interpolationF(const T *x1Data, const T *x2Data, const T *fData, size_t NofElements)
{
	setData(x1Data, x2Data, fData, NofElements);
}

template <typename T>
void interpolationF<T>::setData(const T *x1Data, const T *x2Data, const T *fData, size_t NofElements)
{
	m_variableN = 2;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(x1Data, x1Data + m_dataLength);
	m_data[1] = std::vector<T>(x2Data, x2Data + m_dataLength);
    m_data[2] = std::vector<T>( fData,  fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 3 1D vectors:
template <typename T>
interpolationF<T>::interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &fData)
{
	setData(x1Data, x2Data, fData);
}

template <typename T>
void interpolationF<T>::setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &fData)
{
	m_variableN = 2;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(x1Data.begin(), x1Data.begin() + m_dataLength);
	m_data[1] = std::vector<T>(x2Data.begin(), x2Data.begin() + m_dataLength);
    m_data[2] = std::vector<T>( fData.begin(),  fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 2 1D vectors (grids) and 1 2d vector (function values):
template <typename T>
interpolationF<T>::interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<std::vector<T>> &fData)
{
	setData(x1Data, x2Data, fData);
}

template <typename T>
void interpolationF<T>::setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<std::vector<T>> &fData)
{
	m_variableN = 2;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(x1Data.begin(), x1Data.end());
	m_data[1] = std::vector<T>(x2Data.begin(), x2Data.end());
	for (const auto &row : fData)
		for (const auto &elem : row)
			m_data[2].push_back(elem);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 4 1D arrays:
template <typename T>
interpolationF<T>::interpolationF(const T *x1Data, const T *x2Data, const T *x3Data, const T *fData, size_t NofElements)
{
	setData(x1Data, x2Data, x3Data, fData, NofElements);
}

template <typename T>
void interpolationF<T>::setData(const T *x1Data, const T *x2Data, const T *x3Data, const T *fData, size_t NofElements)
{
	m_variableN = 3;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(x1Data, x1Data + m_dataLength);
	m_data[1] = std::vector<T>(x2Data, x2Data + m_dataLength);
	m_data[2] = std::vector<T>(x3Data, x3Data + m_dataLength);
    m_data[3] = std::vector<T>( fData,  fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 4 1D vectors:
template <typename T>
interpolationF<T>::interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &fData)
{
	setData(x1Data, x2Data, x3Data, fData);
}

template <typename T>
void interpolationF<T>::setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &fData)
{
	m_variableN = 3;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(x1Data.begin(), x1Data.begin() + m_dataLength);
	m_data[1] = std::vector<T>(x2Data.begin(), x2Data.begin() + m_dataLength);
	m_data[2] = std::vector<T>(x3Data.begin(), x3Data.begin() + m_dataLength);
    m_data[3] = std::vector<T>( fData.begin(),  fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 5 1D arrays:
template <typename T>
interpolationF<T>::interpolationF(const T *x1Data, const T *x2Data, const T *x3Data, const T *x4Data, const T *fData, size_t NofElements)
{
	setData(x1Data, x2Data, x3Data, x4Data, fData, NofElements);
}

template <typename T>
void interpolationF<T>::setData(const T *x1Data, const T *x2Data, const T *x3Data, const T *x4Data, const T *fData, size_t NofElements)
{
	m_variableN = 4;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(x1Data, x1Data + m_dataLength);
	m_data[1] = std::vector<T>(x2Data, x2Data + m_dataLength);
	m_data[2] = std::vector<T>(x3Data, x3Data + m_dataLength);
	m_data[3] = std::vector<T>(x4Data, x4Data + m_dataLength);
    m_data[4] = std::vector<T>( fData,  fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 5 1D vectors:
template <typename T>
interpolationF<T>::interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &x4Data, const std::vector<T> &fData)
{
	setData(x1Data, x2Data, x3Data, x4Data, fData);
}

template <typename T>
void interpolationF<T>::setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &x4Data, const std::vector<T> &fData)
{
	m_variableN = 4;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<T>(x1Data.begin(), x1Data.begin() + m_dataLength);
	m_data[1] = std::vector<T>(x2Data.begin(), x2Data.begin() + m_dataLength);
	m_data[2] = std::vector<T>(x3Data.begin(), x3Data.begin() + m_dataLength);
	m_data[3] = std::vector<T>(x4Data.begin(), x4Data.begin() + m_dataLength);
    m_data[4] = std::vector<T>( fData.begin(),  fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//DESTRUCTORS:
template <typename T>
interpolationF<T>::~interpolationF() {}

//INTERPOLATION FUNCTIONS:
//1D interpolation
template <typename T>
T interpolationF<T>::interpolation(T pointValue) const
{
	if (m_variableN > 1) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<T>::quiet_NaN();
	}
	else {
		if (pointValue < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		return interpolation1D(pointValue);
	}
}

//2D interpolation
template <typename T>
T interpolationF<T>::interpolation(T pointValue1, T pointValue2) const
{
	if (m_variableN < 2) {
		std::cerr << "Error: too much points for interpolation." << std::endl;
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (m_variableN > 2) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<T>::quiet_NaN();
	}
	else {
		if (pointValue1 < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue1 > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue2 < m_domain[1][0]) {
			std::cerr << "Error: point value in dimension 2 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue2 > m_domain[1][1]) {
			std::cerr << "Error: point value in dimension 2 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		return interpolation2D(pointValue1, pointValue2);
	}
}

//3D interpolation
template <typename T>
T interpolationF<T>::interpolation(T pointValue1, T pointValue2, T pointValue3) const
{
	if (m_variableN < 3) {
		std::cerr << "Error: too much points for interpolation." << std::endl;
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (m_variableN > 3) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<T>::quiet_NaN();
	}
	else {
		if (pointValue1 < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue1 > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue2 < m_domain[1][0]) {
			std::cerr << "Error: point value in dimension 2 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue2 > m_domain[1][1]) {
			std::cerr << "Error: point value in dimension 2 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue3 < m_domain[2][0]) {
			std::cerr << "Error: point value in dimension 3 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue3 > m_domain[2][1]) {
			std::cerr << "Error: point value in dimension 3 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		return interpolation3D(pointValue1, pointValue2, pointValue3);
	}
	return 0.0;
}

//4D interpolation
template <typename T>
T interpolationF<T>::interpolation(T pointValue1, T pointValue2, T pointValue3, T pointValue4) const
{
	if (m_variableN < 4) {
		std::cerr << "Error: too much points for interpolation." << std::endl;
		return std::numeric_limits<T>::quiet_NaN();
	}
	else if (m_variableN > 4) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<T>::quiet_NaN();
	}
	else {
		if (pointValue1 < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue1 > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue2 < m_domain[1][0]) {
			std::cerr << "Error: point value in dimension 2 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue2 > m_domain[1][1]) {
			std::cerr << "Error: point value in dimension 2 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue3 < m_domain[2][0]) {
			std::cerr << "Error: point value in dimension 3 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue3 > m_domain[2][1]) {
			std::cerr << "Error: point value in dimension 3 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue4 < m_domain[3][0]) {
			std::cerr << "Error: point value in dimension 4 smaller than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		if (pointValue4 > m_domain[3][1]) {
			std::cerr << "Error: point value in dimension 4 larger than domain." << std::endl;
			return std::numeric_limits<T>::quiet_NaN();
		}
		return interpolation4D(pointValue1, pointValue2, pointValue3, pointValue4);
	}
	return 0.0;
}

template <typename T>
const std::vector<std::vector<T>> & interpolationF<T>::domain() const
{
	return m_domain;
}

template <typename T>
const std::vector<T> & interpolationF<T>::codomain() const
{
	return m_codomain;
}

template <typename T>
void interpolationF<T>::createGrids()
{
	for (size_t iv=0; iv<m_variableN; iv++) {
		std::sort(m_data[iv].begin(), m_data[iv].end());
		m_data[iv].erase(std::unique(m_data[iv].begin(), m_data[iv].end()), m_data[iv].end());
		m_domain.push_back({m_data[iv].front(), m_data[iv].back()});
	}
	m_codomain.push_back(*std::min_element(m_data[m_variableN].begin(), m_data[m_variableN].end()));
	m_codomain.push_back(*std::max_element(m_data[m_variableN].begin(), m_data[m_variableN].end()));
}

template <typename T>
void interpolationF<T>::locatePointF(const std::vector<T> &points, std::vector<size_t> &positions) const
{
	positions.resize(points.size(), 0);
    int ju, jm, jl, mm = 1 + 1;
	bool ascnd;
	for (size_t iv=0; iv<m_data.size()-1; iv++)
	{
		jl = 0;
		ju = m_data[iv].size() - 1;
		ascnd = (m_data[iv].back() >= m_data[iv][0]);

		while ((ju - jl) >1) {
			jm = (ju + jl) >> 1;
			if ((points[iv] >= m_data[iv][jm]) == ascnd) {
				jl = jm;
			}
			else {
				ju = jm;
			}
		}
		int n = static_cast<int>(m_data[iv].size());
		positions[iv] = static_cast<size_t>(std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1))));
	}
}

//1D linear interpolation
template <typename T>
T interpolationF<T>::lin1DInterpolation(const T x[2], const T f[2], T xx) const
{
	return (f[0] + (xx - x[0])*(f[1] - f[0]) / (x[1] - x[0]));
}

//1D interpolation (full function)
template <typename T>
T interpolationF<T>::interpolation1D(T pointValue) const
{
	//searching for position
	const std::vector<T> points{pointValue};
	std::vector<size_t> positions;
	locatePointF(points, positions);

	//setting x and Q values
	T x[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	T Q[] = {m_data[1][positions[0]], m_data[1][positions[0] + 1]};

	return lin1DInterpolation(x, Q, pointValue);
}

//2D interpolation
template <typename T>
T interpolationF<T>::interpolation2D(T pointValue1, T pointValue2) const
{
	//searching for position
	const std::vector<T> points{pointValue1, pointValue2};
	std::vector<size_t> positions;
	locatePointF(points, positions);

	T x1[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	T x2[] = {m_data[1][positions[1]], m_data[1][positions[1] + 1]};	

	T Q2[2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			Q2[i1][i2] = m_data[2][(positions[0] + i1)*m_data[1].size() + (positions[1] + i2)];

	T Q1[2];
	for (int i1=0; i1<2; i1++)
		Q1[i1] = lin1DInterpolation(x2, Q2[i1], pointValue2);

	return lin1DInterpolation(x1, Q1, pointValue1);
}

//3D interpolation
template <typename T>
T interpolationF<T>::interpolation3D(T pointValue1, T pointValue2, T pointValue3) const
{
	//searching for position
	const std::vector<T> points{pointValue1, pointValue2, pointValue3};
	std::vector<size_t> positions;
	locatePointF(points, positions);
	
	T x1[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	T x2[] = {m_data[1][positions[1]], m_data[1][positions[1] + 1]};
	T x3[] = {m_data[2][positions[2]], m_data[2][positions[2] + 1]};

	T Q3[2][2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			for (size_t i3=0; i3<2; i3++)
				Q3[i1][i2][i3] = m_data[3][(positions[0] + i1)*m_data[2].size()*m_data[1].size() + 
										   (positions[1] + i2)*m_data[2].size() +
										   (positions[2] + i3)]; 

	T Q2[2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			Q2[i1][i2] = lin1DInterpolation(x3, Q3[i1][i2], pointValue3);

	T Q1[2];
	for (size_t i1=0; i1<2; i1++)
		Q1[i1] = lin1DInterpolation(x2, Q2[i1], pointValue2);

	return lin1DInterpolation(x1, Q1, pointValue1);
}

//4D interpolation
template <typename T>
T interpolationF<T>::interpolation4D(T pointValue1, T pointValue2, T pointValue3, T pointValue4) const
{
	//searching for position
	const std::vector<T> points{pointValue1, pointValue2, pointValue3, pointValue4};
	std::vector<size_t> positions;
	locatePointF(points, positions);

	T x1[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	T x2[] = {m_data[1][positions[1]], m_data[1][positions[1] + 1]};
	T x3[] = {m_data[2][positions[2]], m_data[2][positions[2] + 1]};
	T x4[] = {m_data[3][positions[3]], m_data[3][positions[3] + 1]};

	T Q4[2][2][2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			for (size_t i3=0; i3<2; i3++)
				for (size_t i4=0; i4<2; i4++)
					Q4[i1][i2][i3][i4] = m_data[4][(positions[0] + i1)*m_data[3].size()*m_data[2].size()*m_data[1].size() +
												   (positions[1] + i2)*m_data[3].size()*m_data[2].size() +
												   (positions[2] + i3)*m_data[3].size() +
												   (positions[3] + i4)];

	T Q3[2][2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			for (size_t i3=0; i3<2; i3++)
				Q3[i1][i2][i3] = lin1DInterpolation(x4, Q4[i1][i2][i3], pointValue4);
	
	T Q2[2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			Q2[i1][i2] = lin1DInterpolation(x3, Q3[i1][i2], pointValue3);
	
	
	T Q1[2];
	for (size_t i1=0; i1<2; i1++)
		Q1[i1] = lin1DInterpolation(x2, Q2[i1], pointValue2);

	return lin1DInterpolation(x1, Q1, pointValue1);
}

template class interpolationF<float>;
template class interpolationF<double>;
template class interpolationF<long double>;