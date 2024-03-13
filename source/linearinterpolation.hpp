#ifndef HEADERFILE_LINEARINTERPOLATION
#define HEADERFILE_LINEARINTERPOLATION

#include <vector>

template<typename T>
class interpolationF {
public:
	//CONSTRUCTORS:
	interpolationF();

	//input is 2 1D arrays:
	interpolationF(const T *xData, const T *fData, size_t NofElements);
	void setData(const T *xData, const T *fData, size_t NofElements);

	//input is 2 1D vectors:
	interpolationF(const std::vector<T> &xData, const std::vector<T> &fData);
	void setData(const std::vector<T> &xData, const std::vector<T> &fData);

	//input is 3 1D arrays:
	interpolationF(const T *x1Data, const T *x2Data, const T *fData, size_t NofElements);
	void setData(const T *x1Data, const T *x2Data, const T *fData, size_t NofElements);

	//input is 3 1D vectors:
	interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &fData);

	//input is 2 1D vectors (grids) and 1 2d vector (function values):
	interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<std::vector<T>> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<std::vector<T>> &fData);

	//input is 4 1D arrays:
	interpolationF(const T *x1Data, const T *x2Data, const T *x3Data, const T *fData, size_t NofElements);
	void setData(const T *x1Data, const T *x2Data, const T *x3Data, const T *fData, size_t NofElements);

	//input is 4 1D vectors:
	interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &fData);

	//input is 5 1D arrays:
	interpolationF(const T *x1Data, const T *x2Data, const T *x3Data, const T *x4Data, const T *fData, size_t NofElements);
	void setData(const T *x1Data, const T *x2Data, const T *x3Data, const T *x4Data, const T *fData, size_t NofElements);

	//input is 5 1D vectors:
	interpolationF(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &x4Data, const std::vector<T> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &x4Data, const std::vector<T> &fData);

	//DESTRUCTOR:
	~interpolationF();

	//INTERPOLATION FUNCTIONS:
	//1D interpolation
	T interpolation(T pointValue) const;

	//2D interpolation
	T interpolation(T pointValue1, T pointValue2) const;

	//3D interpolation
	T interpolation(T pointValue1, T pointValue2, T pointValue3) const;

	//4D interpolation
	T interpolation(T pointValue1, T pointValue2, T pointValue3, T pointValue4) const;

	//miscellaneous FUNCTIONS:
	//function that returns domains:
	const std::vector<std::vector<T>> & domain() const;

	//function that returns codomain:
	const std::vector<T> & codomain() const;

private:
	size_t m_dataLength;
	std::vector<std::vector<T>> m_data;
	size_t m_variableN;
	std::vector<size_t> m_gridLengths;
	std::vector<size_t> m_relPosition;
	std::vector<std::vector<T>> m_domain;
	std::vector<T> m_codomain;

	void createGrids();

	//function that locates points
	void locatePointF(const std::vector<T> &points, std::vector<size_t> &positions) const;

	T lin1DInterpolation(const T x[2], const T f[2], T xx) const;

	//1D interpolation (full function)
	T interpolation1D(T pointValue) const;

	//2D interpolation
	T interpolation2D(T pt1, T pt2) const;

	//3D interpolation
	T interpolation3D(T pt1, T pt2, T pt3) const;

	//4D interpolation
	T interpolation4D(T pt1, T pt2, T pt3, T pt4) const;
};

#endif