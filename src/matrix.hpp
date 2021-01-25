#pragma once
#include <vector>
#include <stdexcept>
#include <iosfwd>
#include <iomanip>

/**
	\file matrix.hpp 
	\brief Definiuje klasę \ref matrix do operacji macierzowych
	\author Jacek Wieczorek
*/

/**
	\brief Klasa ułatwiająca operacje macierzowe
*/
template <typename T>
class matrix
{
public:
	matrix() :
		m_w(0),
		m_h(0)
	{}

	/**
		\brief Inicjalizuje macierz o określonych rozmiarach
		\param h wysokość
		\param w szerokość
	*/
	matrix(int h, int w) :
		m_matrix(w * h),
		m_w(w),
		m_h(h)
	{}

	~matrix() = default;

	matrix(const matrix<T> &) = default;
	matrix(matrix<T> &&) noexcept = default;
	matrix<T> &operator=(const matrix<T> &) = default;
	matrix<T> &operator=(matrix<T> &&) noexcept = default;
	
	/**
		\brief Zwraca szerokość macierzy
	*/
	int get_width() const
	{
		return m_w;
	}

	/**
		\brief Zwraca wysokość macierzy
	*/
	int get_height() const
	{
		return m_h;
	}

	/**
		\brief Dostęp do danych w podlegającym macierzy std::vector
	*/
	T *data()
	{
		return m_matrix.data();
	}

	/**
		\brief Dostęp do danych w macierzy
	*/
	T &operator()(int row, int col)
	{
		if (row < 0 || row >= get_height() || col < 0 || col >= get_width())
			throw std::out_of_range("access outside of matrix");

		return m_matrix[row * m_w + col];
	}

	/**
		\brief Dostęp do danych w macierzy
	*/
	const T &operator()(int row, int col) const
	{
		if (row < 0 || row >= get_height() || col < 0 || col >= get_width())
			throw std::out_of_range("access outside of matrix");

		return m_matrix[row * m_w + col];
	}

	/**
		\brief Dostęp do danych w macierzy
	*/
	const T &at(int row, int col) const
	{
		if (row < 0 || row >= m_h || col < 0 || col >= m_w)
			throw std::out_of_range("access outside of matrix");
		
		return this->operator()(row, col);
	}

	/**
		\brief Dostęp do danych w macierzy
	*/
	T &at(int row, int col) 
	{
		return const_cast<T&>(const_cast<const matrix<T>&>(*this).at(row, col));
	}

	/**
		Nadpisuje fragment macierzy inną mniejszą macierzą

		\param row początkowy wiersz
		\param col początkowa kolumna
		\param mat macierz do wpisania
		\throw std::out_of_range jeśli operacja wymagałaby wykroczenia poza macierz
	*/
	void replace(int row, int col, const matrix<T> &mat)
	{
		if (row < 0 
			|| col < 0
			|| row + mat.get_height() - 1 >= get_height()
			|| col + mat.get_width() - 1 >= get_width())
			throw std::out_of_range("matrix<T>::replace() - out of range");

		for (int y = 0; y < mat.get_height(); y++)
			for (int x = 0; x < mat.get_width(); x++)
				at(y + row, x + col) = mat.at(y, x);
	}

	/**
		\brief Zwraca transpozycję macierzy
	*/
	matrix<T> transpose() const
	{
		matrix<T> mat(get_width(), get_height());

		for (int y = 0; y < get_height(); y++)
			for (int x = 0; x < get_width(); x++)
				mat(x, y) = this->at(y, x);
		
		return mat;
	}

	/**
		\brief Możenie macierzy przez skalar
	*/
	matrix<T> &operator*(const T &scalar)
	{
		for (int y = 0; y < get_height(); y++)
			for (int x = 0; x < get_width(); x++)
				at(y, x) *= scalar;
		return *this;
	}

	/**
		\brief Dodawanie skalara do macierzy
	*/
	matrix<T> &operator+(const T &scalar)
	{
		for (int y = 0; y < get_height(); y++)
			for (int x = 0; x < get_width(); x++)
				at(y, x) += scalar;
		return *this;
	}

private:
	std::vector<T> m_matrix; //!< Dane macierzy
	int m_w, m_h;
};

/**
	\brief Złącza macierze w poziomie
*/
template <typename T>
matrix<T> join_matrices_horizontal(const matrix<T> &l, const matrix<T> &r)
{
	if (l.get_height() != r.get_height())
		throw std::runtime_error("cannot horizontally join matrices of different heights");

	matrix<T> mat(l.get_height(), l.get_width() + r.get_width());
	mat.replace(0, 0, l);
	mat.replace(0, l.get_width(), r);
	return mat;
}

/**
	\brief Złącza macierze w pionie
*/
template <typename T>
matrix<T> join_matrices_vertical(const matrix<T> &u, const matrix<T> &d)
{
	if (u.get_width() != d.get_width())
		throw std::runtime_error("cannot vertically join matrices of different widths");

	matrix<T> mat(u.get_height() + d.get_height(), u.get_width());
	mat.replace(0, 0, u);
	mat.replace(u.get_height(), 0, d);
	return mat;
}

/**
	\brief Operator mnożenia macierzowego
*/
template <typename T, typename U, typename V = std::common_type_t<T, U>>
matrix<V> operator*(const matrix<T> &lhs, const matrix<U> &rhs)
{
	if (lhs.get_width() != rhs.get_height())
		throw std::runtime_error("invalid matrix dimensions in multiplication");

	matrix<V> res(lhs.get_height(), rhs.get_width());

	for (int i = 0; i < lhs.get_height(); i++)
	{
		for (int j = 0; j < rhs.get_width(); j++)
		{
			res(i, j) = 0;
			for (int k = 0; k < lhs.get_width(); k++)
				res(i, j) += lhs(i, k) * rhs(k, j);
		}
	}

	return res;
}

/**
	\brief Operator wypisania dla macierzy
*/
template <typename T>
std::ostream &operator<<(std::ostream &s, const matrix<T> &mat)
{
	for (int y = 0; y < mat.get_height(); y++)
	{
		for (int x = 0; x < mat.get_width(); x++)
		{
			s << std::setw(6) << mat.at(y, x) << " ";
		}

		s << std::endl;
	}

	return s;
}
