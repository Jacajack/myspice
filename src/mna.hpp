#pragma once

#include <utility>
#include <complex>
#include <vector>

#include "matrix.hpp"

/**
	\file mna.hpp
	\brief Definicje typów solvera MNA
	\author Jacek Wieczorek
*/

/**
	\brief Solver układów metodą MNA
*/
namespace mna {

/**
	\brief Admitancja międzywęzłowa.
	
	Każdy element pasywny jest do takiej uogólniany.
*/
struct admittance
{
	std::pair<int, int> nodes;
	std::complex<double> Y;
};

/**
	\brief Idealna siła elektromotoryczna

	\note Przyjmujemy, że pierwszy węzeł to "plus"
*/
struct voltage_source
{
	std::pair<int, int> nodes;
	double V;
};

/**
	\brief Idealna siła prądomotoryczna

	\note Przyjmujemy, że pierwszy węzeł to "plus"
*/
struct current_source
{
	std::pair<int, int> nodes;
	double I;
};

/**
	\brief Idealny wzmacniacz operacyjny

	\warning Zakładamy, że wzmacniacz pracuje z ujemnym sprzężeniem zwrotnym. 
	Analiza opiera się na założeniu, że napięcie na wejściach jest równe.
	Dodatkowym skutkiem ubocznym jest fakt, że podłączenie wzmacniacza "na odwrót" niczego nie zmienia.
*/
struct opamp
{
	int pos_input_node;
	int neg_input_node;
	int output_node;
};

/**
	\brief Wynik analizy układu metodą MNA

	Wynik analizy układu. Dostarcza informacje o potencjałach węzłowych i prądach
	płynących przez siły elektromotoryczne.
*/
class mna_solution
{
public:
	mna_solution(const matrix<std::complex<double>> &solution, int node_count, int vs_count);

	std::complex<double> voltage(int pos, int neg = -1) const;
	std::complex<double> voltage_source_current(int id) const;
	std::complex<double> opamp_current(int id) const;

	const matrix<std::complex<double>> get_matrix() const;

private:
	matrix<std::complex<double>> m_solution;
	int m_node_count;
	int m_voltage_source_count;
};


/**
	\brief Układ do rozwiązania metodą MNA

	Układ jest zdegenerowany do listy admitancji międzywęzłowych, sił napięciowych
	i prądowych.

	\note Węzły poniżej 0 to napięcie odniesienia (masa). 

	\note Numeracja węzłów wg. macierzy - użycie wysokich liczb w tej strukturze
	poskutkuje obliczeniami na dużej macierzy.
*/
struct mna_problem
{
	std::vector<admittance> admittances;
	std::vector<voltage_source> voltage_sources;
	std::vector<current_source> current_sources;
	std::vector<opamp> opamps;

	mna_solution solve() const;

private:
	int get_max_node() const;
	matrix<std::complex<double>> compute_matrix_A(int node_count) const;
	matrix<std::complex<double>> compute_matrix_z(int node_count) const;
};

}