#include "mna.hpp"
#include <iostream>
using namespace mna;

/**
	\file mna.cpp
	\brief Implementacja algorytmu MNA z eliminacją Gaussa

	Ten plik implementuje narzędzia do analizy układów na "najniższym poziomie".
	Analizowany układ musi zostać uprzednio zdegenerowany do opisu problemu, na
	który składają się admitancje międzywęzłowe, źródła napięciowe i prądowe
	oraz idealne wzmacniacze operacyjne (pracujące z ujemnym sprzężeniem zwrotnym).

	Na tym poziomie obowiązuje numeracja węzłów od 0. Węzły o numerach ujemnych
	traktowane są jako napięcie odniesienia (masa).

	\see https://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA3.html
	\author Jacek Wieczorek
*/

/**
	\brief Rozwiązuje układ zespolonych równań liniowych metodą eliminacji Gaussa

	\param mat Macierz o rozmiarze Nx(N+1) opisująca układ równań
	\returns Macierz o rozmiarze Nx1 zawierająca rozwiązanie
*/
static matrix<std::complex<double>> gaussian_elimination(matrix<std::complex<double>> mat)
{
	const int N = mat.get_height();

	if (mat.get_width() != N + 1)
		throw std::runtime_error("Invalid equation system dimensions");

	// Funkcja pomocnicza do zamiany wierszy
	auto swap_rows = [&](int a, int b){
		for (int i = 0; i < mat.get_width(); i++)
			std::swap(mat(a, i), mat(b, i));
	};

	// Funkcja pomocnicza do sumowania wierszy
	auto add_rows = [&](int dest, int src){
		for (int i = 0; i < mat.get_width(); i++)
			mat(dest, i) += mat(src, i);
	};

	// Funkcja pomocnicza do mnożenia wiersza przez stałą
	auto multiply_row = [&](int r, std::complex<double> k){
		for (int i = 0; i < mat.get_width(); i++)
			mat(r, i) *= k;
	};

	// Redukcja macierzy do macierzy trójkątnej (row echelon form)
	for (int k = 0; k < N; k++)
	{
		// Wyszukanie wiersza z "największą" wartością w k-tej kolumnie
		int row_max = k;
		auto max = std::abs(mat(k, k));
		for (int i = k + 1; i < N; i++)
		{
			auto x = std::abs(mat(i, k));
			if (x > max)
			{
				max = x;
				row_max = i;
			}
		}

		// Żadne równanie nie korzysta z tej zmiennej - brak jednoznacznego rozwiązania
		if (max == 0.0)
			throw std::runtime_error("Could not solve equation system (Gaussian elimination failed)");

		// Przerzucenie "maksymalnego" wiersza na górę (zamiast k-tego wiersza)
		swap_rows(row_max, k);

		// Wartość "maksymalnego" współczynnika
		auto v = mat(k, k);

		// Redukujemy współczynniki przy tej zmiennej do 0 we wszystkich
		// równaniach poniżej
		for (int i = k + 1; i < N; i++)
			if (mat(i, k) != 0.0)
			{
				multiply_row(i, -v / mat(i, k));
				add_rows(i, k);
			}
	}

	// Macierz zawierająca rozwiązanie
	matrix<std::complex<double>> solution(N, 1);

	// Podstawianie wsteczne (back substitution) celem
	// wyznaczenia wartości zmiennych
	for (int i = N - 1; i >= 0; i--)
	{
		auto sum = mat(i, N); // Wyraz wolny i-tego równania

		// Podstawienie reszty zmiennych
		for (int j = i; j < N; j++)
			sum -= mat(i, j) * solution(j, 0);
	
		// Wartość i-tej zmiennej
		solution(i, 0) = sum / mat(i, i); 
	}

	return solution;
}

/**
	\brief Wyznacza maksymalny numer węzła występujący w rozwiązywanym układzie
*/
int mna_problem::get_max_node() const
{
	int max_node = -1;

	for (const auto &v : admittances)
	{
		max_node = std::max(max_node, v.nodes.first);
		max_node = std::max(max_node, v.nodes.second);
	}

	for (const auto &v : voltage_sources)
	{
		max_node = std::max(max_node, v.nodes.first);
		max_node = std::max(max_node, v.nodes.second);
	}

	for (const auto &v : current_sources)
	{
		max_node = std::max(max_node, v.nodes.first);
		max_node = std::max(max_node, v.nodes.second);
	}

	for (const auto &v : opamps)
	{
		max_node = std::max(max_node, v.pos_input_node);
		max_node = std::max(max_node, v.neg_input_node);
		max_node = std::max(max_node, v.output_node);
	}

	return max_node;
}

/**
	\brief Wyznacza i zwraca rozwiązanie (potencjały węzłowe) układu.
	
	Implementuje macierzowy algorytm opisany [tutaj](https://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA3.html).

	\see https://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA3.html
*/
mna_solution mna_problem::solve() const
{
	// Wyznaczanie maksymalnego numeru węzła
	const int node_count = get_max_node() + 1;
	auto A = compute_matrix_A(node_count);
	auto z = compute_matrix_z(node_count);
	auto system = join_matrices_horizontal(A, z);
	auto x = gaussian_elimination(system);

	// DEBUG
	// std::cout << system << std::endl;
	// std::cout << x << std::endl;

	return mna_solution(x, node_count, voltage_sources.size());
}

/**
	\brief Tworzy i zwraca macierz A potrzebną do wyznaczenia rozwiązania.
*/
matrix<std::complex<double>> mna_problem::compute_matrix_A(int node_count) const
{
	// Liczba węzłów (bez masy) i sił elektromotorycznych
	const auto n = node_count;
	const auto m = opamps.size() + voltage_sources.size();

	// Macierze składowe macierzy A
	matrix<std::complex<double>> G(n, n);
	matrix<std::complex<double>> B(n, m);
	matrix<std::complex<double>> C(m, n);
	matrix<std::complex<double>> D(m, m);

	// Budowa macierzy G na podstawie admitancji międzywęzłowych
	for (const auto &elem : admittances)
	{
		// Każdy element na przekątnej jest sumą admitancji elementów
		// podpiętych do danego węzła
		if (elem.nodes.first >= 0) G(elem.nodes.first, elem.nodes.first) += elem.Y;
		if (elem.nodes.second >= 0) G(elem.nodes.second, elem.nodes.second) += elem.Y;

		// Ujemne konduktancje na elementach odpowiadających parom łączonych przez
		// element węzłów
		if (elem.nodes.first >= 0 && elem.nodes.second >= 0)
		{
			G(elem.nodes.first, elem.nodes.second) -= elem.Y;
			G(elem.nodes.second, elem.nodes.first) -= elem.Y;
		}
	}

	// Budowa macierzy B na podstawie źródeł napięciowych
	for (unsigned int i = 0; i < voltage_sources.size(); i++)
	{
		const auto &vs = voltage_sources[i];
		if (vs.nodes.first >= 0) B(vs.nodes.first, i) = 1;
		if (vs.nodes.second >= 0) B(vs.nodes.second, i) = -1;
	}

	C = B.transpose();

	// Uwzględnienie źródeł napięciowych we wzmacniaczach operacyjnych w macierzy B
	// Źródło napięciowe wzmacniacza jest podłączone do masy i do wyjścia
	int i = voltage_sources.size();
	for (const auto &opa : opamps)
	{
		if (opa.output_node >= 0) B(opa.output_node, i) = 1;
		i++;
	}

	// C jest transpozycją B uwzględniającego tylko niezależne źródła napięciowe
	// z później uwzględnionymi wejściami wzmacniaczy. Wyjścia wzmacniaczy *nie*
	// są uwzględnione w C
	i = voltage_sources.size();
	for (const auto &opa : opamps)
	{
		if (opa.pos_input_node >= 0) C(i, opa.pos_input_node) = 1;
		if (opa.neg_input_node >= 0) C(i, opa.neg_input_node) = -1;
		i++;
	}	

	return join_matrices_vertical(join_matrices_horizontal(G, B), join_matrices_horizontal(C, D));
}

/**
	\brief Buduje i zwraca macierz (wektor) z potrzebny do wyznaczenia rozwiązania
*/
matrix<std::complex<double>> mna_problem::compute_matrix_z(int node_count) const
{
	// Liczba węzłów (bez masy) i sił elektromotorycznych
	const auto n = node_count;
	const auto m = opamps.size() + voltage_sources.size();

	matrix<std::complex<double>> I(n, 1);
	matrix<std::complex<double>> E(m, 1);

	// Uwzględnienie źródeł prądowych
	for (const auto &cs : current_sources)
	{
		if (cs.nodes.first >= 0) I(cs.nodes.first, 0) += cs.I;
		if (cs.nodes.second >= 0) I(cs.nodes.second, 0) += -cs.I;
	}

	// Uwzględnienie źródeł napięciowych
	for (unsigned int i = 0; i < voltage_sources.size(); i++)
		E(i, 0) = voltage_sources[i].V;

	return join_matrices_vertical(I, E);
}

/**
	\brief Tworzy klasę zawierającą rozwiązanie na podstawie wektora napięć i prądów płynących
	przez siły elektromotoryczne.

	\param solution Macierz zawierająca rozwiązanie
	\param node_count Liczba węzłów w układzie
	\param vs_count Liczba SEM w układzie (nie licząc wzmacniaczy operacyjnych)
*/
mna_solution::mna_solution(const matrix<std::complex<double>> &solution, int node_count, int vs_count) : 
	m_solution(solution),
	m_node_count(node_count),
	m_voltage_source_count(vs_count)
{
}

/**
	\brief Zwraca napięcie między dwoma węzłami

	\note Węzły poniżej 0 traktowane są jako masa.

	\param pos numer mierzonego
	\param neg numer węzła odniesienia (domyślnie -1)
*/
std::complex<double> mna_solution::voltage(int pos, int neg) const
{
	if (pos >= m_node_count || neg >= m_node_count)
		throw std::out_of_range("mna_solution::voltage() invalid node ID");

	auto vpos = pos < 0 ? 0.0 : m_solution(pos, 0);
	auto vneg = neg < 0 ? 0.0 : m_solution(neg, 0);
	return vpos - vneg;
}

/**
	\brief Zwraca prąd pobierany ze źródła napięciowego

	\param id numer źródła napięciowego (numeracja od 0)
*/
std::complex<double> mna_solution::voltage_source_current(int id) const
{
	if (id < 0 || id >= m_voltage_source_count)
		throw std::out_of_range("mna_solution::voltage_source_current() invalid source ID");

	return m_solution.at(m_node_count + id, 0);
}

/**
	\brief Zwraca prąd pobierany z wyjścia wzmacniacza operacyjnego

	\param id numer wzmacniacza
*/
std::complex<double> mna_solution::opamp_current(int id) const
{
	if (id < 0 || id >= m_solution.get_height() - m_node_count - m_voltage_source_count)
		throw std::out_of_range("mna_solution::opamp_current() invalid source ID");

	return m_solution.at(m_node_count + m_voltage_source_count + id, 0);
}

/**
	\brief Zwraca macierz zawierająca rozwiązanie
*/
const matrix<std::complex<double>> mna_solution::get_matrix() const
{
	return m_solution;
}