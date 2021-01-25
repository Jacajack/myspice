#pragma once
#include <utility>
#include <map>
#include <memory>
#include <complex>
#include <optional>
#include "mna.hpp"

/**
	\file circuit.hpp
	\brief Definicje typów \ref circuit, \ref circuit_solver i elementów obwodu
	\author Jacek Wieczorek
*/

/**
	\brief Klasa bazowa dla wszystkich komponentów, które mogą się znaleźć w obwodzie
*/
struct circuit_component
{
	virtual ~circuit_component() {}
};

/**
	\brief Klasa bazowa dla elementów z dwoma wyprowadzeniami
*/
struct bipole_component : public circuit_component
{
	bipole_component(const std::pair<int, int> &p) :
		nodes(p)
	{}

	//! Para węzłów połączonych przez komponent
	std::pair<int, int> nodes;
};

/**
	\brief Klasa bazowa dla komponentów posiadających admitancję zależną od częstotliwości
*/
struct passive_component : public bipole_component
{
	using bipole_component::bipole_component;
	virtual std::complex<double> admittance(double omega) const = 0;
};

/**
	\brief Idealna siła elektromotoryczna
*/
struct voltage_source : public bipole_component
{
	voltage_source(const std::pair<int, int> &p, double v, double ac = 0.0) :
		bipole_component(p),
		dcV(v),
		acV(ac)
	{}

	//! Wartość napięcia [V] dla analizy DC
	double dcV;

	//! Wartość napięcie [V] dla analizy AC
	double acV;
};

/**
	\brief Idealna siła prądomotoryczna
*/
struct current_source : public bipole_component
{
	current_source(const std::pair<int, int> &p, double i, double ac = 0.0) :
		bipole_component(p),
		dcI(i),
		acI(ac)
	{}

	//! Wartość prądu [A] dla analizy DC
	double dcI;

	//! Wartość prądu [A] dla analizy AC
	double acI;
};

/**
	\brief Idealny rezystor
*/
struct resistor : public passive_component
{
	resistor(const std::pair<int, int> &p, double r) :
		passive_component(p),
		R(r)
	{}

	std::complex<double> admittance(double omega) const override;

	//! Rezystancja [Ohm]
	double R;
};

/**
	\brief Prawie idealna cewka

	\note Dla omega = 0 przyjmowana jest minimalna rezystancja,
	by nie dopuścić do dzielenia przez 0
*/
struct inductor : public passive_component
{
	inductor(const std::pair<int, int> &p, double l) :
		passive_component(p),
		L(l)
	{}

	std::complex<double> admittance(double omega) const override;

	//! Indukcyjność [H]
	double L;
};

/**
	\brief Idealny kondensator
*/
struct capacitor : public passive_component
{
	capacitor(const std::pair<int, int> &p, double c) :
		passive_component(p),
		C(c)
	{}

	std::complex<double> admittance(double omega) const override;

	//! Pojemność [F]
	double C;
};

/**
	\brief Idealny wzmacniacz operacyjny
*/
struct opamp : public circuit_component
{
	opamp(int pos, int neg, int out) : 
		pos_input_node(pos),
		neg_input_node(neg),
		output_node(out)
	{}

	int pos_input_node; //! Numer węzła wejścia nieodwracającego
	int neg_input_node; //! Numer węzła wejścia odwracającego
	int output_node; //! Numer węzła wyjściowego
};

/**
	\brief Obwód elektryczny - zbiór nazwanych elementów
*/
using circuit = std::map<std::string, std::shared_ptr<circuit_component>>;

/**
	\brief Analizator układów liniowych

	Analizator jest tworzony na podstawie układu (\ref circuit). Pozwala na przeprowadzenie
	analizy metodą MNA wykorzystując \ref mna::mna_problem. Głównym zadaniem tej klasy
	jest wprowadzenie dodatkowej abstrakcji - pozwala ona na tworzenie obwodów poprzez
	proste dodawanie różnych elementów do powiązanej klasy \ref circuit i rozluźnia
	wymagania dot. numeracji węzłów.

	Po wywołaniu \ref solve(), możliwy jest pomiar napięć, prądów i mocy
	w układzie za pomocą \ref voltage(), \ref current() i \ref power().
*/
class circuit_solver
{
public:
	explicit circuit_solver(const circuit &circ);

	void update();	
	void solve(double omega);

	const mna::mna_solution &get_solution() const;
	const std::map<int, int> &get_node_map() const;
	double get_solution_omega() const;

	std::complex<double> voltage(int pos, int neg = 0) const;
	std::complex<double> voltage(const circuit_component &comp) const;
	std::complex<double> current(const circuit_component &comp) const;
	std::complex<double> power(const circuit_component &comp) const;

	std::complex<double> voltage(const std::string &ref) const;
	std::complex<double> current(const std::string &ref) const;
	std::complex<double> power(const std::string &ref) const;

private:
	void update_node_map();

	//! Analizowany obwód
	const circuit *m_circuit;

	//! Mapowanie numerów węzłów do bardziej restrykcyjnej numeracji mna::mna_problem
	std::map<int, int> m_node_map;

	//! Obwód w wersji mna_problem
	mna::mna_problem m_problem;

	//! Rozwiązanie (może być nieobecne)
	std::optional<mna::mna_solution> m_solution;

	//! Pulsacja dla której układ był ostatnio analizowany
	std::optional<double> m_solution_omega;
};