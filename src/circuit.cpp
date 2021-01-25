#include "circuit.hpp"
#include <iostream>
using namespace std::complex_literals;

/**
	\file circuit.cpp
	\brief Implementacja funkcji związanych z \ref circuit i jego elementami
	\author Jacek Wieczorek
*/

/**
	\brief Admitancja rezystancji
*/
std::complex<double> resistor::admittance(double omega) const
{
	(void) omega;
	return 1.0 / R;
}

/**
	\brief Admitancja cewki

	\note Przy analizie DC cewka jest zastępowana minimalną rezystancją
	(wielką admitancją), by nie dopuścić do dzielenia przez 0.
*/
std::complex<double> inductor::admittance(double omega) const
{
	return 1.0 / (omega == 0 ? 1e-9 : (1i * omega * L));
}

/**
	\brief Admitancja kondensatora
*/
std::complex<double> capacitor::admittance(double omega) const
{
	return 1i * omega * C;
}

/**
	\brief Tworzy solver układów
*/
circuit_solver::circuit_solver(const circuit &circ) :
	m_circuit(&circ)
{
	update_node_map();
}

/**
	\brief Zwraca mapowania węzłów (tylko do odczytu)
*/
const std::map<int, int> &circuit_solver::get_node_map() const
{
	return m_node_map;
}

/**
	\brief Zwraca pulsację dla której wyznaczone zostało rozwiązanie
*/
double circuit_solver::get_solution_omega() const
{
	return *m_solution_omega;
}


/**
	\brief Aktualizuje mapę numeracji węzłów i rozwiązanie (jeżeli było gotowe)
*/
void circuit_solver::update()
{
	update_node_map();
	
	if (m_solution.has_value())
		solve(*m_solution_omega);
}

/**
	\brief Poddaje obwód analizie dla zadanej pulsacji

	Przy analizie AC wszystkie źródła DC są pomijane i na odwrót.
	Analiza DC uruchamiana jest przez podanie omega = 0.

	\param omega pulsacja sygnału źródeł AC. 0 oznacza analizę DC.
*/
void circuit_solver::solve(double omega)
{
	// Zapisujemy omegę, dla której była przeprowadzona analiza
	m_solution_omega = omega;

	// Reset obwodu zapisanego w mna_problem
	m_problem.admittances.clear();
	m_problem.voltage_sources.clear();
	m_problem.current_sources.clear();
	m_problem.opamps.clear();
	
	// Mapowanie par węzłów
	auto map_node_pair = [&](std::pair<int, int> p)->std::pair<int, int>{
		return {m_node_map.at(p.first), m_node_map.at(p.second)};
	};

	// Dodawanie elementów do problemu do analizy
	for (const auto &[ref, comp_ptr] : *m_circuit)
	{
		// Elementy pasywne
		if (auto pcomp = std::dynamic_pointer_cast<const passive_component>(comp_ptr))
			m_problem.admittances.emplace_back(map_node_pair(pcomp->nodes), pcomp->admittance(omega));

		// Wzmacniacze operacyjne
		if (auto opa = std::dynamic_pointer_cast<const opamp>(comp_ptr))
			m_problem.opamps.emplace_back(
				m_node_map.at(opa->pos_input_node),
				m_node_map.at(opa->neg_input_node),
				m_node_map.at(opa->output_node));

		// Źródła napięciowe
		if (auto vs = std::dynamic_pointer_cast<const voltage_source>(comp_ptr))
		{
			auto x = (omega == 0) ? vs->dcV : vs->acV;
			m_problem.voltage_sources.emplace_back(map_node_pair(vs->nodes), x);
		}

		// Źródła prądowe
		if (auto cs = std::dynamic_pointer_cast<const current_source>(comp_ptr))
		{
			auto x = (omega == 0) ? cs->dcI : cs->acI;
			m_problem.current_sources.emplace_back(map_node_pair(cs->nodes), x);
		}
	}

	// Analiza
	try
	{
		m_solution = m_problem.solve();
	}
	catch (const std::runtime_error &ex)
	{
		throw std::runtime_error(std::string{"Could not compute operating point - reason: "} + ex.what());
	}
}

/**
	\brief Zwraca rozwiązanie jako mna_solution
*/
const mna::mna_solution &circuit_solver::get_solution() const
{
	return *m_solution;
}

/**
	\brief Pomiar napięcia między węzłami
	\param pos Numer mierzonego węzła
	\param neg Numer węzła odniesienia
*/
std::complex<double> circuit_solver::voltage(int pos, int neg) const
{
	return m_solution->voltage(m_node_map.at(pos), m_node_map.at(neg));
}

/**
	\brief Pomiar napięcia na komponencie
	\note Pomiar napięcia jest możliwy tylko na elementach z dwoma wyprowadzeniami
*/
std::complex<double> circuit_solver::voltage(const circuit_component &comp) const
{
	if (auto bipole = dynamic_cast<const bipole_component*>(&comp))
	{
		return voltage(bipole->nodes.first, bipole->nodes.second);
	}

	// Wzmacniacz operacyjny
	if (auto opa = dynamic_cast<const opamp*>(&comp))
	{
		return voltage(opa->output_node);
	}

	throw std::runtime_error("Cannot measure voltage on component");
}

/**
	\brief Pomiar prądu płynącego przez komponent
	\note Pomiar prądu jest możliwy tylko na elementach z dwoma wyprowadzeniami i na wzmacniaczach operacyjnych (prąd wyjścia)
*/
std::complex<double> circuit_solver::current(const circuit_component &comp) const
{
	// Komponent pasywny
	if (auto passive = dynamic_cast<const passive_component*>(&comp))
	{
		return voltage(comp) * passive->admittance(*m_solution_omega);
	}
	
	// SEM
	if (auto vs = dynamic_cast<const voltage_source*>(&comp))
	{
		// Określamy numer źródła
		int cnt = 0;
		for (const auto &[k, v] : *m_circuit)
		{
			if (auto c = dynamic_cast<const voltage_source*>(v.get()))
			{
				if (c == vs) break;
				cnt++;
			}			
		}

		return m_solution->voltage_source_current(cnt);
	}

	// SPM
	if (auto cs = dynamic_cast<const current_source*>(&comp))
	{
		return *m_solution_omega == 0 ? -cs->dcI : -cs->acI;
	}

	// Wzmacniacz operacyjny
	if (auto opa = dynamic_cast<const opamp*>(&comp))
	{
		// Określamy numer wzmacniacza
		int cnt = 0;
		for (const auto &[k, v] : *m_circuit)
		{
			if (auto c = std::dynamic_pointer_cast<const opamp>(v))
			{
				if (c.get() == opa) break;
				cnt++;
			}			
		}

		return m_solution->opamp_current(cnt);
	}

	throw std::runtime_error("Cannot measure current through component");
}

/**
	\brief Pomiar mocy traconej na komponencie
*/
std::complex<double> circuit_solver::power(const circuit_component &comp) const
{
	return voltage(comp) * current(comp);
}

/**
	\brief Pomiar spadku napięcia na komponencie
*/
std::complex<double> circuit_solver::voltage(const std::string &ref) const
{
	return voltage(*m_circuit->at(ref));
}

/**
	\brief Pomiar prądu płynącego przez komponent
*/
std::complex<double> circuit_solver::current(const std::string &ref) const
{
	return current(*m_circuit->at(ref));
}

/**
	\brief Pomiar mocy traconej na komponencie
*/
std::complex<double> circuit_solver::power(const std::string &ref) const
{
	return power(*m_circuit->at(ref));
}

/**
	\brief Aktualizuje mapowania węzłów

	Węzeł 0 to węzeł odniesienia - mapowany jest na -1 - indeks węzła masy
	w zapisie macierzowym. Wszystkie inne numery węzłów są mapowane na nieujemne
	liczby całkowite.

	Ta funkcja powinna zostać wywołana po wprowadzeniu zmian w układzie
	z którym powiązana jest klasa circuit_solver.

	\see update()
*/
void circuit_solver::update_node_map()
{
	m_node_map.clear();
	m_node_map[0] = -1;
	int cnt = 0;

	auto add_node = [&](int n)
	{
		try
		{
			m_node_map.at(n);
		}
		catch (const std::out_of_range &ex)
		{
			m_node_map[n] = cnt++;
		}
	};

	for (const auto &[ref, comp_ptr] : *m_circuit)
	{
		// Elementy z dwoma wyprowadzeniami
		if (auto comp = dynamic_cast<const bipole_component*>(comp_ptr.get()))
		{
			add_node(comp->nodes.first);
			add_node(comp->nodes.second);
		}

		// Wzmacniacze operacyjne
		if (auto opa = dynamic_cast<const opamp*>(comp_ptr.get()))
		{
			add_node(opa->pos_input_node);
			add_node(opa->neg_input_node);
			add_node(opa->output_node);
		}
	}
}
