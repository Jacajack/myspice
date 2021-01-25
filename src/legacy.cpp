#include "legacy.hpp"
#include "circuit.hpp"
#include <string>
#include <cstdio>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <optional>

/**
	\file legacy.cpp
	\brief Główny plik podstawowej wersji programu
	\author Jacek Wieczorek
*/


/**
	\brief Przekształca prostą netlistę na obwód
	\param netlist Strumień dostarczający netlistę
	\returns Obwód opisany przez netlistę

	Wczytywana netlista musi być zgodna z formatem opisanym w treści zadania.
	Najbardziej bolesne jest to, że netlista, która ma być wczytana przez program,
	nie podaje, który węzeł jest węzłem odniesienia (masą). Wybór węzła odniesienia
	jest jednak niezbędny do przeprowadzenia analizy układu. Z tego powodu zakładam,
	przesuwam numerację węzłów o jeden - w zadaniu zaczyna się od 1, a tutaj od zera.
	
	\warning Brak węzła o numerze 1 może poskutkować błędem symulacji lub otrzymaniem niepoprawnych wyników.

	\see circuit_solver
*/
static circuit read_netlist(std::istream &netlist)
{
	circuit circ;
	int R_count = 0;
	int E_count = 0;
	int I_count = 0;
	
	// Tokenizacja stringów to ból w C++ :(
	int line_count = 0;
	std::string line;
	while (line_count++, std::getline(netlist, line))
	{
		char ref;
		double value;
		std::pair<int, int> nodes;

		// Usuwanie białych znaków z początku linii
		line.erase(line.begin(), std::find_if_not(line.begin(), line.end(), [](char c){return std::isspace(c);}));
		
		// Pomijamy puste linie
		if (line.empty())
			continue;

		std::istringstream line_iss(line);

		/**
			Para węzłów jest wczytywana w odwrotnej kolejności. To przez to, że
			mój solver przyjmuje numerację węzłów w SEM i SPM zgodną ze SPICE,
			a nie z treścią zadania.
		*/
		if (line_iss >> ref >> nodes.second >> nodes.first >> value)
		{
			nodes.first--;
			nodes.second--;

			switch (ref)
			{
				case 'R':
					circ["R" + std::to_string(++R_count)] = std::make_unique<resistor>(nodes, value);
					break;

				case 'I':
					circ["I" + std::to_string(++I_count)] = std::make_unique<current_source>(nodes, value);
					break;

				case 'E':
					circ["E" + std::to_string(++E_count)] = std::make_unique<voltage_source>(nodes, value);
					break;

				default:
					throw std::runtime_error(std::string{"Niepoprawny typ elementu (linia "} + std::to_string(line_count) + ")");
			}
		}
		else
		{
			throw std::runtime_error(std::string{"Niepoprawna netlista (linia "} + std::to_string(line_count) + ")");
		}
	}

	return circ;
}

/**
	\brief Wypisuje wynik symulacji układu do podanego strumienia
	\param circ Obwód
	\param sol Solver z rozwiązaniem
	\param f Strumień wyjściowy
*/
static void print_solution(const circuit &circ, const circuit_solver &sol, std::ostream &f)
{
	// Potencjały węzłowe
	f << "Potencjaly wezlowe:" << std::endl;
	for (const auto &[k, v] : sol.get_node_map())
		f << "\tV(" << k + 1 << ") = " << sol.voltage(k).real() << " V" << std::endl;

	f << std::endl;

	// Spadki, prądy i moce
	for (const auto &[ref, comp_ptr] : circ)
	{
		const auto &comp = *comp_ptr;
		const auto &bp = dynamic_cast<const bipole_component&>(comp);
		
		f << ref << " - [" << bp.nodes.second + 1 << ", " << bp.nodes.first + 1 << "]:" << std::endl;
		
		try
		{
			auto V = sol.voltage(comp).real();
			auto I = sol.current(comp).real();
			auto P = sol.power(comp).real();

			f << "\tV(" << ref << ") = " << V << " V" << std::endl;
			f << "\tI(" << ref << ") = " << I << " A" << std::endl;
			f << "\tP(" << ref << ") = " << P << " W" << std::endl;
		}
		catch (const std::exception &ex)
		{
		}

		f << std::endl;
	}

	// Bilans mocy
	double P = 0.0;
	for (const auto &[ref, comp_ptr] : circ)
		if (auto pc = std::dynamic_pointer_cast<const passive_component>(comp_ptr))
			P += sol.power(*pc).real();

	f << "Moc calkowita: " << P << " W." << std::endl;
}

/**
	\brief Rozwiązuje układ zgodny ze specyfikacją zadania
	\param netlist Strumień dostarczający netliste
	\param output Strumień wyjściowy gdzie ma być wypisany wynik
*/
static void solve_legacy(std::istream &netlist, std::ostream &output)
{
	try
	{
		auto circ = read_netlist(netlist);
		try
		{
			circuit_solver solver(circ);
			solver.solve(0);
			print_solution(circ, solver, output);
		}
		catch (const std::exception &ex)
		{
			std::cerr << "Analiza ukladu nie powiodla sie..." << std::endl;
			return; 
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Wczytywanie netlisty nie powiodlo sie...\nPowod: " << ex.what() << std::endl;
		return;
	}
}

/**
	\brief Wypisuje informację pomocy
*/
static void help()
{
	std::cerr << 
		"Sposob uzycia: myspice NETLISTA [WYNIK]\n"
		"\t NETLISTA - plik z netlistą\n"
		"\t WYNIK - plik wynikowy (opcjonalny)\n"
		"\nAutor: Jacek Wieczorek, 2020r."
		<< std::endl;
}

/**
	\brief Główna funkcja uproszczonej wersji programu
*/
int main_legacy(int argc, char *argv[])
{
	if (argc < 2 || argc > 3)
	{
		help();
		return 0;
	}

	// Plik wej.
	std::ifstream fin(argv[1]);
	if (!fin)
	{
		std::cerr << "Nie mozna otworzyc pliku '" << argv[1] << "'!" << std::endl;
		return 1;
	}

	// Plik wyj.
	std::optional<std::ofstream> fout;
	if (argc == 3)
	{
		fout.emplace(argv[2]);
		if (!*fout)
		{
			std::cerr << "Nie mozna otworzyc pliku '" << argv[2] << "'!" << std::endl;
			return 1;
		}
	}

	solve_legacy(fin, fout.has_value() ? fout.value() : std::cout);
	return 0;
}
