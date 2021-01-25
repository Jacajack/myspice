#include "extended.hpp"
#include <iterator>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <memory>
#include <numeric>
#include <optional>
#include <cmath>
#include <regex>
#include "circuit.hpp"

using namespace std::string_literals;

/**
	\file extended.cpp
	\brief Rozszerzona część programu
	\author Jacek Wieczorek
*/

/**
	\brief Parametry do analizy AC
*/
struct ac_analysis_params
{
	double start;    //!< Dolna częstotliwość [Hz]
	double stop;     //!< Górna częstotliwość [Hz]
	double exponent; //!< Wykładnik sweep'a. 0 to sweep liniowy.
	int steps;       //!< Liczba punktów na exponent-krotną zmianę częstotliwości/łącznie
};

/**
	\brief Metoda pomiaru zespolonych wielkości fizycznych
*/
enum class complex_probing_method
{
	DEFAULT,
	MAGNITUDE,
	PHASE,
	REAL,
	IMAGINARY
};

/**
	\brief Zwraca określony komponent zespolonej wielkości fizycznej

	W trybie pomiaru DEFAULT zwraca część rzeczywistą dla pulsacji równej 0,
	a w trybie AC zwraca moduł.
*/
double probe_complex(std::complex<double> c, complex_probing_method method, double omega = 0)
{
	switch (method)
	{
		case complex_probing_method::DEFAULT:
			if (omega == 0.0)
				return c.real();
			else
				return std::abs(c);
			break;

		case complex_probing_method::MAGNITUDE:
			return std::abs(c);
			break;

		case complex_probing_method::PHASE:
			return std::arg(c);
			break;

		case complex_probing_method::REAL:
			return c.real();
			break;

		case complex_probing_method::IMAGINARY:
			return c.imag();
			break;

		default:
			return 0.0;
			break;
	}
}

/**
	\brief Zwraca suffix dodawany do oznaczenia mierzonej wartości zespolonej
*/
std::string probing_method_suffix(complex_probing_method method)
{
	switch (method)
	{
		case complex_probing_method::DEFAULT:
			return "";
			break;

		case complex_probing_method::MAGNITUDE:
			return "mag";
			break;

		case complex_probing_method::PHASE:
			return "p";
			break;

		case complex_probing_method::REAL:
			return "re";
			break;

		case complex_probing_method::IMAGINARY:
			return "im";
			break; 

		default:
			return "";
	}
}

/**
	\brief Uogólnienie mierników
*/
class probe
{
public:
	std::string get_name() const
	{
		return m_name;
	}

	virtual double get_value(const circuit_solver &solver) const = 0;

protected:
	std::string m_name;
};

/**
	\brief Miernik napięcia międzywęzłowego/na elemencie
*/
class voltage_probe : public probe
{
public:
	voltage_probe(const circuit &circ, const std::string &ref, complex_probing_method pm) :
		m_probing_method(pm)
	{
		try
		{
			auto bp = dynamic_cast<const bipole_component&>(*circ.at(ref));
			m_nodes = bp.nodes;
		}
		catch (const std::bad_cast &ex)
		{
			throw std::runtime_error("Cannot probe voltage on non-bipole component");
		}

		m_name = std::string{"V"} + probing_method_suffix(pm) + "(" + ref + ")";
	}

	voltage_probe(int pos, int neg, complex_probing_method pm) :
		m_nodes(pos, neg),
		m_probing_method(pm)
	{
		if (neg != 0)
		{
			m_name = std::string{"V"} + probing_method_suffix(pm) + "(" 
				+ std::to_string(m_nodes.first) + ", " + std::to_string(m_nodes.second) + ")";
		}
		else
			m_name = std::string{"V"} + probing_method_suffix(pm) + "(" + std::to_string(m_nodes.first) + ")";
	}

	double get_value(const circuit_solver &solver) const override
	{
		try
		{
			return probe_complex(solver.voltage(m_nodes.first, m_nodes.second), m_probing_method, solver.get_solution_omega());
		}
		catch (const std::exception &ex)
		{
			throw std::runtime_error("Probing '"s + m_name + "' failed");
		}
	}

private:
	std::pair<int, int> m_nodes;
	complex_probing_method m_probing_method;
};

/**
	\brief Miernik prądu płynącego przez element
*/
class current_probe : public probe
{
public:
	current_probe(const circuit &circ, const std::string &ref, complex_probing_method pm) :
		m_ref(ref),
		m_probing_method(pm)
	{
		m_name = std::string{"I"} + probing_method_suffix(pm) + "(" + ref + ")";
	}

	double get_value(const circuit_solver &solver) const override
	{
		try
		{
			return probe_complex(solver.current(m_ref), m_probing_method, solver.get_solution_omega());
		}
		catch (const std::exception &ex)
		{
			throw std::runtime_error("Probing '"s + m_name + "' failed");
		}
	}

protected:
	std::string m_ref;
	complex_probing_method m_probing_method;
};

/**
	\brief Miernik mocy wydzielanej na elemencie
*/
class power_probe : public current_probe
{
public:
	power_probe(const circuit &circ, const std::string &ref, complex_probing_method pm) :
		current_probe(circ, ref, pm)
	{
		m_name = std::string{"P"} + probing_method_suffix(pm) + "(" + ref + ")";
	}

	double get_value(const circuit_solver &solver) const override
	{
		try
		{
			return probe_complex(solver.power(m_ref), m_probing_method, solver.get_solution_omega());
		}
		catch (const std::exception &ex)
		{
			throw std::runtime_error("Probing '"s + m_name + "' failed");
		}
	}
};

/**
	\brief Symulacja układu
*/
struct circuit_simulation
{
	std::string title;
	circuit circ;
	std::optional<ac_analysis_params> ac;
	std::vector<std::shared_ptr<probe>> probes;
};

/**
	\brief Zamienia litery w stringu na małe
*/
static std::string tolower(std::string s)
{
	for (char &c : s)
		c = std::tolower(c);
	return s;
}

/**
	\brief Zwraca wektor fragmentów tekstu rozdzielonych znakami białymi
*/
static std::vector<std::string> tokenize_string(const std::string &s, std::function<bool(char c)> predicate)
{
	std::vector<std::string> tokens;

	auto bit = s.begin();
	auto eit = bit;

	while (bit != s.end() && eit != s.end())
	{
		eit = std::find_if(bit, s.end(), predicate);
		if (bit != eit)
			tokens.emplace_back(bit, eit);		
		bit = std::next(eit);
	}

	return tokens;
}

/**
	\brief Zamienia tekst będący liczbą z przedrostkiem SI na wartość
*/
static double si_string_to_double(const std::string &s)
{
	const static std::map<std::string, double> muls = {
		{"p", 1e-12},
		{"n", 1e-9},
		{"u", 1e-6},
		{"m", 1e-3},
		{"k", 1e3},
		{"Meg", 1e6},
		{"G", 1e9},
	};

	std::istringstream ss(s);
	std::string prefix;
	double val;
	ss >> val >> prefix;

	if (!prefix.empty())
	{
		try
		{
			return val * muls.at(prefix);
		}
		catch (const std::out_of_range &ex)
		{
			throw std::runtime_error("Invalid SI prefix");
		}
	}
	else
		return val;
}


/**
	\brief Tworzy komponent na podstawie "stokenizowanej" linii pliku SPICE

	\note Akceptuje nazwy węzłów typu '2z'. Nie jest to piękne, ale też na razie nie ma 
	potrzeby żeby to na siłę naprawiać.
*/
static std::shared_ptr<circuit_component> create_component(const std::vector<std::string> &tokens)
{
	auto ref = tokens.at(0);
	std::string ref_type{ref.begin(), std::find_if(ref.begin(), ref.end(), [](char c){return std::isdigit(c);})};

	// Interpretuje trójkę wartości, która opisuje wszystkie komponenty "bipolowe"
	auto parse_component_triplet = [](const std::vector<std::string> &tokens, std::pair<int, int> &nodes, double &value){
		try
		{
			nodes.first = std::stoi(tokens.at(1));
			nodes.second = std::stoi(tokens.at(2));
			value = si_string_to_double(tokens.at(3));
		}
		catch (const std::out_of_range &ex)
		{
			throw std::runtime_error("Missing arguments (or invalid value)");
		}
		catch (const std::invalid_argument &ex)
		{
			throw std::runtime_error("Invalid value");
		}
	};

	std::pair<int, int> nodes;
	double value;

	if (ref_type == "R")
	{
		parse_component_triplet(tokens, nodes, value);
		return std::make_shared<resistor>(nodes, value);
	}

	if (ref_type == "L")
	{
		parse_component_triplet(tokens, nodes, value);
		return std::make_shared<inductor>(nodes, value);
	}

	if (ref_type == "C")
	{
		parse_component_triplet(tokens, nodes, value);
		return std::make_shared<capacitor>(nodes, value);
	}

	if (ref_type == "E" || ref_type == "V")
	{
		parse_component_triplet(tokens, nodes, value);
		double ac = 0.0;

		if (tokens.size() >= 6 && tolower(tokens[4]) == "ac")
			ac = std::stoi(tokens[5]);

		return std::make_shared<voltage_source>(nodes, value, ac);
	}

	if (ref_type == "I")
	{
		parse_component_triplet(tokens, nodes, value);
		double ac = 0.0;

		if (tokens.size() >= 6 && tolower(tokens[4]) == "ac")
			ac = std::stoi(tokens[5]);

		return std::make_shared<current_source>(nodes, value, ac);
	}

	if (ref_type == "OPA")
	{
		if (tokens.size() < 4)
			throw std::runtime_error("Missing nodes!");

		int pos = std::stoi(tokens[1]);
		int neg = std::stoi(tokens[2]);
		int out = std::stoi(tokens[3]);
		return std::make_shared<opamp>(pos, neg, out);
	}

	throw std::runtime_error("Invalid component type");
}

/**
	\brief Tworzy symulację na podstawie pliku częściowo kompatybilnego z formatem SPICE
*/
static circuit_simulation read_spice_file(std::istream &netlist)
{
	circuit_simulation sim;

	// Pierwsza linia zawiera tytuł
	std::getline(netlist, sim.title);

	// Wszystkie napotkane polecenia
	std::vector<std::string> commands;

	int line_number = 1;
	std::string line;
	while (line_number++, std::getline(netlist, line))
	{
		auto tokens = tokenize_string(line, [](char c){return isspace(c);});
		if (!tokens.size()) continue;

		// Polecenie SPICE do obsłużenia później
		if (tokens[0][0] == '.')
			commands.push_back(line);
		else // Element układu
		{
			auto ref = tokens[0];
			if (sim.circ.find(ref) == sim.circ.end())
			{
				try
				{
					sim.circ[ref] = create_component(tokens);
				}
				catch (const std::exception &ex)
				{
					throw std::runtime_error("Could not parse component in line "s + std::to_string(line_number) + " - reason: " + ex.what());
				}
			}
			else
				throw std::runtime_error("Duplicate components found! (line "s + std::to_string(line_number) + ")");
		}
	}

	// Interpretacja poleceń SPICE
	for (const auto &cmd : commands)
	{
		auto tokens = tokenize_string(cmd, [](char c){return isspace(c);});
		auto lowercase_command = tolower(tokens[0]);

		if (lowercase_command == ".ac")
		{
			if (tokens.size() != 5)
				throw std::runtime_error("Invalid use of .ac command!");

			// Typ sweepa
			auto sweep_type = tolower(tokens[1]);
			double exponent;
			if (sweep_type == "lin")
				exponent = 0;
			else if (sweep_type == "dec")
				exponent = 10;
			else if (sweep_type == "oct")
				exponent = 2;
			else
				throw std::runtime_error("Invalid .ac sweep type!");

			// Parametry liczbowe
			int n;
			double fstart;
			double fstop;
			try
			{
				n = std::stoi(tokens[2]);
				fstart = si_string_to_double(tokens[3]);
				fstop = si_string_to_double(tokens[4]);

				if (fstart <= 0 || fstop <= fstart || n <= 0)
					throw std::runtime_error("Invalid .ac command parameter value");
			}
			catch (const std::exception &ex)
			{
				throw std::runtime_error("Malformed .ac command parameter");
			}

			sim.ac = ac_analysis_params{.start = fstart, .stop = fstop, .exponent = exponent, .steps = n};
		}
		else if (lowercase_command == ".print")
		{
			// Nazwy na różne typy interpretacji zespolonych wielkości fizycznych
			const static std::map<std::string, complex_probing_method> complex_probing_names =
			{
				{"", complex_probing_method::DEFAULT},
				{"re", complex_probing_method::REAL},
				{"im", complex_probing_method::IMAGINARY},
				{"mag", complex_probing_method::MAGNITUDE},
				{"ph", complex_probing_method::PHASE},
			};

			const std::regex probe_regex("([VPI])(re|im|mag|ph)?\\(\\s*([^\\s,]*)(\\s*,\\s*([^\\s,]*))?\\s*\\)",
				std::regex_constants::icase);

			std::smatch match;
			for (std::string str = cmd; std::regex_search(str, match, probe_regex); str = match.suffix())
			{
				auto probe_type = tolower(match[1]);
				auto probing_method_name = tolower(match[2]);
			
				// Metoda pomiaru wielkości zespolonej
				auto probing_method = complex_probing_method::DEFAULT;
				if (auto it = complex_probing_names.find(probing_method_name); it != complex_probing_names.end())
					probing_method = it->second;
				else
					throw std::runtime_error("Invalid probing method '"s + probing_method_name + "'");

				try
				{ 
					// Pomiar napięcia
					if (probe_type == "v")
					{
						int pos, neg = 0;
						try
						{
							pos = std::stoi(match[3]);
							if (!match[4].str().empty())
								neg = std::stoi(match[5]);
						}
						catch (const std::exception &ex)
						{
							throw std::runtime_error("Invalid node numbers in '"s + line + "'");
						}

						sim.probes.push_back(std::make_shared<voltage_probe>(pos, neg, probing_method));
					}

					// Pomiar prądu
					if (probe_type == "i") 
						sim.probes.push_back(std::make_shared<current_probe>(sim.circ, match[3], probing_method));

					// Pomiar mocy
					if (probe_type == "p") 
						sim.probes.push_back(std::make_shared<power_probe>(sim.circ, match[3], probing_method));
				}
				catch (const std::exception &ex)
				{
					throw std::runtime_error("Could not probe '"s + line + "' - reason: "s + ex.what());
				}
			}
		}
		else
		{
			std::cerr << "Ignoring command '" << lowercase_command << "'..." << std::endl;
		}
	}

	return sim;
}

/**
	\brief Funkcja main dla rozszerzonej wersji programu
*/
int main_extended(int argc, char *argv[])
{
	auto &fin = std::cin;
	auto &fout = std::cout;

	try
	{
		auto sim = read_spice_file(fin);

		try
		{
			circuit_solver solver(sim.circ);
	
			if (sim.ac)
			{
				auto &params = *sim.ac;
				auto start_omega = 2 * M_PI * params.start;
				auto stop_omega = 2 * M_PI * params.stop;
				int steps = params.steps;
				bool linear = params.exponent == 0 || params.exponent == 1;

				// Sweep nieliniowy - określona liczba kroków
				// na exponent-krotną zmianę częstotliwości
				if (!linear)
					steps = std::floor(params.steps * std::log(params.stop / params.start) / std::log(params.exponent));

				// Wypisanie nagłówków
				fout << "step\tfrequency\t";
				for (const auto &p : sim.probes)
					fout << p->get_name() << "\t";
				fout << std::endl;

				for (int i = 0; i < steps; i++)
				{
					// Pulsacja dla tego kroku
					double omega;
					if (linear)
						omega = start_omega + (stop_omega - start_omega) * i / (steps - 1);
					else
					{
						auto s = std::log(start_omega) / std::log(params.exponent);
						auto e = std::log(stop_omega) / std::log(params.exponent);
						omega = std::pow(params.exponent, s + (e - s) * i / (steps - 1));
					}

					try
					{
						solver.solve(omega);
					}
					catch (const std::exception &ex)
					{
						throw std::runtime_error("Could not perform "s + std::to_string(i) 
							+ " step of small signal AC analysis - reason: "s + ex.what());
					}

					// Wypisanie mierzonych wartości
					try
					{
						fout << i << "\t" << omega / 2.0 / M_PI << "\t"; 
						for (const auto &p : sim.probes)
						{
							auto val = p->get_value(solver);
							fout << val << "\t";
						}
						fout << std::endl;
					}
					catch (const std::exception &ex)
					{
						throw std::runtime_error("AC probing failed - reason: "s + ex.what());
					}
				}
			}
			else
			{
				try
				{
					solver.solve(0);
				}
				catch (const std::exception &ex)
				{
					throw std::runtime_error(ex.what());
				}

				try
				{
					for (const auto &p : sim.probes)
					{
						auto val = p->get_value(solver);
						fout << p->get_name() << " = " << val << std::endl;
					}
				}
				catch (const std::exception &ex)
				{
					throw std::runtime_error("DC probing failed - reason: "s + ex.what());
				}
			}
		}
		catch (const std::exception &ex)
		{
			std::cerr << "Simulation failed...\nReason: " << ex.what() << std::endl;
			return 1;
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Could not parse SPICE file...\nReason: " << ex.what() << std::endl;
		return 1;
	}

	return 0;
}
