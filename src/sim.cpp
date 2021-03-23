#include "sim.h"

Sim::Sim(char* config_file) {
	libconfig::Config cfg;

	try {
		cfg.readFile(config_file);
	} catch(const libconfig::FileIOException &fioex) {
		std::cerr << "I/O error while reading configuration file." << std::endl;
		throw;
	} catch(const libconfig::ParseException &pex) {
		std::cerr << "Error reading configuration file. Parse error at " << pex.getFile() 
			<< ":" << pex.getLine() << " - " << pex.getError() << std::endl;
		throw;
	}

	const libconfig::Setting& root = cfg.getRoot();
	const libconfig::Setting& parameters = root["parameters"];

	read_param<std::string>(root, "logfile", logfile);
	read_param<int>(parameters, "nbrhood", nbrhood);
	read_param<float>(parameters, "diff_rate", diff_rate);
	read_param<float>(parameters, "ox_supply_level", ox_supply_level);
	read_param<float>(parameters, "ox_supply_rate", ox_supply_rate);
	read_param<float>(parameters, "healthy_ox_rate", healthy_ox_rate);
	read_param<float>(parameters, "immune_ox_rate", immune_ox_rate);
	read_param<float>(parameters, "tumor_ox_rate", tumor_ox_rate);

	//std::cout << "nbrhood: " << nbrhood << std::endl;
	//std::cout << "diff_rate: " << diff_rate << std::endl;
	//std::cout << "ox_supply_level: " << ox_supply_level << std::endl;
	//std::cout << "ox_supply_rate: " << ox_supply_rate << std::endl;
	
	for(size_t i = 0; i < size; i++) {
		for(size_t j = 0; j < size; j++) {
			cells[i][j] = Cell::Empty;
		}
	}

	for(size_t i = 40; i < 60; i++) {
		for(size_t j = 40; j < 60; j++) {
			cells[i][j] = Cell::Healthy;
		}
	}
}

Sim::~Sim() {}

template<class T>
void Sim::read_param(const libconfig::Setting& setting, const char* name, T& var) {
	try {
		var = T(setting.lookup(name));
	} catch(const libconfig::SettingNotFoundException &snfex) {
		std::cerr << "Setting '" << name << "' not found in configuration file." << std::endl;
		throw;
	}
}

void* Sim::get_cells_ptr() {
	return (void *)&cells;
}

void* Sim::get_oxygen_ptr() {
	return (void *)&oxygen;
}

void* Sim::get_toxin_ptr() {
	return (void *)&toxin;
}

void Sim::diffuse() {
	switch(nbrhood) {
		case 4:
			diffuse_4(oxygen);
			diffuse_4(toxin);
			break;
		case 8:
			diffuse_8(oxygen);
			diffuse_8(toxin);
			break;
		case 6:
			diffuse_6_my(oxygen);
			diffuse_6_my(toxin);
			break;
		case 7:
			diffuse_6_ext(oxygen);
			diffuse_6_ext(toxin);
			break;
		default:
			break;
	}
}
void Sim::diffuse_4(float subst[size][size]) {
	float d;

	std::memcpy(temp_float, subst, size * size * sizeof(float));
	for(size_t i = 1; i < size-1; i++) {
		for(size_t j = 1; j < size-1; j++) {
			d = temp_float[i][j] - temp_float[i-1][j];
			subst[i][j] -= d * diff_rate;
			subst[i-1][j] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i][j-1];
			subst[i][j] -= d * diff_rate;
			subst[i][j-1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i][j+1];
			subst[i][j] -= d * diff_rate;
			subst[i][j+1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i+1][j];
			subst[i][j] -= d * diff_rate;
			subst[i+1][j] += d * diff_rate;
		}
	}
}
void Sim::diffuse_8(float subst[size][size]) {
	float d;

	std::memcpy(temp_float, subst, size * size * sizeof(float));
	for(size_t i = 1; i < size-1; i++) {
		for(size_t j = 1; j < size-1; j++) {
			d = temp_float[i][j] - temp_float[i-1][j-1];
			subst[i][j] -= d * diff_rate;
			subst[i-1][j-1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i-1][j];
			subst[i][j] -= d * diff_rate;
			subst[i-1][j] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i-1][j+1];
			subst[i][j] -= d * diff_rate;
			subst[i-1][j+1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i][j-1];
			subst[i][j] -= d * diff_rate;
			subst[i][j-1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i][j+1];
			subst[i][j] -= d * diff_rate;
			subst[i][j+1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i+1][j-1];
			subst[i][j] -= d * diff_rate;
			subst[i+1][j-1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i+1][j];
			subst[i][j] -= d * diff_rate;
			subst[i+1][j] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i+1][j+1];
			subst[i][j] -= d * diff_rate;
			subst[i+1][j+1] += d * diff_rate;
		}
	}
}
void Sim::diffuse_6_my(float subst[size][size]) {
	float d;

	std::memcpy(temp_float, subst, size * size * sizeof(float));
	for(size_t i = 1; i < size-1; i++) {
		for(size_t j = 1; j < size-1; j++) {
			if(j % 2) {
				d = temp_float[i][j] - temp_float[i-1][j];
				subst[i][j] -= d * diff_rate;
				subst[i-1][j] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i][j-1];
				subst[i][j] -= d * diff_rate;
				subst[i][j-1] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i][j+1];
				subst[i][j] -= d * diff_rate;
				subst[i][j+1] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i+1][j-1];
				subst[i][j] -= d * diff_rate;
				subst[i+1][j-1] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i+1][j];
				subst[i][j] -= d * diff_rate;
				subst[i+1][j] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i+1][j+1];
				subst[i][j] -= d * diff_rate;
				subst[i+1][j+1] += d * diff_rate;
			} else {
				d = temp_float[i][j] - temp_float[i-1][j-1];
				subst[i][j] -= d * diff_rate;
				subst[i-1][j-1] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i-1][j];
				subst[i][j] -= d * diff_rate;
				subst[i-1][j] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i-1][j+1];
				subst[i][j] -= d * diff_rate;
				subst[i-1][j+1] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i][j-1];
				subst[i][j] -= d * diff_rate;
				subst[i][j-1] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i][j+1];
				subst[i][j] -= d * diff_rate;
				subst[i][j+1] += d * diff_rate;

				d = temp_float[i][j] - temp_float[i+1][j];
				subst[i][j] -= d * diff_rate;
				subst[i+1][j] += d * diff_rate;
			}
		}
	}
}
void Sim::diffuse_6_ext(float subst[size][size]) {
	float d;

	std::memcpy(temp_float, subst, size * size * sizeof(float));
	for(size_t i = 1; i < size-1; i++) {
		for(size_t j = 1; j < size-1; j++) {
			d = temp_float[i][j] - temp_float[i-1][j-1];
			subst[i][j] -= d * diff_rate;
			subst[i-1][j-1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i-1][j];
			subst[i][j] -= d * diff_rate;
			subst[i-1][j] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i][j-1];
			subst[i][j] -= d * diff_rate;
			subst[i][j-1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i][j+1];
			subst[i][j] -= d * diff_rate;
			subst[i][j+1] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i+1][j];
			subst[i][j] -= d * diff_rate;
			subst[i+1][j] += d * diff_rate;

			d = temp_float[i][j] - temp_float[i+1][j+1];
			subst[i][j] -= d * diff_rate;
			subst[i+1][j+1] += d * diff_rate;
		}
	}
}

void Sim::oxygenate() {
	float d;

	for(size_t i = 0; i < size; i++) {
		for(size_t j = 0; j < size; j++) {
			d = ox_supply_level - oxygen[i][j];
			oxygen[i][j] += d * ox_supply_rate;
		}
	}
}

void Sim::uptake_ox() {
	float new_val;
	
	for(size_t i = 0; i < size; i++) {
		for(size_t j = 0; j < size; j++) {
			switch(cells[i][j]) {
				case Cell::Healthy:
					new_val = oxygen[i][j] - healthy_ox_rate;
					break;
				case Cell::Immune:
					new_val = oxygen[i][j] - immune_ox_rate;
					break;
				case Cell::Tumor:
					new_val = oxygen[i][j] - tumor_ox_rate;
					break;
				default:
					new_val = oxygen[i][j];
					break;
			}
			
			oxygen[i][j] = (new_val > 0.0f) ? new_val : 0.0f;
		}
	}
}
