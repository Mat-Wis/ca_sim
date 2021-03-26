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
	read_param<float>(parameters, "toxin_secrete_rate", toxin_secrete_rate);
	read_param<float>(parameters, "init_immune_ratio", init_immune_ratio);

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			cells[i][j] = Cell::Empty;
		}
	}

	for(size_t i = 20; i < 40; ++i) {
		for(size_t j = 20; j < 40; ++j) {
			cells[i][j] = Cell::Healthy;
		}
	}
	for(size_t i = 20; i < 40; ++i) {
		for(size_t j = 60; j < 80; ++j) {
			cells[i][j] = Cell::Tumor;
		}
	}

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			oxygen[i][j] = 0.9;
		}
	}

	std::default_random_engine gen;
	std::uniform_real_distribution<float> dist;
	float num;

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			num = dist(gen);
			if(num < init_immune_ratio) {
				immune[i][j] = Cell::Immune;
			}
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
	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
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
	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
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
	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
			if(j % 2) { // i % 2
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
	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
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

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			d = ox_supply_level - oxygen[i][j];
			oxygen[i][j] += d * ox_supply_rate;
		}
	}
}

void Sim::secrete_toxin() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(cells[i][j] == Cell::Tumor) {
				toxin[i][j] += toxin_secrete_rate;
			}
		}
	}
}

void Sim::uptake_ox() {
	float new_val;
	
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			switch(cells[i][j]) {
				case Cell::Healthy:
					new_val = oxygen[i][j] - healthy_ox_rate;
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

void Sim::move_immune() {
	std::default_random_engine gen;
	std::uniform_int_distribution<int> dist(-1, 1);
	int x, y;
	
	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
			if(immune[i][j] == Cell::Immune) {
				x = dist(gen);
				y = dist(gen);
				
				if(immune[i+x][j+y] == Cell::Empty) {
					immune[i][j] = Cell::Empty;
					immune[i+x][j+y] = Cell::Immune;
				}
			}
		}
	}
}

void Sim::hypoxia() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(oxygen[i][j] < ox_surv_thr) {
				cells[i][j] = Cell::Empty;
				immune[i][j] = Cell::Empty;
				prolif_cnt[i][j] = 0;
			}
		}
	}
}
