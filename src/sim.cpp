#include "sim.h"

Sim::Sim(char* config_file) : 
	gen(std::chrono::system_clock::now().time_since_epoch().count()), 
	dist_1(-1, 1)
{
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
	read_param<float>(parameters, "diff_rate", diff_rate);
	read_param<float>(parameters, "ox_supply_level", ox_supply_level);
	read_param<float>(parameters, "ox_supply_rate", ox_supply_rate);
	read_param<float>(parameters, "healthy_ox_rate", healthy_ox_rate);
	read_param<float>(parameters, "immune_ox_rate", immune_ox_rate);
	read_param<float>(parameters, "tumor_ox_rate", tumor_ox_rate);
	read_param<float>(parameters, "ox_surv_thr", ox_surv_thr);
	read_param<float>(parameters, "ox_prolif_thr", ox_prolif_thr);
	read_param<float>(parameters, "toxin_secrete_rate", toxin_secrete_rate);
	read_param<float>(parameters, "init_immune_ratio", init_immune_ratio);
	read_param<int>(parameters, "t_cycle", t_cycle);
	read_param<int>(parameters, "kill_limit", kill_limit);

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			cells[i][j] = Cell::Healthy;
			//cells[i][j] = Cell::Empty;
		}
	}

	for(size_t i = 40; i <= 60; ++i) {
		for(size_t j = 40; j <= 60; ++j) {
			if((i-50)*(i-50) + (j-50)*(j-50) <= 100) {
				cells[i][j] = Cell::Tumor;
			}
		}
	}

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			oxygen[i][j] = 0.9;
			prolif_cnt[i][j] = 0;
			kill_cnt[i][j] = 0;
		}
	}

	std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float num;

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			num = dist(gen);
			if(num <= init_immune_ratio) {
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
	diffuse_(oxygen);
	diffuse_(toxin);
}

void Sim::diffuse_(float subst[size][size]) {
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
	int x, y;
	
	memcpy(temp_cell, immune, size * size * sizeof(Cell));
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(temp_cell[i][j] == Cell::Immune) {
				x = dist_1(gen);
				y = dist_1(gen);

				if(i+x < size && i+x >= 0 && j+x < size && j+x >= 0 && temp_cell[i+x][j+y] == Cell::Empty && immune[i+x][j+y] == Cell::Empty) {
					immune[i][j] = Cell::Empty;
					immune[i+x][j+y] = Cell::Immune;
					kill_cnt[i+x][j+y] = kill_cnt[i][j];
					kill_cnt[i][j] = 0;
				}
			}
		}
	}
}

inline void Sim::cell_die(size_t i, size_t j) {
	cells[i][j] = Cell::Empty;
	prolif_cnt[i][j] = 0;
}

inline void Sim::immune_die(size_t i, size_t j) {
	immune[i][j] = Cell::Empty;
	kill_cnt[i][j] = 0;
}

void Sim::kill_tumor() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(immune[i][j] == Cell::Immune && cells[i][j] == Cell::Tumor) {
				cell_die(i, j);
				++kill_cnt[i][j];
			}
		}
	}
}

void Sim::kill_immune() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(immune[i][j] == Cell::Immune && kill_cnt[i][j] >= kill_limit) {
				immune_die(i, j);
			}
		}
	}
}

void Sim::proliferate() {
	int x, y;
	int n;
	std::vector<int> i_vec;
	std::vector<int> j_vec;

	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
			if(cells[i][j] == Cell::Tumor && prolif_cnt[i][j] > ox_prolif_thr) {
				++prolif_cnt[i][j];
				if(prolif_cnt[i][j] >= t_cycle) {
					i_vec.clear();
					j_vec.clear();
					
					for(int n = 0; n < nbrhood; ++n) {
						//if(j % 2) {
							//x = nbr_even[n][0];
							//y = nbr_even[n][1];
						//} else {
							//x = nbr_odd[n][0];
							//y = nbr_odd[n][1];
						//}
						x = nbr[n][0];
						y = nbr[n][1];

						if(cells[i+x][j+y] == Cell::Empty) {
							i_vec.push_back(i+x);
							j_vec.push_back(j+y);
						}
					}

					if(!i_vec.empty()) {
						std::uniform_int_distribution<int> dist_n(0, i_vec.size()-1);
						n = dist_n(gen);
						x = i_vec[n];
						y = j_vec[n];

						cells[x][y] = Cell::Tumor;
					}
				}
			}
		}
	}
}

void Sim::hypoxia() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(oxygen[i][j] < ox_surv_thr) {
				cell_die(i, j);
				//immune_die(i, j);
			}
		}
	}
}
