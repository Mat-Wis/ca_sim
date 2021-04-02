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
	read_param<float>(parameters, "toxin_thr", toxin_thr);
	read_param<float>(parameters, "init_immune_ratio", init_immune_ratio);
	read_param<int>(parameters, "t_cycle", t_cycle);
	read_param<int>(parameters, "kill_limit", kill_limit);
	read_param<int>(parameters, "life_limit", life_limit);
	read_param<int>(parameters, "n_steps", n_steps);

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			cells[i][j] = Cell::Healthy;
			//cells[i][j] = Cell::Empty;
		}
	}

	for(size_t i = size/2-10; i <= size/2+10; ++i) {
		for(size_t j = size/2-10; j <= size/2+10; ++j) {
			if((i-size/2)*(i-size/2) + (j-size/2)*(j-size/2) <= 100) {
				cells[i][j] = Cell::Tumor;
			}
		}
	}

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			oxygen[i][j] = 0.9;
			prolif_cnt[i][j] = 0;
			kill_cnt[i][j] = 0;
			life_cnt[i][j] = 0;
		}
	}

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
			d = temp_float[i-1][j] + temp_float[i+1][j] + temp_float[i][j-1] + temp_float[i][j+1] - 4 * temp_float[i][j];
			subst[i][j] += d * diff_rate;
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
	size_t i, j;
	std::vector<Coord> immune_cells;
	
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(immune[i][j] == Cell::Immune) {
				immune_cells.push_back({i, j});
			}
		}
	}
	std::random_shuffle(immune_cells.begin(), immune_cells.end());

	for(auto const& c : immune_cells) {
		i = c.x;
		j = c.y;

		x = dist_1(gen);
		y = dist_1(gen);

		if(i+x < size && i+x >= 0 && j+x < size && j+x >= 0 && immune[i+x][j+y] == Cell::Empty) {
			immune[i][j] = Cell::Empty;
			immune[i+x][j+y] = Cell::Immune;
			kill_cnt[i+x][j+y] = kill_cnt[i][j];
			kill_cnt[i][j] = 0;
			life_cnt[i+x][j+y] = life_cnt[i][j];
			life_cnt[i][j] = 0;
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
	life_cnt[i][j] = 0;
}

void Sim::kill_tumor() {
	std::uniform_int_distribution<int> dist_5(-5, 5);
	std::uniform_real_distribution<float> dist_f(0.0f, 1.0f);
	int x, y;

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(immune[i][j] == Cell::Immune && cells[i][j] == Cell::Tumor) {
				cell_die(i, j);
				++kill_cnt[i][j];

				x = dist_5(gen);
				y = dist_5(gen);
				if(immune[i+x][j+y] == Cell::Empty && dist_f(gen) < 0.2) {
					immune[i+x][j+y] = Cell::Immune;
				}
			}
		}
	}
}

void Sim::kill_immune() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(immune[i][j] == Cell::Immune) {
				++life_cnt[i][j];
				if((kill_cnt[i][j] >= kill_limit || life_cnt[i][j] >= life_limit)) {
					immune_die(i, j);
				}
			}
		}
	}
}

void Sim::kill_healthy() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(cells[i][j] == Cell::Healthy && toxin[i][j] >= toxin_thr) {
				cell_die(i, j);
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

void Sim::proliferate() {
	int x, y;
	int n;
	size_t i, j;
	std::vector<int> i_vec;
	std::vector<int> j_vec;
	std::vector<Coord> tumor_cells;

	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
			if(cells[i][j] == Cell::Tumor && oxygen[i][j] > ox_prolif_thr) {
				tumor_cells.push_back({i, j});
			}
		}
	}
	std::random_shuffle(tumor_cells.begin(), tumor_cells.end());

	for(auto const& c : tumor_cells) {
		i = c.x;
		j = c.y;
		
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

				if(0 <= i+x && i+x < size && 0 <= j+x && j+x < size && cells[i+x][j+y] == Cell::Empty) {
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

void Sim::recruit_immune() {
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float num;

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			num = dist(gen);
			if(num <= init_immune_ratio / life_limit && immune[i][j] == Cell::Empty) {
				immune[i][j] = Cell::Immune;
			}
		}
	}	
}

void Sim::count_cells() {
	num_healthy = 0;
	num_tumor = 0;
	num_immune = 0;

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(immune[i][j] == Cell::Immune) {
				++num_immune;
			}
			switch(cells[i][j]) {
				case Cell::Healthy:
					++num_healthy;
					break;
				case Cell::Tumor:
					++num_tumor;
					break;
				default:
					break;
			}
		}
	}
}
