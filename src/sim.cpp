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
	read_param<float>(parameters, "sim_time", sim_time);
	read_param<float>(parameters, "dt", dt);
	read_param<int>(parameters, "log_step", log_step);
	read_param<int>(parameters, "t_cycle", t_cycle);
	read_param<int>(parameters, "kill_limit", kill_limit);
	read_param<int>(parameters, "life_limit", life_limit);

	n_steps = static_cast<int>(sim_time / dt * 60.0f);

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

	for(size_t j = 0; j < size; ++j) {
		oxygen[0][j] = 1.0;
		oxygen[size-1][j] = 1.0;
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
	} catch(const libconfig::SettingTypeException &stex) {
		std::cerr << "Wrong type in setting '" << name << std::endl;
		throw;
	}
}

void Sim::diffuse() {
	float d, c;
	float diff_dt = 0.1;
	float diff, max_diff;
	
	for(int n = 0; n < 1000; ++n) {
		max_diff = 0.0f;
		std::memcpy(temp_float, oxygen, size * size * sizeof(float));
		
		for(size_t i = 1; i < size-1; ++i) {
			/* diffuse top row (boundary conditions) */
			d = temp_float[i-1][0] + temp_float[i+1][0] + temp_float[i][size-1] + temp_float[i][1] - 4 * temp_float[i][0];

			if(cells[i][0] == Cell::Healthy) {
				c = healthy_ox_rate;
			} else if(cells[i][0] == Cell::Tumor) {
				c = tumor_ox_rate;
			} else {
				c = 0.0f;
			}

			oxygen[i][0] += (d * diff_rate - c) * diff_dt;
			if(oxygen[i][0] < 0) {
				oxygen[i][0] = 0;
			}

			diff = std::abs(oxygen[i][0] - temp_float[i][0]);
			if(diff > max_diff) {
				max_diff = diff;
			}

			/* diffuse bottom row (boundary conditions) */
			d = temp_float[i-1][size-1] + temp_float[i+1][size-1] + temp_float[i][size-2] + temp_float[i][0] - 4 * temp_float[i][size-1];

			if(cells[i][size-1] == Cell::Healthy) {
				c = healthy_ox_rate;
			} else if(cells[i][size-1] == Cell::Tumor) {
				c = tumor_ox_rate;
			} else {
				c = 0.0f;
			}

			oxygen[i][size-1] += (d * diff_rate - c) * diff_dt;
			if(oxygen[i][size-1] < 0) {
				oxygen[i][size-1] = 0;
			}

			diff = std::abs(oxygen[i][size-1] - temp_float[i][size-1]);
			if(diff > max_diff) {
				max_diff = diff;
			}
			
			for(size_t j = 1; j < size-1; ++j) {
				d = temp_float[i-1][j] + temp_float[i+1][j] + temp_float[i][j-1] + temp_float[i][j+1] - 4 * temp_float[i][j];

				if(cells[i][j] == Cell::Healthy) {
					c = healthy_ox_rate;
				} else if(cells[i][j] == Cell::Tumor) {
					c = tumor_ox_rate;
				} else {
					c = 0.0f;
				}

				oxygen[i][j] += (d * diff_rate - c) * diff_dt;
				if(oxygen[i][j] < 0) {
					oxygen[i][j] = 0;
				}

				diff = std::abs(oxygen[i][j] - temp_float[i][j]);
				if(diff > max_diff) {
					max_diff = diff;
				}
			}
		}
		
		if(max_diff < 0.00001) {
			break;
		}
	}
}

void Sim::damage_ecm() {
	int x, y;
	int n;
	std::vector<int> i_vec;
	std::vector<int> j_vec;
	std::uniform_real_distribution<float> dist_f(0.0f, 0.1f);
	
	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
			if(cells[i][j] == Cell::Tumor) {
				i_vec.clear();
				j_vec.clear();
				
				for(int n = 0; n < nbrhood; ++n) {
					x = nbr[n][0];
					y = nbr[n][1];
					
					if(0 <= i+x && i+x < size && 0 <= j+x && j+x < size && cells[i+x][j+y] == Cell::Healthy) {
						i_vec.push_back(x);
						j_vec.push_back(y);
					}
				}
				
				if(!i_vec.empty()) {
					std::uniform_int_distribution<int> dist_n(0, i_vec.size()-1);
					n = dist_n(gen);
					x = i_vec[n];
					y = j_vec[n];
					
					ecm_stress[i+x][j+y] += dist_f(gen);
				}
			}
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

		x = dist_1(gen)+1;
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

inline void Sim::healthy_die(size_t i, size_t j) {
	cells[i][j] = Cell::Empty;
	ecm_stress[i][j] = 0.0f;
}

inline void Sim::tumor_apoptosis(size_t i, size_t j) {
	cells[i][j] = Cell::Empty;
	prolif_cnt[i][j] = 0;
}

inline void Sim::tumor_necrosis(size_t i, size_t j) {
	cells[i][j] = Cell::DeadTumor;
	prolif_cnt[i][j] = 0;
}

inline void Sim::immune_die(size_t i, size_t j) {
	immune[i][j] = Cell::Empty;
	kill_cnt[i][j] = 0;
	life_cnt[i][j] = 0;
}

void Sim::kill_tumor() {
	std::uniform_real_distribution<float> dist_f(0.0f, 1.0f);
	std::uniform_int_distribution<int> dist_5(-5, 5);
	int x, y;

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(cells[i][j] == Cell::Tumor) {
				if(immune[i][j] == Cell::Immune) {
					tumor_apoptosis(i, j);
					++kill_cnt[i][j];

					x = dist_5(gen);
					y = dist_5(gen);
					if(i+x >= 0 && i+x < size && j+y >= 0 && j+y < size && immune[i+x][j+y] == Cell::Empty && dist_f(gen) < 0.2) {
						immune[i+x][j+y] = Cell::Immune;
					}
				}
				if(oxygen[i][j] < ox_surv_thr) {
					tumor_necrosis(i, j);
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
			if(cells[i][j] == Cell::Healthy && (ecm_stress[i][j] >= 3.0f || oxygen[i][j] < ox_surv_thr)) {
				healthy_die(i, j);
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
				x = nbr[n][0];
				y = nbr[n][1];

				if(0 <= i+x && i+x < size && 0 <= j+x && j+x < size && cells[i+x][j+y] == Cell::Empty) {
					i_vec.push_back(x);
					j_vec.push_back(y);
				}
			}

			if(!i_vec.empty()) {
				std::uniform_int_distribution<int> dist_n(0, i_vec.size()-1);
				n = dist_n(gen);
				x = i_vec[n];
				y = j_vec[n];

				cells[i+x][j+y] = Cell::Tumor;
			}
		}
	}
}

void Sim::recruit_immune() {
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float num;

	for(size_t j = 0; j < size; ++j) {
		num = dist(gen);
		if(num <= init_immune_ratio / life_limit * size && immune[0][j] == Cell::Empty) {
			immune[0][j] = Cell::Immune;
		}

		num = dist(gen);
		if(num <= init_immune_ratio / life_limit * size && immune[size-1][j] == Cell::Empty) {
			immune[size-1][j] = Cell::Immune;
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
