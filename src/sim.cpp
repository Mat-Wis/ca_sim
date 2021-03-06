#include "sim.h"

Sim::Sim(char* config_file) : 
	gen(std::chrono::system_clock::now().time_since_epoch().count())
{
	libconfig::Config cfg;

	/* read configuration file */
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

	/* read all parameters */
	read_param<std::string>(root, "mat_file", mat_file);
	read_param<std::string>(root, "num_file", num_file);
	read_param<float>(parameters, "alpha2", alpha2);
	read_param<float>(parameters, "lambda", lambda);
	read_param<float>(parameters, "beta2", beta2);
	read_param<float>(parameters, "nutr_surv_thr", nutr_surv_thr);
	read_param<float>(parameters, "nutr_prolif_thr", nutr_prolif_thr);
	read_param<float>(parameters, "stress_thr", stress_thr);
	read_param<float>(parameters, "imm_rnd", imm_rnd);
	read_param<float>(parameters, "init_immune_ratio", init_immune_ratio);
	read_param<float>(parameters, "sim_time", sim_time);
	read_param<float>(parameters, "dt", dt);
	read_param<int>(parameters, "log_step", log_step);
	read_param<bool>(parameters, "log_mat", log_mat);
	read_param<int>(parameters, "t_cycle", t_cycle);
	read_param<int>(parameters, "kill_limit", kill_limit);
	read_param<int>(parameters, "life_limit", life_limit);
	read_param<bool>(parameters, "vessels_on_borders", vessels_on_borders);
	read_param<float>(parameters, "vessel_num", vessel_num);
	
	t_steps = static_cast<int>(t_cycle * 60.0f / dt);
	n_steps = static_cast<int>(sim_time * 60.0f / dt);
	life_steps = static_cast<int>(life_limit * 24 * 60.0f / dt);
	
	/* fill simulation area with healthy cells */
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			cells[i][j] = Cell::Healthy;
		}
	}

	/* initialize all layers */
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			nutrient[i][j] = 0.9;
			prolif_cnt[i][j] = 0;
			kill_cnt[i][j] = 0;
			life_cnt[i][j] = 0;
			attr[i][j] = 0.0;
		}
	}
	
	/* add tumor cells */
	std::uniform_int_distribution<int> dist_prolif(static_cast<int>(-2 * 60.0f / dt), t_steps);
	int init_tumor = parameters["tumor_x"].getLength();
	int x, y;
	for(int i = 0; i < init_tumor; ++i) {
		try {
			x = parameters["tumor_x"][i];
			y = parameters["tumor_y"][i];

			cells[x][y] = Cell::Tumor;
			prolif_cnt[x][y] = dist_prolif(gen);
		} catch(const libconfig::SettingTypeException &stex) {
			std::cerr << "Wrong type in tumor_x or tumor_y." << std::endl;
			throw;
		}
	}

	/* add blood vessels */
	if(vessels_on_borders) {
		for(size_t j = 0; j < size; ++j) {
			cells[0][j] = Cell::Vessel;
			nutrient[0][j] = 1.0f;
			cells[size-1][j] = Cell::Vessel;
			nutrient[size-1][j] = 1.0f;
		}
	} else {
		std::uniform_real_distribution<float> dist_v(0.0f, 1.0f);
		float thr = vessel_num * static_cast<float>(size*size) / (size*size - init_tumor);
		int R = 10;
		bool clear = true;
		float val;

		for(size_t i = 0; i < size; ++i) {
			for(size_t j = 0; j < size; ++j) {
				val = dist_v(gen);
				if((j <= 85 && val < thr) || (j >= 85 && val < 3*thr)) {
					clear = true;
					for(int m = -R; m <= R; ++m) {
						for(int n = -R; n <= R; ++n) {
							if((m*m + n*n <= R*R) && (i+m >= 0) && (i+m < size) && (j+n >= 0) && (j+n < size) && cells[i+m][j+n] == Cell::Tumor) {
								clear = false;
							}
						}
					}
					if(clear) {
						cells[i][j] = Cell::Vessel;
						nutrient[i][j] = 1.0f;
					}
				}
			}
		}
	}

	/* add immune cells */
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	std::uniform_int_distribution<int> dist_int(0, life_steps);
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(dist(gen) < init_immune_ratio) {
				immune[i][j] = Cell::Immune;
				life_cnt[i][j] = dist_int(gen);
			}
		}
	}

	count_cells();
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
		std::cerr << "Wrong type in setting '" << name << "'." << std::endl;
		throw;
	}
}

inline void Sim::diffuse_nutr(size_t i, size_t j, float& max_diff) {
	if(cells[i][j] == Cell::Vessel) {
		return;
	}

	float A = alpha2 * (static_cast<float>(cells[i][j] == Cell::Healthy) +
					   	static_cast<float>(immune[i][j] == Cell::Immune) + 
						lambda * static_cast<float>(cells[i][j] == Cell::Tumor)); 

	float nabla2 = temp_float[i-1][j] + temp_float[i+1][j] + temp_float[i][j-1] + temp_float[i][j+1] - 4 * temp_float[i][j];

	float r = temp_float[i][j];
	float k1 = nabla2 - A*r;
	float k2 = nabla2 - A*(r + diff_dt*k1/2);
	float k3 = nabla2 - A*(r + diff_dt*k2/2);
	float k4 = nabla2 - A*(r + diff_dt*k3);

	nutrient[i][j] = r + 1.0f/6.0f * diff_dt * (k1 + 2*k2 + 2*k3 + k4);;

	if(nutrient[i][j] < 0.0f) {
		nutrient[i][j] = 0.0f;
	}

	float diff = std::abs(nutrient[i][j] - temp_float[i][j]);
	if(diff > max_diff) {
		max_diff = diff;
	}
}

inline void Sim::diffuse_nutr(size_t i, size_t j, size_t i_m, size_t i_p, size_t j_m, size_t j_p, float& max_diff) {
	if(cells[i][j] == Cell::Vessel) {
		return;
	}

	float A = alpha2 * (static_cast<float>(cells[i][j] == Cell::Healthy) +
					   	static_cast<float>(immune[i][j] == Cell::Immune) + 
						lambda * static_cast<float>(cells[i][j] == Cell::Tumor)); 

	float nabla2 = temp_float[i_m][j] + temp_float[i_p][j] + temp_float[i][j_m] + temp_float[i][j_p] - 4 * temp_float[i][j];

	float r = temp_float[i][j];
	float k1 = nabla2 - A*r;
	float k2 = nabla2 - A*(r + diff_dt*k1/2);
	float k3 = nabla2 - A*(r + diff_dt*k2/2);
	float k4 = nabla2 - A*(r + diff_dt*k3);

	nutrient[i][j] = r + 1.0f/6.0f * diff_dt * (k1 + 2*k2 + 2*k3 + k4);;

	if(nutrient[i][j] < 0.0f) {
		nutrient[i][j] = 0.0f;
	}

	float diff = std::abs(nutrient[i][j] - temp_float[i][j]);
	if(diff > max_diff) {
		max_diff = diff;
	}
}

inline void Sim::diffuse_attr(size_t i, size_t j, float& max_diff) {
	float A = beta2 * static_cast<float>(cells[i][j] == Cell::Tumor || cells[i][j] == Cell::DeadTumor);
	float nabla2 = temp_float[i-1][j] + temp_float[i+1][j] + temp_float[i][j-1] + temp_float[i][j+1] - 4 * temp_float[i][j];
	
	attr[i][j] = temp_float[i][j] + diff_dt * (nabla2 + A);
	
	if(attr[i][j] < 0) {
		attr[i][j] = 0;
	}
	
	float diff = std::abs(attr[i][j] - temp_float[i][j]);
	if(diff > max_diff) {
		max_diff = diff;
	}
}

void Sim::diffuse() {
	float max_diff;
	
	/* Diffuse nutrient */
	size_t n;
	for(n = 0; n < 100000; ++n) {
		std::memcpy(temp_float, nutrient, size * size * sizeof(float));
		max_diff = 0.0f;
		
		/* diffuse left column (boundary conditions) */
		for(size_t j = 0; j < size; ++j) {
			diffuse_nutr(0, j, size-1, 1, j-1, j+1, max_diff);
		}

		/* diffuse right column (boundary conditions) */
		for(size_t j = 0; j < size; ++j) {
			diffuse_nutr(size-1, j, size-2, 0, j-1, j+1, max_diff);
		}

		for(size_t i = 1; i < size-1; ++i) {
			/* diffuse top row (boundary conditions) */
			diffuse_nutr(i, 0, i-1, i+1, size-1, 1, max_diff);

			/* diffuse bottom row (boundary conditions) */
			diffuse_nutr(i, size-1, i-1, i+1, size-2, 0, max_diff);
			
			for(size_t j = 1; j < size-1; ++j) {
				diffuse_nutr(i, j, max_diff);
			}
		}

		if(max_diff < 0.0000003) {
			break;
		}
	}

	/* Diffuse attractant */
	for(size_t n = 0; n < 100000; ++n) {
		std::memcpy(temp_float, attr, size * size * sizeof(float));
		max_diff = 0.0f;
		
		for(size_t i = 1; i < size-1; ++i) {
			for(size_t j = 1; j < size-1; ++j) {
				diffuse_attr(i, j, max_diff);
			}
		}

		if(max_diff < 0.0000003) {
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
					
					ecm_stress[i+x][j+y] += dist_f(gen) * (dt / 20.0f);
				}
			}
		}
	}
}

void Sim::move_immune() {
	size_t i, j, x, y;
	std::vector<Coord> immune_cells;
	float gx, gy, rndx, rndy, rnd_norm, vecx, vecy, angle;
	int n_angle;
	std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
	
	for(size_t n = 0; n < size; ++n) {
		for(size_t m = 0; m < size; ++m) {
			if(immune[n][m] == Cell::Immune) {
				immune_cells.push_back({n, m});
			}
		}
	}
	std::random_shuffle(immune_cells.begin(), immune_cells.end());

	for(auto const& c : immune_cells) {
		i = c.x;
		j = c.y;

		if(i == 0) {
			gx = attr[i+1][j] - attr[i][j];
		} else if(i == size-1) {
			gx = attr[i][j] - attr[i-1][j];
		} else {
			gx = (attr[i+1][j] - attr[i-1][j]) / 2;
		}

		if(j == 0) {
			gy = attr[i][j+1] - attr[i][j];
		} else if(j == size-1) {
			gy = attr[i][j] - attr[i][j-1];
		} else {
			gy = (attr[i][j+1] - attr[i][j-1]) / 2;
		}
		
		rndx = dist(gen);
		rndy = dist(gen);
		rnd_norm = std::sqrt(rndx*rndx + rndy*rndy) / imm_rnd;
		rndx /= rnd_norm;
		rndy /= rnd_norm;
		
		vecx = gx + rndx;
		vecy = gy + rndy;

		angle = std::atan2(vecx, vecy);
		n_angle = std::floor(angle / (M_PI / 8.0f));
		
		switch(n_angle) {
			case 7:
			case -8:
				x = i; y = j-1; break;
			case -7:
			case -6:
				x = i-1; y = j-1; break;
			case -5:
			case -4:
				x = i-1; y = j; break;
			case -3:
			case -2:
				x = i-1; y = j+1; break;
			case -1:
			case 0:
				x = i; y = j+1; break;
			case 1:
			case 2:
				x = i+1; y = j+1; break;
			case 3:
			case 4:
				x = i+1; y = j; break;
			case 5:
			case 6:
				x = i+1; y = j-1; break;
			default:
				x = i; y = j; break;
		}

		if(x >= 0 && x < size && y >= 0 && y < size && immune[x][y] == Cell::Empty) {
			immune[i][j] = Cell::Empty;
			immune[x][y] = Cell::Immune;
			kill_cnt[x][y] = kill_cnt[i][j];
			kill_cnt[i][j] = 0;
			life_cnt[x][y] = life_cnt[i][j];
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
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(cells[i][j] == Cell::Tumor) {
				if(immune[i][j] == Cell::Immune) {
					tumor_apoptosis(i, j);
					++kill_cnt[i][j];
				} else if(nutrient[i][j] < nutr_surv_thr) {
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
				if(kill_cnt[i][j] >= kill_limit || life_cnt[i][j] >= life_steps || nutrient[i][j] < nutr_surv_thr) {
					immune_die(i, j);
				}
			}
		}
	}
}

void Sim::kill_healthy() {
	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(cells[i][j] == Cell::Healthy && (ecm_stress[i][j] >= stress_thr || nutrient[i][j] < nutr_surv_thr)) {
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
	std::uniform_int_distribution<int> dist_prolif(static_cast<int>(-2 * 60.0f / dt), static_cast<int>(2 * 60.0f / dt));

	for(size_t i = 1; i < size-1; ++i) {
		for(size_t j = 1; j < size-1; ++j) {
			if(cells[i][j] == Cell::Tumor && nutrient[i][j] > nutr_prolif_thr) {
				tumor_cells.push_back({i, j});
			}
		}
	}
	std::random_shuffle(tumor_cells.begin(), tumor_cells.end());

	for(auto const& c : tumor_cells) {
		i = c.x;
		j = c.y;
		
		++prolif_cnt[i][j];
		if(prolif_cnt[i][j] >= t_steps) {
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
				prolif_cnt[i+x][j+y] = dist_prolif(gen);
				prolif_cnt[i][j] = dist_prolif(gen);
			}
		}
	}
}

void Sim::recruit_immune() {
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);
	float num, thr;
	std::vector<Coord> vessels;
	float ves_n;

	vessels.clear();

	for(size_t i = 0; i < size; ++i) {
		for(size_t j = 0; j < size; ++j) {
			if(cells[i][j] == Cell::Vessel) {
				vessels.push_back({i, j});
			}
		}
	}
	
	ves_n = vessels.size();
	thr = (init_immune_ratio * size * size - num_immune) / ves_n;

	for(auto v : vessels) {
		num = dist(gen);
		if((num <= thr) && (immune[v.x][v.y] == Cell::Empty)) {
			immune[v.x][v.y] = Cell::Immune;
		}
	}
}

void Sim::count_cells() {
	num_healthy = 0;
	num_tumor = 0;
	num_deadtumor = 0;
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
					if(j < 70) {
						++num_tumor;
					}
					break;
				case Cell::DeadTumor:
					++num_deadtumor;
					break;
				default:
					break;
			}
		}
	}
}

bool Sim::tumor_killed() {
	return num_tumor == 0;
}
