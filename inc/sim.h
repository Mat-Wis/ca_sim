#pragma once

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <libconfig.h++>

enum class Cell : int { 
	Empty		= 0,
	Healthy		= 10, 
	Tumor		= 20,
	DeadTumor	= 30,
	Immune		= 40
};

struct Coord {
	size_t x;
	size_t y;
};

class Sim {
	public:
		Sim(char* config_file);
		~Sim();

		void diffuse();
		void oxygenate();
		void damage_ecm();
		void uptake_ox();
		void move_immune();
		void kill_tumor();
		void kill_immune();
		void kill_healthy();
		void proliferate();
		void recruit_immune();
		void count_cells();

		int n_steps;
		int log_step;

	private:
		static constexpr size_t size = 100;
		static constexpr int nbrhood = 8;
		
		/* Matrices */
		Cell cells[size][size];
		Cell immune[size][size];
		int prolif_cnt[size][size];
		int kill_cnt[size][size];
		int life_cnt[size][size];
		float oxygen[size][size];
		float ecm_stress[size][size];

		Cell temp_cell[size][size];
		int temp_int[size][size];
		float temp_float[size][size];
		
		/* Parameters */
		float sim_time;
		float dt;
		float diff_rate;
		float ox_supply_level;
		float ox_supply_rate;
		float healthy_ox_rate;
		float immune_ox_rate;
		float tumor_ox_rate;
		float ox_surv_thr;
		float ox_prolif_thr;
		float toxin_secrete_rate;
		float toxin_thr;
		float init_immune_ratio;
		int t_cycle;
		int kill_limit;
		int life_limit;

		/* Neighbourhood*/
		static constexpr int nbr[][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};
		
		/* Cell numbers */
		int num_healthy;
		int num_tumor;
		int num_immune;

		/* Random numbers */
		std::default_random_engine gen;
		std::uniform_int_distribution<int> dist_1;

		/* Functions */
		template<class T>
		void read_param(const libconfig::Setting& setting, const char* name, T& var);

		void healthy_die(size_t i, size_t j);
		void tumor_apoptosis(size_t i, size_t j);
		void tumor_necrosis(size_t i, size_t j);
		void immune_die(size_t i, size_t j);

		friend class Logger;
		std::string logfile;
};
