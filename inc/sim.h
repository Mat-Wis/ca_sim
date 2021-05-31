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
	Vessel		= 40,
	Immune		= 50
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
		void damage_ecm();
		void move_immune();
		void kill_tumor();
		void kill_immune();
		void kill_healthy();
		void proliferate();
		void recruit_immune();
		void count_cells();
		bool tumor_killed();

		int n_steps;
		int log_step;

	private:
		static constexpr size_t size = 200;
		static constexpr int nbrhood = 8;
		
		/* Matrices */
		Cell cells[size][size];
		Cell immune[size][size];
		int prolif_cnt[size][size];
		int kill_cnt[size][size];
		int life_cnt[size][size];
		float nutrient[size][size];
		float attr[size][size];
		float ecm_stress[size][size];

		float temp_float[size][size];
		
		/* Parameters */
		float sim_time;
		float dt;
		float alpha2;
		float lambda;
		float beta2;
		float nutr_surv_thr;
		float nutr_prolif_thr;
		float stress_thr;
		float imm_rnd;
		float init_immune_ratio;
		int t_cycle;
		int t_steps;
		int kill_limit;
		int life_limit;
		int life_steps;
		bool vessels_on_borders;
		float vessel_num;

		/* Neighbourhood*/
		static constexpr int nbr[][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};
		
		/* Cell numbers */
		int num_healthy;
		int num_tumor;
		int num_deadtumor;
		int num_immune;

		/* Random numbers */
		std::default_random_engine gen;

		/* Functions */
		template<class T>
		void read_param(const libconfig::Setting& setting, const char* name, T& var);

		void healthy_die(size_t i, size_t j);
		void tumor_apoptosis(size_t i, size_t j);
		void tumor_necrosis(size_t i, size_t j);
		void immune_die(size_t i, size_t j);
		void diffuse_nutr(size_t i, size_t j, float& max_diff);
		void diffuse_nutr(size_t i, size_t j, size_t i_m, size_t i_p, size_t j_m, size_t j_p, float& max_diff);
		void diffuse_attr(size_t i, size_t j, float& max_diff);

		friend class Logger;
		std::string mat_file;
		std::string num_file;
};
