#pragma once

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <libconfig.h++>

enum class Cell : int { 
	Empty	= 0,
	Healthy	= 10, 
	Tumor	= 20,
	Immune	= 30
};

class Sim {
	public:
		Sim(char* config_file);
		~Sim();

		void diffuse();
		void oxygenate();
		void secrete_toxin();
		void uptake_ox();
		void move_immune();
		void kill_tumor();
		void hypoxia();
		void proliferate();

	private:
		static constexpr size_t size = 100;
		
		/* MATRICES */
		Cell cells[size][size];
		Cell immune[size][size];
		int prolif_cnt[size][size];
		float oxygen[size][size];
		float toxin[size][size];

		Cell temp_cell[size][size];
		int temp_int[size][size];
		float temp_float[size][size];
		
		/* Parameters */
		int nbrhood;
		float diff_rate;
		float ox_supply_level;
		float ox_supply_rate;
		float healthy_ox_rate;
		float immune_ox_rate;
		float tumor_ox_rate;
		float ox_surv_thr;
		float ox_prolif_thr;
		float toxin_secrete_rate;
		float init_immune_ratio;
		int t_cycle;

		/* Random numbers */
		std::default_random_engine gen;
		std::uniform_int_distribution<int> dist_1;

		/* Functions */
		template<class T>
		void read_param(const libconfig::Setting& setting, const char* name, T& var);

		void diffuse_4(float subst[size][size]);
		void diffuse_8(float subst[size][size]);
		void diffuse_6_my(float subst[size][size]);

		void cell_die(size_t i, size_t j);
		void immune_die(size_t i, size_t j);

		friend class Logger;
		std::string logfile;
};
