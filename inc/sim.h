#pragma once

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <random>
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

		void* get_cells_ptr();
		void* get_oxygen_ptr();
		void* get_toxin_ptr();

		void diffuse();
		void oxygenate();
		void uptake_ox();

		static constexpr size_t size = 100;

	private:
		/* MATRICES */
		Cell cells[size][size];
		float oxygen[size][size];
		float toxin[size][size];

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

		/* Functions */
		template<class T>
		void read_param(const libconfig::Setting& setting, const char* name, T& var);

		void diffuse_4(float subst[size][size]);
		void diffuse_8(float subst[size][size]);
		void diffuse_6_my(float subst[size][size]);
		void diffuse_6_ext(float subst[size][size]);

		friend class Logger;
		std::string logfile;
};
