#pragma once

#include "matio.h"
#include "sim.h"

class Logger {
	public:
		Logger(Sim& sim);
		~Logger();

		void log();

	private:
		mat_t* file;
		matvar_t* cells_var;
		matvar_t* immune_var;
		matvar_t* oxygen_var;
		matvar_t* toxin_var;

		void saveParam(int* var, const char* name);
		void saveParam(float* var, const char* name);
};
