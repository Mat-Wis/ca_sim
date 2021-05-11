#pragma once

#include "matio.h"
#include "sim.h"

class Logger {
	public:
		Logger(Sim& sim);
		~Logger();

		void log_num();
		void log_mat();

	private:
		mat_t* file;
		matvar_t* cells_var;
		matvar_t* immune_var;
		matvar_t* nutrient_var;
		matvar_t* attr_var;
		matvar_t* ecm_var;
		matvar_t* num_healthy_var;
		matvar_t* num_tumor_var;
		matvar_t* num_immune_var;

		void saveParam(int* var, const char* name);
		void saveParam(float* var, const char* name);
};
