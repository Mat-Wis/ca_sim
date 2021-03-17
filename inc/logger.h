#pragma once

#include "matio.h"

class Logger {
	public:
		Logger(void* cells, void* oxygen, void* nutrient, void* toxin, size_t size, char* filename);
		~Logger();

		void log();

	private:
		mat_t* file;
		matvar_t* cells_var;
		matvar_t* oxygen_var;
		matvar_t* nutrient_var;
		matvar_t* toxin_var;

};
