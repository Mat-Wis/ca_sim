#include "matio.h"
#include "logger.h"

Logger::Logger(void* cells, void* oxygen, void* nutrient, void* toxin, size_t size, char* filename) {
	file = Mat_CreateVer(filename, NULL, MAT_FT_MAT73);

	size_t dims[3] = {size, size, 1};

	cells_var = Mat_VarCreate("cells", MAT_C_INT32, MAT_T_INT32, 3, dims, cells, MAT_F_DONT_COPY_DATA);
	oxygen_var = Mat_VarCreate("oxygen", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, oxygen, MAT_F_DONT_COPY_DATA);
	nutrient_var = Mat_VarCreate("nutrient", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, nutrient, MAT_F_DONT_COPY_DATA);
	toxin_var = Mat_VarCreate("toxin", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, toxin, MAT_F_DONT_COPY_DATA);
}

Logger::~Logger() {
	Mat_VarFree(cells_var);
	Mat_VarFree(oxygen_var);
	Mat_VarFree(nutrient_var);
	Mat_VarFree(toxin_var);

	Mat_Close(file);
}

void Logger::log() {
	Mat_VarWriteAppend(file, cells_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(file, oxygen_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(file, nutrient_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(file, toxin_var, MAT_COMPRESSION_NONE, 3);
}
