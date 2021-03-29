#include "logger.h"
#include "matio.h"
#include "sim.h"

Logger::Logger(Sim& sim) {
	file = Mat_CreateVer(sim.logfile.c_str(), NULL, MAT_FT_MAT73);
	
	size_t dims[3] = {sim.size, sim.size, 1};

	cells_var = Mat_VarCreate("cells", MAT_C_INT32, MAT_T_INT32, 3, dims, &(sim.cells), MAT_F_DONT_COPY_DATA);
	immune_var = Mat_VarCreate("immune", MAT_C_INT32, MAT_T_INT32, 3, dims, &(sim.immune), MAT_F_DONT_COPY_DATA);
	oxygen_var = Mat_VarCreate("oxygen", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, &(sim.oxygen), MAT_F_DONT_COPY_DATA);
	toxin_var = Mat_VarCreate("toxin", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, &(sim.toxin), MAT_F_DONT_COPY_DATA);

	saveParam(&(sim.diff_rate), "diff_rate");
	saveParam(&(sim.ox_supply_level), "ox_supply_level");
	saveParam(&(sim.ox_supply_rate), "ox_supply_rate");
	saveParam(&(sim.healthy_ox_rate), "healthy_ox_rate");
	saveParam(&(sim.immune_ox_rate), "immune_ox_rate");
	saveParam(&(sim.tumor_ox_rate), "tumor_ox_rate");
	saveParam(&(sim.toxin_secrete_rate), "toxin_secrete_rate");
	saveParam(&(sim.toxin_thr), "toxin_thr");
	saveParam(&(sim.init_immune_ratio), "init_immune_ratio");
	saveParam(&(sim.t_cycle), "t_cycle");
	saveParam(&(sim.kill_limit), "kill_limit");
}

Logger::~Logger() {
	Mat_VarFree(cells_var);
	Mat_VarFree(oxygen_var);
	Mat_VarFree(toxin_var);

	Mat_Close(file);
}

void Logger::saveParam(int* var, const char* name) {
	size_t param_dim[2] = {1, 1};
	matvar_t* param_var = Mat_VarCreate(name, MAT_C_INT32, MAT_T_INT32, 2, param_dim, var, 0);
	Mat_VarWrite(file, param_var, MAT_COMPRESSION_NONE);
	Mat_VarFree(param_var);
}

void Logger::saveParam(float* var, const char* name) {
	size_t param_dim[2] = {1, 1};
	matvar_t* param_var = Mat_VarCreate(name, MAT_C_SINGLE, MAT_T_SINGLE, 2, param_dim, var, 0);
	Mat_VarWrite(file, param_var, MAT_COMPRESSION_NONE);
	Mat_VarFree(param_var);
}

void Logger::log() {
	Mat_VarWriteAppend(file, cells_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(file, immune_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(file, oxygen_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(file, toxin_var, MAT_COMPRESSION_NONE, 3);
}
