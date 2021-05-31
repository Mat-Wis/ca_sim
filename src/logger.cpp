#include "logger.h"
#include "matio.h"
#include "sim.h"

Logger::Logger(Sim& sim) {
	mat_file = Mat_CreateVer(sim.mat_file.c_str(), NULL, MAT_FT_MAT73);
	num_file = Mat_CreateVer(sim.num_file.c_str(), NULL, MAT_FT_MAT73);
	
	size_t dims[3] = {sim.size, sim.size, 1};
	size_t dims_1[2] = {1, 1};

	cells_var = Mat_VarCreate("cells", MAT_C_INT32, MAT_T_INT32, 3, dims, &(sim.cells), MAT_F_DONT_COPY_DATA);
	immune_var = Mat_VarCreate("immune", MAT_C_INT32, MAT_T_INT32, 3, dims, &(sim.immune), MAT_F_DONT_COPY_DATA);
	nutrient_var = Mat_VarCreate("nutrient", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, &(sim.nutrient), MAT_F_DONT_COPY_DATA);
	attr_var = Mat_VarCreate("attr", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, &(sim.attr), MAT_F_DONT_COPY_DATA);
	ecm_var = Mat_VarCreate("ecm_stress", MAT_C_SINGLE, MAT_T_SINGLE, 3, dims, &(sim.ecm_stress), MAT_F_DONT_COPY_DATA);
	num_healthy_var = Mat_VarCreate("num_healthy", MAT_C_INT32, MAT_T_INT32, 2, dims_1, &(sim.num_healthy), MAT_F_DONT_COPY_DATA);
	num_tumor_var = Mat_VarCreate("num_tumor", MAT_C_INT32, MAT_T_INT32, 2, dims_1, &(sim.num_tumor), MAT_F_DONT_COPY_DATA);
	num_deadtumor_var = Mat_VarCreate("num_deadtumor", MAT_C_INT32, MAT_T_INT32, 2, dims_1, &(sim.num_deadtumor), MAT_F_DONT_COPY_DATA);
	num_immune_var = Mat_VarCreate("num_immune", MAT_C_INT32, MAT_T_INT32, 2, dims_1, &(sim.num_immune), MAT_F_DONT_COPY_DATA);

	saveParam(&(sim.sim_time), "sim_time");
	saveParam(&(sim.dt), "dt");
	saveParam(&(sim.log_step), "log_step");

	saveParam(&(sim.alpha2), "alpha2");
	saveParam(&(sim.lambda), "lambda");
	saveParam(&(sim.beta2), "beta2");
	saveParam(&(sim.nutr_surv_thr), "nutr_surv_thr");
	saveParam(&(sim.nutr_prolif_thr), "nutr_prolif_thr");
	saveParam(&(sim.stress_thr), "stress_thr");
	saveParam(&(sim.imm_rnd), "imm_rnd");
	saveParam(&(sim.vessel_num), "vessel_num");

	saveParam(&(sim.t_cycle), "t_cycle");
	saveParam(&(sim.init_immune_ratio), "init_immune_ratio");
	saveParam(&(sim.kill_limit), "kill_limit");
	saveParam(&(sim.life_limit), "life_limit");
}

Logger::~Logger() {
	Mat_VarFree(cells_var);
	Mat_VarFree(nutrient_var);
	Mat_VarFree(attr_var);
	Mat_VarFree(ecm_var);
	Mat_VarFree(num_healthy_var);
	Mat_VarFree(num_tumor_var);
	Mat_VarFree(num_deadtumor_var);
	Mat_VarFree(num_immune_var);

	Mat_Close(mat_file);
	Mat_Close(num_file);
}

void Logger::saveParam(int* var, const char* name) {
	size_t param_dim[2] = {1, 1};
	matvar_t* param_var = Mat_VarCreate(name, MAT_C_INT32, MAT_T_INT32, 2, param_dim, var, 0);
	Mat_VarWrite(mat_file, param_var, MAT_COMPRESSION_NONE);
	Mat_VarWrite(num_file, param_var, MAT_COMPRESSION_NONE);
	Mat_VarFree(param_var);
}

void Logger::saveParam(float* var, const char* name) {
	size_t param_dim[2] = {1, 1};
	matvar_t* param_var = Mat_VarCreate(name, MAT_C_SINGLE, MAT_T_SINGLE, 2, param_dim, var, 0);
	Mat_VarWrite(mat_file, param_var, MAT_COMPRESSION_NONE);
	Mat_VarWrite(num_file, param_var, MAT_COMPRESSION_NONE);
	Mat_VarFree(param_var);
}

void Logger::log_num() {
	Mat_VarWriteAppend(num_file, num_healthy_var, MAT_COMPRESSION_NONE, 2);
	Mat_VarWriteAppend(num_file, num_tumor_var, MAT_COMPRESSION_NONE, 2);
	Mat_VarWriteAppend(num_file, num_deadtumor_var, MAT_COMPRESSION_NONE, 2);
	Mat_VarWriteAppend(num_file, num_immune_var, MAT_COMPRESSION_NONE, 2);
}

void Logger::log_mat() {
	Mat_VarWriteAppend(mat_file, cells_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(mat_file, immune_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(mat_file, nutrient_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(mat_file, attr_var, MAT_COMPRESSION_NONE, 3);
	Mat_VarWriteAppend(mat_file, ecm_var, MAT_COMPRESSION_NONE, 3);
}
