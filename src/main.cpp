#include <iostream>
#include "logger.h"
#include "sim.h"

char config_file[] = "../config.cfg";

int main(int argc, char** argv)
{
	Sim sim(config_file);

	Logger logger(sim);

	for(int n = 0; n < sim.n_steps; ++n) {
		sim.damage_ecm();
		sim.diffuse();

		sim.move_immune();
		sim.recruit_immune();
		sim.kill_tumor();
		sim.kill_immune();
		sim.kill_healthy();
		sim.proliferate();
		sim.count_cells();

		logger.log_num();

		if(sim.log_mat && n % sim.log_step == 0) {
			logger.log_mat();
		}

		if(n % 100 == 0) {
			std::cout << "n = " << n << std::endl;
		}

		if(sim.tumor_killed()) {
			break;
		}
	}

	return 0;
}
