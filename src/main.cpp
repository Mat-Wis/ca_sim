#include <iostream>
#include "logger.h"
#include "sim.h"

char config_file[] = "../config.cfg";

int main()
{
	Sim sim(config_file);
	Logger logger(sim);

	int n_steps = static_cast<int>(sim.sim_time / sim.dt * 3600.0f);
	int log_step = static_cast<int>(sim.log_step * 3600.0f);

	for(int n = 0; n < n_steps; ++n) {
		sim.diffuse();
		sim.uptake_ox();
		sim.oxygenate();
		sim.secrete_toxin();
		sim.move_immune();
		sim.recruit_immune();
		sim.kill_tumor();
		sim.kill_immune();
		sim.kill_healthy();
		sim.hypoxia();
		sim.proliferate();
		sim.count_cells();

		logger.log_num();

		if(n % log_step == 0) {
			logger.log_mat();
		}

		if(n % 1000 == 0) {
			std::cout << "n = " << n << std::endl;
		}
	}

	return 0;
}
