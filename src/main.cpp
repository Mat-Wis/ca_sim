#include <iostream>
#include "logger.h"
#include "sim.h"

char config_file[] = "../config.cfg";

int main()
{
	Sim sim(config_file);
	Logger logger(sim);

	for(int n = 0; n < sim.n_steps; ++n) {
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
		logger.log();

		if(n % 1000 == 0) {
			std::cout << "n = " << n << std::endl;
		}
	}

	return 0;
}
