#include <iostream>
#include "logger.h"
#include "sim.h"

char config_file[] = "../config.cfg";

int main()
{
	Sim sim(config_file);
	Logger logger(sim);

	for(int n = 0; n < sim.n_steps; ++n) {
		sim.uptake_ox();
		sim.oxygenate();
		sim.diffuse();
		sim.secrete_toxin();
		sim.move_immune();
		sim.recruitImmune();
		sim.kill_tumor();
		sim.kill_immune();
		sim.kill_healthy();
		sim.hypoxia();
		sim.proliferate();
		logger.log();
	}

	return 0;
}
