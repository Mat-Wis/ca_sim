#include <iostream>
#include "logger.h"
#include "sim.h"

char config_file[] = "../config.cfg";

int main()
{
	Sim sim(config_file);
	Logger logger(sim);

	for(int n = 0; n < 1000; ++n) {
		sim.oxygenate();
		sim.diffuse();
		sim.secrete_toxin();
		sim.uptake_ox();
		sim.move_immune();
		sim.kill_tumor();
		logger.log();
	}

	return 0;
}
