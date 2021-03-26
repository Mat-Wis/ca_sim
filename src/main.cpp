#include <iostream>
#include "logger.h"
#include "sim.h"

char config_file[] = "../config.cfg";

int main()
{
	Sim sim(config_file);
	Logger logger(sim);

	for(int n = 0; n < 1000; ++n) {
		sim.diffuse();
		sim.oxygenate();
		sim.secrete_toxin();
		sim.uptake_ox();
		logger.log();
	}

	return 0;
}
