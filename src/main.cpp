#include <iostream>
#include <cstring>
#include <random>
#include "logger.h"


constexpr size_t size = 100;
char filename[] = "../vis/data.mat";

int cells[size][size];
float oxygen[size][size];
float nutrient[size][size];
float toxin[size][size];
float temp[size][size];

int main()
{
	Logger logger((void *)&cells, (void *)&oxygen, (void *)&nutrient, (void *)&toxin, size, filename);

	return 0;
}
