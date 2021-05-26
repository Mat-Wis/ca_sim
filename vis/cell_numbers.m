close all
clear variables
load('data.mat');

figure('units', 'normalized', 'outerposition', [0 0 1 1]);

time = (0:size(cells, 3)-1) * dt / 60 / 24;

subplot(3, 1, 1);
plot(time, num_healthy);
title('Healthy cells');
xlim([time(1), time(end)]);

subplot(3, 1, 2);
plot(time, num_tumor);
title('Tumor cells');
xlim([time(1), time(end)]);

subplot(3, 1, 3);
plot(time, num_immune);
title('Immune cells');
xlim([time(1), time(end)]);