close all
% clear variables
load('data_num.mat');

figure('units', 'normalized', 'outerposition', [0 0 1 1]);

subplot(3, 1, 1);
plot(num_healthy);
title('Healthy cells');

subplot(3, 1, 2);
plot(num_tumor);
title('Tumor cells');

subplot(3, 1, 3);
plot(num_immune);
title('Immune cells');