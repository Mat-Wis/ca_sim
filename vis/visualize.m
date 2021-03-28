close all
clear variables
load('data.mat');

n = 4;
dims = [1, 1;
        1, 2;
        1, 3;
        2, 2;
        2, 3;
        2, 3; ];

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
colormap('hot');

subplot(dims(n, 1), dims(n, 2), 1);
f_cells = imagesc(cells(:, : ,1), [0, 30]);
axis('equal');

subplot(dims(n, 1), dims(n, 2), 2);
f_immune = imagesc(immune(:, : ,1), [0, 30]);
axis('equal');

subplot(dims(n, 1), dims(n, 2), 3);
f_ox = imagesc(oxygen(:, :, 1), [0, 1]);
axis('equal');

subplot(dims(n, 1), dims(n, 2), 4);
f_tox = imagesc(toxin(:, :, 1), [0, 1]);
axis('equal');

waitforbuttonpress;

for i = 1:size(cells, 3)
    pause(0.001);
    set(f_cells, 'CData', cells(:, :, i));
    set(f_immune, 'CData', immune(:, :, i));
    set(f_ox, 'CData', oxygen(:, :, i));
    set(f_tox, 'CData', toxin(:, :, i));
end