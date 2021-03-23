close all
clear variables
load('data.mat');

colormap('hot');

f = imagesc(oxygen(:,:,1), [0, 1]);
axis('equal');

for i = 1:1000
    pause(0.001);
    set(f, 'CData', oxygen(:, :, i));
end