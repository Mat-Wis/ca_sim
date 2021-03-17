close all
clear variables
load('data.mat');

colormap('hot');

f = imagesc(oxygen(:,:,1), [0, max(oxygen(:))]);
axis('equal');

for i = 1:1000
    pause(0.01);
    set(f, 'CData', oxygen(:, :, i));
end