close all
clear variables
load('data.mat');

n = 5;
dims = [1, 1;
        1, 2;
        1, 3;
        2, 2;
        2, 3;
        2, 3; ];

f = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
colormap('hot');
sz = size(cells, 1);

subplot(dims(n, 1), dims(n, 2), 1);
f_cells = imagesc(cells(:, : ,1), [0, 30]);
axis('equal');
title('Cells');
xlim([0, sz]);

subplot(dims(n, 1), dims(n, 2), 2);
f_immune = imagesc(immune(:, : ,1), [0, 30]);
axis('equal');
title('Immune cells');
xlim([0, sz]);

subplot(dims(n, 1), dims(n, 2), 3);
f_ox = imagesc(nutrient(:, :, 1), [0, 1]);
axis('equal');
title('Nutrient');
xlim([0, sz]);

subplot(dims(n, 1), dims(n, 2), 4);
f_attr = imagesc(attr(:, :, 1), [0, 10]);
axis('equal');
title('Immune attractant');
xlim([0, sz]);

subplot(dims(n, 1), dims(n, 2), 5);
f_ecm = imagesc(ecm_stress(:, :, 1), [0, 1]);
axis('equal');
title('ECM stress');
xlim([0, sz]);

sgtitle(print_time(0));

waitforbuttonpress;

for i = 1:size(cells, 3)
%     pause(0.001);
    drawnow;
    set(f_cells, 'CData', cells(:, :, i));
    set(f_immune, 'CData', immune(:, :, i));
    set(f_ox, 'CData', nutrient(:, :, i));
    set(f_attr, 'CData', attr(:, :, i));
    set(f_ecm, 'CData', ecm_stress(:, :, i));
    
    d = floor(i / 24 / 3);
    h = floor((i - d * 24 * 3) / 3);
    m = (i - d * 24 * 3 - h * 3) * 20;
    sgtitle(print_time(i));
    
%     frame = getframe(f);
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if i == 1
<<<<<<< HEAD
%         imwrite(imind, cm, 'imm_cell0.gif','gif', 'Loopcount', inf, 'DelayTime', 0.1); 
%     else 
%         imwrite(imind, cm, 'imm_cell0.gif','gif','WriteMode','append', 'DelayTime', 0.1); 
=======
%         imwrite(imind,cm, 'sepdiff_imm','gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm, 'sepdiff_imm','gif','WriteMode','append'); 
>>>>>>> parent of 76bb31a (Changed immune cell recruitment process to maintain constant immune cell level)
%     end 
end

function str = print_time(i)
    d = floor(i / 24 / 3);
    h = floor((i - d * 24 * 3) / 3);
    m = (i - d * 24 * 3 - h * 3) * 20;
    
    str = [num2str(d, '%02d'), ' days ', num2str(h, '%02d'), ' hours ', num2str(m, '%02d'), ' minutes'];
end
