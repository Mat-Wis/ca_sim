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

subplot(dims(n, 1), dims(n, 2), 1);
f_cells = imagesc(cells(:, : ,1), [0, 30]);
axis('equal');
title('Cells');

subplot(dims(n, 1), dims(n, 2), 2);
f_immune = imagesc(immune(:, : ,1), [0, 30]);
axis('equal');
title('Immune cells');

subplot(dims(n, 1), dims(n, 2), 3);
f_ox = imagesc(oxygen(:, :, 1), [0, 1]);
axis('equal');
title('Oxygen');

subplot(dims(n, 1), dims(n, 2), 4);
f_attr = imagesc(attr(:, :, 1), [0, 1]);
axis('equal');
title('Immune attractant');

subplot(dims(n, 1), dims(n, 2), 5);
f_ecm = imagesc(ecm_stress(:, :, 1), [0, 1]);
axis('equal');
title('ECM stress');

sgtitle('n = 0');

waitforbuttonpress;

for i = 1:size(cells, 3)
    pause(0.001);
    set(f_cells, 'CData', cells(:, :, i));
    set(f_immune, 'CData', immune(:, :, i));
    set(f_ox, 'CData', oxygen(:, :, i));
    set(f_attr, 'CData', attr(:, :, i));
    set(f_ecm, 'CData', ecm_stress(:, :, i));
    sgtitle(['n = ', num2str(i)]);
    
%     frame = getframe(f);
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if i == 1
%         imwrite(imind,cm, 'sepdiff_imm','gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm, 'sepdiff_imm','gif','WriteMode','append'); 
%     end 
end
