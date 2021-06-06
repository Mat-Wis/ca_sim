close all
clear variables

N = 100;

final_size = zeros(1, N);

for i = 1:N
    i
    [status, cmdout] = system('../build/ca_sim');
    
    if ~status
        load('data_num.mat');
        final_size(i) = num_tumor(end);
    end
end

histogram(final_size);
xlabel('Number of tumor cells after 200 hours');
ylabel('Number of simulations');

m = matfile('immune_speed_i10.mat', 'Writable', true);
m.dt15 = final_size;