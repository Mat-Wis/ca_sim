close all
clear variables

N = 100;

final_size = zeros(1, N);

m = matfile('asym_rand.mat', 'Writable', true);

for i = 1:N
	disp([num2str(n), '.', num2str(i)]);
	[status, cmdout] = system(['../build/ca_sim ', num2str(n)]);

	if ~status
		load('data_num.mat');
		final_size(i) = num_tumor(end);
	end
end

m.d70 = final_size;

figure;
histogram(final_size, 0:5:150)
xmax = 150;
ymax = 100;
xlim([0, xmax]);
ylim([0, ymax]);
xlabel('Number of tumor cells after 200 hours');
ylabel('Number of simulations');
