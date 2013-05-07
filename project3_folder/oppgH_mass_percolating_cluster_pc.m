clear all;
close all;
pc = 0.59275;
num_exp_per_L =3000;
k = [4 5 6 7 8 9 10];
%k=[5,6, 7];
L_list = 2.^k;
%all_bins = zeros(length(L_list), L_list(length(L_list)));
Mass = zeros(length(k),1);
% M = load('last_mass.mat', 'Mass');
% Mass = M.Mass;
for i=1:length(k),
    i
    L = L_list(i);
    bin_index = unique(round(logspace(0, log10(L*L),100)));
    l = length(bin_index);
    bins = zeros(l-1,1);
    hist_Areas = zeros(L*L, 1);
    for j=1:num_exp_per_L,
        r = rand(L,L);
        z = r<pc;
        [lw, num] = bwlabel(z,4);
        spanning = find_span_cl_numbers(lw);
        stats = regionprops(lw, 'Area');
        areas = cat(1,stats.Area);
        Mass(i) = Mass(i) +  sum(areas(spanning));
    end
    %name = sprintf('bins_%d.mat',i);
    %name_2 = sprintf('bin_index_%d.mat', i);
    %save(name, 'bins');
    %save(name_2, 'bin_index');
    name_mass = sprintf('last_mass.mat');
    save(name_mass, 'Mass');    
end
figure
plot(log10(L_list),log10(Mass), '.');
Pol = polyfit(log10(L_list), log10(Mass)', 1);
hold all 
plot(log10(L_list), Pol(1)*log10(L_list) + Pol(2));
legend('Experiment', sprintf('Linear approx: D = %g', Pol(1)));
title('Mass of the percolating cluster at pc');
xlabel('Log10(L)');
ylabel('Log10(Mass)');

