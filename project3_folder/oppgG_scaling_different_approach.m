clear all;
close all;
num_p = 20;
L = 500;
tau = 1.87;
pc = 0.59275;
dp = 0.1;
p_f_above = linspace(pc, pc+dp, num_p);
%p_f_above = logspace(log10(pc), log10(pc+dp), num_p);
%p_f_below = logspace(log10(pc-dp), log10(pc), num_p);
p_f_below = linspace(pc-dp, pc, num_p);
p_tot = [p_f_below, p_f_above];
sig = 0.45;
bin_index = unique(round(logspace(0, log10(L*L),100)));
l = length(bin_index);
s = bin_index(1:length(bin_index)-1);
leg = {};
counter = 1;
for(i=1:num_p-1),
    name = sprintf('bins_%d.mat', i);
    name2 = sprintf('bins2_%d.mat', i);
    B_b = load(name, 'bins');
    B_a = load(name2, 'bins2');
    n_sp_b = B_b.bins/L^2;
    n_sp_a = B_a.bins2/L^2;
    n_sp = [n_sp_b, n_sp_a];
    p = p_tot(i);
    x = (p - pc)^(1/sig).*s;
    y = s.^tau.*n_sp_b';
    plot(log10(x(y~=0)),log10(y(y~=0)))
    hold all
    leg{counter} = sprintf('p=%g',p);
    counter  = counter+1;
end
legend(leg);
title_string = sprintf('${F(s(p-pc)}^{1/\\sigma})$ for $\\sigma$ = %g',sig); 
title(title_string, 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$s|p-pc|^{1/\sigma}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$s^\tau (s,p)$', 'Interpreter', 'latex', 'FontSize', 12);


    

