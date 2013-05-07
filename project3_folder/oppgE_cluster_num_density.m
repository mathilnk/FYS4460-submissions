close all;
clear all;
pc = 0.59275;
num_exp_per_p = 2000;
num_p = 20;
L = 500;
dp = 0.1;
p_f_above = linspace(pc, pc+dp, num_p);
%p_f_above = logspace(log10(pc), log10(pc+dp), num_p);
%p_f_below = logspace(log10(pc-dp), log10(pc), num_p);
p_f_below = linspace(pc-dp, pc, num_p);
have_plot = false;
legends = {};
%j=19;
j=1;
all_bins = [];
%hist_Areas = zeros(L*L,1);
bin_index = unique(round(logspace(0, log10(L*L),100)));
l = length(bin_index);
bins = zeros(l-1,1);
bins2 = zeros(l-1,1);
j=12;
while j<=num_p,
    hist_Areas = zeros(L*L,1);
    hist_Areas2 = zeros(L*L,1);
    p = p_f_below(j);
    p2 = p_f_above(j);
    for k=1:num_exp_per_p,
        %clear lw;
        r = rand(L,L);
        %[z,lw,num] = gen_bin_array(r,p);
        z = r<p;
        z2 = r<p2;
        
        [lw,num] = bwlabel(z,4);
        [lw2, num2] = bwlabel(z2, 4);
        spanning = find_span_cl_numbers(lw);
        spanning2 = find_span_cl_numbers(lw2);
%         for i=1:length(spanning),
%             [x,y] = find(lw == spanning(i));
%             lw(x,y) = 0;
%         end
        stats = regionprops(lw,'Area');
        stats2 = regionprops(lw2, 'Area');
        areas = cat(1, stats.Area);
        areas2 = cat(1, stats2.Area);
        areas(spanning) = [];
        areas2(spanning2) = [];

        max_a = max(areas);
        min_a = min(areas);
        max_a2 = max(areas2);
        min_a2 = min(areas2);
        [hist_thing,h] = hist(areas,max_a-min_a+1);
        [hist_thing2, h2] = hist(areas2,max_a2-min_a2+1);
        hist_Areas(min_a:max_a) = hist_Areas(min_a:max_a) + hist_thing';
        hist_Areas2(min_a2:max_a2) = hist_Areas2(min_a2:max_a2) + hist_thing2';
    end
    %Areas(j,:) =hist_Areas; 
    if(length(hist_Areas~=0))
        %bin_index = 2.^(0:11)

        for i=1:l-1,
            ds = bin_index(i+1) - bin_index(i);
            %bins(i) = sum((Areas > bin_index(i)).*(Areas < bin_index(i+1)))/ds;
            bins(i) = sum(hist_Areas(bin_index(i):bin_index(i+1)-1))/ds;
            bins2(i) = sum(hist_Areas2(bin_index(i):bin_index(i+1)-1))/ds;
        end
        name = sprintf('bins_%d.mat',j);
        name2 = sprintf('bins2_%d.mat',j);
        save(name, 'bins');
        save(name2, 'bins2');
        
            
% 
%         %all_bins(j,:) = bins;
%         x = log10(bin_index(1:length(bin_index)-1)/L^2)
%         y = log10(bins)'
%         
%         plot(x,y);
%         legends{j} = sprintf('p = %f', p);
%         hold all
%         x2 = x(y~=-Inf);
%         y2 = y(y~=-Inf);
%         Pol = polyfit(x2, y2,1)
%         %plot(x2, Pol(1)*x2 + Pol(2));
%         hold all
%         
%         have_plot = true;
        j= j+1
    else
        'tries again'
    end
    %sum(H)/length(Areas)
    %[H,s] = hist((Areas), max(Areas));
    
    
            
end

for i=1:num_p,
    name = sprintf('bins_%d.mat', i);
    %name_2 = sprintf('bins2_%d.mat', i);
    B = load(name, 'bins');
    %B2 = load(name_2, 'bins2');
    %L = L_list(i);
    x = log10(bin_index(1:length(bin_index)-1));
    y = log10(B.bins/L^2)';
    %y2 = log10(B2.bins/L^2)';

    plot(x,y);
    legends{i} = sprintf('p = %f', p_f_below(i));
    hold all
    x2 = x(y~=-Inf);
    y12 = y(y~=-Inf);
    Pol = polyfit(x2, y12,1)
end
legend(legends);
title('loglog(n(s,p)) p->pc from below');
xlabel('log10(s)');
ylabel('log10(n(s,p))');

figure
for i=1:num_p,
    %name = sprintf('bins_%d.mat', i);
    name_2 = sprintf('bins2_%d.mat', i);
    %B = load(name, 'bins');
    B2 = load(name_2, 'bins2');
    %L = L_list(i);
    x = log10(bin_index(1:length(bin_index)-1));
    %y = log10(B.bins/L^2)';
    y2 = log10(B2.bins2/L^2)';
    plot(x,y2);
    legends{i} = sprintf('p = %f', p_f_above(i));
    hold all
    x2 = x(y~=-Inf);
    y22 = y(y~=-Inf);
    Pol = polyfit(x2, y22,1)
end
legend(legends);
title('loglog(n(s,p)) p->pc from above');
xlabel('log10(s)');
ylabel('log10(n(s,p))'); 




