clear all;
close all;
pc = 0.59275;
num_exp_per_L =5000;
k = [4 5 6 7 8 9];
%k=[5,6];
L_list = 2.^k;
%all_bins = zeros(length(L_list), L_list(length(L_list)));

% for i=1:length(k),
%     i
%     L = L_list(i);
%     bin_index = unique(round(logspace(0, log10(L*L),100)));
%     l = length(bin_index);
%     bins = zeros(l-1,1);
%     hist_Areas = zeros(L*L, 1);
%     for j=1:num_exp_per_L,
%         r = rand(L,L);
%         z = r<pc;
%         [lw, num] = bwlabel(z,4);
%         spanning = find_span_cl_numbers(lw);
%         stats = regionprops(lw, 'Area');
%         areas = cat(1,stats.Area);
%         areas(spanning) = [];
%         max_a = max(areas);
%         min_a = min(areas);
%         [hist_thing,h] = hist(areas,max_a-min_a+1);
%         size(hist_Areas(min_a:max_a));
%         size(hist_thing');
%         if(length(hist_thing~=0))     
%             hist_Areas(min_a:max_a) = hist_Areas(min_a:max_a) + hist_thing';
%         end
% 
%     end
%    
%     for p=1:l-1,
%         ds = bin_index(p+1) - bin_index(p);
%         %bins(i) = sum((Areas > bin_index(i)).*(Areas < bin_index(i+1)))/ds;
%         bins(p) = sum(hist_Areas(bin_index(p):bin_index(p+1)-1))/ds;
%     end
%     %all_bins(i,:) = bins;
%     name = sprintf('bins_%d.mat',i);
%     name_2 = sprintf('bin_index_%d.mat', i);
%     save(name, 'bins');
%     save(name_2, 'bin_index');
%     
% end
legends = {};

for i=1:length(k),
    name = sprintf('bins_%d.mat', i);
    name_2 = sprintf('bin_index_%d.mat', i);
    Bins = load(name, 'bins');
    Bin_index = load(name_2, 'bin_index');
    L = L_list(i);
    x = log10(Bin_index.bin_index(1:length(Bin_index.bin_index)-1));
    y = log10(Bins.bins/L^2)';

    plot(x,y);
    legends{i} = sprintf('L = %f', L_list(i));
    hold all
    x2 = x(y~=-Inf);
    y2 = y(y~=-Inf);
    Pol = polyfit(x2, y2,1)
end
legend(legends);
title('loglog(n(s,pc)) for various L');
xlabel('log10(s)');
ylabel('log10(n(s,pc))');