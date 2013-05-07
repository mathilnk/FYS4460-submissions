function [H, s] = hist_log(data, num_of_bins)
max_data = max(data);
min_data = min(data);
X_axis = logspace(log10(min_data), log10(max_data), num_of_bins);
[H_log, x_log] = hist(log10(data), num_of_bins);
for i=1:num_of_bins-1,
    dx = X_axis(i+1) - X_axis(i);
    H_log(i) = H_log(i)/dx;
end
H_log = H_log/length(data); 
start = 1;
stop = length(H_log) -1;
H = H_log(start:stop);
s = (x_log(1:length(x_log)-1));