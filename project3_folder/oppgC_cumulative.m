close all;
z = rand(1e6,1).^(-2);
% %hist(z)
% 
% 
% 
% figure
% %make a logarithmic binning
max_z = max(z);
min_z = 1;
N = 1001;
% z_log =log10(z);
% hist(z_log, N)
% figure
X = logspace(log10(min_z), log10(max_z),N)
% 
% [H, x_log] = hist(z_log, N);
% for i=1:length(X)-1,
%     dz = X(i+1) - X(i);
%     H(i) = H(i)/dz;
% end
% plot(x_log(1:length(x_log)-1),H(1:length(H)-1)/length(z));
% figure
% [H_test, x_test] = hist_log(z, N);
% plot(x_test, H_test);
% title('test histlog');

%version 2:
%find the cumulative distribution with logarithmic binning?
Y = linspace(min_z, 100, N);
cumulative = zeros(N,1);
for i=1:N,
    num = length(find(z<X(i)));
    cumulative(i) = cumulative(i) + num;
end
plot(log10(X),log10(1-cumulative/length(z)));
title('log of the cumulative distribution with logarithmic binning');
xlabel('log10(z) values');
ylabel('log10(P(Z>z))');
size(X)
size(cumulative)
log10(1-cumulative/length(z))
poly = polyfit(transpose(log10(X)), log10(1-cumulative/length(z)),1);
alph = poly(1) -1;
hold on
plot(log10(X), log10(X)*poly(1) + poly(2), 'r');
legend('log10(P(x))', sprintf('Linear approx, alpha = %f', alph));


max_z
figure

%find f by finding the derivative
f = zeros(N-1, 1);
%dz = (X(length(X)) - X(1))/length(X);
for i=1:N-1,
    dz = X(i+1) - X(i);
    f(i) = (cumulative(i+1) - cumulative(i))/dz;
end
plot(log10(X(1:length(X)-1)), (f/length(z)));
title('f_Z(z) with logarithmic binning')
xlabel('log10(z) values')
ylabel('f_Z(z)');
fz = -(alph +1).*X.^(alph);
hold on 
plot(log10(X), fz, 'r');
legend('f found by dP/dz', 'f = -(alpha+1)z^{alpha}')








