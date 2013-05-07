% pc
% dp = (1-pc)*0.2;
% index = find(and(p_arr>(pc+0.01),p_arr<(pc+dp)))
% short_p_arr = p_arr(index);
% x = (short_p_arr - ones(size(short_p_arr))*pc);
% P_arr_m = P_matrix(1,:)
% short_P_matrix = P_arr_m(index);
% figure
% plot(log(x), log(short_P_matrix),'*b')
% title(sprintf('Log P(p,L) for L = %d', L));
% xlabel('log(p-pc)')
% ylabel('log(P(p,L))');
% hold on
% %figure
% 
% p = polyfit(log(x),log(short_P_matrix),1)
% plot(log(x), p(1)*log(x) + p(2), 'r')
% legend('log(P(p,L))', sprintf('linear regression, beta = %f', p(1)));
% 
% x = logspace(0,4);
% y = x;
% figure
% plot(log10(x),log10(y),'.r')
figure

j=4;
y = log10(all_bins(j,:));
x2 = x(y~=-Inf);
y2 = y(y~=-Inf);
plot(x2,y2)
hold all
Pol = polyfit(x2(1:length(x2)-20), y2(1:length(x2)-20),1);
plot(x2, Pol(1)*x2 + Pol(2));
legend(legends{j}, sprintf('slope = %f', Pol(1)));
title('log10(n(s,p)) with linear approximation');
xlabel('log10(s)');
ylabel('log10(n(s,p))');