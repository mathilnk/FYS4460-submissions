clear all;
close all;
num_p = 30;
L = 300;
tau = 1.87;
pc = 0.59275;
bin_index = unique(round(logspace(0, log10(L*L),100)));
l = length(bin_index);
s = bin_index(1:length(bin_index)-1);
n_s_pc =s.^(-tau);
%tolerance = 0.0001;
s_ksi1 = zeros(num_p,1);
s_ksi2 = zeros(num_p,1);
dp = 0.1;
dp = 0.1;
p_f_above = linspace(pc, pc+dp, num_p);
%p_f_above = logspace(log10(pc), log10(pc+dp), num_p);
%p_f_below = logspace(log10(pc-dp), log10(pc), num_p);
p_f_below = linspace(pc-dp, pc, num_p);

%figure
for i=1:30,
    name = sprintf('bins_%d.mat',i);
    name2 = sprintf('bins2_%d.mat',i);
    B_below = load(name, 'bins');
    B_above = load(name2, 'bins2');
    n_sp_below = B_below.bins/L^2;
    n_sp_above = B_above.bins2/L^2;
    
    f_below = n_sp_below./n_s_pc';
    f_above = n_sp_above./n_s_pc';
    [max_f, max_f_i] = max(f_below(5:length(f_below)));
    [max_f_above, max_f_ai] = max(f_above(5:length(f_above)));
    
    half_b = ones(size(f_below))*max_f*0.5;
    half_a = ones(size(f_above))*max_f_above*0.5;
    [X_b, Y_b] = intersections(s, f_below, s, half_b)
    s_xi_array = X_b(X_b>s(max_f_i));
    s_ksi1(i) = s_xi_array(1)
    
    [X_a, Y_a] = intersections(s, f_above, s, half_a)
    s_xi_array_a = X_a(X_a>s(max_f_ai));
    s_ksi2(i) = s_xi_array_a(1)
    plot(log10(s), log10(f_below));
    hold all
    plot(log10(s), log10(half_b));
end
figure
p_tot = [p_f_below, p_f_above];
s_xi = [s_ksi1', s_ksi2'];
plot(p_tot, s_xi);
x = log10(abs(p_tot - pc));
y = log10(s_xi);
figure
plot(x,y, '.');
x_a = log10(abs(p_f_above-ones(size(p_f_above))*pc));
y_a = log10(s_ksi2);
figure
plot(x_a, y_a)
hold on
x_b = log10(abs(p_f_below-ones(size(p_f_below))*pc));
y_b = log10(s_ksi1);
plot(x_b, y_b, '.');
x_bb = x_b(x_b~= -Inf);
y_bb = y_b(x_b~=-Inf);
Pb = polyfit(x_bb,y_bb',1)
x_aa = x_a(x_a~= -Inf);
y_aa = y_a(x_a~=-Inf);
Pa = polyfit(x_aa, y_aa', 1)
sig = -1/((Pb(1) + Pa(1))*0.5)
%sig = 0.5;
figure
plot(p_tot, s_xi, '.')
hold all
%figure
plot(p_tot, abs(p_tot- pc).^(-1/sig));
title('Estimation of s_{xi}')
xlabel('p');
ylabel('s_{xi}');
legend('Experiment', sprintf('Estimated, sigma = %f', sig));

% p_tot = [p_f_below, p_f_above];%     s_ksi1(i) = s(R_below);
%     s_ksi2(i) = s(R_above);
%     plot(log10(s), log10(ones(length(s),1)*0.5*max_f));
%     hold all
%     plot(log10(s),log(f_below));
%     
%     %loglog(s,f_below);
%     
% end
% %figure
% p_tot = [p_f_below, p_f_above];
% s_ksi = [s_ksi1', s_ksi2'];
% %plot(p_tot, s_ksi, '.');
% %plot((-pc*ones(size(p_f_below))+p_f_below), s_ksi1, '.');
% %hold on
% %plot((p_f_above- pc*ones(size(p_f_above))),s_ksi2, '.');
% x_a0 = log10(abs(-p_f_above + pc*ones(size(p_f_above))));
% x_b0 = log10(abs(p_f_below - pc*ones(size(p_f_below))));
% y_b0 = log10(s_ksi1);
% y_a0 = log10(s_ksi2);
% figure
% plot(x_b0,y_b0, '.');
% x_b = x_b0(x_b0~=-Inf);%(x_b0>-1.5);
% y_b = y_b0(x_b0~=-Inf)';%(x_b0>-1.5)';
% Pol_below = polyfit(x_b,y_b,1)
% hold all 
% plot(x_b, Pol_below(1)*x_b + Pol_below(2));
% title('log10(s_{ksi}) for p < p_c')
% xlabel('log10(|p-pc|)')
% ylabel('log10(s_{ksi})');
% sig = -1.0/Pol_below(1);
% legend('s_{ksi}', sprintf('linear app: sig = %f', sig));
% figure
% %plot(x_a0, y_a0, '.');
% x_a = x_a0(x_a0~=-Inf);
% y_a = y_a0(x_a0~=-Inf)';
% plot(x_a, y_a, '.');
% Pol_abv = polyfit(x_a, y_a, 1)
% hold all 
% plot(x_a, Pol_abv(1)*x_a + Pol_abv(2))
% title('log10(s_{ksi}) for p > p_c')
% xlabel('log10(|p-pc|)')
% ylabel('log10(s_{ksi})');
% sig = -1.0/Pol_abv(1);
% sig = 0.4
% legend('s_{ksi}', sprintf('linear app: sig = %f', sig));
% figure
% plot(p_tot, s_ksi, '.');
% hold all
% s_ksi_theory = abs(p_tot - pc*ones(size(p_tot))).^(-1.0/sig);
% plot(p_tot, s_ksi_theory);
% title('s_{ksi}');
% xlabel('p');
% ylabel('s_{ksi}');
% legend('Experimental', 'Theoretical');

% s_ksi = [s_ksi1', s_ksi2'];
% %plot(p_tot, s_ksi, '.');
% %plot((-pc*ones(size(p_f_below))+p_f_below), s_ksi1, '.');
% %hold on
% %plot((p_f_above- pc*ones(size(p_f_above))),s_ksi2, '.');
% x_a0 = log10(abs(-p_f_above + pc*ones(size(p_f_above))));
% x_b0 = log10(abs(p_f_below - pc*ones(size(p_f_below))));
% y_b0 = log10(s_ksi1);
% y_a0 = log10(s_ksi2);
% figure
% plot(x_b0,y_b0, '.');
% x_b = x_b0(x_b0~=-Inf);%(x_b0>-1.5);
% y_b = y_b0(x_b0~=-Inf)';%(x_b0>-1.5)';
% Pol_below = polyfit(x_b,y_b,1)
% hold all 
% plot(x_b, Pol_below(1)*x_b + Pol_below(2));
% title('log10(s_{ksi}) for p < p_c')
% xlabel('log10(|p-pc|)')
% ylabel('log10(s_{ksi})');
% sig = -1.0/Pol_below(1);
% legend('s_{ksi}', sprintf('linear app: sig = %f', sig));
% figure
% %plot(x_a0, y_a0, '.');
% x_a = x_a0(x_a0~=-Inf);
% y_a = y_a0(x_a0~=-Inf)';
% plot(x_a, y_a, '.');
% Pol_abv = polyfit(x_a, y_a, 1)
% hold all 
% plot(x_a, Pol_abv(1)*x_a + Pol_abv(2))
% title('log10(s_{ksi}) for p > p_c')
% xlabel('log10(|p-pc|)')
% ylabel('log10(s_{ksi})');
% sig = -1.0/Pol_abv(1);
% sig = 0.4
% legend('s_{ksi}', sprintf('linear app: sig = %f', sig));
% figure
% plot(p_tot, s_ksi, '.');
% hold all
% s_ksi_theory = abs(p_tot - pc*ones(size(p_tot))).^(-1.0/sig);
% plot(p_tot, s_ksi_theory);
% title('s_{ksi}');
% xlabel('p');
% ylabel('s_{ksi}');
% legend('Experimental', 'Theoretical');

