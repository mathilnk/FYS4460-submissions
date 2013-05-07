% L = 200;
% lx = L;
% ly = L;
% pc = 0.59274;
% num_p = 5;
% dp = 0.1;
% p_list = linspace(pc, pc+dp, num_p);
% num_exp_per_N = 500;
% %k=1:30;
% %N_list =2.^k; 
% num_N = 20;
% N_list = unique(round(logspace(log10(1), log10(1e6), num_N)));
% %N_list = linspace(1, 1e4, num_N);
% for i=1:num_p,
%     p = p_list(i);
%     i
%     num_p
%     R_sq = zeros(size(N_list));
%     for j=1:length(N_list),
%         nstep = N_list(j);
%         nnstep = nstep+1;
%         for exp=1:num_exp_per_N,
%             ncount = 0;
%             perc = [];
%             while (size(perc ,1)==0)
%                     ncount = ncount + 1;
%                     if (ncount >1000)
%                         return
%                     end
%                     z=rand(lx,ly)<p;
%                     [lw,num]=bwlabel(z,4);
%                     perc_x = intersect(lw(1,:),lw(lx ,:));
%                     perc = perc_x(find(perc_x >0));
%             end
%             zz = lw == perc(1);
%             rz = 1.0*zz;
%             n = 1;
%             while (n<=1)
%                 r = rand(nnstep ,1);
%                 [w,n] = percwalk(rz,r,0);
%             end
%             x = w(1,:);
%             y = w(2,:);
%             
%             R_sq(j) = R_sq(j) + (x(1)-x(length(x)))^2 + (y(1)-y(length(y)))^2;
%         end
%     end
%     R_sq = R_sq./num_exp_per_N;
%     save(sprintf('R_square_%d.mat', i), 'R_sq');
% end

close all;
legends = {};
num_p
start = 1;
stop = num_p;
counter = 1;
list = [num_p];
%list =  start:stop;
for i=list,
    R_load = load(sprintf('R_square_%d.mat', i), 'R_sq');
    plot(log10(N_list), log10(R_load.R_sq));
    hold all
    %Pol = polyfit(log10(N_list), log10(R_load.R_sq), 1);
    %plot(log10(N_list), log10(N_list)*Pol(1) + Pol(2));
    legends{counter} = sprintf('p = %g', p_list(i));
    counter = counter +1;
    dw = 2/Pol(1);
    %legends{counter} = sprintf('Linear approx: dw = %g', dw);
    %counter = counter +1;
end
legend(legends)
title('<R^2>')
xlabel('N');
ylabel('<R^2>');
figure
legends = {};
num_p
start = 1;
stop = num_p;
counter = 1;
list = [num_p];
list =  start:stop;
v = 1.35;
B = 0.2;
dw = 2.4;
mu = (dw-2)*v+B;
for i=list,
    R_load = load(sprintf('R_square_%d.mat', i), 'R_sq');
    p = p_list(i);
    x = N_list./(p-pc)^-(2*v+mu-B);
    y = R_load.R_sq;
    plot(x, y);
    hold all
    %Pol = polyfit(log10(N_list), log10(R_load.R_sq), 1);
    %plot(log10(N_list), log10(N_list)*Pol(1) + Pol(2));
    legends{counter} = sprintf('p = %g', p_list(i));
    counter = counter +1;
    %dw = 2/Pol(1);
    %legends{counter} = sprintf('Linear approx: dw = %g', dw);
    %counter = counter +1;
end
legend(legends)
title('<R^2>')
xlabel('N/(p-pc)^{-(2v+\mu-\beta)}');
ylabel('<R^2>');
    