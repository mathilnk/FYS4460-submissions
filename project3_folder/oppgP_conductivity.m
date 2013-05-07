%
clear all; clf; close all;
% First , find the backbone
% Generate spanning cluster (l-r spanning)
L_list = [200 400];
num_exp_per_p = 500;
pc = 0.5927;
num_p = 100;
dp = 0.2;
p_list = linspace(pc, pc +dp, num_p);
%L_list = 30;
% for i=1:length(L_list),
%     L = L_list(i)
%     lx = L;
%     ly = L;
%     conductivity = zeros(size(p_list));
%     for j=1:length(p_list),
%             prob = p_list(j);
%          
%         for exp=1:num_exp_per_p,
%             ncount = 0;
%             perc = [];
%             while (size(perc ,1)==0)
%                 ncount = ncount + 1;
%                 if (ncount >1000)
%                     'hei'
%                     return
%                 end
%                 z=rand(lx,ly)<prob;
%                 [lw,num]=bwlabel(z,4);
%                 perc_x = intersect(lw(1,:),lw(lx ,:));
%                 perc = perc_x(find(perc_x >0));
%             end
%     %         s = regionprops(lw,'Area');
%     %         clusterareas = cat(1,s.Area);
%     %         maxarea = max(clusterareas);
%     %         i = find(clusterareas==maxarea);
%             zz = lw == perc(1);
%             % zz now contains the spanning cluster
%             % Transpose
%             zzz = zz';
%             % Generate bond lattice from this
%             g = sitetobond(zzz);
%             % Generate conductivity matrix
%             [p c_eff] = FIND_COND(g,lx,ly);
%             % Transform this onto a nx x ny lattice
%             x = coltomat(full(p),lx,ly);
%             P = x.*zzz;
%             g1 = g(:,1);
%             g2 = g(:,2);
%             z1 = coltomat(g1,lx,ly);
%             z2 = coltomat(g2,lx,ly);
%             % Plotting
% 
%             f2 = zeros(lx,ly);
%             for iy = 1:ly-1
%                 f2(:,iy) = (P(:,iy) - P(:,iy+1)).*z2(:,iy);
%             end
%             f1 = zeros(lx,ly);
%             for ix = 1:lx-1
%                 f1(ix ,:) = (P(ix ,:) - P(ix+1,:)).*z1(ix ,:);
%             end
%             % Find the sum of absolute fluxes into each site
%             fn = zeros(lx,ly);
%             fn = fn + abs(f1);
%             fn = fn + abs(f2);
%             fn(:,2:ly) = fn(:,2:ly) + abs(f2(:,1:ly -1));
%             fn(:,1) = fn(:,1) + abs((P(:,1) - 1.0).*(zzz(:,1)));
%             fn(2:lx ,:) = fn(2:lx ,:) + abs(f1(1:lx -1,:));
% 
%             sum_flux_inn = sum(fn(:,1));
%             f_size = size(fn);
%             len = f_size(2);
%             sum_flux_ut = sum(fn(:,len));
%             system_flux = sum_flux_inn + sum_flux_ut;
%             conductivity(j) = conductivity(j) + system_flux;
%         end
%         
% 
%     end
%     conductivity = conductivity/num_exp_per_p;
%     save(sprintf('conductivity_%d.mat', i), 'conductivity');
%     
% end

legends = {};
counter = 1;
for i=1:length(L_list),
    C = load(sprintf('conductivity_%d.mat', i), 'conductivity');
    x = (p_list - pc);
    y = (C.conductivity);
    plot(x,y, '.');
    legends{counter} = sprintf('Experiment, L = %g', L_list(i));
    counter = counter +1;
    hold all
    Pol = polyfit(x,y,1);
    %plot(x, x*Pol(1) + Pol(2));
    %legends{counter} = sprintf('Linear approx, \\mu = %g', Pol(1));
    %counter = counter +1;
    %hold all
end
legend(legends);
title('Conductivity');
xlabel('(p-p_c)');
ylabel('sigma');





