%
clear all; clf; close all;
% First , find the backbone
% Generate spanning cluster (l-r spanning)
%L_list = [200 400];% 800];
%
L_list = 300:25:1000
% num_exp_per_L = 500;
% pc = 0.5927;
% %L_list = 30;
% prob = pc;
% conductivity = zeros(size(L_list)); 
% for j=1:length(L_list),
%     L = L_list(j)
%     lx = L;
%     ly = L;
%        
%     for exp=1:num_exp_per_L,
%         ncount = 0;
%         perc = [];
%         while (size(perc ,1)==0)
%             ncount = ncount + 1;
%             if (ncount >1000)
%                 'hei'
%                 return
%             end
%             z=rand(lx,ly)<prob;
%             [lw,num]=bwlabel(z,4);
%             perc_x = intersect(lw(1,:),lw(lx ,:));
%             perc = perc_x(find(perc_x >0));
%         end
% %         s = regionprops(lw,'Area');
% %         clusterareas = cat(1,s.Area);
% %         maxarea = max(clusterareas);
% %         i = find(clusterareas==maxarea);
%         zz = lw == perc(1);
%         % zz now contains the spanning cluster
%         % Transpose
%         zzz = zz';
%         % Generate bond lattice from this
%         g = sitetobond(zzz);
%         % Generate conductivity matrix
%         [p c_eff] = FIND_COND(g,lx,ly);
%         % Transform this onto a nx x ny lattice
%         x = coltomat(full(p),lx,ly);
%         P = x.*zzz;
%         g1 = g(:,1);
%         g2 = g(:,2);
%         z1 = coltomat(g1,lx,ly);
%         z2 = coltomat(g2,lx,ly);
%         % Plotting
% 
%         f2 = zeros(lx,ly);
%         for iy = 1:ly-1
%             f2(:,iy) = (P(:,iy) - P(:,iy+1)).*z2(:,iy);
%         end
%         f1 = zeros(lx,ly);
%         for ix = 1:lx-1
%             f1(ix ,:) = (P(ix ,:) - P(ix+1,:)).*z1(ix ,:);
%         end
%         % Find the sum of absolute fluxes into each site
%         fn = zeros(lx,ly);
%         fn = fn + abs(f1);
%         fn = fn + abs(f2);
%         fn(:,2:ly) = fn(:,2:ly) + abs(f2(:,1:ly -1));
%         fn(:,1) = fn(:,1) + abs((P(:,1) - 1.0).*(zzz(:,1)));
%         fn(2:lx ,:) = fn(2:lx ,:) + abs(f1(1:lx -1,:));
% 
%         sum_flux_inn = sum(fn(:,1));
%         f_size = size(fn);
%         len = f_size(2);
%         sum_flux_ut = sum(fn(:,len));
%         system_flux = sum_flux_inn + sum_flux_ut;
%         conductivity(j) = conductivity(j) + system_flux;
%     end
%         
% 
%     end
%     conductivity = conductivity/num_exp_per_L;
%     save(sprintf('conductivity_L.mat'), 'conductivity');
    


close all;
legends = {};
counter = 1;
for i=1:1,
    C = load(sprintf('conductivity_L.mat'), 'conductivity');
    x = log10(L_list);
    y = log10(C.conductivity);
    plot(x,y, '.');
    legends{counter} = sprintf('Experiment');%, L = %g', L_list(i));
    counter = counter +1;
    hold all
    Pol = polyfit(x,y,1);
    plot(x, x*Pol(1) + Pol(2));
    legends{counter} = sprintf('Linear approx, \\zeta = %g', -Pol(1));
    counter = counter +1;
    hold all
end
legend(legends);
title('Conductivity - loglog');
xlabel('log10(L)');
ylabel('log10(sigma)');


figure
legends = {};
counter = 1;
for i=1:1,
    C = load(sprintf('conductivity_L.mat'), 'conductivity');
    x = (L_list);
    y = (C.conductivity);
    plot(x,y, '.');
    legends{counter} = sprintf('Experiment');%, L = %g', L_list(i));
    counter = counter +1;
    hold all
    %Pol = polyfit(x,y,1);
    y_new = x.^(Pol(1));
    y_new = y_new*y(1)/y_new(1);
    plot(x, y_new);
    legends{counter} = sprintf('sigma =A L^{-zeta}');
    counter = counter +1;
    hold all
end
legend(legends);
title('Conductivity');
xlabel('(L)');
ylabel('(sigma)');

