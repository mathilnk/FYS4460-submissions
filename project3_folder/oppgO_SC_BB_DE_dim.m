%
clear all; clf; close all;
% First , find the backbone
% Generate spanning cluster (l-r spanning)
L_list = 50:50:800;
num_exp_per_L = 1000;
pc = 0.5927;
num_p = 10;
dp = 0.1;
p_list = linspace(pc, pc +dp, num_p);
%L_list = 30;
% for j=1:length(p_list),
%     prob = p_list(j);
%     mass_DE = zeros(size(L_list));
%     mass_SC = zeros(size(L_list));
%     mass_BB = zeros(size(L_list));
%     for i=1:length(L_list),
%             L = L_list(i)
%             lx = L;
%             ly = L;
%         for exp=1:num_exp_per_L,
%    
%             ncount = 0;
%             perc = [];
%             while (size(perc ,1)==0)
%                 ncount = ncount + 1;
%                 if (ncount >1000)
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
%             SC = (fn>=system_flux);
%             conductivity = system_flux;
%             limit = 0.00001;
%             zfn = fn>limit;
%             zbb = (zzz + 2*zfn);
%             zbb = zbb/max(max(zbb));
%             eps = 0.00001;
%             mass_BB(i) = mass_BB(i) + sum(sum(abs(zbb-ones(size(zbb)))<eps),2);
%             mass_SC(i) = mass_SC(i) + sum(sum(abs(SC-ones(size(zbb)))<eps),2);
%             mass_DE(i) = mass_DE(i) + sum(sum(abs(zbb-ones(size(zbb)))>eps),2);
%         end
%         
% 
%     end
%     mass_BB = mass_BB/num_exp_per_L;
%     mass_SC = mass_SC/num_exp_per_L;
%     mass_DE = mass_DE/num_exp_per_L;
%     save(sprintf('mass_BB_%d.mat', j), 'mass_BB');
%     save(sprintf('mass_SC_%d.mat', j), 'mass_SC');
%     save(sprintf('mass_DE_%d.mat', j), 'mass_DE');
%     
% end





legends = {};
counter = 1;
close all;
figure(1)
for j=1:length(p_list),
    hold all
    p = p_list(j);
    BB = load(sprintf('mass_BB_%d.mat', j), 'mass_BB');
    y_BB = log10(BB.mass_BB);
    x = log10(L_list);
    legends{counter} = sprintf('Experiment p =%g', p);
    counter = counter +1;
    plot(x,y_BB, '.');
    
    PBB = polyfit(x, y_BB, 1);
    PBB_l(j) = PBB(1);
    legends{counter} = sprintf('Linear approx, D_{BB} = %g', PBB(1));
    counter = counter +1;
    plot(x, x*PBB(1) + PBB(2));
    %hold all
    
    
end

legend(legends);
title('loglog M_{BB}');
xlabel('log10 L')
ylabel('log10 M_{BB}');



legends = {};
counter = 1;
figure(3)
for j=1:length(p_list),
    hold all
    p = p_list(j);
    DE = load(sprintf('mass_DE_%d.mat', j), 'mass_DE');
    y_DE = log10(DE.mass_DE);
    x = log10(L_list);
    legends{counter} = sprintf('Experiment p =%g', p);
    counter = counter +1;
    plot(x,y_DE, '.');
    
    PDE = polyfit(x, y_DE, 1);
    PDE_l(j) = PDE(1);
    legends{counter} = sprintf('Linear approx, D_{DE} = %g', PDE(1));
    counter = counter +1;
    plot(x, x*PDE(1) + PDE(2));
    %hold all
    
    
end

legend(legends);
title('loglog M_{DE}');
xlabel('log10 L')
ylabel('log10 M_{DE}');



legends = {};
counter = 1;
figure(2)
for j=1:length(p_list),
    hold all
    p = p_list(j);
    SC = load(sprintf('mass_SC_%d.mat', j), 'mass_SC');
    y_SC = log10(SC.mass_SC);
    x = log10(L_list);
    legends{counter} = sprintf('Experiment p =%g', p);
    counter = counter +1;
    plot(x,y_SC, '.');
    
    PSC = polyfit(x, y_SC, 1);
    
    legends{counter} = sprintf('Linear approx, D_{SC} = %g', PSC(1));
    counter = counter +1;
    plot(x, x*PSC(1) + PSC(2));
    %hold all
    
    
end

legend(legends);
title('loglog M_{SC}');
xlabel('log10 L')
ylabel('log10 M_{SC}');





