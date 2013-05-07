% L_min = 100;
% L_max = 600;
% L_num = 100;
% dL = round((L_max -L_min)/L_num) +1;
% L_max = dL*L_num;
% L_list = L_min:dL:L_max
% Mass = zeros(size(L_list));
% p = 0.59275;
% 
% 
% num_exp_per_L = 1000;
% for i=1:length(L_list),
%     L = L_list(i)
%     for exp=1:num_exp_per_L,
%         perc = [];
%         ncount = 0;
%         while (size(perc ,1)==0)
%             ncount = ncount + 1;
%             if (ncount >1000)
%                 return
%             end
%             z=rand(L,L)<p;
%             [lw,num]=bwlabel(z,4);
%             perc_x = intersect(lw(1,:),lw(L ,:));
%             perc = perc_x(find(perc_x >0));
%         end
%         %s = regionprops(lw,'Area');
%         %clusterareas = cat(1,s.Area);
%         %maxarea = max(clusterareas);
%         %i = find(clusterareas==maxarea);
%         zz = zeros(size(lw));
%         for j=1:length(perc),
%             zz = zz + lw==perc(j);
%         end
%         %zz = lw == i;
%         % zz now contains the spanning cluster
%         %imagesc(zz);
%         % Display spanning cluster
%         % Run walk on this cluster
%         [l,r] = walk(zz);
%         zzz = l.*r;
%         Mass(i) = Mass(i) + length(find(zzz~=0));
%     end
%     %Mass(i) = Mass(i)/num_exp_per_L;
% end

close all;
figure
x = log10(L_list);
y = log10(Mass);
plot(x,y, '.');
Pol = polyfit(x,y,1)
hold all
plot(x, Pol(1)*x + Pol(2));
title('LogLog of the mass of the singly connected bonds');
xlabel('log10(L)');
ylabel('log10($M_{SC}$)', 'Interpreter', 'latex');
legend('Experiment', sprintf('Linear approx: $D_{SC}$ = %g', Pol(1)));
figure
plot(L_list, Mass(1)*L_list.^Pol(1)/L_list(1)^Pol(1));
hold all
plot(L_list, Mass, '.');
title('Mass of the singly connected bonds');
xlabel('L');
ylabel('$M_{SC}$', 'Interpreter', 'latex');
legend('$M_{SC}\propto L^{D_{SC}}$','Experiment', 'Interpreter', 'latex');