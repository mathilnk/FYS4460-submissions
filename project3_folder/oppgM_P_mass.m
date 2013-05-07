% L_list = [50 100 200 400 800];
% %L_list = [10];
% dp = 0.2;
% p_c = 0.59275;
% num_p = 30;
% p_list = linspace(p_c, p_c+dp, num_p);
% num_exp_per_p = 100;
% d = 2;
% 
% for l=1:length(L_list),
%     L = L_list(l)
%     prob = zeros(size(p_list));
%     for i=1:length(p_list),
%         p = p_list(i);
%         for exp=1:num_exp_per_p,
%             perc = [];
%             ncount = 0;
%             while (size(perc ,1)==0)
%                 ncount = ncount + 1;
%                 if (ncount >1000)
%                     'her'
%                     return
%                 end
%                 z=rand(L,L)<p;
%                 [lw,num]=bwlabel(z,4);
%                 perc_x = intersect(lw(1,:),lw(L ,:));
%                 perc = perc_x(find(perc_x >0));
%             end
%             zz = zeros(size(lw));
%             for j=1:length(perc),
%                 zz = zz + lw==perc(j);
%             end
%             [l,r] = walk(zz);
%             zzz = l.*r;
%             mass = length(find(zzz~=0));
%             prob(i) = prob(i) + mass;
%         end
% 
%     end
%     prob = prob/(L^d*num_exp_per_p);
%     save(sprintf('prob_mass_sc_L%d.mat', L), 'prob')
%     
% end
close all;
legends = {};
for i=1:length(L_list),
    L=L_list(i);
    P = load(sprintf('prob_mass_sc_L%d.mat', L), 'prob');
    plot(p_list-p_c, P.prob);
    legends{i} = sprintf('L = %d', L);
    hold all
end
title('P_{SC}')%, 'Interpreter', 'latex');
legend(legends);
xlabel('p-p_c')%, 'Interpreter', 'latex');
ylabel('P_{SC} = M_{SC}/L^2')%, 'Interpreter', 'latex');
