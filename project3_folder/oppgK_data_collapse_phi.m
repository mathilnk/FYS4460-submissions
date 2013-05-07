v = 1.35;
pc = 0.59275;
num_exp_per_p = 3000;
L_list = [25 50 100 200 400 800];
%L_list = [25];
p_list = linspace(pc-0.1, pc+0.1, 100);
for i=1:length(L_list),
    L = L_list(i)
    Pi = zeros(size(p_list));
    for j=1:length(p_list),
        p = p_list(j);
        count_success = 0;
        for k=1:num_exp_per_p,
            r = rand(L,L);
            z = r<p;
            [lw,num] = bwlabel(z,4);
            spanning = find_span_cl_numbers(lw);
            if(length(spanning)>=1)
                count_success = count_success +1;
            end
        end
        Pi(j) = count_success/num_exp_per_p;
    end
    save(sprintf('Pi_%d.mat', L), 'Pi');
end



%plot pi(p,L) to find phi(u)
figure
legends = {};
for i=1:length(L_list),
    L = L_list(i);
    Pii = load(sprintf('Pi_%d.mat', L), 'Pi');
    x = (p_list - pc)*L^(1/v);
    y = Pii.Pi;
    plot(x,y);
    hold all
    legends{i} = sprintf('L = %d', L);
end
title('$\Phi (u)$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$(p-p_c)L^{1/\nu}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Phi(u)$', 'Interpreter', 'latex', 'FontSize', 14);
legend(legends);