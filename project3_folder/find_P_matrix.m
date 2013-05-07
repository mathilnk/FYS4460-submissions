function [P_matrix, p_arr] = find_P_matrix(L_arr, N)
P_matrix = zeros(length(L_arr), N);
pc = 0.59275;
p_arr = logspace(log10(pc-0.1), log10(1), N);
%p_arr = linspace(pc-0.1,1, N); 
P_arr = zeros(1,N);
num_experiment = 50;
%z= gen_bin_array(L);
for i=1:length(L_arr)
    sprintf('Starting with L=%d',L_arr(i))
    L = L_arr(i);
    for k=1:num_experiment
        r = rand(L,L);      
        for j=1:N,
            p = p_arr(j);
            [z, lw, num] = gen_bin_array(r,p);
            common = find_span_cl_numbers(lw);
            P_arr(j) = P_arr(j) + find_P(common, lw)./num_experiment;
            
        end
        
    end
    P_matrix(i,:) = P_arr;
end

end
