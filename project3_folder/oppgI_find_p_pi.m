close all;
clear all;
pc = 0.59275;
num_exp_per_p = 2000;
L_list = [25 50 100 200 400 800];
%L_list = [25];

error = 100;
tol = 0.002;
Pi_goal = 0.3;
p_pi_x = zeros(size(L_list));

tol_p = 0.00006;
for i=1:length(L_list)
    L = L_list(i)
    p_large = 1;
    p_small = 0;
    i
    error = 100;
    counter = 0;
    p_old = 100;
    while abs(error)>tol,
        counter = counter +1;
        hist_Areas = zeros(L*L,1);
        hist_Areas2 = zeros(L*L,1);
        %p = p_start;
        p = p_small + (p_large - p_small)/2
        count_success = 0;
        for k=1:num_exp_per_p,
            r = rand(L,L);
            z = r<p;
            [lw,num] = bwlabel(z,4);     
            spanning = find_span_cl_numbers(lw);
            if(length(spanning)>=1)
                count_success = count_success + 1;
            end
        end
        Pi_test = count_success/num_exp_per_p;
        error = Pi_goal - Pi_test;
        Pi_test
        if(abs(p-p_old)<tol_p)
            break
        end
        p_old =p;
        if(error<0)
            p_large = p;
        else
            p_small = p;
        end
    
    end
    p_pi_x(i) = p;
end
save(sprintf('p_pi_x_0_3.mat'), 'p_pi_x');


pc = 0.59275;
num_exp_per_p = 3000;
L_list = [25 50 100 200 400 800];
%L_list = [25];

error = 100;
tol = 0.002;
Pi_goal = 0.8;
p_pi_x = zeros(size(L_list));

for i=1:length(L_list)
    L = L_list(i)
    p_large = 1;
    p_small = 0;
    i
    error = 100;
    counter = 0;
    p_old = 100;
    while abs(error)>=tol,
        counter = counter +1;
        hist_Areas = zeros(L*L,1);
        hist_Areas2 = zeros(L*L,1);
        %p = p_start;
        p = p_small + (p_large - p_small)/2
        count_success = 0;
        for k=1:num_exp_per_p,
            r = rand(L,L);
            z = r<p;
            [lw,num] = bwlabel(z,4);     
            spanning = find_span_cl_numbers(lw);
            if(length(spanning)>=1)
                count_success = count_success + 1;
            end
        end
        Pi_test = count_success/num_exp_per_p;
        error = Pi_goal - Pi_test;
        Pi_test
        if(abs(p-p_old)<tol_p)
            break
        end
        p_old =p;
        if(error<0)
            p_large = p;
        else
            p_small = p;
        end
        if(counter >20)
            error = tol-1;
        end
    end
    p_pi_x(i) = p;
end
save(sprintf('p_pi_x_0_8.mat'), 'p_pi_x');



%j
P3 = load('p_pi_x_0_3.mat', 'p_pi_x');
P8 = load('p_pi_x_0_8.mat', 'p_pi_x');
dp = (P8.p_pi_x - P3.p_pi_x);
x = log10(L_list);
y = log10(dp);
figure
plot(x,y);
hold all
Pol = polyfit(x,y,1)
plot(x, Pol(1)*x + Pol(2));
v = -1/Pol(1)










