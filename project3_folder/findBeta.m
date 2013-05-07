function beta  = findBeta(P_array, p_array, pc)
index = find(p_array>pc);
new_p = p_array(index)
new_P = P_array(index);
x = new_p-ones(size(new_p))*pc;
%y = x.^4;
y = new_P;
%loglog(x,y);
plot(log(x),log(y));
title('loglog');
beta = 1;
end