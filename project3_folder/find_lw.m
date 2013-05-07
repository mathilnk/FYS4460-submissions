function [lw,num] = find_lw(z)
[lw,num] = bwlabel(z,4);
end