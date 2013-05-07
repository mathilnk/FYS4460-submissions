function [bin_array, lw, num] = gen_bin_array(r, p)
z = r<p;
[lw,num] = bwlabel(z,4);
bin_array = z;
end
