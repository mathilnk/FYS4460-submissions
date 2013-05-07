L=100;
r = rand(L,L);
p=0.6;
z = r<p;
[lw,num] = bwlabel(z,4);
num
img=label2rgb(lw,'jet', 'k', 'shuffle');
image(img);
s=regionprops(lw,'Area');
s.Area;
area = cat(1,s.Area);
s2 = regionprops(lw, 'BoundingBox');
bbox = cat(1,s2.BoundingBox);
