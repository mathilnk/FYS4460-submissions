% Example of use of the walk routine
% Generate spanning cluster (l-r spanning)
figure
lx =64;
ly = 64;
p = 0.585;
ncount = 0;
perc = [];
while (size(perc ,1)==0)
    ncount = ncount + 1;
    if (ncount >1000)
        return
    end
    z=rand(lx,ly)<p;
    [lw,num]=bwlabel(z,4);
    perc_x = intersect(lw(1,:),lw(lx ,:));
    perc = perc_x(find(perc_x >0));
end
%s = regionprops(lw,'Area');
%clusterareas = cat(1,s.Area);
%maxarea = max(clusterareas);
%i = find(clusterareas==maxarea);
zz = zeros(size(lw));
for j=1:length(perc),
    zz = zz + lw==perc(j);
end
%zz = lw == i;
% zz now contains the spanning cluster
imagesc(zz);
% Display spanning cluster
% Run walk on this cluster
[l,r] = walk(zz);
zzz = l.*r;
% Find points where both l and r are non-zero
zadd = zz + zzz;
subplot(2,2,1), imagesc(zz);
subplot(2,2,2), imagesc(zadd);
subplot(2,2,3), imagesc(zzz >0);
subplot(2,2,4), imagesc(l+r>0);

