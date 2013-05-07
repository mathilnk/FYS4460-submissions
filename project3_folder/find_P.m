function P = find_P(common, lw)
num = 0;
sizes = size(lw);
x_size = sizes(2);
y_size = sizes(1);
for i = 1:length(common),
    num = num + length(find(lw==common(i)));
end
% for row=1:y_size,
%     for col=1:x_size,
%         value = lw(row, col);
%         num = num + length(find(common==value));
%     end
% end
P = num/(y_size*x_size);
end