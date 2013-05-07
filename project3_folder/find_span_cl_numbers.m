function common = find_span_cl_numbers(lw)
sizes = size(lw);
y_size = sizes(1);
x_size = sizes(2);
%sjekker hvor mange tall som er felles i den øverste og nederste raden
common_top_bottom = intersect(lw(1,:),lw(y_size,:));
%sjekker hvor mange tall som er felles i den venstre og høyre kolonna
common_left_right = intersect(lw(:,1),lw(:,x_size));
%tar unionen i tilfelle et spanning cluster går fra for eks høyre til top;
%passer på dobbelttelling
common_w_zero = union(common_top_bottom, common_left_right);
%fjerner null
indices = find(common_w_zero);
common = common_w_zero(indices);
end
