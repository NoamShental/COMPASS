% This function returns the lengths of elements in a cell array
function l = length_cell(c)
[m n] = size(c); l=zeros(m,n);
for i=1:m
    for j=1:n
        l(i,j)= length(c{i,j});
    end
end
