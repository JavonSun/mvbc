% show distribution of a discrete vector

function tbl = tbl(x) 
    u = unique(x);
    tbl = zeros(2, length(u));
    for i = 1: length(u)
        tbl(1, i) = u(i);
        tbl(2, i) = length(find(x == u(i)));
    end;
end