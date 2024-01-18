p_11 = 3;
for q = 1:10
    if q > p_11-1
        [q mod(q-p_11,m)+1]
    end
end