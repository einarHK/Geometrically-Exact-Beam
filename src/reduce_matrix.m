
% reduces the number of columns and rows, for a given matrix.
% n1 - start index, n2 - end index. 
function S2 = reduce_matrix(S, n1, n2)
    [n_rows, n_cols] = size(S); 
    dn = n2 - n1+1; 
    S2 = zeros(n_rows - dn, n_cols - dn);
    disp(size(S2));
    for i1=1:n_rows
        for i2=1:n_cols
            if ((n1 <= i1 <= n2) && (n1 <= i2 <= n2))
                
            else
                S2(i1,i2) = S(i1,i2);
            end
        end
    end
end

