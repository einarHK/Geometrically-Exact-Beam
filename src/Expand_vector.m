
function v = Expand_vector(v1, length2)

length1 = length(v1);
v = zeros(length1 + length2,1); 
v(1:length1) = v1; 

for i=1:length2
    v(length1 + i) = v(length1) + i; 
end


end