% This function calculates the adjacency matrix (Adj) needed 
% for the Harandi's method.
% For within class points the value of adjacency matrix is 1
% and for outer class points this is -1.
%
function Adj = Harandi_adjmat( labels )

Adj = zeros(numel(labels));

for i = 1:numel(labels)
    for j = i:numel(labels)
        if( labels(i) == labels(j) )
            Adj(i,j) = 1;
        else
            Adj(i,j) = -1;
        end
        Adj(j,i) = Adj(i,j);
    end  
end
