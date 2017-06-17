function [ out ] = vect( X )
        
    sample = size(X,3);
    N = size(X, 1);
    out = zeros(sample, N * (N + 1) / 2);

    for i = 1:sample 
        out(i, 1:N) = diag(X(:,:,i), 0);% * sqrt(2);
        sum = N;
        for k = 1:N-1       
            out(i, sum+1:sum+N-k) = diag(X(:,:,i), k);
            sum = sum + N-k;
        end
    end
end

