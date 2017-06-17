function [ out ] = unvect( X )
        
    sample = size(X,1);
    N = floor(sqrt(size(X, 2) * 2));
    out = zeros(N, N, sample);

    for i = 1:sample 
        
        sum = N;
        for k = 1:N-1  
            out(:,:,i) = out(:,:,i) + diag(X(i,sum+1:sum+N-k), k);
            sum = sum + N-k;
        end
        out(:,:,i) = out(:,:,i) + out(:,:,i)';
        out(:,:,i) = out(:,:,i) + diag(X(i,1:N), 0);
    end
end
