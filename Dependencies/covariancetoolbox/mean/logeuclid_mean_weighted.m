
function A = logeuclid_mean_weighted(B, W)

K = size(B,3); % Nombre de matrices
% Initialisation avec la matrice moyenne arythmetique
fc = zeros(size(B,1));
for i=1:K
    fc = fc + W(i) * logm(B(:,:,i));
end
A = expm( (1/K)*fc / sum(W));

end

