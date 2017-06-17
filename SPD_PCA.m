% This is an unsupervised dimensionality reduction algorithm proposed by
% Horev et al. [1].
% 
% PARAMS:
% 	data: An m*m*n dimensional matrix which contains n SPD matrix of size m*m
% 	dim : The dimensionality of SPD matrices after DM
%	distType: The distance type. Either `riemann` or `ld`
%	verbose	: More verbosity if this is 1
%
% RETUTNS:
%	U	: The calculated transformation matrix which can be used as 
%		  below to transform an m*m SPD matrix `x` to an dim*dim 
%		  SPD matrix y:
%				y = U'*x*U
%	obj	: Final value of the objective function
%	M	: The riemannian mean of data
%
% References:
% 			[1] Horev, Inbal, Florian Yger, and Masashi Sugiyama. 
%				"Geometry-aware principal component analysis for symmetric positive 
%				definite matrices." Asian Conference on Machine Learning. 2016.
%
function [U, obj, M] = SPD_PCA( data, dim, distType, varargin)

verbose = 1;
if(~isempty(varargin))
    verbose = varargin{1};
end

[M, ~, ~] = riemann_mean(data);

MM = M^-0.5;
for k = 1:size(data,3)
    data(:,:,k) = MM * data(:,:,k) * MM;
end

I = eye(dim);

u0 = randn(size(data, 1), dim);
u0 = orth(u0);

opts.record = verbose;
opts.mxitr  = 1000;
opts.xtol = .1;
opts.gtol = .1;
opts.ftol = .1;
%obj.tau = 1e-3;
%opts.nt = 1;

t = tic;
[U, obj]= OptStiefelGBB(u0, @objfunc, opts);
tsolve = toc(t);

if(verbose)
    [f,~] = objfunc(u0);
    [flast,~] = objfunc(U);
    
    fprintf('Elapsed time: %f\n', tsolve);
    fprintf('Stoped: %s , Init f: %f, Last f: %f \n',obj.msg, f, flast);
end


    function [F, G] = objfunc(u)
        
        F = 0;
        G = 0;
        for i = 1:size(data, 3)
            if(strcmp(distType,'ld'))
                F = F + distance_ld(u' * data(:,:,i) * u, I);
                
                G = G + 2 * inv(u')  - u ...
                    - data(:,:,i) * u * inv(u' * data(:,:,i) * u);
                
            elseif(strcmp(distType,'riemann'))
                F = F + distance_riemann(u' * data(:,:,i) * u, I);
                
                G = G + 4 * (data(:,:,i) * u * inv(u'*data(:,:,i)*u) - u) ...
                    * logm(u'*data(:,:,i)*u);
            end
            
        end
    end

end