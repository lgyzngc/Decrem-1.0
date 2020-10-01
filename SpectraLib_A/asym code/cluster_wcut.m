function asig=cluster_wcut (A, T, k, normalize) 

% asig=cluster_wcut (A, T, k) 
% Implements BestWCut algorithm as described in Marina Meila and William
% Pentney 'Clustering by weighted cuts in directed graphs'
%
% A - n-by-n asymmetric affinity matrix with _non-negative_ elements
% T - column vector of _positive_ weights of length n
% k - number of clusters
% normalize - if 1 (default value) then the rows of eigen vectors (n-by-k
% matrix) are normalized to unit length before kmeans clustering (it was
% noticed to provide better clustering in practice) 
%
% returns _row_ vector of cluster assignments (clusters are marked 1:k)
%
% EXAMPLE: 
%         webkb_readdata;
%         iifill = zeros( size( Do));
%         iifill( ido0 ) = 1;
%         Dofill = Do+1*iifill; %use 1 for outdegrees where outdegrees is 0
%         T=Dofill.^1.5;
%         asig=cluster_wcut(A,T,k);
%         disp(clustering_error( asig, asig_true ));
%         disp(VI( asig, asig_true )); % compute variation of information
%
%
% $Authors: Marina Meila, Tatiana Maravina
% $Part of SpectraLib_A
% $Last revision: 06-June-2007
%

if nargin<4, normalize=1; end

n=size(A,1);
opts.disp = 0;
opts.issim = 1;
opts.p = 2*k;

sqrtTinv = 1./sqrt( T );
Do=sum(A,2);
B = (diag(Do)-A).* repmat( sqrtTinv, [1, n] ).*repmat( sqrtTinv', [n, 1] );
H = (B + B')/2;

% Eigenproblems
[yy, ll] = eigs( H, k+1, 'sa',  opts );
vv = yy(:,1:k).*repmat( sqrtTinv, [1, k ] );

% Cluster
[ center, asig, distort] = cluster_normalized_kmeans( vv, k, normalize);


