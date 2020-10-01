function [ S, diffS ] = genBlockSumSymMatrix( sizes, Phat, convergence, degs, S0)

% function [ S, diffS ] = genBlockSumSymMatrix( sizes, Phat, convergence )
% 
% Generates a square symmetric matrix S such that matix P=D^(-1)*S 
% (where D=diag(degs)) has PCE with the block sums ('hat') matrix given by Phat.
%
% size( S ) = sum( sizes )
% The size of block (i,j) is given by (sizes (i),sizes(j))
%
% CONSISTENCY CONDITION:
% diag(untelescope(degs', sizes))*Phat=H is symmetric matrix
%  untelescope computes sums of degs over the subsets given by sizes
%
% degs is a _column_ vector (default value = random (uniform) normalized 
% to satisfy the consistency condition with H=diag(sizes)*Phat. To assure
% that H is symmetric in this case, please form Phat as shown in the
% example below)
%
%
% convergence: convergence(1) is the S-difference tolerance (default value = 10^(-6))
%              convergence(2) is the maximum number of iterations (default value = 1000) 
%
% S0 is a starting matrix (default=random(uniform)), not necessarily
% symmetric
%
% diffS is max S-difference at the last iteration
% 
% 
% REMARKS: 
% 1) If each row in Phat sums to 1, then row sums of S are equal to degs.
%    Hence, in this case we generate a square matrix S such that 
%    P=D^(-1)*S (where D is a diagonal matrix of row sums of S) 
%    is block-stochastic with the transition _probability_ matrix given by
%    Phat.
%
% 2) If the consistency condition doesn't hold the algorithm in 
%    general still converges but the transition ('hat') matrix is not equal to Phat! 
%
% 3) The resulting P is not block-stochastic when H is asymmetric with some of
%    the (off-diagonal) elements equal to 0
%    This is the only case of which I know when algorithm doesn't converge
%
%   EXAMPLE:
%     sizes = [ 10 20 30 20 20 ];  
%     %Create a symmetric (k,k) matrix that sums to sizes
%     H = [ 5 1 2 1 1; 1 9 5 3 2; 2 5 14 4 5; 1 3 4 8 4; 1 2 5 4 8 ];
%     Phat = diag(sizes)\H;
%     S = genBlockSumSymMatrix( sizes, Phat );
%
%     % Check results:
%     % 1) S is symmetric (S-S' should be equal to zero)
%     imagesc(S-S'); 
%     % 2) P is block-stochastic with the transition matrix equal to Phat
%     P=S./repmat(sum(S,2),[1,size(S,2)]); 
%     untelescope(P,sizes)
%     imagesc( untelescope( P, sizes )); 
%    
% $Authors: Marina Meila, Tatiana Maravina
% $Part of SpectraLib_A
% $Last revision: 06-June-2007


if nargin<3,  convergence=[]; end
if isempty(convergence), convergence=10^(-6); end
if length(convergence)<2, convergence(2)=1000; end
if nargin<4, degs=[]; end

%preliminaries
nsizes = length( sizes );
ji = [0 cumsum(sizes) ];
nn = sum( sizes );
if nargin<5
    S = rand( nn );
else 
    S=S0;
end

if isempty(degs) % degs not specified
    degs = rand(nn,1);
    degs=(degs./telescope(untelescope(degs', sizes)', sizes)).*telescope(sizes', sizes);
end

% check the consitency condition
H=diag(untelescope(degs', sizes))*Phat; 
H=abs(H-H');
if max(max(H))>1e-12 % zero is not exactly zero...
    warning('genBlockSumSymMatrix:consistency', 'H is not symmetric. The transition matrix might not be equal to Phat!')
end

strengths = diag( degs )*telescope( Phat, sizes ); % The row sums in block (i,j) equal to strengths(i,j).

% loop
loop=0;
done=0;
S = (S+S')/2;
while (~done)
    % for the stopping criteria
    loop=loop+1;
    Sold = S;
    % normalize
    for isc = 1:nsizes;
        ii = ji(isc)+1:ji(isc+1);
        sumS = sum( S( :, ii), 2 );
        sumS(sumS==0)=1; % to avoid division by zero below
        S( :, ii ) = diag( strengths(:, isc) ./ sumS ) * S( :, ii );
    end;
    
    % symmetrize
    S = (S+S')/2;
    
    % for the stopping criteria
    diffS = max(max(abs(Sold-S)));
    done = (diffS<convergence(1))|(loop>=convergence(2));
end; %loop

