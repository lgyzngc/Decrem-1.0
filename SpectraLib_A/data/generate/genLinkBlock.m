function [ aa, dd ] = genLinkBlock( m, n, tresh, dd, istriu )

% function [ aa, dd ] = genLinkBlock( m, n, tresh,[] ,istriu )
%
% generates a random {0,1} matrix of size m x n, for which
%	the marginal probability of 1 = thresh 
%       the probability of 1 in row i = dd( i )
%		where dd( i ) is a non-uniform distribution generated
% 		by this function
% istriu = 1 ==> an upper triangular matrix is generated
% 0 < thresh < 1 
% dd(m,1) = distribution of row weigths
%	    if dd = [] then dd is generated internally
%
% This function is used in genBlockACitationMatrix()
%
%
% $Authors: Marina Meila
% $Part of SpectraLib_A
% $Last revision: 06-June-2007
%


mfrac = 30;  % what fraction of rows are dense

if isempty( dd )
   dd = rand( m, 1 );
   dd(1:round(m/mfrac))=dd(1:round(m/mfrac))*10;
   dd = dd./sum( dd );
   dd = dd( randperm( m ));
end;

aa = rand( m, n );
if istriu
   aa = triu( aa, 1 );
end;
aa = aa.*repmat( dd, [1,n] );
aa = aa/sum(sum( aa ));
dum = sort( reshape( aa, 1, m*n ) );
dum = dum( find( dum ));   % eliminate 0's if istriu
tresh2 = dum( round((1-tresh)*length( dum )) );
aa = (aa > tresh2 );


if 0 
    pmarg = sum( aa, 2 );
    pmarg = pmarg/sum(pmarg);
    if istriu
        plot( 1:m, pmarg, ':or', 1:m, dd.*((m:-1:1)'+n-m)/n, '-b') % !!!??
    else
        plot( 1:m, pmarg, ':or', 1:m, dd, '-b')
    end;
end;
