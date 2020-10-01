function [ aa ] = genBlockBandMatrix( bands, sizes, strengths )

% function [ aa ] = genBlockBandMatrix( bands, sizes, strengths )
%
% generates a block diagonal matrix. sizes is the vector containing
% the size of each block. block i is equal to a band matrix of size
% sizes( i ) with bandwidth bands( i ) and elements strengths( i )
%
% $Authors: Marina Meila, Tatiana Maravina
% $Part of SpectraLib_A
% $Last revision: 06-June-2007

nbloc = length( sizes );
aa = zeros( sum( sizes ));
ioff = 0;

for ii = 1:nbloc;
  bb = toeplitz( [ ones( 1, bands( ii )) zeros( 1, sizes( ii ) - bands( ii )) ] );
  aa( (ioff+1):(ioff+sizes(ii)), (ioff+1):(ioff+sizes(ii))) = strengths(ii)*bb;
  ioff = ioff + sizes( ii );
end;

