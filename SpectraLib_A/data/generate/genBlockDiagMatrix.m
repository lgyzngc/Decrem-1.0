function [ aa ] = genBlockDiagMatrix( sizes, strengths )

% function [ aa ] = genBlockDiagMatrix( sizes, strengths )
%
% generates a block diagonal matrix. sizes is the vector containing
% the size of each block. block i is equal to ones( sizes( i )) *
% strengths( i )
% 
% $Authors: Marina Meila
% $Part of SpectraLib_A
% $Last revision: 06-June-2007

nbloc = length( sizes );
aa = zeros( sum( sizes ));
ioff = 0;

for ii = 1:nbloc;
  aa( (ioff+1):(ioff+sizes(ii)), (ioff+1):(ioff+sizes(ii))) = strengths(ii)*ones(sizes(ii));
  ioff = ioff + sizes( ii );
end;

