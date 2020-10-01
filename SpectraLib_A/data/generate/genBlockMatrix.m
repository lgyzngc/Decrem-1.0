function [ aa ] = genBlockMatrix( sizes, strengths )

% function [ aa ] = genBlockMatrix( sizes, strengths )
%
% generates a matrix of constant blocks. sizes is the vector containing
% the size of each block. block i,j is equal to 
%  ones( sizes( i ), sizes( j )) * strengths( i, j )
%
% $Authors: Marina Meila
% $Part of SpectraLib_A
% $Last revision: 06-June-2007

nbloc = length( sizes );
aa = zeros( sum( sizes ));
ioffi = 0;

for ii = 1:nbloc;
   ioffj = 0;
   for jj = 1: nbloc;
      aa( (ioffi+1):(ioffi+sizes(ii)), (ioffj+1):(ioffj+sizes(jj))) = strengths(ii,jj)*ones( sizes(ii), sizes( jj ));
      ioffj = ioffj + sizes( jj );
  end;
  ioffi = ioffi + sizes( ii );
end;

