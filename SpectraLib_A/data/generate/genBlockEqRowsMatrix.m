function [ aa ] = genBlockEqRowsMatrix( sizes, strengths )

% function [ aa ] = genBlockEqRowsMatrix( sizes, strengths )
%
% generates a block matrix with >0 elements that has groups of equal rows.
% sizes is the vector containing the size of each block. 
% a row block contains equal rows. sizes also defines a partitioning
% of the columns. 
% strengths( i,j ) defines the row sum of the elements in block i,j
%
% $Authors: Marina Meila
% $Part of SpectraLib_A
% $Last revision: 06-June-2007

nbloc = length( sizes );
jj = [ 0 cumsum( sizes ) ];
n1 = jj( nbloc+1 );
aa = zeros( n1 );

bb = rand( nbloc, n1 );


for ii = 1:nbloc; % generate rows with given sums
  ss = sum( bb( :, jj( ii )+1:jj(ii+1) ), 2 );
  bb( :, jj( ii )+1:jj(ii+1) ) = diag( 1./ss )*diag( strengths( :, ii ))*bb( :, jj( ii )+1:jj(ii+1) );
end;

for ii = 1:nbloc;
  aa( jj( ii )+1:jj(ii+1), : ) = repmat( bb( ii, : ), [ sizes( ii ), 1 ] );
end;


