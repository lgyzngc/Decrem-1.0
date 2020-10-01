function [ aa ] = genBlockSumMatrix( sizes, strenghts )

%function [ aa ] = genBlockSumMatrix( sizes, strenghts )
%
% generates a matrix aa partitioned into blocks as given by the vector sizes.
% the row sums in block (i,j)  is equal to strenghts(i,j).
% if strengths is a stochastic matrix so will be aa.
% size( strenghts ) = sum( sizes ) for both rows and columns
%
% $Authors: Marina Meila
% $Part of SpectraLib_A
% $Last revision: 06-June-2007
%

aa = [];
nsizes = length( sizes );

for isc = 1:nsizes;
  aatt = [];
  for isr = 1:nsizes;
    aat = rand( sizes( isr ), sizes( isc) );
    if( isc == isr ) 
      aat = aat + aat';
    end;
    saat = 1./sum( aat, 2 );
    aat = diag( saat ) * aat * strenghts( isr, isc );
    aatt = [ aatt; aat ];
  end;
  aa = [ aa aatt ];
end;  

