function [ aa ] = genBlockACitationMatrix( sizes, strengths, isdiag, isdistrib )

% function [ aa ] = genBlockACitationMatrix( sizes, strengths )
%
% generates an asymmetric block matrix. 
% "Citation" means that aa(i,j)= 0 for i>j
% sizes is the vector containing the size of each block.
%	in a diagonal block, aa(i,j) = 1 w.p strengths( block) 
%	in a supra-diagonal block aa(i,j) = 1 w.p strengths( block) 
%	in a sub-diagonal block aa(i,j) = 0
% strengths is a square matrix containing the above probabilities. its lower
% triangular part is ignored
%
% isdiag = 1  (default value) ==> put 1 on the diagonal of aa
% isdistrib = 1 (default value) ==> make the outdegrees be nonuniformly distributed
%
% EXAMPLE: 
%   sizes=[20 30 20 10];
%   strengths=[0.8 0.1 0 0.1; 0 0.6 0.3 0.1; 0 0 0.7 0.3; 0 0 0 0.5];
%   S=genBlockACitationMatrix(sizes, strengths);
%   % check probabilities: compare with strengths
%   disp(untelescope(untelescope(S',sizes)', sizes)./(sizes'*sizes)); 
%
% $Authors: Marina Meila
% $Part of SpectraLib_A
% $Last revision: 06-June-2007
%

if nargin<3,  isdiag=[]; end
if nargin<4,  isdistrib=[]; end
if isempty(isdiag), isdiag=1; end
if isempty(isdistrib), isdistrib=1; end

nbloc = length( sizes );
aa = zeros( sum( sizes ));
ioff = 0;

for ii = 1:nbloc;
  % Make diagonal block
  nn = sizes( ii );
  if ~isdistrib
    adum = rand( nn );  
    adum = triu(adum < strengths( ii, ii ), 1);
  else
     [adum, dd ] = genLinkBlock( nn, nn, strengths( ii, ii ), [], 1 );
  end;
  if( isdiag )
      adum = adum + eye( nn );
  end;
  aa( (ioff+1):(ioff+nn), (ioff+1):(ioff+nn)) = adum;

  % Make supra diagonal blocks
  joff = ioff+nn;
  for jj = ii+1:nbloc;
      if isdistrib
	     adum = genLinkBlock( nn, sizes( jj ), strengths( ii, jj ), dd, 0 );
      else
	    adum = rand( nn, sizes( jj ));  
        adum = (adum < strengths( ii, jj ));
      end;
      aa( (ioff+1):(ioff+sizes(ii)), (joff+1):(joff+sizes(jj))) = adum;
      joff = joff + sizes( jj );
  end;
  ioff = ioff + sizes( ii );
end;

