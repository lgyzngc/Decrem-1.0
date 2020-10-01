function [AA S]=genPerfectH (sizes, degs, HH)

% ---work in progress---
% return S only for testing... remove later
% add cc as a parameter later
%
% generates asymmetric affinity matrix AA which is 'perfect' for WCut
% asymmetric clustering (as follows from the asymmetric multicut lemma 
% from Marina Meila and William Pentney 
% 'Clustering by weighted cuts in directed graphs')
% 
% algorithm is described in ....
% 
% 
% CONSTISTENCY CONDITIONS:
% 1) sum of all elements in HH is zero; its off-diagonal elements are <0
%
% 2) diag(untelescope(degs', sizes))*Phat=HH (used internally to compute
% Phat)
%
% REMARKS: 
% 1) Troubles with large n (linear program works terribly long - actually
% reaches the MaxIter and stops uncoverged)
%
% 2) Tricky to specify initial parameters to get reasonable results (so
% that by plotting A we could easily see the clusters)
%    for example, should consider same scale (for cc which determines the
%           diagonal elements of S and non-diagonal elements of S...)
%    meaning of HH??? how does it affect the final AA?
%    What the parameters should be for the exactly k first eigen vectors to
%    be PCE...
%
%
% $Authors: Marina Meila, Tatiana Maravina
% $Part of SpectraLib_A
% $Last revision: 06-June-2007
%


%%%%%%% PART I. GENERATE SPECIAL S %%%%%%%
nsizes = length( sizes );
nn = sum( sizes );

Phat=diag(untelescope(degs',sizes))\HH;

% make special initial S0
S0 = rand( nn );
S0 = (S0+S0')/2; % all elements positive
S0=S0-diag(diag(S0)); % zero out diagonal elements
S0=-S0; % all off-diagonal elements are negative, all diagonal are zeros
cc=rand( [nsizes 1])+max(0,diag(Phat));
Phat_new=Phat-diag(cc);

disp('Part I: compute special S');
S = genBlockSumSymMatrix( sizes, Phat_new, [], degs, S0);
% alter diagonal elements
S(logical(eye(nn)))=telescope(cc, sizes).*degs; 


%%%%%%%% PART II. GENERATE H %%%%%%%%%%%%%%

disp('Part II: get AA from S')


%    Linear program to solve
%   
%    variables [ rr xi ]
%
%    min xi
%    s.t. A*rr = D
%        rr_ij-ss_ij-xi <= 0
%        -rr_ij-ss_ij-xi <= 0
%        -xi <= 0


D = sum(S,2);
dum = -triu(S,1);
ss = dum( 1, 2:nn );
for ii = 2:nn-1;
    ss = [ ss  dum( ii, ii+1:nn )];
end

dimr = nn*(nn-1)/2;  % number rr variables

% constraints for rr_ij-ss_ij-xi <= 0
A1 = [ speye( dimr ) -ones( dimr, 1 ) ];

% constraints for -rr_ij-ss_ij-xi <= 0
A2 = [ -speye( dimr ) -ones( dimr, 1 ) ];

% constraint for -xi <= 0
A3=sparse(1, dimr+1,-1,1, dimr+1);

% constraints for Ae*rr = D;
jdum=zeros(2*dimr,1);
idum=jdum;
sdum=jdum;
last=0;
for ii=1:nn-1;
   start=((nn-1)+(nn-ii+1))*(ii-1)/2+1;
   jdum((last+1):(last+2*(nn-ii)),1)=[(start:(start+nn-ii-1))';(start:(start+nn-ii-1))'];
   idum((last+1):(last+2*(nn-ii)),1)=[ii*ones(nn-ii,1); ((ii+1):nn)'];
   sdum((last+1):(last+2*(nn-ii)),1)=[ones(nn-ii,1); -ones(nn-ii, 1)] ;
   last=last+2*(nn-ii);
end
A0=sparse( idum, jdum, sdum, nn, dimr+1 );

% objective vector for sum_ij
c = [ zeros( 1, dimr ) 1 ];

% Linear program
disp('solving linear optimization program...');
xsol = linprog( c, [A1; A2; A3], [ss'; ss'; 0], A0, D );
rr = xsol( 1:dimr );
%xi = xsol( dimr+1 );

% Construct AA
disp('constructing AA...');
rr = rr';
nind = nn;
R = [0 rr( 1:nn-1 )];
for ii = 2:nn-1;
    R = [R; [zeros( 1, ii ) rr( nind:nind+nn-ii-1)]];
    nind = nind+nn-ii;
end;
R = [ R; zeros( 1, nn )];
R = R-R';
AA = R-(S-diag(diag(S)));





