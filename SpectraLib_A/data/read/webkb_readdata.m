% reads asymmetric Web Graph Data stored in webkbdata.mat file. 
% See Marina Meila, William Pentney 'Clustering by weighted cuts in directed graphs'.
%
% Please set up the parameters in the section between EDIT and END EDIT
% before running the script!
%
% OUTPUT: 
% A - asymmetric affinity matrix Aij=S' (indegrees are actually
% important...)
% k - number of clusters
% asig_true - row vector of true assignments
% Do, Di, ido0, idi0
%
% Authors: Marina Meila, Tatiana Maravina
% Part of SpectraLib_A
% Last revision: 06-June-2007
%


%%% -------------EDIT--------------%%%

% Initialization
istranspose=1; % if 1 then transpose the affinity matrix (meaning that in-links are actually of importance)
is01=1; % if 1 then substitute all non-zero elements with 1

%add path to the directory with the WebKB data
addpath /costila/speclust/Data/WebKB; 
load webkbdata.mat

%%% -----------END EDIT------------%%%

[dumi, dumj, dums ] = find(S);  % use in-links, 1/0 weights
A = S;
if istranspose, A = A'; end 
if is01, A(A>0)=1; end
A = A-diag( diag( A )); 
Do = sum( A, 2 );  % out degrees
Di = sum( A, 1 );  % in degrees

% remove nodes with no connections
idkeep = find( Di+Do' > 0 );
A = A( idkeep, idkeep );

% find connected components in A and keep only the largest one
idkeep2=getlargestcomponent(A);
A = A( idkeep2, idkeep2 );
Di = Di( idkeep( idkeep2 ) );  
Do = Do( idkeep( idkeep2 ) );
ido0 = find( Do == 0 );% nodes with out degrees equal 0
idi0 = find( Di == 0 );% nodes with in degrees equal 0
asig_true = cluniv( idkeep( idkeep2 ))'; % need row vector
k = length(unique(asig_true));  % nr clusters (or classes) 


% clear temporary variables
clear idkeep idkeep2;


