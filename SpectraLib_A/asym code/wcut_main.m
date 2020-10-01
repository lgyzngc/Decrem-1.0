% Runs wcut with many different T matrices, and selects the result by 
% distortion
% 
% Please, set up the parameters in the sections between EDIT and END EDIT
% before running the script
%
% $Authors: Marina Meila, Tatiana Maravina
% $Part of SpectraLib_A
% $Last revision: 06-June-2007
%

%%% -------------EDIT--------------%%%
%  Initializations
isviewdata = 0;   % display various statistics
isreaddata = 1;
fontsize= 10;
eps_diag = 1;

if isreaddata 
   webkb_readdata;
   % webkb_readdata creates:
   % A - affinity matrix
   % Do, Di - out, in-degrees
   % ido0, idi0 - nodes with zero out/in degrees
   % asig_true - true classes assignment
   % k - number of classes
end;
n=size(A,1);


%  What should T be? Define it by strings to be evaluated 
%  why am I doing this? to have 1 loop try all the T's

%Dpowers = [ -1.5 -1 -0.85 -0.75 -0.5 -0.35 -0.25 -0.1 0 0.1 0.25 0.35 0.5 0.75 0.85 1 1.5 2 2.5 3 ];
%Dpowers = [  0 0.1 0.25 0.35 0.5 0.75 0.85 1 1.1 1.2 1.3  1.4 1.5 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5];
Dpowers= [0 1 1.5];
npowers = length( Dpowers );

%%% -----------END EDIT------------%%%


Tcommand_all = {};
for ipow = 1:npowers;
    ttt = { ['T = Dofill.^(' num2str( Dpowers( ipow )) ');']};
    Tcommand_all = [ Tcommand_all, ttt ];
end;
ttt = { 'T = min( Dofill, 5 );'};
Tcommand_all = [ Tcommand_all,  ttt ];
nalgs = length( Tcommand_all );


%%% -------------EDIT--------------%%%
onlyalgs = ones( 1, nalgs );  % run only selected algorithms
%onlyalgs = [0 0 0 0 0 0 0 1 0 0 0 0  ];
%%% -----------END EDIT------------%%%


ialgs = find( onlyalgs );
if max( ialgs ) > nalgs 
   disp( '*******  Error ***** inexistent algorithm ' );
end;


vi_all = zeros( 1, nalgs );
ce_all = zeros( 1, nalgs );
distort_all = zeros( 1, nalgs );
separation_all = zeros( 1, nalgs );
wcut_all =  zeros( 1, nalgs );
gaps_all =  zeros( 1, nalgs );
gapegap_all =  zeros( 1, nalgs );

iifill = zeros( size( Do));
iifill( ido0 ) = 1;
Dofill = Do+eps_diag*iifill;

opts.disp = 0;
opts.issim = 1;
opts.p = 2*k;

kv = k;   % use e-vectors 1:kv in clustering, kv <= k


% Loop over selected algorithms (T's)
global_options; % global constants from Spectrul Clustering library
for ia = 1:nalgs
    if onlyalgs( ia )
       Tcommand = char(Tcommand_all( ia ));
       eval( Tcommand );
       sqrtTinv = 1./sqrt( T );
       B = (diag(Do)-A).* repmat( sqrtTinv, [1, n] ).*repmat( sqrtTinv', [n, 1] );
       H = (B + B')/2;

       % Eigenproblems
       disp( [ num2str( ia ) '...now computing eigenvectors of H...'] );
       [ yy, ll ] = eigs( H, k+1, 'sa',  opts );
       ll = diag( ll );
       vv = yy(:,1:kv).*repmat( sqrtTinv, [1, kv ] );

       % Cluster
       [ center, asig, distort_all( ia ) ] = cluster_normalized_kmeans( vv, k);
       centerc = center - center*ones(k,k)/k; % "centered centers"
       separation_all(ia) = distort_all(ia)/n/sum( diag( centerc'*centerc))*k;
       ce_all( ia ) = clustering_error( asig, asig_true );
       vi_all( ia ) = VI( asig, asig_true );
       yasig = zeros( n, k );
       for ik = 1:k;
           iik = find( asig == ik );
           yasig( iik, ik ) = 1;
           yasig( :, ik ) = yasig( :, ik )./norm(yasig( :,ik ));
       end;
       wcut_all( ia ) = trace( yasig'*B*yasig );
       gaps_all( ia ) = wcut_all( ia ) - sum( ll( 1:k));
       gapegap_all( ia ) = gaps_all( ia )/(ll(k+1)-ll(k));
    end; %if onlyalgs
end; % for ia


% Display results for all algorithms
disp( Tcommand_all( ialgs ));
disp( ['VI:     ', num2str( vi_all( ialgs )) ]);
disp( ['CE:     ', num2str( ce_all( ialgs )) ]);
disp( ['distort:', num2str( distort_all( ialgs ))]);

disp('Now choose minimum distortion:');
[ distmin, algmin ] = min( distort_all );
disp( [char(Tcommand_all( algmin )), ' VI= ', num2str( vi_all( algmin )), '  CE= ', num2str( ce_all( algmin ))]);

[ cemin0, algmin0 ] = min( ce_all );
if ce_all( algmin0 ) < ce_all( algmin )
   disp( [ '* * * imperfect selection: best CE = ' num2str( cemin0, 3 )] );
end; 

Tcommandmin = Tcommand_all( algmin );
cemin = ce_all( algmin );
vimin = vi_all( algmin );

clf;
hp = plot( Dpowers, distort_all( 1:npowers )/max(distort_all)*max(ce_all), ':o', Dpowers, ce_all( 1:npowers ), 'r-x', Dpowers, gapegap_all( 1:npowers )/max(gapegap_all(gapegap_all<Inf))*max(ce_all), 'k:s', Dpowers, separation_all(1:npowers)/max(separation_all)*max(ce_all), 'm:v', Dpowers, gaps_all(1:npowers)/max(gaps_all)*max(ce_all), 'c:^' );
hl = legend( 'distortion (scaled)', 'CE', 'gap/eigengap (scaled)', 'separation (scaled)', 'gap (scaled)' );
set( gca, 'FontSize', fontsize );


% TO DO: Clear intermediate variables


