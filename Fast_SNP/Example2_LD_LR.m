% Example2_LD_LR.m: Performs loop detection in iJO1366 Fast-SNP
% Pedro Saa 2016 UQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load model
load('ecoli_core');

% Create a random flux distribution
queryFlux = rand(87,1);

% (1) Run LR using any solver available
[vloopless,modelMIQP] = fastLoopRemoval(model,queryFlux);

% (2) Check whether the solutions is loopless
feasibleFlux = fastLoopDetection(modelMIQP,vloopless);
if feasibleFlux; disp('The flux vector is loopless'); end;