function modelLD = loopDetection_CPLEX(model,queryFlux)
% Formulates Loop Detection (LD) problem (cplex)
% Inputs:  model structure, queryPoint, solver Name
% Outputs: modelLD  (LP formulation)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Check inputs
if nargin<1
    disp('Not enough inputs.');
end

% Check size
[m,n] = size(model.Nint);

% Model structure definition
modelLD.f     = zeros(1,n);
modelLD.Aeq   = full(model.Nint); % Loop condition
modelLD.Aineq = [];
modelLD.bineq = [];
modelLD.lb    = -1e3*(1+sign(queryFlux(model.internal)))-sign(queryFlux(model.internal));
modelLD.ub    =  1e3*(1-sign(queryFlux(model.internal)))-sign(queryFlux(model.internal));

% Right-hand side definition
modelLD.beq = zeros(m,1);

% Initial point
modelLD.x0 = [];