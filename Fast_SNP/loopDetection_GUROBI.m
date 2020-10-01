function modelLD = loopDetection_GUROBI(model,queryFlux)
% Formulates Loop Detection (LD) problem (gurobi)
% Inputs:  model structure, queryPoint, solver Name
% Outputs: modelLD  (LP formulation)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs
if nargin<1
    disp('Not enough inputs.');
end

% Check size
[m,n] = size(model.Nint);

% Model structuredefinition
modelLD.obj        = zeros(n,1);
modelLD.modelsense = 'min';
modelLD.A          = sparse(model.Nint); % Loop condition
modelLD.lb         = -1e3*(1+sign(queryFlux(model.internal)))-sign(queryFlux(model.internal));
modelLD.ub         =  1e3*(1-sign(queryFlux(model.internal)))-sign(queryFlux(model.internal));

% Sign assignation (stays the same)
modelLD.sense = '=';

% Assignation of variable types (stays the same)
modelLD.vtype = 'C';

% Right-hand side definition
modelLD.rhs = zeros(m,1);