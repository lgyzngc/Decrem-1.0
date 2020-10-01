function modelSNP = nullSparseBasisStructure_CPLEX(modelSNP)
% Defines model structure for Fast-SNP
% Inputs:     model structure
%
% Outputs:    model structure for Fast-SNP (cplex)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter initialization
[m,n] = size(modelSNP.S);

% Definition of model structure
modelSNP.Aeq   = full([modelSNP.S,zeros(m,n)]);
modelSNP.Aineq = full([eye(n),-eye(n);-eye(n),-eye(n);zeros(1,2*n)]);

% Redefine LP problem:
% RHS
modelSNP.beq   = zeros(m,1);
modelSNP.bineq = zeros(2*n+1,1);

% Bounds
modelSNP.lb = [modelSNP.lb;-1e3*ones(n,1)];
modelSNP.ub = [modelSNP.ub;1e3*ones(n,1)];

% Obj fxn
modelSNP.f = [zeros(1,n),ones(1,n)];

% Initial point
modelSNP.x0 = [];