function modelMILP = looplessStructureMILP_CPLEX(model,Nint)
% Builds loopless structure problem
% Inputs:  model structure, Nint
% Outputs: loopless model structure
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l,n] = size(model.S);

% Reorganize the stoichiometric matrix as [ìnternalRxns|exchangeRxns]
modelMILP.numRxns  = n;
modelMILP.internal = model.internal;
modelMILP.exchange = model.exchange;
modelMILP.S        = model.S;
modelMILP.rxns     = model.rxns;
modelMILP.rxnNames = model.rxnNames;
modelMILP.mets     = model.mets;
modelMILP.metNames = model.metNames;
modelMILP.rev      = (model.lb<0).*(model.ub>0);
modelMILP.description = model.description;

% Calculate null basis of Sint (if necessary)
modelMILP.Nint = Nint';
[p,m]          = size(modelMILP.Nint);

% Determine constrained internal reactions
M = eye(n);
M(modelMILP.exchange,:) = [];

% Define binary constant, model objective and sense
K                 = 1e3;
modelMILP.f       = zeros(1,n+2*m);
modelMILP.sostype = [];              % Additional parameters
modelMILP.sosind  = [];
modelMILP.soswr   = [];
modelMILP.x0      = [];

% Constraints definition
modelMILP.Aeq   = full([model.S,zeros(l,2*m);zeros(p,n),modelMILP.Nint,zeros(p,m)]);
modelMILP.Aineq = full([zeros(m,n),eye(m),(1+K)*eye(m);zeros(m,n),-eye(m),-(1+K)*eye(m);...
                        M,zeros(m),-diag(model.ub(modelMILP.internal));...
                        -M,zeros(m),-diag(model.lb(modelMILP.internal))]);

% RHS model
modelMILP.beq   = zeros(l+p,1);
modelMILP.bineq = [K*ones(m,1);-ones(m,1);zeros(m,1);-model.lb(modelMILP.internal)];


% Assignation of variable types
modelMILP.ctype = blanks(n+2*m);
for i = 1:n+2*m
    if i <= n+m
        modelMILP.ctype(i) = 'C';
    else
        modelMILP.ctype(i) = 'B';
    end
end

% Bounds assignation 
modelMILP.lb  = [model.lb(1:n);-K*ones(m,1);zeros(m,1)];
modelMILP.ub  = [model.ub(1:n);K*ones(m,1);ones(m,1)];