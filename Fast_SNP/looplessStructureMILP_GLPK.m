function modelMILP = looplessStructureMILP_GLPK(model,Nint)
% Builds loopless structure problem (glpk)
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
K                    = 1e3;
modelMILP.c          = zeros(n+2*m,1);
modelMILP.sense      = -1;               % Minimization

% Constraints definition
modelMILP.A = full([model.S,zeros(l,2*m);...
    zeros(p,n),modelMILP.Nint,zeros(p,m);...
    zeros(m,n),eye(m),(1+K)*eye(m);...
    zeros(m,n),-eye(m),-(1+K)*eye(m);...
    M,zeros(m),-diag(model.ub(modelMILP.internal));...
    -M,zeros(m),-diag(model.lb(modelMILP.internal))]);

% RHS model
modelMILP.b = [zeros(l+p,1);K*ones(m,1);-ones(m,1);zeros(m,1);-model.lb(modelMILP.internal)];

% Sign assignation
modelMILP.ctype = blanks(l+p+4*m);
for i = 1:l+p+4*m
    if i <= l+p
        modelMILP.ctype(i) = 'S';
    else
        modelMILP.ctype(i) = 'U';
    end
end

% Assignation of variable types
modelMILP.vartype = blanks(n+2*m);
for i = 1:n+2*m
    if i <= n+m
        modelMILP.vartype(i) = 'C';
    else
        modelMILP.vartype(i) = 'I';
    end
end

% Bounds assignation 
modelMILP.lb  = [model.lb(1:n);-K*ones(m,1);zeros(m,1)];
modelMILP.ub  = [model.ub(1:n);K*ones(m,1);ones(m,1)];