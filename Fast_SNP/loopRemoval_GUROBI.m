function modelMIQP = loopRemoval_GUROBI(model,queryPoint)
% Loop Removal (LR): finds closest feasible point to the query point
% Inputs:  modelMILP (MILP structure)
% Outputs: modelMIQ  (MIQP formulation)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign previous structure
modelMIQP = model;
n = model.numRxns;
m = model.internal(end);
K = 1e3;

% Add distances for flux vectors
modelMIQP.A   = sparse([model.A,zeros(size(model.A,1),n);...
                eye(n),zeros(n,2*m),-eye(n)]);
modelMIQP.obj = zeros(size(modelMIQP.A,2),1);
modelMIQP.Q   = sparse([zeros(n+2*m,2*n+2*m);zeros(n,n+2*m),eye(n)]);
modelMIQP.rhs = [model.rhs;queryPoint];

% Update remaining fields
modelMIQP.lb         = [modelMIQP.lb;-K*ones(n,1)];
modelMIQP.ub         = [modelMIQP.ub;K*ones(n,1)];
modelMIQP.vtype      = [modelMIQP.vtype,repmat('C',1,n)];
modelMIQP.sense      = [modelMIQP.sense,repmat('=',1,n)];
modelMIQP.modelsense = 'min';