function modelMIQP = loopRemoval_CPLEX(modelMILP,queryPoint)
% Loop Removal (LR): finds closest feasible point to the query point
% Inputs:  modelMILP (MILP structure)
% Outputs: modelMIQ  (MIQP formulation)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
% Assign previous structure
modelMIQP = modelMILP;
n = modelMILP.numRxns;
m = modelMILP.internal(end);
K = 1e3;

% Add distances for flux vectors
modelMIQP.Aeq   = full([modelMILP.Aeq,zeros(size(modelMILP.Aeq,1),n);...
                  eye(n),zeros(n,2*m),-eye(n)]);
modelMIQP.beq   = [modelMILP.beq;queryPoint];
modelMIQP.Aineq = [modelMILP.Aineq,zeros(size(modelMILP.Aineq,1),n)];
modelMIQP.bineq = [modelMILP.bineq];
modelMIQP.f     = zeros(1,size(modelMIQP.Aeq,2));
modelMIQP.H     = sparse([zeros(n+2*m,2*n+2*m);zeros(n,n+2*m),eye(n)]);

% Update remaining fields
modelMIQP.lb         = [modelMIQP.lb;-K*ones(n,1)];
modelMIQP.ub         = [modelMIQP.ub;K*ones(n,1)];
modelMIQP.ctype      = [modelMIQP.ctype,repmat('C',1,n)];

% Extra params
modelMIQP.sostype = [];    % Additional parameters
modelMIQP.sosind  = [];
modelMIQP.soswr   = [];
modelMIQP.x0      = [];