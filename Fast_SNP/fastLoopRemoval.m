function [vloopless,modelMIQP] = fastLoopRemoval(model,queryFlux,solverName)
% Finds closest feasible point to the query point (not supported in glpk)
% Inputs:  model structure, queryPoint (possibly unfeasible flux vector),
%          solver Name
% Outputs: closest feasible flux vector, model structure (MIQP)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Check inputs
if nargin<2
    disp('Not enough input arguments.'); return;
elseif nargin<3
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
       disp('No suitable solver found'); return; 
    end
end

% Assign parameters to the appropiate solver
if strcmp(solverName,'gurobi')
    options.outputflag      = 0;
    options.OptimalityTol   = 1e-9;
    options.FeasibilityTol  = 1e-9;
    options.IntFeasTol      = 1e-6;
    options.MIPGapAbs       = 1e-6;
    options.MIPGap          = 1e-6;
    %     options.TimeLimit       = 1*60*60;    % 1 h time limit (uncomment if desired)
elseif strcmp(solverName,'cplex')
    options.Display         = 'off';
    options.TolXInteger     = 1e-4;
    %     options.MaxTime         = 1*60*60;
end
tol = 1e-7;

% Determine whether Nint is present
if ~isfield(model,'Nint')
    model.Nint = fastSNP(model,solverName);
    model      = parseInternalRxns(model);
end

% Define solver and solve MIQP
t0  = cputime;
if strcmp(solverName,'gurobi')
    modelMILP = looplessStructureMILP_GUROBI(model,model.Nint);
    modelMIQP = loopRemoval_GUROBI(modelMILP,queryFlux);
    sol       = gurobi(modelMIQP,options);
    if strcmp('OPTIMAL',sol.status)
        vloopless = sol.x(1:modelMIQP.numRxns).*(abs(sol.x(1:modelMIQP.numRxns))>tol);
    else
        disp('No optimal value found');
        vloopless = [];
    end
elseif strcmp(solverName,'cplex')
    modelMILP = looplessStructureMILP_CPLEX(model,model.Nint);
    modelMIQP = loopRemoval_CPLEX(modelMILP,queryFlux);
    modelMIQP.options = options;
    [xopt,~,exitFlag] = cplexmiqp(modelMIQP);
    if (exitFlag==1)
        vloopless = xopt(1:modelMIQP.numRxns).*(abs(xopt(1:modelMIQP.numRxns))>tol);
    else
        disp('No optimal value found');
        vloopless = [];
    end
end

% Assing final time
modelMIQP.removalTime = (cputime-t0);
