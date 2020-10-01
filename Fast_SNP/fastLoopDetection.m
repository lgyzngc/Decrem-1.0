function [feasibleFlux,modelLP] = fastLoopDetection(model,queryFlux,solverName)
% Determines whether the query flux is feasible (in the loopless sense)
% Inputs:  model structure, queryPoint (possibly unfeasible flux vector),
%          solver Name
% Outputs: feasibleFlux (boolean, 1 flux is feasible, 0 unfeasible)
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
    options.outputflag     = 0;
    options.OptimalityTol  = 1e-9;
    options.FeasibilityTol = 1e-9;
    options.IntFeasTol     = 1e-9;
elseif strcmp(solverName,'cplex')
    options.Display        = 'off';
    options.TolFun         = 1e-9;
    options.TolRLPFun      = 1e-6;
elseif strcmp(solverName,'glpk')
    options.tolobj         = 1e-9;
    options.tolbnd         = 1e-8;
    options.toldj          = 1e-8;
end

% Determine whether Nint is present
if ~isfield(model,'Nint')
    model.Nint = fastSNP(model,solverName)';
    model      = parseInternalRxns(model);
end

% Define solver and solve LD
feasibleFlux = 0;
t0  = cputime;
if strcmp(solverName,'gurobi')
    modelLP = loopDetection_GUROBI(model,queryFlux);
    sol     = gurobi(modelLP,options);
    if strcmp('OPTIMAL',sol.status)
        feasibleFlux = 1;    
    end
elseif strcmp(solverName,'cplex')
    modelLP = loopDetection_CPLEX(model,queryFlux);
    modelLP.options = options;
    [~,~,exitFlag]  = cplexlp(modelLP);
    if (exitFlag==1)
        feasibleFlux = 1;
    end
elseif strcmp(solverName,'glpk')
    modelLP = loopDetection_GLPK(model,queryFlux);
    [~,~,exitFlag] = glpk(modelLP.c,modelLP.A,modelLP.b,modelLP.lb,modelLP.ub,modelLP.ctype,modelLP.vartype,modelLP.sense,options);
    if (exitFlag==5)
        feasibleFlux = 1;
    end
end

% Assing final time
modelLP.detectionTime = (cputime-t0);
