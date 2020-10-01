function [model,fvaRange] = fastLooplessFVA(model,obj_vector,alpha,method,solverName)
% Performs conventional or fast ll-FVA. In the latter case, Fast-SNP is
% employed to find a suitable basis for Nint
% Inputs:  model structure, obj_vector (optional), alpha (% to optimality, optional)
% Outputs: model structure with new bounds and directionalities, ll-FVA
%          ranges
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs
if nargin<1
    disp('Not enough input arguments.'); return;
elseif nargin<2
    obj_vector = [];           % No obj fxn defined
    alpha  = 0;                % %-optimality
    method = 'fast';
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
elseif nargin<3
    alpha  = 0;                % %-optimality
    method = 'fast';
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
elseif nargin<4
    method = 'fast';
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
elseif nargin<5
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
end

% Assign parameters to the appropiate solver
if strcmp(solverName,'gurobi')
    options.outputflag      = 0;
    options.OptimalityTol   = 1e-9;
    options.FeasibilityTol  = 1e-6;
    options.IntFeasTol      = 1e-6;
    options.MIPGapAbs       = 1e-4;
    options.MIPGap          = 1e-3;
    %     options.TimeLimit       = 10*60;    % 10 min time limit (uncomment if desired)
elseif strcmp(solverName,'cplex')
    options.Display         = 'off';
    options.TolFun          = 1e-9;
    options.TolRLPFun       = 1e-6;
    options.TolXInteger     = 1e-4;
    %     options.MaxTime         = 10*60;
else
    options.tolint          = 1e-4;
    options.tolobj          = 1e-6;
    options.tolbnd          = 1e-6;
    options.toldj           = 1e-6;
    %     options.tmlim           = 10*60;
end
tol = 1e-8;

% Main script
t0 = cputime;

% Determine sparse basis for Nint
if strcmp('fast',method)
    disp('Performing ll-FVA using Fast-SNP...');
    Nint  = fastSNP(model,solverName);
    [model,rxnOrder] = parseInternalRxns(model);
else
    disp('Performing traditional ll-FVA...');
    [model,rxnOrder] = parseInternalRxns(model);
    Nint  = null(full(model.S(:,model.internal)),'r');
    Nint(abs(Nint)<tol) = 0; % Remove elements below numerical tolerance
end

%Initialize FVA ranges
fvaRange = zeros(size(model.S,2),2);

%% Set-up problem (GUROBI)
if strcmp(solverName,'gurobi')
    model = looplessStructureMILP_GUROBI(model,Nint);
    fobj  = zeros(size(model.obj));
    
    % (1) Run initial FBA problem if indicated
    if any(obj_vector)
        model.obj(obj_vector(rxnOrder)~=0) = abs(obj_vector(obj_vector~=0));
        if sum(obj_vector)>0
            model.modelsense = 'max';
        else
            model.modelsense = 'min';
        end
        sol = gurobi(model,options);
        
        % Fix the value to alpha% of the optimal
        if sum(obj_vector)>0
            model.lb(model.obj~=0) = sol.objval*alpha;
            model.ub(model.obj~=0) = sol.objval;
        else
            model.lb(model.obj~=0) = sol.objval*(2-alpha); % This is in the case the optimal val is negative
            model.ub(model.obj~=0) = sol.objval;
        end
    end
    model.obj = fobj;
    
    % (2) ll-FVA optimization loop
    for i = 1:model.numRxns
        model.obj(i) = 1;
        
        % Minimization problem
        model.modelsense = 'min';
        sol           = gurobi(model,options);
        fvaRange(i,1) = sol.objval*(abs(sol.objval)>tol);
        
        % Maximization problem
        model.modelsense = 'max';
        sol           = gurobi(model,options);
        fvaRange(i,2) = sol.objval*(abs(sol.objval)>tol);
        
        % Restart objective
        model.obj = fobj;
    end
    
    %% Set-up problem (CPLEX)
elseif strcmp(solverName,'cplex')
    model = looplessStructureMILP_CPLEX(model,Nint);
    fobj  = zeros(size(model.f));
    
    % Run initial FBA problem if indicated
    if any(obj_vector)
        model.f(obj_vector(rxnOrder)~=0) = -obj_vector(obj_vector~=0); % deals with either max or min
        [~,fval] = cplexmilp(model);

        % Fix the value to alpha% of the optimal
        if sum(obj_vector)>0
            model.lb(model.f~=0) = -fval*alpha;
            model.ub(model.f~=0) = -fval;
        else
            model.lb(model.f~=0) = -fval*(2-alpha);
            model.ub(model.f~=0) = -fval;
        end
    end
    model.f = fobj;
    
    % (2) ll-FVA optimization loop
    for i = 1:model.numRxns
        
        % Minimization problem
        model.f(i)    = 1;
        [~,fval]      = cplexmilp(model);
        fvaRange(i,1) = fval*(abs(fval)>tol);
        
        % Maximization problem
        model.f(i)    = -1;
        [~,fval]      = cplexmilp(model);
        fvaRange(i,2) = -fval*(abs(fval)>tol);
        
        % Restart objective
        model.f = fobj;
    end
    
    %% Set-up problem (GLPK)
elseif strcmp(solverName,'glpk')
    model = looplessStructureMILP_GLPK(model,Nint);
    fobj  = zeros(size(model.c));
    
    % Run initial FBA problem if indicated
    if any(obj_vector)
        model.c(obj_vector(rxnOrder)~=0) = abs(obj_vector(obj_vector~=0)); % deals with either max or min
        [~,fval] = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);

        % Fix the value to alpha% of the optimal
        if sum(obj_vector)>0
            model.lb(model.c~=0) = fval*alpha;
            model.ub(model.c~=0) = fval;
        else
            model.lb(model.c~=0) = fval*(2-alpha);
            model.ub(model.c~=0) = fval;
        end
    end
    model.c = fobj;
    
    % (2) ll-FVA optimization loop
    for i = 1:model.numRxns
        model.c(i)    = 1;
        
        % Minimization problem
        model.sense   = 1;
        [~,fval]      = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);
        fvaRange(i,1) = fval*(abs(fval)>tol);
        
        % Maximization problem
        model.sense   = -1;
        [~,fval]      = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);
        fvaRange(i,2) = fval*(abs(fval)>tol);
        
        % Restart objective
        model.c = fobj;
    end
end

% Save ll-FVA time
elapsedTime = (cputime-t0);

% Find zero reactions
zeroRxns = find(abs(fvaRange(:,1))+abs(fvaRange(:,2))==0);
if isempty(zeroRxns)
    
    % Reset objective function
    model.obj = fobj;
    
    % Redefine model boundaries
    model.lb(1:model.numRxns) = fvaRange(:,1);
    model.ub(1:model.numRxns) = fvaRange(:,2);    
else
    % Remove zero rxns
    modelTemp.S = model.S;
    modelTemp.S(:,zeroRxns) = [];
    modelTemp.c = zeros(size(modelTemp.S,2),1);
    Nint        = Nint(~ismember(model.internal,zeroRxns),:);
    
    % Find orphan metabolites
    orphanMets = find(sum(full(abs(modelTemp.S)),2)==0);
    modelTemp.S(orphanMets,:) = [];
    
    % Re-define other quantities
    modelTemp.b                    = zeros(size(modelTemp.S,1),1);
    modelTemp.lb                   = fvaRange(:,1);
    modelTemp.lb(zeroRxns)         = [];
    modelTemp.ub                   = fvaRange(:,2);
    modelTemp.ub(zeroRxns)         = [];
    modelTemp.rxns                 = model.rxns;
    modelTemp.rxns(zeroRxns)       = [];
    modelTemp.rxnNames             = model.rxnNames;
    modelTemp.rxnNames(zeroRxns)   = [];
    modelTemp.mets                 = model.mets;
    modelTemp.mets(orphanMets)     = [];
    modelTemp.metNames             = model.metNames;
    modelTemp.metNames(orphanMets) = [];
    modelTemp.description          = model.description;
    
    % Re-format the problem
    if strcmp(solverName,'gurobi')
        model = looplessStructureMILP_GUROBI(parseInternalRxns(modelTemp),Nint);
    elseif strcmp(solverName,'cplex')
        model = looplessStructureMILP_CPLEX(parseInternalRxns(modelTemp),Nint);
    elseif strcmp(solverName,'glpk')
        model = looplessStructureMILP_GLPK(parseInternalRxns(modelTemp),Nint);
    end
end

% Assign reversibilities to thermodynamically allowable rxns
model.rev         = (model.lb(1:model.numRxns)<0).*(model.ub(1:model.numRxns)>0);
model.prepTime    = elapsedTime;
disp(['Done. Reactions removed ',num2str(length(zeroRxns))]);