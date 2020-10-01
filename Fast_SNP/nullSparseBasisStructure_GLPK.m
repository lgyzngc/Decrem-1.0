function modelSNP = nullSparseBasisStructure_GLPK(modelSNP,x)
% Defines model structure for Fast-SNP
% Inputs:     model structure
%
% Outputs:    model structure for Fast-SNP (glpk)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter initialization
delta=1.0e-7;
maxFlux=1e0;
scale_para=1e3;
[m,n] = size(modelSNP.S);

% Definition of model structure
modelSNP.A = full([modelSNP.S,zeros(m,2*n);...
                -eye(n),eye(n),zeros(n);...
                eye(n),eye(n),zeros(n);...
                eye(n).*maxFlux,zeros(n),-1*eye(n);...
                -1*eye(n).*maxFlux,zeros(n),eye(n);...
                zeros(size(modelSNP.binaryMatrix,1),2*n),modelSNP.binaryMatrix*-1;...
                x,zeros(1,2*n)]);

% Redefine LP problem:
% (1) RHS
modelSNP.b = zeros(size(modelSNP.A,1),1);
% modelSNP.b(m+2*n+1) = 1e-3; %% add by gaoyang li
modelSNP.b(m+2*n+1:m+3*n) = -1+delta;
modelSNP.b(m+4*n+1:end-1) = -1;
% (2) Constraints
modelSNP.ctype = blanks(size(modelSNP.A,1));
for ix = 1:size(modelSNP.A,1)
    if ix<=m
        modelSNP.ctype(ix) = 'S';
    else
        modelSNP.ctype(ix) = 'L';
    end
end

% (3) Bounds
% modelSNP.lb = [modelSNP.lb;-1e3*ones(n,1);zeros(n,1)];
modelSNP.lb = [modelSNP.lb;zeros(n,1);zeros(n,1)];
modelSNP.ub = [modelSNP.ub/scale_para;1e3*ones(n,1)/scale_para;ones(n,1)];

% (4) Model sense
modelSNP.sense = 1;       % Minimize

% (5) Variables
modelSNP.vartype = blanks(3*n);
for ix = 1:3*n
    if ix<=2*n
        modelSNP.vartype(ix) = 'C';
    else
        modelSNP.vartype(ix) = 'B';
    end
end

% (6) Obj fxn
modelSNP.c = [zeros(n,1);ones(n,1);zeros(n,1)];

% Remove conflicting information
modelSNP = rmfield(modelSNP,{'obj','rhs','modelsense','vtype'});