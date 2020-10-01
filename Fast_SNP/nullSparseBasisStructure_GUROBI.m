function modelSNP = nullSparseBasisStructure_GUROBI(modelSNP,x)
% Defines model structure for Fast-SNP
% Inputs:     model structure
%
% Outputs:    model structure for Fast-SNP (gurobi)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter initialization
[m,n] = size(modelSNP.S);

% Definition of model structure
modelSNP.A = sparse([modelSNP.S,zeros(m,n);...
    -eye(n),eye(n);...
    eye(n),eye(n);...
    x,zeros(1,n)]);

% Redefine LP problem:
% (1) RHS
modelSNP.rhs = zeros(m+2*n+1,1);

% (2) Sense
modelSNP.sense = blanks(m+2*n+1);
for ix = 1:m+2*n+1
    if ix<=m
        modelSNP.sense(ix) = '=';
    else
        modelSNP.sense(ix) = '>';
    end
end
% add quadric constrint
% QUC=modelSNP.binaryMatrix;
% for i=1:size(QUC,1)
%     index=find(QUC(i,:)==1);
%     temp = eye(2*n)*0.00000100001;
% %     temp(index(1),:)=0;
% %     temp(index(2),:)=0;
%     temp(index(1),index(2))=100;
% %      temp(index(2),index(1))=1/2;
%     
%     modelSNP.quadcon(i).Qc = sparse(temp);
%     modelSNP.quadcon(i).q = zeros(2*n,1);
%     modelSNP.quadcon(i).rhs = 0.0;
%     
% %     modelSNP.quadcon(2*i-1).Qc = sparse(1000*temp);
% %     modelSNP.quadcon(2*i-1).q = zeros(2*n,1);
% %     modelSNP.quadcon(2*i-1).rhs = 0.0;
% %     
% %     modelSNP.quadcon(2*i).Qc = sparse(-1000*temp);
% %     modelSNP.quadcon(2*i).q = zeros(2*n,1);
% %     modelSNP.quadcon(2*i).rhs = 0.0;
% end

% (3) Bounds
modelSNP.lb = [modelSNP.lb;-1e3*ones(n,1)];
% modelSNP.lb = [modelSNP.lb;zeros(n,1)];
modelSNP.ub = [modelSNP.ub;1e3*ones(n,1)];

% (4) Model sense
modelSNP.modelsense = 'min';

% (5) Variables
modelSNP.vtype = 'C';

% (6) Obj fxn
modelSNP.obj = [zeros(n,1);ones(n,1)];