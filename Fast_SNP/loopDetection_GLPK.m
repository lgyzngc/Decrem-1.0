function modelLD = loopDetection_GLPK(model,queryFlux)
% Formulates Loop Detection (LD) problem (glpk)
% Inputs:  model structure, queryPoint, solver Name
% Outputs: modelLD  (LP formulation)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs
if nargin<1
    disp('Not enough inputs.');
end

% Check size
[m,n] = size(model.Nint);

% Model structuredefinition
modelLD.c     = zeros(n,1);
modelLD.sense = 1;                  % Minimization
modelLD.A     = full(model.Nint);   % Loop condition
modelLD.lb    = -1e3*(1+sign(queryFlux(model.internal)))-sign(queryFlux(model.internal));
modelLD.ub    =  1e3*(1-sign(queryFlux(model.internal)))-sign(queryFlux(model.internal));

% Sign assignation (stays the same)
modelLD.ctype = blanks(m);
for ix = 1:m
    modelLD.ctype(ix) = 'S';
end

% Assignation of variable types (stays the same)
modelLD.vartype = blanks(n);
for ix = 1:n
    modelLD.vartype(ix) = 'C';
end

% Right-hand side definition
modelLD.b = zeros(m,1);