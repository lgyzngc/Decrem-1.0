function [model,rxnOrder] = parseInternalRxns(model)
% Re-arrange internal and exchange reactions
% Inputs:  model (COBRA structure)
% Outputs: model (optim structure)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
internalRxns = []; exchangeRxns = [];
for ix = 1:size(model.S,2)
    if sum(full(model.S(:,ix)~=0)) > 0
%     if sum(full(model.S(:,ix)~=0)) > 1
        internalRxns = [internalRxns,ix];
    else
        exchangeRxns = [exchangeRxns,ix];
    end
end

% Reorganize reaction order
model.S        = model.S(:,[internalRxns,exchangeRxns]);
model.rxns     = model.rxns([internalRxns,exchangeRxns]);
model.rxnNames = model.rxnNames([internalRxns,exchangeRxns]);
model.lb       = model.lb([internalRxns,exchangeRxns]);
model.ub       = model.ub([internalRxns,exchangeRxns]);
model.internal = 1:length(internalRxns);
model.exchange = length(internalRxns)+1:length(internalRxns)+length(exchangeRxns);
% model.binaryMatrix = sparse(model.binaryMatrix(:,[internalRxns,exchangeRxns]));
rxnOrder       = [internalRxns,exchangeRxns]';