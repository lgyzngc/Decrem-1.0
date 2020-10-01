% Example1_ll_FVA.m: Performs ll-FVA in Ecoli core using Fast-SNP
% Pedro Saa 2016 UQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load model
load('ecoli_core.mat');

% (1) Run fast ll-FVA using Fast-SNP: No objFxn, alpha = 0, use any solver
% available
[model_fast1,fvaRange_fast1] = fastLooplessFVA(model);

% (2) Run fast ll-FVA without using Fast-SNP: No objFxn, alpha = 0, use any solver
% available
[model_norm1,fvaRange_norm1] = fastLooplessFVA(model,[],0,'normal');

% (3) Run fast ll-FVA using Fast-SNP: Growth max, alpha = 0.8, use any solver
% available
[model_fast2,fvaRange_fast2] = fastLooplessFVA(model,model.c,.8);

% (4) Run fast ll-FVA without using Fast-SNP: Growth max, alpha = 0.8, use any solver
% available
[model_norm2,fvaRange_norm2] = fastLooplessFVA(model,model.c,.8,'normal');

% (5) Compare final results
figure (1)
subplot(1,2,1)
plot(fvaRange_norm1(:,2)-fvaRange_norm1(:,1),fvaRange_fast1(:,2)-fvaRange_fast1(:,1),'+b');
xlabel('Flux range normal approach')
ylabel('Flux range with Fast-SNP')
title(['\alpha = 0, Time with pre-proc: ',num2str(model_fast1.prepTime),', Time without pre-proc: ',num2str(model_norm1.prepTime)])
subplot(1,2,2)
plot(fvaRange_norm2(:,2)-fvaRange_norm2(:,1),fvaRange_fast2(:,2)-fvaRange_fast2(:,1),'+b');
xlabel('Flux range normal approach')
ylabel('Flux range with Fast-SNP')
title(['\alpha = 0.8, Time with pre-proc: ',num2str(model_fast2.prepTime),', Time without pre-proc: ',num2str(model_norm2.prepTime)])