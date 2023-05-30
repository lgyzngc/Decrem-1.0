%% loading data
File_geneExpession = fopen('.\kinetic_data\PEMC_geneExpression.txt');
geneExpession = [];
strainName = {};
geneName={};
index=0;
while(~feof(File_geneExpession))  
    line = fgetl(File_geneExpession);
    str = regexp(line,'\t','split');
    index = index +1;
    if index== 1
        for i=2:length(str)
            strainName{i-1} = str{i};
        end
    else 
        geneName{index-1} = str{1};
        temp = [];
        for i=2:length(str)
            temp = [temp,str2num(str{i})];
        end
        geneExpession = [geneExpession;temp];
    end
end
fclose(File_geneExpession);

File_metaboliteExpession = fopen('.\kinetic_data\\PEMC_metaboliteConcentration.txt');
metaboliteExpession = [];
metaboliteName={};
index=0;
while(~feof(File_metaboliteExpession))
    line = fgetl(File_metaboliteExpession);
    str = regexp(line,'\t','split');
    index = index +1;
    if index== 1
        continue;
    else 
        metaboliteName{index} = str{1};
        temp = [];
        for i=2:length(str)
            temp = [temp,str2num(str{i})];
        end
        metaboliteExpession = [metaboliteExpession;temp];
    end
end
fclose(File_metaboliteExpession);

File_metabolicFlux = fopen('.\kinetic_data\PEMC_metabolicFlux.txt');
metabolicFlux = [];
fluxName={};
index=0;
while(~feof(File_metabolicFlux))
    line = fgetl(File_metabolicFlux);
    str = regexp(line,'\t','split');
    index = index +1;
    if index== 1
        continue;
    else 
        fluxName{index} = str{1};
        temp = [];
        for i=2:length(str)
            temp = [temp,str2num(str{i})];
        end
        metabolicFlux = [metabolicFlux;temp];
    end
end
fclose(File_metabolicFlux);
metabolicFlux(1:49,:) = metabolicFlux(1:49,:).*(ones(49,1)*metabolicFlux(50,:));


File_externalMetabolite = fopen('.\kinetic_data\PEMC_externalMetabolite.txt');
externalMetabolite = [];
externalMetaboliteName={};
index=0;
while(~feof(File_externalMetabolite))
    line = fgetl(File_externalMetabolite);
    str = regexp(line,'\t','split');
    index = index +1;
    if index== 1
        continue;
    else 
        externalMetaboliteName{index} = str{1};
        temp = [];
        for i=2:length(str)
            temp = [temp,str2num(str{i})];
        end
        externalMetabolite = [externalMetabolite;temp];
    end
end
fclose(File_externalMetabolite);

%% data analysis 
%% the metabolite indexes are identified from the heapmap clustering ayalysis
selectedSampleCoeff = [];
[coeff,score,latent,tsquared,explained,mu] = pca(metaboliteExpession([2,6:7,11,21,26:29,31,33:34,40,42,44:45],:)')
temp1 = corr(score,externalMetabolite')
index=length(metabolicFlux)-2:length(metabolicFlux);
temp2 = corr(score,metabolicFlux(index,:)')
selectedSampleCoeff = [selectedSampleCoeff;temp2(1,:),temp1(1,:)];

[coeff,score,latent,tsquared,explained,mu] = pca(metaboliteExpession([1,3:5,15,17,18,20,23,24,30,36,37,43],:)')
temp1 = corr(score,externalMetabolite')
index=length(metabolicFlux)-2:length(metabolicFlux);
temp2 = corr(score,metabolicFlux(index,:)')
selectedSampleCoeff = [selectedSampleCoeff;temp2(1,:),temp1(1,:)];
% random smapling 100 times for File_metabolicFlux and AMP cluster
% and computing the corr(random samples, File_metabolicFlux(index,:)'),
% built the distribution of correlation of random samples and present the p
% value of the selected grobal growth related gene set.
sample_matrix = [];
for i=1:10000
    sample_matrix(i,:) = randperm(size(metaboliteExpession,1),16);
end
% sample_matrix = randint(100,16,[1,size(metaboliteExpession,1)]);
randomSampleCoeff = [];
for i=1:size(sample_matrix,1)
    [coeff,score,latent,tsquared,explained,mu] = pca(metaboliteExpession(sample_matrix(i,:),:)');
    temp1 = corr(score(:,1),metabolicFlux(index,:)');
    temp2 = corr(score(:,1),externalMetabolite');
    randomSampleCoeff(i,:) = [temp1,temp2];
end
[density,freq] = ksdensity(randomSampleCoeff(:,1));
figure(1)
plot(freq,density, 'bo')
figure(2)
[mu,sigma] = normfit(randomSampleCoeff(:,1))
randx=sort(randn(1,500)*sqrt(sigma)+mu);
d=pdf('norm',randx,mu,sigma);
para1=plot(randx,d,'r',[0.8134,0.8134],[0,0.2],'g',[-0.1328,-0.1328],[0,0.2],'b')
set(gca,'Fontsize',15,'Fontname','Timesnewroman','FontWeight','bold','linew',2);
set(para1,'LineWidth',2);
p_value = cdf('norm',0.8134,mu,sigma,'upper')
p_value = cdf('norm',-0.1328,mu,sigma,'upper')

pls_response_coeff = [];
pls_loading_coeff = [];
pls_bata = [];
for i=1:length(geneName)
%     [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(metaboliteExpession([2,6:7,11,21,26:29,31,33:34,40,42,44:45],:)',geneExpession(i,:)',16);
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(metaboliteExpession([2,6:7,11,23:45],:)',geneExpession(i,:)',10);
%      [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(metaboliteExpession([11,23:26,28:31,33:44],:)',geneExpession(i,:)',16);
%    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(metaboliteExpession([1,3:5,15,17,18,20,23,24,30,36,37,43],:)',geneExpession(i,:)',10);

    pls_response_coeff = [pls_response_coeff;PCTVAR(2,:)];
    pls_loading_coeff = [pls_loading_coeff;PCTVAR(1,:)];
    pls_bata = [pls_bata';BETA']';
end
index = (find(pls_response_coeff(:,1)>=0.38))'
index_total = find(sum(pls_response_coeff')>=0.84);
%index = (find(pls_response_coeff(:,2)>=0.33))'
%index_total = find(sum(pls_response_coeff')>=0.7);
index = intersect(index,index_total);
index_rev = setdiff(1:85,index);
selected_metabolicFlux = [4,9,10,11,21,24];

[coeff,score,latent,tsquared,explained,mu] = pca(metaboliteExpession([2,6:7,11,23:45],:)');
temp11 = corr(score,geneExpession(index,:)');
[coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(geneExpession(index,:)');
temp12 = corr(score,score1);
temp13 = corr(score,geneExpession')
[coeff2,score2,latent2,tsquared2,explained2,mu2] = pca(geneExpession');

[coeff21,score21,latent21,tsquared21,explained21,mu21] = pca(metaboliteExpession');
temp21 = corr(score21,geneExpession(index,:)');
[coeff31,score31,latent31,tsquared31,explained31,mu31] = pca(metaboliteExpession([1:10,12:22,27,32,45],:)');
temp31 = corr(score31,geneExpession(index,:)');

% canonical correlation analysis
% [A,B,r,U,V] = canoncorr(score(:,1:7),score2(:,1:7));
 [A,B,r,U,V] = canoncorr(score21(:,1:2),score2(:,1:7));
% [A,B,r,U,V] = canoncorr(score31(:,1:3),score2(:,1:7));
 [A,B,r,U,V] = canoncorr(score(:,1:7),score1(:,1:10));
% [A,B,r,U,V] = canoncorr(score21(:,1:3),score1(:,1:7));


t = tiledlayout(2,2);
title(t,'Canonical Scores of  metabolite vs gene expression')
xlabel(t,'Canonical Variables of metabolite')
ylabel(t,'Canonical Variables of gene expression')
t.TileSpacing = 'compact';

nexttile
plot(U(:,1),V(:,1),'.')
xlabel('u1')
ylabel('v1')

nexttile
plot(U(:,2),V(:,1),'.')
xlabel('u2')
ylabel('v1')

nexttile
plot(U(:,1),V(:,2),'.')
xlabel('u1')
ylabel('v2')

nexttile
plot(U(:,2),V(:,2),'.')
xlabel('u2')
ylabel('v2')

[XL_combin,YL_combin,XS_combin,YS_combin,BETA_combin,PCTVAR_combin,MSE_combin,stats_combin] = plsregress(metaboliteExpession([2,6:7,11,23:45],:)',geneExpession(index,:)',10);
predicted_gene_Expression = [];
for i=1:length(index)
    yfit1 = [ones(size(metaboliteExpession([2,6:7,11,23:45],:)',1),1),metaboliteExpession([2,6:7,11,23:45],:)']*pls_bata(:,index(i));
    yfit2 = [ones(size(metaboliteExpession([2,6:7,11,23:45],:)',1),1),metaboliteExpession([2,6:7,11,23:45],:)']*BETA_combin(:,i);
    corr(yfit1,yfit2,'type','Spearman')
    corr(geneExpession(index(i),:)',yfit1,'type','Spearman')
    corr(geneExpession(index(i),:)',yfit2,'type','Spearman')
    figure(2+i)
    plot(geneExpession(index(i),:),yfit1','r+',geneExpession(index(i),:),yfit2','b*');
    predicted_gene_Expression = [predicted_gene_Expression;yfit2'];
end
for i=1:length(index)
    disp(geneName{index(i)})
end

related_reaction_index = [1,10,4,6,11,7,12,9,24,25,18,21,29,8];
% extra_flux=[0.00E+00	0.00E+00	0.00E+00	2.35E-02	0.00E+00	9.00E-03	0.00E+00	1.50E-03	1.80E-02	1.65E-02	0.00E+00	0.00E+00	0.00E+00	0.00E+00	9.50E-02	0.00E+00	9.00E-03	3.24E-01	0.00E+00	7.25E-02	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	8.00E-03	4.90E-01	2.48E-01
% 6.49E-03	3.71E-03	5.90E-03	4.17E-03	3.22E-03	5.34E-03	3.50E-03	4.64E-03	3.69E-03	0.00E+00	0.00E+00	8.94E-03	4.95E-03	6.10E-03	0.00E+00	1.44E-02	5.06E-03	8.50E-03	5.80E-03	0.00E+00	9.50E-03	0.00E+00	4.90E-03	5.32E-03	4.79E-03	2.42E-02	1.64E+00	3.12E+00
% 0.00E+00	2.13E-04	1.42E-04	0.00E+00	0.00E+00	0.00E+00	4.98E-04	9.24E-04	4.98E-04	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	3.20E-03	0.00E+00	1.05E-02	0.00E+00	0.00E+00	0.00E+00	0.00E+00	3.44E-03	0.00E+00	0.00E+00	0.00E+00	0.00E+00	6.68E-03
% 0.00E+00	1.36E-03	0.00E+00	2.15E-03	9.87E-04	1.29E-03	7.89E-04	1.29E-03	2.51E-03	2.87E-04	4.30E-03	2.75E-03	1.08E-03	1.87E-03	6.34E-03	6.82E-03	1.15E-03	3.37E-02	3.09E-02	1.14E-02	2.15E-03	0.00E+00	4.30E-03	1.55E-03	1.00E-03	1.00E-03	0.00E+00	1.79E-04
% ];
extra_flux=[5.99E-03	7.43E-03	8.96E-03	1.27E-01	9.83E-03	1.04E-02	7.24E-03	8.79E-02	4.65E-03	4.94E-03	1.22E-02	1.12E-02	1.97E-03	0.00E+00	9.59E-03	6.95E-03	1.51E-02	3.45E-03	9.30E-03	9.01E-03	5.70E-03	8.53E-03	9.92E-03	9.54E-03	3.31E-03	4.55E-03	4.07E-03	6.79E-01
0.00E+00	0.00E+00	0.00E+00	2.35E-02	0.00E+00	9.00E-03	0.00E+00	1.50E-03	1.80E-02	1.65E-02	0.00E+00	0.00E+00	0.00E+00	0.00E+00	9.50E-02	0.00E+00	9.00E-03	3.24E-01	0.00E+00	7.25E-02	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	8.00E-03	4.90E-01	2.48E-01
6.49E-03	3.71E-03	5.90E-03	4.17E-03	3.22E-03	5.34E-03	3.50E-03	4.64E-03	3.69E-03	0.00E+00	0.00E+00	8.94E-03	4.95E-03	6.10E-03	0.00E+00	1.44E-02	5.06E-03	8.50E-03	5.80E-03	0.00E+00	9.50E-03	0.00E+00	4.90E-03	5.32E-03	4.79E-03	2.42E-02	1.64E+00	3.12E+00
0.00E+00	2.13E-04	1.42E-04	0.00E+00	0.00E+00	0.00E+00	4.98E-04	9.24E-04	4.98E-04	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	3.20E-03	0.00E+00	1.05E-02	0.00E+00	0.00E+00	0.00E+00	0.00E+00	3.44E-03	0.00E+00	0.00E+00	0.00E+00	0.00E+00	6.68E-03
0.00E+00	1.36E-03	0.00E+00	2.15E-03	9.87E-04	1.29E-03	7.89E-04	1.29E-03	2.51E-03	2.87E-04	4.30E-03	2.75E-03	1.08E-03	1.87E-03	6.34E-03	6.82E-03	1.15E-03	3.37E-02	3.09E-02	1.14E-02	2.15E-03	0.00E+00	4.30E-03	1.55E-03	1.00E-03	1.00E-03	0.00E+00	1.79E-04
2.27E-02	5.85E-03	2.92E-02	7.31E-03	8.77E-03	8.14E-02	0.00E+00	0.00E+00	2.44E-03	1.61E-02	0.00E+00	1.46E-02	4.30E-03	0.00E+00	0.00E+00	0.00E+00	1.95E-03	2.27E-02	0.00E+00	3.98E-02	0.00E+00	0.00E+00	1.15E-02	1.04E-02	6.82E-03	1.61E-02	0.00E+00	0.00E+00
6.28E-03	4.95E-03	2.73E-03	2.92E-03	1.78E-03	1.93E-03	4.06E-03	2.92E-03	6.59E-03	4.06E-03	2.06E-03	3.39E-03	1.32E-03	4.05E-03	4.29E-03	0.00E+00	1.88E-03	8.67E-03	7.60E-03	9.50E-03	6.19E-03	4.21E-03	0.00E+00	1.65E-03	1.29E-03	2.18E-03	3.67E-03	9.02E-03
0.00E+00	1.95E-02	0.00E+00	1.68E-02	2.42E-04	0.00E+00	0.00E+00	0.00E+00	5.23E-02	0.00E+00	0.00E+00	8.31E-03	2.82E-04	0.00E+00	0.00E+00	1.93E-02	2.59E-02	3.62E-02	8.04E-03	1.94E-03	5.73E-03	4.36E-03	4.84E-04	1.21E-03	0.00E+00	1.31E-02	4.22E-01	1.29E-01
];
corr(metabolicFlux(related_reaction_index,:)',geneExpession(index,:)')
extra_index=[10,11,12,14,15,22,26,27,31,32];
corr(extra_flux',geneExpession(index(extra_index),:)')
[XL_flux,YL_flux,XS_flux,YS_flux,BETA_flux,PCTVAR_flux,MSE_flux,stats_flux] = plsregress(geneExpession(index(extra_index),:)',extra_flux',8);
%% the Michaelis-Menten equation fit the metabolite data
Accoa_metabolite = [1.85E-01	2.79E-02	2.87E-02	5.70E-01	6.90E-02	1.27E-01	1.60E-01	1.23E-01	1.62E-01	8.38E-02	4.93E-01	3.48E-02	6.27E-02	2.42E-01	7.10E-02	5.43E-02	5.43E-01	2.45E-01	1.19E-01	2.32E-01	2.89E-01	2.24E-01	5.61E-01	3.71E-01	2.47E-01	6.47E-02	1.53E-01	3.84E-01];
pyruate = [0.00E+00	1.35E-01	2.58E-01	0.00E+00	0.00E+00	2.54E-01	1.53E-01	1.13E-01	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	2.70E-01	4.78E-01	2.23E-01	1.91E-01	2.54E-01	1.30E-01	2.59E-01	2.66E-01	0.00E+00	0.00E+00	0.00E+00	2.69E-01];
% metaboliteExpession = [metaboliteExpession;Accoa_metabolite];

aceE = predicted_gene_Expression(10,:);
adhE = predicted_gene_Expression(26,:);
ackA = predicted_gene_Expression(27,:);
ldha = predicted_gene_Expression(22,:);
% aceE = geneExpession(index(10),:);
% adhE = geneExpession(index(26),:);
% ackA = geneExpession(index(12),:);
% ldha = geneExpession(index(11),:);
AMP = metaboliteExpession(15,:);
ATP = metaboliteExpession(14,:);

% nonzero_index1=find(extra_flux(3,:)~=0);
nonzero_index=find(extra_flux(3,:)~=0);
nonzero_index=nonzero_index(1:end);
% nonzero_index = intersect(nonzero_index,nonzero_index1);
nonzero_acetate = extra_flux(3,nonzero_index);
Accoa_metabolite = Accoa_metabolite(nonzero_index);
aceE = aceE(nonzero_index);
adhE = adhE(nonzero_index);
ackA = ackA(nonzero_index);
ldha = ldha(nonzero_index);
AMP = AMP(nonzero_index);
ATP = ATP(nonzero_index);
pyruate = pyruate(nonzero_index);
% the michaelis menten equation of acetate i.e log(extra_flux(ace)) = log(aceE) +
% log(Accoa_metabolite) - log(km+Accoa_metabolite) - ki*log(1+ATP/AMP) +
% bata; km=18.26mM

train_data_label=log(nonzero_acetate)-log(aceE)-log(Accoa_metabolite)+log(mean(Accoa_metabolite)+Accoa_metabolite);
train_data=[ones(1,length(AMP));log(1./(1+pyruate));log(AMP./(AMP+ATP))];
[b, bint,r,rint,stats]=regress(train_data_label',train_data');
corr(exp(train_data_label)',exp(b'*train_data+log(aceE)+log(Accoa_metabolite)-log(mean(Accoa_metabolite)+Accoa_metabolite))')
figure(40)
plot(train_data_label,(b'*train_data+log(aceE)+log(Accoa_metabolite)-log(mean(Accoa_metabolite)+Accoa_metabolite)),'r*')

train_data_label=log(nonzero_acetate);
% train_data=[ones(1,length(AMP));log(1./(1+pyruate));log(AMP./(AMP+ATP));log(aceE);log(adhE);log(ackA);log(ldha);log(Accoa_metabolite);log(1.0./(mean(Accoa_metabolite)+Accoa_metabolite))];
train_data=[ones(1,length(AMP));log(aceE);log(adhE);log(ackA);log(ldha);log(Accoa_metabolite)];
[b, bint,r,rint,stats]=regress(train_data_label',train_data')
% [XL_flux,YL_flux,XS_flux,YS_flux,BETA_flux,PCTVAR_flux,MSE_flux,stats_flux] = plsregress(train_data',train_data_label',8);
corr(exp(train_data_label)',exp(train_data'*b))
figure(42)
plot(train_data_label',(train_data'*b),'r*',-7:0.5:2,-7:0.5:2,'b--')

% Ethonal predicting with adhE
pyruate = [0.00E+00	1.35E-01	2.58E-01	0.00E+00	0.00E+00	2.54E-01	1.53E-01	1.13E-01	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	0.00E+00	2.70E-01	4.78E-01	2.23E-01	1.91E-01	2.54E-01	1.30E-01	2.59E-01	2.66E-01	0.00E+00	0.00E+00	0.00E+00	2.69E-01];
Accoa_metabolite = [1.85E-01	2.79E-02	2.87E-02	5.70E-01	6.90E-02	1.27E-01	1.60E-01	1.23E-01	1.62E-01	8.38E-02	4.93E-01	3.48E-02	6.27E-02	2.42E-01	7.10E-02	5.43E-02	5.43E-01	2.45E-01	1.19E-01	2.32E-01	2.89E-01	2.24E-01	5.61E-01	3.71E-01	2.47E-01	6.47E-02	1.53E-01	3.84E-01];

adhE = predicted_gene_Expression(26,:);
ackA = predicted_gene_Expression(27,:);
aceE = predicted_gene_Expression(10,:);
AMP = metaboliteExpession(15,:);
ATP = metaboliteExpession(14,:);

nonzero_index=find(pyruate~=0);
nonzero_ethonal = extra_flux(2,nonzero_index);
nonzero_ethonal=1+nonzero_ethonal;
Accoa_metabolite = Accoa_metabolite(nonzero_index);
adhE = adhE(nonzero_index);
ackA = ackA(nonzero_index);
aceE = aceE(nonzero_index);
AMP = AMP(nonzero_index);
ATP = ATP(nonzero_index);
pyruate = pyruate(nonzero_index);

train_data_label=log(nonzero_ethonal);
train_data=[ones(1,length(AMP));log(1./(1+pyruate));log(AMP./(AMP+ATP));log(adhE);log(ackA);log(aceE);log(Accoa_metabolite);log(1.0./(mean(Accoa_metabolite)+Accoa_metabolite))];
[b, bint,r,rint,stats]=regress(train_data_label',train_data')
corr(exp(train_data_label)',exp(train_data'*b))
figure(43)
plot(train_data_label',(train_data'*b),'r*',-1:0.5:1,-1:0.5:1,'b--')

% The Glucose flux predicting
glucose = extra_flux(1,:);
glu_index=find(glucose~=0);
glucose = glucose(glu_index);
g6p = metaboliteExpession(10,glu_index);
AMP = metaboliteExpession(15,glu_index);
ATP = metaboliteExpession(14,glu_index);
ptsH = predicted_gene_Expression(31,glu_index);
ptsL = predicted_gene_Expression(32,glu_index);
add1 = predicted_gene_Expression(1,glu_index);
add2 = predicted_gene_Expression(2,glu_index);
add3 = predicted_gene_Expression(29,glu_index);
add4 = predicted_gene_Expression(30,glu_index);
glu_flux = metabolicFlux(1,glu_index);

% 
train_data_label=log(glu_flux);
train_data=[ones(1,length(AMP));log(AMP./(AMP+ATP));log(ptsH);log(ptsL);log(add1);log(add2);log(add3);log(add4);log(glucose);log(g6p./glucose);log(1.0./(glucose./mean(glucose)+g6p./mean(g6p)))];
[b, bint,r,rint,stats]=regress(train_data_label',train_data')
corr(exp(train_data_label)',exp(train_data'*b))
figure(44)
plot(train_data_label',(train_data'*b),'r*',4.5:0.5:7.5,4.5:0.5:7.5,'b--')

