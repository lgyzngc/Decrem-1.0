function linear_model=decoupledModelConstruct(model,cofactor_path,secrated_path,nutrient_path,general_IO_path)
	% input:
	%	model: metabolic models with mat format
	%	cofactor_path: the file path of cofactors
	%	secrated_path: the file path of secrated reactions
	%	nutrient_path: the file path of nutrient reactions
	%	general_IO_path: the file path of general IO reactions
	%output:
	%	the reconstructed decoupled metabolic model
	usr_path_flag = 0;
	if nargin == 5
		usr_path_flag = 1
	end
	%% read model and cofactor
	load('..\bacillus\iYO844.mat','-mat');
	Ecoli_Model = iYO844;
	if usr_path_flag == 1
		Ecoli_Model = model
	end

	% Ecoli_Model.lb=Ecoli_Model.lb/10000;
	% Ecoli_Model.ub=Ecoli_Model.ub/10000;
	Ecoli_Model.lb(find(Ecoli_Model.lb<=-1000))=-1000;
	Ecoli_Model.ub(find(Ecoli_Model.ub>=1000))=1000;
	for i=1:length(Ecoli_Model.lb)
		if Ecoli_Model.lb(i)==Ecoli_Model.ub(i)
			Ecoli_Model.ub(i)=1000;
			if Ecoli_Model.rev(i)==1
				Ecoli_Model.lb(i)=-1000;
			else
				Ecoli_Model.lb(i)=0;
			end
		end
	end
	cofactor_hash_temp = java.util.Hashtable;
	cofactor_file=fopen('..\bacillus\cofactor.txt');
	if usr_path_flag == 1
		cofactor_file=fopen(cofactor_path);
	end
	while(~feof(cofactor_file))
		cofactor_hash_temp.put(fgetl(cofactor_file),1);
	end
	fclose(cofactor_file);
	cofactor_hash = java.util.Hashtable;
	for i=1:length(Ecoli_Model.mets)
		temp=regexp(Ecoli_Model.mets{i},'[a-zA-Z0-9]+','match');
		if ~isempty(temp{1})
			if cofactor_hash_temp.containsKey(temp{1})
				cofactor_hash.put(Ecoli_Model.mets{i},i);
			end
		end
	end
	%% reaction constraint
	exchange_reaction_index=[];
	input_nutrient={};
	index=1;

	#secrated_metabolite_set={'SUCCt2r','CITt2r'};
	secrated_nutrient_file=fopen('..\bacillus\secrated_bacillus.txt');
	if usr_path_flag == 1
		secrated_nutrient_file=fopen(secrated_path);
	end
	secrated_metabolite_set={};
	while(~feof(secrated_nutrient_file))
		line=fgetl(secrated_nutrient_file);
		if(~isempty(line))
			secrated_metabolite_set{index}=line;
			index=index+1;
		end
	end

	% secrated_metabolite_positive={'ACt2r'};
	% secrated_metabolite_set={};
	secrated_metabolite_positive={};
	input_nutrient_file=fopen('..\bacillus\nutrient_bacillus.txt');
	if usr_path_flag == 1
		input_nutrient_file=fopen(nutrient_path);
	end
	while(~feof(input_nutrient_file))
		line=fgetl(input_nutrient_file);
		if(~isempty(line))
			input_nutrient{index}=line;
			index=index+1;
		end
	end
	
	input_output_nutrient={};
	index=1;
	input_nutrient_file=fopen('..\bacillus\general_IO_bacillus.txt');
	if usr_path_flag == 1
		input_nutrient_file=fopen(general_IO_path);
	end
	while(~feof(input_nutrient_file))
		line=fgetl(input_nutrient_file);
		if(~isempty(line))
			input_output_nutrient{index}=line;
			index=index+1;
		end
	end

	for i=1:length(Ecoli_Model.rxns)
		% shold be modified to adapt the specific model
	%     if ~(isempty(regexp(Ecoli_Model.rxns{i},'\w+tex$','match')) && isempty(regexp(Ecoli_Model.rxns{i},'\w+texi$','match')))
		  if ~(isempty(regexp(Ecoli_Model.rxns{i},'.+t2r$','match')) && isempty(regexp(Ecoli_Model.rxns{i},'.+t$','match')) && isempty(regexp(Ecoli_Model.rxns{i},'.+ti$','match')))

			flag=0;
			exchange_reaction_index=[exchange_reaction_index,i];
			for j=1:length(input_nutrient)
				if strcmp(Ecoli_Model.rxns{i},input_nutrient{j})
	%                 Ecoli_Model.rev(i)=0;
					disp(Ecoli_Model.rxns{i})
					Ecoli_Model.lb(i)=-1;
					flag=1;
					break;
				end
			end
			if flag==1
				continue;
			end
			for j=1:length(input_output_nutrient)
				if strcmp(Ecoli_Model.rxns{i},input_output_nutrient{j})
					disp(Ecoli_Model.rxns{i})
					flag=1;
					break;
				end
			end
			if flag==1
				continue;
			end
			for j=1:length(secrated_metabolite_positive)
				if strcmp(Ecoli_Model.rxns{i},secrated_metabolite_positive{j})
	%                 Ecoli_Model.lb(i)=-640-10;
	%                 Ecoli_Model.ub(i)=-640+10;
					if strcmp(secrated_metabolite_positive{j},'ACTNabc1')
						continue;
					end
					Ecoli_Model.ub(i)=5;
					flag=1;
					break;
				end
			end
			if flag==1
				continue;
			end
			 for j=1:length(secrated_metabolite_set)
				if strcmp(Ecoli_Model.rxns{i},secrated_metabolite_set{j})
					Ecoli_Model.lb(i)=-1;
					Ecoli_Model.ub(i)=1;
					flag=1;
					break;
				end
			end
			if flag==1
				continue;
			end
	%         Ecoli_Model.rev(i)=0;
			Ecoli_Model.ub(i)=1;
		end
	end

	%% construct condition specific network based on nutrient
	lowFlux_reaction=[];
	initCobraToolbox();
	specific_reaction=1:length(Ecoli_Model.rxns);
	for i=1:length(specific_reaction)
		Ecoli_Model_specific.rxns{i}=Ecoli_Model.rxns{specific_reaction(i)};
	end
	Ecoli_Model_specific.mets=Ecoli_Model.mets;
	Ecoli_Model_specific.S=Ecoli_Model.S(:,specific_reaction);
	Ecoli_Model_specific.rev=Ecoli_Model.rev(specific_reaction);
	Ecoli_Model_specific.lb=Ecoli_Model.lb(specific_reaction);
	Ecoli_Model_specific.ub=Ecoli_Model.ub(specific_reaction);
	Ecoli_Model_specific.c=Ecoli_Model.c(specific_reaction);
	if (isfield(Ecoli_Model,'rules'))
		for i=1:length(specific_reaction)
			Ecoli_Model_specific.rules{i}=Ecoli_Model.rules{specific_reaction(i)};
		end
	end
	if (isfield(Ecoli_Model,'rxnGeneMat'))
		Ecoli_Model_specific.rxnGeneMat=Ecoli_Model.rxnGeneMat(specific_reaction,:);
	end

	Ecoli_Model_specific.genes=Ecoli_Model.genes;
	for i=1:length(specific_reaction)
		Ecoli_Model_specific.grRules{i}=Ecoli_Model.grRules{specific_reaction(i)};
	end
	for i=1:length(specific_reaction)
		Ecoli_Model_specific.subSystems{i}=Ecoli_Model.subSystems{specific_reaction(i)};
	end
	for i=1:length(specific_reaction)
		Ecoli_Model_specific.rxnNames{i}=Ecoli_Model.rxnNames{specific_reaction(i)};
	end
	Ecoli_Model_specific.metNames=Ecoli_Model.metNames;
	Ecoli_Model_specific.b=Ecoli_Model.b;

	% Ecoli_Model_specific=Ecoli_Model;
	%% construct connect graph
	network=ConnectGraphConstruct(Ecoli_Model_specific.S,Ecoli_Model_specific.rxns,Ecoli_Model_specific.mets,cofactor_hash,Ecoli_Model_specific.rev);
	%% calculate the simple cycle
	% flag_java=simplecycle_java();
	similarity_matrix=load('.\bacillus\similarity_matrix_5len_rec4.txt');
	IO_reaction_index=[];
	for i=1:size(Ecoli_Model_specific.S,2)
		if length(find(Ecoli_Model_specific.S(:,i)~=0))==1
			IO_reaction_index=[IO_reaction_index,i];
		end
	end
	similarity_matrix(IO_reaction_index,:)=0;
	similarity_matrix(:,IO_reaction_index)=0;
	%% cluster analysis
	 nonzeroReactionSet=[];
	 noClusterReactionSet=[];
	 for i=1:size(similarity_matrix,1)
		 if ~(sum(similarity_matrix(:,i))==0 && sum(similarity_matrix(i,:))==0)
			 nonzeroReactionSet=[nonzeroReactionSet,i];
		 else
			 noClusterReactionSet=[noClusterReactionSet,i];
		 end
	 end
	 
	  similarityMatrix_nozero=zeros(length(nonzeroReactionSet));
	 for i=1:length(nonzeroReactionSet)
		 for j=1:length(nonzeroReactionSet)
			 similarityMatrix_nozero(i,j)=similarity_matrix(nonzeroReactionSet(i),nonzeroReactionSet(j));
		 end
	 end

	  cluster_num=1;
	%   clusterStru=cluster_similarityMatrix(similarityMatrix_nozero,cluster_num);
	  cluster_set={};
	  clusterStru=ones(length(nonzeroReactionSet),1);
	  for i=1:cluster_num
		  temp_index=find(clusterStru==i);
		  cluster_set{i}=nonzeroReactionSet(temp_index);
	  end
	  %% construct sub metabolic model for each sub cluster
	  Stichimetic_Submatrix_set={};
	  reaction_set=1:size(Ecoli_Model_specific.S,2);
	  for i=1:length(cluster_set)
		  temp=Ecoli_Model_specific.S(:,cluster_set{i});
		  zero_rowIndex=[];
		  rest_reaction_set=setdiff(reaction_set,cluster_set{i});
		  for j=1:size(temp,1)
			  if sum(abs(temp(j,:)))==0
				  zero_rowIndex=[zero_rowIndex,j];
			  end
		  end
		  metabolite_inset=1:size(temp,1);
		  metabolite_set=setdiff(metabolite_inset,zero_rowIndex);
		  temp=temp(metabolite_set,:);
		  addtional_colmon=[];
		  index=1;
		  for j=1:length(metabolite_set)
			  input_temp=find(Ecoli_Model_specific.S(metabolite_set(j),rest_reaction_set) > 0);
			  input_test=find(temp(j,:) < 0);
			  if ~(isempty(input_temp) || isempty(input_test))
				temp_colmon=zeros(1,length(metabolite_set));
				temp_colmon(j)=1;
				addtional_colmon=[addtional_colmon;temp_colmon];
			  end
			  
			  output_temp=find(Ecoli_Model_specific.S(metabolite_set(j),rest_reaction_set) < 0);
			  output_test=find(temp(j,:) > 0);
			  if ~(isempty(output_temp) || isempty(output_test))
				temp_colmon=zeros(1,length(metabolite_set));
				temp_colmon(j)=-1;
				addtional_colmon=[addtional_colmon;temp_colmon];
			  end
		  end
		  temp=[temp';addtional_colmon];
		  
		  reversibleReac_index=find(Ecoli_Model_specific.rev(cluster_set{i})==1);
		  temp1=cluster_set{i};
		  temp1=temp1(reversibleReac_index);
		  temp2=Ecoli_Model_specific.S(metabolite_set,temp1)*-1;
		  temp=[temp2';temp];
		  temp=temp';
		  Stichimetic_Submatrix_set{i,1}=temp;
		  Stichimetic_Submatrix_set{i,2}=metabolite_set;
		  Stichimetic_Submatrix_set{i,3}=cluster_set{i};
		  Stichimetic_Submatrix_set{i,4}=temp1;
		  reversible_pair_index=[];
		  reac_index=cluster_set{i};
		  for j=1:length(reversibleReac_index)
			  reversible_pair_index(j,1)=j;
			  reversible_pair_index(j,2)=length(reversibleReac_index)+reversibleReac_index(j);
		  end
		  Stichimetic_Submatrix_set{i,5}=reversible_pair_index;
	  end
	  %% solving the sparse basis vector of sub cluster model based on the L1 approprate method

	  for i=1:size(Stichimetic_Submatrix_set,1)
		  sub_model(i).numRxns=size(Stichimetic_Submatrix_set{i,1},2);
		  sub_model(i).obj=zeros(sub_model(i).numRxns,1);
		  sub_model(i).S=sparse(Stichimetic_Submatrix_set{i,1});
		  sub_model(i).A=sparse(Stichimetic_Submatrix_set{i,1});
		  sub_model(i).rhs=zeros(size(sub_model(i).S,1),1);
		  sub_model(i).InternalReacNum= length(Stichimetic_Submatrix_set{i,3})+ size(Stichimetic_Submatrix_set{i,5},1);
		  sub_model(i).IOReacNum=sub_model(i).numRxns-sub_model(i).InternalReacNum;
		  
		  sub_model(i).lb(1:length(Stichimetic_Submatrix_set{i,4}))=0;
		  sub_model(i).lb(1+length(Stichimetic_Submatrix_set{i,4}):length(Stichimetic_Submatrix_set{i,3})+length(Stichimetic_Submatrix_set{i,4}))=Ecoli_Model.lb(Stichimetic_Submatrix_set{i,3});
		  sub_model(i).lb(length(sub_model(i).lb)+1:sub_model(i).numRxns)=0;
		  sub_model(i).lb(find(sub_model(i).lb<0))=0;
		  sub_model(i).lb=sub_model(i).lb';
		  
		  sub_model(i).ub(1:length(Stichimetic_Submatrix_set{i,4}))=abs(Ecoli_Model.lb(Stichimetic_Submatrix_set{i,4}));
		  sub_model(i).ub(1+length(Stichimetic_Submatrix_set{i,4}):length(Stichimetic_Submatrix_set{i,3})+length(Stichimetic_Submatrix_set{i,4}))=Ecoli_Model.ub(Stichimetic_Submatrix_set{i,3});
		  sub_model(i).ub(length(sub_model(i).ub)+1:sub_model(i).numRxns)=max(sub_model(i).ub);
		  sub_model(i).ub=sub_model(i).ub';
		  
		  temp_pair_reac=Stichimetic_Submatrix_set{i,5};
		  temp_binaryMatrix=zeros(size(temp_pair_reac,1),sub_model(i).numRxns);
		  for j=1:size(temp_pair_reac,1)
			  temp_binaryMatrix(j,temp_pair_reac(j,1))=1;
			  temp_binaryMatrix(j,temp_pair_reac(j,2))=1;
		  end
		  sub_model(i).binaryMatrix=temp_binaryMatrix;
		  
		  
		  sub_model(i).vtype='C';
		  sub_model(i).modelsense='max';
		  sub_model(i).sense='=';
		  
		  temp_index=Stichimetic_Submatrix_set{i,4};
		  for j=1:length(temp_index)
			  sub_model(i).rxns{j}=strcat('rev_',Ecoli_Model.rxns{temp_index(j)});
			  sub_model(i).rxnNames{j}=strcat('rev_',Ecoli_Model.rxnNames{temp_index(j)});
		  end
		  temp_index=Stichimetic_Submatrix_set{i,3};
		  for j=1:length(temp_index)
			  sub_model(i).rxns{j+length(Stichimetic_Submatrix_set{i,4})}=Ecoli_Model.rxns{temp_index(j)};
			  sub_model(i).rxnNames{j+length(Stichimetic_Submatrix_set{i,4})}=Ecoli_Model.rxnNames{temp_index(j)};
		  end
		  temp_index_meta=Stichimetic_Submatrix_set{i,2};
		  for j=length(sub_model(i).rxns)+1:sub_model(i).numRxns
			  metaName=find(sub_model(i).S(:,j)~=0);
			  if length(metaName)==1
				  if sub_model(i).S(metaName(1),j)>0
					  sub_model(i).rxns{j}=strcat('Input_',Ecoli_Model.mets{temp_index_meta(metaName(1))});
					  sub_model(i).rxnNames{j}=strcat('Input_',Ecoli_Model.metNames{temp_index_meta(metaName(1))});
				  else 
					  sub_model(i).rxns{j}=strcat('Output_',Ecoli_Model.mets{temp_index_meta(metaName(1))});
					  sub_model(i).rxnNames{j}=strcat('Output_',Ecoli_Model.metNames{temp_index_meta(metaName(1))});
				  end
			  end
		  end
		  temp_index=Stichimetic_Submatrix_set{i,2};
		  for j=1:length(temp_index)
			  sub_model(i).mets{j}=Ecoli_Model.mets{temp_index(j)};
			  sub_model(i).metNames{j}=Ecoli_Model.metNames{temp_index(j)};
		  end
		  sub_model(i).description=strcat('subnetwork_',num2str(i));
		  sub_model(i).rev=zeros(sub_model(i).numRxns,1);

		  sub_model(i).c(1:length(Stichimetic_Submatrix_set{i,4}))=0;
		  sub_model(i).c(1+length(Stichimetic_Submatrix_set{i,4}):length(Stichimetic_Submatrix_set{i,3})+length(Stichimetic_Submatrix_set{i,4}))=Ecoli_Model.c(Stichimetic_Submatrix_set{i,3});
		  for j=length(sub_model(i).c)+1:sub_model(i).numRxns
			  sub_model(i).c(j)=0;
		  end
	  end
	   % solving the sparse absis vector using the SNP techenique
	  for i=1:length(sub_model)
		Nsnp{i} = fastSNP(sub_model(i),'gurobi');
	  end
	  %% reconstruct new stoichimetic matrix for the sparse basis vector and the rest linear reactions in orginal network
	  reaction_set_basisVector={};
	  reaction_set_basisVector_lb = [];
	  reaction_set_basisVector_ub = [];
	  basisVector_set=[];
	  reaction_basisVector_index={};
	  reaction_basisVector_coff={};
	  repeat_hash=java.util.Hashtable;
	  index=1;
	  for i=1:length(Nsnp)
		  current_basis_vectors=Nsnp{i};
		  
		  metabolite_index=Stichimetic_Submatrix_set{i,2};
		  temp_reactionset=[];
		  internal_reaction_num=sub_model(i).InternalReacNum;
		  
		  reaction_set_index=[];
		  temp_rev=Stichimetic_Submatrix_set{i,5};
		  temp_rev=temp_rev(:,2);
		  temp_rev=temp_rev-length(temp_rev);
		  temp_rev_1=Stichimetic_Submatrix_set{i,3};
		  temp_rev=temp_rev_1(temp_rev);
		  reaction_set_index=[temp_rev,temp_rev_1];
		  reaction_set_rev=zeros(1,length(reaction_set_index));
		  reaction_set_rev(1:length(temp_rev))=1;
		  
		  for j=1:size(current_basis_vectors,2)
			  temp_basisVector=current_basis_vectors(:,j);
			  temp_basisVector=temp_basisVector(1:internal_reaction_num);
			  temp_metabolites=sub_model(i).S(:,1:internal_reaction_num)*temp_basisVector;
			  if sum(temp_metabolites(find(temp_metabolites ~= 0)))<1.0e-8
				  single_reaction=find(temp_basisVector ~= 0);
				  temp_reaction=[];
				  for k=1:length(single_reaction)
					  if repeat_hash.isEmpty || ~repeat_hash.containsKey(reaction_set_index(single_reaction(k)))
						  if single_reaction(k)<=length(temp_rev)
							  repeat_hash.put(reaction_set_index(single_reaction(k)),-1);
							  temp_reaction=[temp_reaction,single_reaction(k)];
						  else
							  repeat_hash.put(reaction_set_index(single_reaction(k)),1);
							  temp_reaction=[temp_reaction,single_reaction(k)];
						  end
					  else
						  if single_reaction(k)<=length(temp_rev) && repeat_hash.get(reaction_set_index(single_reaction(k)))==1
							  temp_reaction=[temp_reaction,single_reaction(k)];
							  repeat_hash.put(reaction_set_index(single_reaction(k)),0);
						  elseif single_reaction(k)> length(temp_rev) && repeat_hash.get(reaction_set_index(single_reaction(k)))==-1
							  temp_reaction=[temp_reaction,single_reaction(k)];
							  repeat_hash.put(reaction_set_index(single_reaction(k)),0);
						  end
					  end
				  end
				  if isempty(temp_reaction)
					  continue;
				  else
					  single_reaction=temp_reaction;
				  end
				  temp_lb=sub_model(i).lb(single_reaction);
				  temp_ub=sub_model(i).ub(single_reaction);
				  reaction_set_basisVector_lb=[reaction_set_basisVector_lb,temp_lb'];
				  reaction_set_basisVector_ub=[reaction_set_basisVector_ub,temp_ub'];
				  
				  single_rev=find(single_reaction<=length(temp_rev));
				  if ~isempty(single_rev)
					single_rev=single_reaction(single_rev);
					rev_reaction=Ecoli_Model_specific.S(:,reaction_set_index(single_rev))*-1;
					single_reaction=setdiff(single_reaction,single_rev);
					temp_reactionset=[temp_reactionset;rev_reaction'];
					for k=1:length(single_rev)
					  reaction_basisVector_index{index}=reaction_set_index(single_rev(k));
					  reaction_basisVector_coff{index}=temp_basisVector(single_rev(k))*-1/abs(temp_basisVector(single_rev(k)));
					  index=index+1;
					end
				  end
				  for_reaction=Ecoli_Model_specific.S(:,reaction_set_index(single_reaction));
				  temp_reactionset=[temp_reactionset;for_reaction'];
				  for k=1:length(single_reaction)
					  reaction_basisVector_index{index}=reaction_set_index(single_reaction(k));
					  reaction_basisVector_coff{index}=temp_basisVector(single_reaction(k))/abs(temp_basisVector(single_reaction(k)));
					  index=index+1;
				  end
				  continue;
			  end
			  temp_reaction_index=find(temp_basisVector ~= 0);
			  temp_reaction=[];
			  if length(temp_reaction_index)==1
				  for k=1:length(temp_reaction_index)
					  if repeat_hash.isEmpty || ~repeat_hash.containsKey(reaction_set_index(temp_reaction_index(k)))
						  if temp_reaction_index(k)<=length(temp_rev)
							  repeat_hash.put(reaction_set_index(temp_reaction_index(k)),-1);
							  temp_reaction=[temp_reaction,temp_reaction_index(k)];
						  else
							  repeat_hash.put(reaction_set_index(temp_reaction_index(k)),1);
							  temp_reaction=[temp_reaction,temp_reaction_index(k)];
						  end
					  else
						  if temp_reaction_index(k)<=length(temp_rev) && repeat_hash.get(reaction_set_index(temp_reaction_index(k)))==1
							  temp_reaction=[temp_reaction,temp_reaction_index(k)];
							  repeat_hash.put(reaction_set_index(temp_reaction_index(k)),0);
						  elseif temp_reaction_index(k)> length(temp_rev) && repeat_hash.get(reaction_set_index(temp_reaction_index(k)))==-1
							  temp_reaction=[temp_reaction,temp_reaction_index(k)];
							  repeat_hash.put(reaction_set_index(temp_reaction_index(k)),0);
						  end
					  end
				  end
				  if isempty(temp_reaction)
					  continue;
				  end
			  end
			  original_temp_reaction_index=temp_reaction_index;
			  temp_reaction_index=reaction_set_index(temp_reaction_index);
			  reaction_basisVector_index{index}=temp_reaction_index;
			  
			  single_rev=find(original_temp_reaction_index<=length(temp_rev));
			  if ~isempty(single_rev)
				  single_rev=original_temp_reaction_index(single_rev);
				  temp_basisVector(single_rev)=temp_basisVector(single_rev)*-1;
			  end
			  reaction_basisVector_coff{index}=temp_basisVector(original_temp_reaction_index);
			  index=index+1;
			  
			  temp_reaction=zeros(size(Ecoli_Model_specific.S,1),1);
			  temp_reaction(metabolite_index)=temp_metabolites;
			  temp_reactionset=[temp_reactionset;temp_reaction'];
			  temp_index=find(temp_basisVector ~= 0);
			  temp_lb=max(sub_model(i).lb(temp_index)./abs(temp_basisVector(temp_index)));
			  temp_ub=min(sub_model(i).ub(temp_index)./abs(temp_basisVector(temp_index)));
			  reaction_set_basisVector_lb=[reaction_set_basisVector_lb,temp_lb];
			  reaction_set_basisVector_ub=[reaction_set_basisVector_ub,temp_ub];
		  end
		  reaction_set_basisVector{i}=temp_reactionset';
	  end
	  
	  linear_matrix=Ecoli_Model_specific.S(:,noClusterReactionSet);
	  linear_matrix=linear_matrix';
	  for i=1:length(reaction_set_basisVector)
		  linear_matrix=[linear_matrix;1*reaction_set_basisVector{i}'];
	  end
	  linear_matrix=linear_matrix';
	  
	  zero_metabolite=[];
	  for i=1:size(linear_matrix,1)
		  if sum(abs(linear_matrix(i,:)))==0
			  zero_metabolite=[zero_metabolite,i];
		  end
	  end
	  metabolite_index=1:size(linear_matrix,1);
	  nonzero_metabolite_index=setdiff(metabolite_index,zero_metabolite);
	  linear_matrix=linear_matrix(nonzero_metabolite_index,:);
	  zero_metaboliteSet={};
	  for i=1:length(zero_metabolite)
		  zero_metaboliteSet{i}=Ecoli_Model_specific.mets{zero_metabolite(i)};
	  end
	  
	  % reconstruct new model
	  for i=1:length(noClusterReactionSet)
		  linear_model.rxns{i}=Ecoli_Model_specific.rxns{noClusterReactionSet(i)};
		  linear_model.rxnNames{i}=Ecoli_Model_specific.rxnNames{noClusterReactionSet(i)};
	  end
	  for i=length(noClusterReactionSet)+1:size(linear_matrix,2)
		  linear_model.rxns{i}=strcat('basis_reaction_',num2str(i));
		  linear_model.rxnNames{i}=strcat('basis_reaction_',num2str(i));
	  end
	  for i=1:length(nonzero_metabolite_index)
		  linear_model.mets{i}=Ecoli_Model_specific.mets{nonzero_metabolite_index(i)};
		  linear_model.metNames{i}=Ecoli_Model_specific.metNames{nonzero_metabolite_index(i)};
	  end
	  linear_model.S=linear_matrix;
	  linear_model.basisVectorIndex=reaction_basisVector_index;
	  linear_model.reaction_basisVector_coff=reaction_basisVector_coff;
	  linear_model.rev=Ecoli_Model_specific.rev(noClusterReactionSet);
	  linear_model.rev(length(noClusterReactionSet)+1:size(linear_matrix,2))=0;
	  linear_model.lb=Ecoli_Model_specific.lb(noClusterReactionSet);
	  linear_model.lb(length(noClusterReactionSet)+1:size(linear_matrix,2))=reaction_set_basisVector_lb;
	  linear_model.ub=Ecoli_Model_specific.ub(noClusterReactionSet);
	  linear_model.ub(length(noClusterReactionSet)+1:size(linear_matrix,2))=reaction_set_basisVector_ub;
	  linear_model.c=Ecoli_Model_specific.c(noClusterReactionSet);
	  linear_model.c(length(noClusterReactionSet)+1:size(linear_matrix,2))=0;
	  biomass_index=0;
	  for i=1:length(linear_model.rxnNames)
		  if ~isempty(regexp(linear_model.rxnNames{i},'biomass|Biomass|BIOMASS','match'))
			  biomass_index=i;
			  break;
		  end
	  end
	  linear_model.c(biomass_index)=1;
	  linear_model.minPathwayStart=length(noClusterReactionSet);
	  linear_model.reserveReac=noClusterReactionSet;
	  linear_model.minPathwayNum=length(linear_model.basisVectorIndex);
	  linear_model.description='linear model';
	  linear_model.b=Ecoli_Model_specific.b(nonzero_metabolite_index);
end