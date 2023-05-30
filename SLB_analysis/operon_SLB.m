load('.\iML1515.mat','-mat');
load('.\linear_model_iML1515.mat','-mat');
genes_file = fopen('.\NC_000913.ptt');
genes_hash=java.util.Hashtable;
genes_name_hash=java.util.Hashtable;
while(~feof(genes_file))
    line=fgetl(genes_file);
    S = regexp(line, '\t', 'split');
    if length(S) < 9
        continue;
    end
    if isempty(regexp(S{1},'Location','match'))
        genes_hash.put(S{4},S{6});
        genes_name_hash.put(S{5},S{6});
    end
end
gene_operon = java.util.Hashtable;
operon = fopen('.\NC_000913.opr');
while(~feof(operon))
    line=fgetl(operon);
    S = regexp(line, '\t', 'split');
    if length(S)>1
        gene_temp = regexp(strtrim(S{2}), ' ', 'split');
        for i=1:length(gene_temp)
            if genes_hash.containsKey(gene_temp{i})
                gene_operon.put(genes_hash.get(gene_temp{i}),S{1});
            end
        end
    end
end
gene_TF = java.util.Hashtable;
TF_file = fopen('network_tf_gene.txt');
while(~feof(TF_file))
    line=fgetl(TF_file);
    if line(1) == '#'
        continue;
    end
    S = regexp(line, '\t', 'split');
    if length(S)>1
        if isempty(strtrim(S{4}))
            continue;
        end
        gene_temp = regexp(strtrim(S{4}), ',', 'split');
        for i=1:length(gene_temp)
            if genes_name_hash.containsKey(strtrim(gene_temp{i}))
                gene_TF.put(genes_name_hash.get(strtrim(gene_temp{i})),S{2});
            end
        end
    end
end
reaction_operon={};
reaction_index=1;
TF_index=1;
sparse_linear_basis={};
unique_TFs={};
for i=1:length(linear_model.basisVectorIndex)
    temp = linear_model.basisVectorIndex{i};
    temp_operon = {};
    index = 1;
    if length(temp)>1
        for j=1:length(temp)
            genes = iML1515.grRules{temp(j)};
            if isempty(genes)
                continue;
            end
            genes = regexp(genes,'b\d+','match');
            for k=1:length(genes)
                %if gene_operon.containsKey(genes{k})
                if gene_TF.containsKey(genes{k})
                    temp_operon{index,1} = genes{k};
                    %temp_operon{index,2} = gene_operon.get(genes{k});
                    temp_operon{index,2} = gene_TF.get(genes{k});
                    unique_TFs{TF_index}=gene_TF.get(genes{k});
                    index=index+1;
                    TF_index=TF_index+1;
                end
            end
        end
        LSB_rec_temp=iML1515.rxns(temp);
        LSB_temp='';
        for i=1:length(LSB_rec_temp)
            if i<length(LSB_rec_temp)
                LSB_temp=[LSB_temp,LSB_rec_temp{i},','];
            else
                LSB_temp=[LSB_temp,LSB_rec_temp{i}];
            end
        end
        sparse_linear_basis{reaction_index}=LSB_temp;
        reaction_operon{reaction_index}=temp_operon;
        reaction_index=reaction_index+1;
    end
end
unique_TFs=unique(unique_TFs);
unique_TFs_hash = java.util.Hashtable;
for i=1:length(unique_TFs)
    unique_TFs_hash.put(unique_TFs{i},i);
end

sparse_TF={};
TF_SBL=zeros(length(sparse_linear_basis),length(unique_TFs));
for i=1:length(reaction_operon)
    temp_tf = reaction_operon{i};
    if isempty(temp_tf)
        continue;
    end
    temp_unique = unique(temp_tf(:,2));
    if length(temp_unique) == 1
        TF_SBL(i,unique_TFs_hash.get(temp_unique))=length(temp_tf(:,2));
        continue;
    end
    for j=1:length(temp_unique)
        temp_num = length(find(strcmp(temp_tf(:,2),temp_unique{j})));
        TF_SBL(i,unique_TFs_hash.get(temp_unique{j}))=temp_num;
    end
end
table_tf = array2table(TF_SBL,'VariableNames',unique_TFs,'RowNames',sparse_linear_basis);
writetable(table_tf,'TF_SLB.csv','WriteRowNames',true);