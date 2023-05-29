function network=ConnectGraphConstruct(stochiemstry,reactions,metabolites,metabolite_cofactors,reversible)
    S=[];
    metabolites_unCofactor={};
    index=1;
    for i=1:length(metabolites)
        if ~metabolite_cofactors.containsKey(metabolites{i})
            S(index,:)=stochiemstry(i,:);
            metabolites_unCofactor{index}=metabolites{i};
            index=index+1;
        end
    end
    metabolite_reaction_connect={};
    index=length(metabolites_unCofactor);
    index_reac=length(reactions);
    for i=1:length(metabolites_unCofactor)
        temp=find(S(i,:)>0);
        temp_rever=find(reversible(temp)==1);
        if ~isempty(temp_rever)
            temp_rever=temp_rever+index_reac;
            temp=[temp temp_rever'];
        end
        temp=temp+index;
        metabolite_reaction_connect{i}=temp;
    end

    for i=1:length(reactions)
        metabolite_reaction_connect{i+index}=find(S(:,i)<0);
        if reversible(i)==1
            metabolite_reaction_connect{i+index+index_reac}=find(S(:,i)>0);
        end
    end
%     metabolite_reaction_connect=sparse(metabolite_reaction_connect);

    network.metabolite_reaction_connect=metabolite_reaction_connect;
    network.metabolite=metabolites_unCofactor;
    network.reaction=reactions;
    adject_list=fopen('metabolites_reaction_network.txt','w');
    for i=1:length(metabolite_reaction_connect)
        temp=metabolite_reaction_connect{i};
        temp_index=find(temp~=0);
        temp=temp(temp_index);
        fprintf(adject_list,'%d\t',temp);
        fprintf(adject_list,'\n');
    end
    fclose(adject_list);
    
    adject_list=fopen('metabolites_reaction_num.txt','w');
    fprintf(adject_list,'%d\t%d\n',length(metabolites_unCofactor),length(reactions));
    fclose(adject_list);
end