work_path = 'E:\fluxmode\flux analysis\bacillus_yeast\'
cofactor_name='cofactor.txt';
input_nutrient_name='nutrient_bacillus.txt';
general_IO_name='general_IO_bacillus.txt';
secretion_name='NULL';
extraflux_name='NULL';
intraflux_name='intracellularflux_bacillus.txt';
knockout_name='NULL';
CBM_name='iYO844.mat';

load(strcat(work_path,CBM_name),'-mat');
CBM_Model = iYO844;
reaction_hash=java.util.Hashtable;
for i=1:length(CBM_Model.rxns)
    reaction_hash.put(CBM_Model.rxns{i},i);
end

cofactors = {};
cofactor_file=fopen(strcat(work_path,cofactor_name));
index = 1;
while(~feof(cofactor_file))
    line = fgetl(cofactor_file)
    str = regexp(line,'\t','split');
    cofactors{index} = str{1};
    index = index + 1;
end
fclose(cofactor_file);

input_nutrient = {};
index = 1;
secretion_file=fopen(strcat(work_path,input_nutrient_name));
while(~feof(secretion_file))
    line = fgetl(secretion_file)
    str = regexp(line,'\t','split');
    input_nutrient{index} = str{1};
    index = index + 1;
end
fclose(secretion_file);

general_IO = {};
index = 1;
general_IO_file=fopen(strcat(work_path,general_IO_name));
while(~feof(general_IO_file))
    line = fgetl(general_IO_file)
    str = regexp(line,'\t','split');
    general_IO{index} = str{1};
    index = index + 1;
end
fclose(general_IO_file);

secretion = {};
index = 1;
if ~strcmp(secretion_name,'NULL')
    secretion_file=fopen(strcat(work_path,secretion_name));
    while(~feof(secretion_file))
        line = fgetl(secretion_file)
        str = regexp(line,'\t','split');
        input_nutrient{index} = str{1};
        index = index + 1;
    end
    fclose(secretion_file);
end

extraflux = [];
index = 1;
if ~strcmp(extraflux_name,'NULL')
    extraflux_file=fopen(strcat(work_path,extraflux_name));
    while(~feof(extraflux_file))
        line = fgetl(extraflux_file)
        str = regexp(line,'\t','split');
        if strcmp(str{1},'reactionname') 
             continue;
        end
        sub_str=str{1};
        %intraflux{index} = str{1};
        if length(sub_str)<2
             continue;
        end
        if reaction_hash.containsKey(sub_str)
             %extraflux(index,1)=reaction_hash.get(sub_str);
             temp=[reaction_hash.get(sub_str)];
             for i=2:length(str)
                 if ~isempty(str{i})
                     temp=[temp,str2num(str{i})];
                 else
                     temp=[temp,0];
                 end
             end
             %extraflux(index,2:(length(temp)+1))=temp;
             extraflux=[extraflux;temp];
             %index=index+1;
         end
    end
    fclose(extraflux_file);
end


intraflux = [];
index = 1;
if ~strcmp(intraflux_name,'NULL')
    intraflux_file=fopen(strcat(work_path,intraflux_name));
    while(~feof(intraflux_file))
        line = fgetl(intraflux_file)
        str = regexp(line,'\t','split');
        if strcmp(str{1},'reactionname') 
             continue;
        end
        sub_str=str{1};
        %intraflux{index} = str{1};
        if length(sub_str)<2
             continue;
        end
        if reaction_hash.containsKey(sub_str)
             %temp=[];
             %intraflux(index,1)=reaction_hash.get(sub_str);
             temp=[reaction_hash.get(sub_str)];
             for i=2:length(str)
                 if ~isempty(str{i})
                     temp=[temp,str2num(str{i})];
                 else
                     temp=[temp,0];
                 end
             end
             %intraflux(index,2:(length(temp)+1))=temp;
             intraflux=[intraflux;temp];
             index=index+1;
         end
        %index = index + 1;
    end
    fclose(intraflux_file);
end

knockout_gene = {};
index = 1;
if ~strcmp(knockout_name,'NULL')
    knockout_file=fopen(strcat(work_path,knockout_name));
    while(~feof(knockout_file))
        line = fgetl(knockout_file)
        str = regexp(line,'\t','split');
        knockout_gene{index} = str{1};
%         for i=1:length(str)
%             knockout_gene{index,1:length(str)} = str;
%         end
%         index = index + 1;
    end
    fclose(knockout_file);
end

if ~isempty(intraflux) && ~isempty(extraflux) && ~isempty(knockout_gene)
    if (size(intraflux,2)-1) == length(knockout_gene) && (size(extraflux,2)-1) == length(knockout_gene)
        for i=2:size(intraflux,2)
            [Decrem_solution,Decrem_model,CBM_model] = Decrem(CBM_Model,cofactors,input_nutrient,secretion,general_IO,1,...
                extraflux(:,[1,i]),intraflux(:,[1,i]),knockout_gene{i});
        end
    end
elseif ~isempty(intraflux) && ~isempty(extraflux)
    if size(intraflux,2) == size(extraflux,2)
        for i=2:size(intraflux,2)
            [Decrem_solution,Decrem_model,CBM_model] = Decrem(CBM_Model,cofactors,input_nutrient,secretion,general_IO,1,...
                extraflux(:,[1,i]),intraflux(:,[1,i]),{});
        end
    end
elseif  ~isempty(intraflux)
    for i=3:size(intraflux,2)
        [Decrem_solution,Decrem_model,CBM_model] = Decrem(CBM_Model,cofactors,input_nutrient,secretion,general_IO,1,...
                [],intraflux(:,[1,i]),{});
    end
else
    [Decrem_solution,Decrem_model,CBM_model] = Decrem(CBM_Model,cofactors,input_nutrient,secretion,general_IO,1,[],[],{});
end