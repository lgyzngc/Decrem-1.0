function isCyc = SimpleCycWithinKStep(source,target,adjacency_graph,minLen,maxLen,reaction_num,metabolite_num)
    global seeked_nodes;
    global CycMatrix;
    isCyc = 0;
    
    if target <= (metabolite_num+reaction_num)
        seeked_nodes = [seeked_nodes,target];
    else
        seeked_nodes = [seeked_nodes,target-reaction_num];
    end
    if length(seeked_nodes) > 2*maxLen
        seeked_nodes=seeked_nodes(1:end-1);
        return
    end
    
    
    edge_nodes = adjacency_graph{target};
    for i=1:length(edge_nodes)
        if edge_nodes(i) == source && length(seeked_nodes)>= 2*minLen
            isCyc = 1;
            CycReac_temp = [];
            for j=1:length(seeked_nodes)
                if seeked_nodes(j) > metabolite_num+reaction_num
                    CycReac_temp = [CycReac_temp,seeked_nodes(j)-(metabolite_num+reaction_num)];
                elseif seeked_nodes(j) > metabolite_num
                    CycReac_temp = [CycReac_temp,seeked_nodes(j)-metabolite_num];
                end
            end
            for j=1:length(CycReac_temp)
                for k=j+1:length(CycReac_temp)
                    CycMatrix(CycReac_temp(j),CycReac_temp(k)) = CycMatrix(CycReac_temp(j),CycReac_temp(k))+1;
                end
            end
        elseif ~(ismember(edge_nodes(i),seeked_nodes) || ismember(edge_nodes(i)-reaction_num,seeked_nodes))
            isCyc = SimpleCycWithinKStep(source,edge_nodes(i),adjacency_graph,minLen,maxLen,reaction_num,metabolite_num);
        end
    end
    seeked_nodes=seeked_nodes(1:end-1);
end