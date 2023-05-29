function CycMatrix = findSimpleCycWithinKStep(adjacency_graph_structure,minLen,maxLen)
    reaction_num = length(adjacency_graph_structure.reaction);
    metabolite_num = length(adjacency_graph_structure.metabolite);
    global CycMatrix
    CycMatrix = zeros(reaction_num);
    adjacency_graph = adjacency_graph_structure.metabolite_reaction_connect;
    
    removed_set = [];
    for i=1:metabolite_num
        global seeked_nodes
        seeked_nodes = [];
        isCyc = SimpleCycWithinKStep(i,i,adjacency_graph,minLen,maxLen,reaction_num,metabolite_num);
        removed_set = [removed_set,i];
        i
    end
end