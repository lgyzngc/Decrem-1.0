function clusterStru = cluster_similarityMatrix(similarity_matrix,clusternum)
    k=clusternum;
    t=ones(size(similarity_matrix,1),1);
    for i=1:size(similarity_matrix,1)
        temp=similarity_matrix(i,:);
%         temp(find(temp>0))=1;
        similarity_matrix(i,:)=temp;
    end
    clusterStru=cluster_wcut(similarity_matrix,t,k,1);
end