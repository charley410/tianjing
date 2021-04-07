function [a,a1,solved]=unification_(idx,U,cluster,k)
% we only get patients with same cluster in all three view, eg:
% idx=[2,2,2]->a=2
% k: number of views
% solved : the index of points that are clustered to same cluster in all views
% a1:[2,2,2]->2,the rest use knn where k=5 to choose its cluster
num = size(idx,1);
a = zeros(num,1);
a1 = zeros(num,1);
a2 = zeros(num,1);
b = zeros(num,cluster);
for i=1:num
    for j=1:cluster
        b(i,j)=length(find(idx(i,:)==j));
    end 
    cc = find(b(i,:)==k);%the point is clustered into same cluseter in all views.
    if length(cc)==1
        a(i)=cc;
        a1(i)=cc;
        a2(i)=cc;
    end
end
unsolved = find(a1==0);
solved = find(a1>0);
cluster_idx = {};
for i=1:cluster
    cluster_idx{i} = find(a1==i);
end
for i=1:length(unsolved)
    for k=1:cluster
        DIST = zeros(length(cluster_idx{k}),1);
        for j=1:length(cluster_idx{k})
            DIST(j) = norm(U(cluster_idx{k}(j),:)-U(unsolved(i),:))+norm(U(cluster_idx{k}(j)+num,:)-U(unsolved(i)+num,:))+...
            norm(U(cluster_idx{k}(j)+2*num,:)-U(unsolved(i)+2*num,:));
        end
        tmp = sort(DIST(:),'ascend');
        dist_with_cluster(k) = mean(tmp(:));
    end
    a1(unsolved(i)) = find(dist_with_cluster==min(dist_with_cluster));
end

a(a==0)=[];



    