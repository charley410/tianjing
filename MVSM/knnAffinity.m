function [W]=knnAffinity(dist,K,cc)
% dist: n*n dist_matrix
% K: for each i, W(i,j)=1,if j belongs to i's K nearest neighbor
% cc: if cc='S', k-nearest neighbor graph; if cc='T', W(i,j)=0.5 if only
% one point belongs to the other point's K nearest neighbor

num = size(dist,1);
W = zeros(num,num);
for i=1:num
    dist(i,i)=inf;
    t = sort(dist(i,:));
    idx = find(dist(i,:)<t(K+1),K);
    W(i,idx)=1;
%    id{i}=idx;
end
% W1 = W+W';
if cc=='S'
    W1 = W+W';
    W(W1==1)=1;
elseif cc=='T'
    W1 = (W+W')/2;
    W = W1;
end
% elseif cc=='M'
%     W(W1==1)=0;
%     index = find(sum(W)==0);
%     for i=1:length(index)
%         j = index(i);
%         k = find(dist(j,:)==min(dist(j,:)));
%         W(j,k)=1;
%         W(k,j)=1;
%     end
% end
    
    