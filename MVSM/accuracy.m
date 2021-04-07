function acc = accuracy(idx,true_idx,cluster)
% Input: idx(n*k): the clustering result idx for all k views;
%        true_idx:
%        cluster: cluster's number;
% Output: (TP+TN)/(TP+TN+FP+FN);

k = size(idx,2);
k1 = size(true_idx,2);
n = size(idx,1);
F = [1:n];
if k==1
    for i=1:k1
        idx_do(:,i)=idx;
    end
    idx = idx_do;
end
for i=1:cluster
    xx = find(true_idx(:,1)==i);
    for j=2:k1
        yy = find(true_idx(:,j)==i);
        xx = intersect(xx,yy);
    end
    F = setdiff(F,xx);
    T{i}=xx;
end
for i=1:cluster
        aa = find(idx(:,1)==i);
        aa1 = aa;
        aa3 =[];
        aa4 = aa;
        for j=2:k1
            bb = find(idx(:,j)==i);
            bb1 = find(idx(:,j)~=i);
            aa3 = union(aa3,bb1);
            aa = intersect(aa,bb);
            clear bb;
            aa1 = intersect(aa1,bb1);
        end
        A1{i}=intersect(aa3,aa4);
        clear aa3;
        clear aa4;
        A2{i}=aa;
        clear aa;
%        A1{i}=aa1;
end
for i=1:cluster
    num1(i) = length(intersect(A1{i},F));
    for j=1:cluster
        num(i,j) = length(intersect(A2{i},T{j}));
    end
end
TP = sum(max(num,[],2));    
TN = sum(num1);
acc = (TP+TN)/n;


     






