function [idx] = commonidx(idxmat,cluster)

index = zeros(size(idxmat,1),1);
rowvec = ones(1,size(idxmat,2));
for i = 1:size(idxmat,1)
    for j =1:cluster
        if idxmat(i,:)==j*rowvec
            index(i)=1;
        end
    end
end
idx = find(index==1);
            