function [dx,U] = SC(A,cluster,beta)
% Input
%       A : N x N x k affinity matrix 
% cluster : desired number of clusters
% Output
%      dx : clustering result    
      %%% compute Laplacian matrix %%%
k=size(A,3);
sizeA=size(A,1);
B_within=zeros(sizeA*k);
B_across=zeros(sizeA*k);

for i=1:k
    D=diag(sum(A(:,:,i),1));
    B_within(((i-1)*sizeA+1):i*sizeA,((i-1)*sizeA+1):i*sizeA)=D-A(:,:,i);
end

%  B_across=repmat(-eye(sizeA),k,k);
%  B_across=B_across-diag(diag(B_across));
%  B_across=beta*B_across;  

B_across1 = mat2cell(B_across,repmat(sizeA,1,k),repmat(sizeA,1,k));
for i=1:k-1
    for j=i+1:k
        B_across1{i,j} = -beta*eye(sizeA);
    end 
end
B_across2 = cell2mat(B_across1);
B_across = B_across2 + B_across2';

% for x=1:(sizeA*k-sizeA)
%       B_across(x,x+sizeA)=-beta;
%       B_across(x+sizeA,x)=-beta;
% end
B=B_within+B_across;%%% eigen decomposition %%%

% [U, eigenvalue] = eigs(B, cluster1, eps);
% [a,b] = sort(diag(eigenvalue),'ascend');
% eigenvalue = eigenvalue(:,b);
% U = U(:,b);
% eigengap = abs(diff(diag(eigenvalue)));
% U = U(:,1:cluster1);
    OPTS.disp = 0;
    [U, ~] = eigs(B , cluster, 'SA', OPTS);%generalized eigenproblem
   % [U, ~] = eigs((B+B'fn)/2, cluster, eps);
    dx = kmeans(U,cluster,'EmptyAction','drop','Replicates',50);
   clear D_;
    clear D;
    clear L;
  
