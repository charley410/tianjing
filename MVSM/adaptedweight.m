function [f,U,iter,F,R,R_0,XXXX,B] = adaptedweight(A,cluster,beta)
%function [U,alpha,iteration,fu,fval,err_fu,err_fval] = adaptedweight(A,cluster)
%Input: 
%       A:n*n*k affinity matrix
%       cluster: number of clusters
% [x,fval]=quadprog(H,g,A,b,Aeq,beq,LB,UB,X0,options);

k=size(A,3);
sizeA=size(A,1);


B_within=zeros(sizeA*k);
B_across=zeros(sizeA*k);
idx=repmat(sizeA,k,1);

for i=1:k
    D=diag(sum(A(:,:,i),1));
    B_within(((i-1)*sizeA+1):i*sizeA,((i-1)*sizeA+1):i*sizeA)=D-A(:,:,i);
end
% B_across=repmat(-eye(sizeA),k,k);
% B_across=B_across-diag(diag(B_across));
B_across1 = mat2cell(B_across,repmat(sizeA,1,k),repmat(sizeA,1,k));
for i=1:k-1
    for j=i+1:k
        B_across1{i,j} = -beta*eye(sizeA);
    end 
end
B_across2 = cell2mat(B_across1);
B_across = B_across2 + B_across2';

%  for x=1:(sizeA*k-sizeA)
%      B_across(x,x+sizeA)=-beta;
%      B_across(x+sizeA,x)=-beta;
%  end
B=B_within+B_across;
[~,V]=eig(B);
%aa = 1.5*max(max(abs(V)));
aa = min(min(abs(V)))+0.1;
B_ = B+aa*eye(size(B,1));
[f,U,iter,F,R,R_0,XXXX]=multiviewcluster_gradient1(B_,idx,cluster);



 
 


