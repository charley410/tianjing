function [A,idx]=rand3module(n,m1,m2,m3,Pmatrix) 

group1=1:m1;
group2=(m1+1):(m1+m2);
group3=(m1+m2+1):(m1+m2+m3);
gs1=m1;
gs2=m2;
gs3=m3;

A=zeros(n,n);
A(group1,group1)= rand(gs1)<Pmatrix(1,1);
A(group2,group2)= rand(gs2)<Pmatrix(2,2);
A(group3,group3)= rand(gs3)<Pmatrix(3,3);
A=triu(A);

A(group1,group2)= rand(gs1,gs2)<Pmatrix(1,2);
A(group1,group3)= rand(gs1,gs3)<Pmatrix(1,3);
A(group2, group3)= rand(gs2,gs3)<Pmatrix(2,3);

%A(group2,group1)= rand(m2,m1)<Pmatrix(2,1);
%A(group3,group1)= rand(m3,m1)<Pmatrix(3,1);
%A(group3, group2)= rand(m3,m2)<Pmatrix(3,2);

%A = triu(A);
A(A<1)=0;
A = A+A';
%A(A>=0.5)=1;
%A(A<0.5)=0;
A=A-diag(diag(A));
%AB = A;
%A = A + A';
idx=[ones(1,gs1),2*ones(1,gs2),3*ones(1,gs3)];
%idx=idx(sum(A)>0);
%A=A(sum(A)>0,sum(A)>0);
idx=idx';
