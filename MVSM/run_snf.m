clear all;
k=3;% number of vies
cluster=3;
n=150;
K = 5; % k for knn
T = 20; % snf iterates for T times
m = 30;
P(:,:,1)=(1/n)*[16,0,0;0,18,0;0,0,17];   
P(:,:,2)=(1/n)*[16,0.4,0.6;0.4,18,0.55;0.6,0.55,17];
P(:,:,3)=(1/n)*[16,0.8,1.2;0.8,18,1.1;1.2,1.1,17];
P(:,:,4)=(1/n)*[16,1.2,1.8;1.2,18,1.65;1.8,1.65,17]; 
classify = [50,50,50;30,90,30;40,60,50];
%classify = [50,50,50;50,50,50;50,50,50];
%classify = [333,333,334;333,333,334;333,333,334;333,333,334;333,333,334;333,333,334];
%classify = [300,300,400;300,300,400;400,300,300;300,350,350;300,400,300;450,250,300;];
    
for y =1:size(P,3)
    for x = 1:m
        sizeA=n;
        AA=cell(1); 
        ddxx=cell(1); 
        A = cell(1);
        for j = 1:size(classify,1)
            [AA{j},ddxx{j}]=rand3module(n,classify(j,1),classify(j,2),classify(j,3),P);
            sizeA_temp = size(AA{j},1);
            sizeA = min(sizeA,sizeA_temp);
        end
        for j = 1:size(classify,1)
            A{j}=AA{j}(1:sizeA,1:sizeA);
            ddx(:,j)=ddxx{j}(1:sizeA);
        end
        W = SNF(A,K,T);
        idx = SpectralClustering_SNF(W,cluster);
        score(x,y)= accuracy(idx,ddx,cluster);
    end
end
std = std(score);
score1 = mean(score,1);
