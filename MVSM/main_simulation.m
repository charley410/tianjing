clear all;
k=6;% number of views
cluster=3;
n=1000;% number of points
m=1;% run experiment for m times
P(:,:,1)=(1/n)*[16,0,0;0,18,0;0,0,17];   
P(:,:,2)=(1/n)*[16,0.4,0.6;0.4,18,0.55;0.6,0.55,17];
P(:,:,3)=(1/n)*[16,0.8,1.2;0.8,18,1.1;1.2,1.1,17];
P(:,:,4)=(1/n)*[16,1.2,1.8;1.2,18,1.65;1.8,1.65,17]; 
%classify = [333,333,334;333,333,334;333,333,334;333,333,334;333,333,334;...
 %           333,333,334;333,333,334;333,333,334;333,333,334;333,333,334;];
%classify = [200,400,400;250,300,450;300,300,400;300,400,300;400,300,300;...
 %           450,300,250;400,350,250;350,400,250;350,350,400;350,300,350];
%classify = [300,300,400;300,300,400;400,300,300;410,290,300;300,370,330;...
%           300,350,350;280,420,300;300,400,300;450,250,300;450,250,300];
%classify = [50,50,50;30,90,30;40,60,50];
%classify = [50,50,50;50,50,50;50,50,50];
%classify = [333,333,334;333,333,334;333,333,334;333,333,334;333,333,334;333,333,334];
classify = [300,300,400;300,300,400;400,300,300;300,350,350;300,400,300;450,250,300;];
        
        
for y =1:size(P,3)
    for x = 1:m
        sizeA=n;
        AA=cell(1); % record k affinity matrices
        ddxx=cell(1); % record k indicator vectors correspongding to classify
        for j = 1:size(classify,1)
            [AA{j},ddxx{j}]=rand3module(n,classify(j,1),classify(j,2),classify(j,3),P(:,:,y));% generate simulation data
            sizeA_temp = size(AA{j},1);
            sizeA = min(sizeA,sizeA_temp);
        end
        for j = 1:size(classify,1) % transform to the form of N*N*k
            A(:,:,j)=AA{j}(1:sizeA,1:sizeA);
            ddx(:,j)=ddxx{j}(1:sizeA);
        end
        
        [~,U1,~,~,~,~,~]=adaptedweight(A,cluster,1); % proposed method
        dx_total_ = kmeans(U1,cluster,'EmptyAction','drop','Replicates',100);
        num = length(dx_total_)/k;
        for i =1:k
            dx_total(:,i) = dx_total_(1+(i-1)*num:i*num);
        end  
        score_our(x,y) = accuracy(dx_total,ddx,cluster);
        [dx_total_MVSC_,U2]=MVSC(A,cluster,1); 
        %dx_total_SC_ = kmeans(U2(:,1:cluster),cluster,'EmptyAction','drop','Replicates',100);
        num1 = length(dx_total_MVSC_)/k;
        for i =1:k
            dx_total_MVSC(:,i) = dx_total_MVSC_(1+(i-1)*num1:i*num1);
        end   
        score_MVSC(x,y) = accuracy(dx_total_MVSC,ddx,cluster);
%         profile viewer
         [dx_aasc,~]=AASC(A,cluster);
         score_AASC(x,y) = accuracy(dx_aasc,ddx,cluster);
        
         
         clear AA;
         clear ddxx;
         clear sizeA;
         clear sizeA_temp;
         clear A;
         clear ddx;
         clear dx;
         clear dx_total;
         clear score_vec;
         clear dx_aasc;
         clear U1;
         clear U2;
         clear dx_total_SC;
    end
end
std_our = std(score_our);
std_MVSC = std(score_MVSC);
std_AASC = std(score_AASC);
score_our_avg = mean(score_our,1);
score_MVSC_avg = mean(score_MVSC,1);
score_AASC_avg = mean(score_AASC,1);
% save('score3.mat','score');