function [idx_new,score_new ] = IDX( cluster,A,A_knn,num,m)
%IDX 此处显示有关此函数的摘要
%   此处显示详细说明
score_ = zeros(3,m); % 
% score_snf_ = zeros(3,m);
% score_aasc_ = zeros(3,m);
% score_mvsc_ = zeros(3,m);
score1_ = zeros(3,m);
% score_mvsc1_ = zeros(3,m);
for K=1:m   
%     W1 = SNF({W_gene,W_DNA,W_miRNA},20,30);
%     idx_snf = SpectralClustering_SNF(W1,cluster);%snf method
%     [idx_aasc,~]=AASC(A_knn,cluster);% AASC method

    [~,U2,~,~,~,~,~,~]=adaptedweight(A_knn,cluster,1); % proposed method
%     [idx_mvsc_,U_mvsc] = MVSC(A_knn,cluster,1);
   idx2 = kmeans(U2,cluster,'EmptyAction','drop','Replicates',100);
    idx = zeros(num,3);
%     idx_mvsc = zeros(num,3);
    for i=1:3
        idx(:,i)=idx2((1+(i-1)*num):i*num);
%         idx_mvsc(:,i)=idx_mvsc_((1+(i-1)*num):i*num);
    end
    [idx1,idx_new,solved] = unification_(idx,U2,cluster,3);
%     [idx_mvsc1,idx_mvsc_new,solved_mvsc] = unification_(idx_mvsc,U_mvsc,cluster,3);
    %idx_mvsc_other = unification(idx_mvsc,cluster);
    %idx_other = unification(idx,cluster);

    for i=1:3
        score_(i,K)=silhouette(A(:,:,i),idx(:,i),cluster);
%         score_snf_(i,K) = silhouette(A(:,:,i),idx_snf,cluster);
%         score_aasc_(i,K) = silhouette(A(:,:,i),idx_aasc,cluster);
%         score_mvsc_(i,K) = silhouette(A(:,:,i),idx_mvsc(:,i),cluster);
        score1_(i,K) = silhouette(A(solved,solved,i),idx1,cluster);
%         score_mvsc1_(i,K) = silhouette(A(solved_mvsc,solved_mvsc,i),idx_mvsc1,cluster);
        score_new_(i,K) = silhouette(A(:,:,i),idx_new,cluster);
%         score_mvsc_new_(i,K) = silhouette(A(:,:,i),idx_mvsc_new,cluster);

    end
end  



score = mean(mean(score_,1));
% score_snf = mean(mean(score_snf_,1));
% score_aasc = mean(mean(score_aasc_,1));
% score_mvsc = mean(mean(score_mvsc_,1));
% scoce1 = mean(mean(score1_,1));
% score_mvsc1 = mean(mean(score_mvsc1_,1));
score_new = mean(mean(score_new_,1));
% score_mvsc_new = mean(mean(score_mvsc_new_,1));

end

