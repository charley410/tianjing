%% load data
COAD_gene=SARGeneExpression
COAD_DNA= SARMethyExpression
COAD_miRNA=SARMirnaExpression
%% construct the graphs for each data type
options = [];
options.NeighborMode = 'KNN';
options.k = 5; 
options.WeightMode = 'HeatKernel';
options.t = 100;

W=constructW(COAD_gene',  options)

W_GE   = full(constructW(COAD_gene',  options));
W_MIR  = full(constructW(COAD_miRNA', options));
W_METH = full(constructW(COAD_DNA',options));

%% Integrative clustering 
n = size(W_GE, 1);

W_COMBINE = zeros(n, n, 3);
W_COMBINE(:,:,1) = W_GE(:,:);
W_COMBINE(:,:,2) = W_METH(:,:);
W_COMBINE(:,:,3) = W_MIR(:,:);

%% compute distance martix
dist_gene = dist2(COAD_gene',COAD_gene');
dist_DNA = dist2(COAD_DNA',COAD_DNA');
dist_miRNA = dist2(COAD_miRNA',COAD_miRNA');

%% construct similarity matrix
W_gene = affinityMatrix(dist_gene,20 ,0.41);
W_DNA = affinityMatrix(dist_DNA,20,0.41);
W_miRNA = affinityMatrix(dist_miRNA,20,0.41);

W_gene_knn = knnAffinity(dist_gene,5,'T');
W_DNA_knn = knnAffinity(dist_DNA,5,'T');
W_miRNA_knn = knnAffinity(dist_miRNA,5,'T');

A(:,:,1)=W_gene;
A(:,:,2)=W_DNA;
A(:,:,3)=W_miRNA;


A_knn(:,:,1)=W_COMBINE(:,:,1);
A_knn(:,:,2)=W_COMBINE(:,:,2);
A_knn(:,:,3)=W_COMBINE(:,:,3);
 

%% compute indicators
m=3;% we will run for m times and count the average as the resul
num = size(A,1);
cluster=3:5
Idex_new=zeros(num,size(cluster,2))
Score_new=zeros(1,size(cluster,2))
for i=1:size(cluster,2)
    [idex_new,score_new]=IDX(cluster(i),A,A_knn,num,m)
    Idex_new(:,i)=idex_new
    Score_new(i)=score_new
end
score_new=max(Score_new)
cluster=cluster(find(Score_new==max(Score_new)))
idx=Idex_new(:,find(Score_new==max(Score_new)))
% score_ = zeros(3,m); % 
% % score_snf_ = zeros(3,m);
% % score_aasc_ = zeros(3,m);
% % score_mvsc_ = zeros(3,m);
% score1_ = zeros(3,m);
% % score_mvsc1_ = zeros(3,m);
% for K=1:m   
% %     W1 = SNF({W_gene,W_DNA,W_miRNA},20,30);
% %     idx_snf = SpectralClustering_SNF(W1,cluster);%snf method
% %     [idx_aasc,~]=AASC(A_knn,cluster);% AASC method
% 
%     [~,U2,~,~,~,~,~,~]=adaptedweight(A_knn,cluster,1); % proposed method
% %     [idx_mvsc_,U_mvsc] = MVSC(A_knn,cluster,1);
%    idx2 = kmeans(U2,cluster,'EmptyAction','drop','Replicates',100);
%     idx = zeros(num,3);
% %     idx_mvsc = zeros(num,3);
%     for i=1:3
%         idx(:,i)=idx2((1+(i-1)*num):i*num);
% %         idx_mvsc(:,i)=idx_mvsc_((1+(i-1)*num):i*num);
%     end
%     [idx1,idx_new,solved] = unification_(idx,U2,cluster,3);
% %     [idx_mvsc1,idx_mvsc_new,solved_mvsc] = unification_(idx_mvsc,U_mvsc,cluster,3);
%     %idx_mvsc_other = unification(idx_mvsc,cluster);
%     %idx_other = unification(idx,cluster);
% 
%     for i=1:3
%         score_(i,K)=silhouette(A(:,:,i),idx(:,i),cluster);
% %         score_snf_(i,K) = silhouette(A(:,:,i),idx_snf,cluster);
% %         score_aasc_(i,K) = silhouette(A(:,:,i),idx_aasc,cluster);
% %         score_mvsc_(i,K) = silhouette(A(:,:,i),idx_mvsc(:,i),cluster);
%         score1_(i,K) = silhouette(A(solved,solved,i),idx1,cluster);
% %         score_mvsc1_(i,K) = silhouette(A(solved_mvsc,solved_mvsc,i),idx_mvsc1,cluster);
%         score_new_(i,K) = silhouette(A(:,:,i),idx_new,cluster);
% %         score_mvsc_new_(i,K) = silhouette(A(:,:,i),idx_mvsc_new,cluster);
% 
%     end
% end  
% 
% 
% 
% score = mean(mean(score_,1));
% % score_snf = mean(mean(score_snf_,1));
% % score_aasc = mean(mean(score_aasc_,1));
% % score_mvsc = mean(mean(score_mvsc_,1));
% % scoce1 = mean(mean(score1_,1));
% % score_mvsc1 = mean(mean(score_mvsc1_,1));
% score_new = mean(mean(score_new_,1));
% % score_mvsc_new = mean(mean(score_mvsc_new_,1));






        