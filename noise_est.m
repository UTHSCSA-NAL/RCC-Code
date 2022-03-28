
% clc;
% clear;

addpath(genpath('nmfv1_4'))


a = dataset('XLSFile','merged_MCI.csv');
a = dataset2table(a);
start_col = 6; end_col = 63;%157;%

sbjFea = table2array(a(1:end-1, start_col:end_col)); %, a(:, 266:523)double(brain);

sumSbjFea = repmat(sum(sbjFea, 2), 1, size(sbjFea,2));
nSbjFea = sbjFea ./ sumSbjFea;
nSbjFea = nSbjFea-min(nSbjFea(:));
nSbjFea = nSbjFea./max(nSbjFea(:));

rind_vec1 = []; rind_vec2 = [];
[numSbj, numFea] = size(nSbjFea);
% [h, w] = size(nSbjFea);
ornSbjFea = nSbjFea;
NS = round( 1 * numSbj );
expec = 0;
lambda_rng = [-15:0.5:-5];
Kcv = 3;  %4 :5 number of subject sub-cluster
Krv = 9; %7:21 number of feature sub-cluster
for i = 1:length(Kcv)
    for j = 1:length(Krv)
        for rep = 2:4
            for loglambda = -15:0.5:-5
                
                clustermatr1 = zeros(NS, numFea);
                clustermatr2 = zeros(NS, numSbj);
                lambda = 2^loglambda;
                
                for kk = 1:NS
                    
                    randid = randperm(numSbj);
                    nSbjFea = ornSbjFea(randid, :);
                    randid2 = randperm(numFea);
                    nSbjFea = nSbjFea(:, randid2);
                    
                    close all;
                    Kc = Kcv(i); Kr = Krv(j);
                    
                    option.orthogonal = [1,1];
                    option.iter = 10000;
                    option.dis = 0;
                    option.residual = 1e-5;
                    option.tof = 1e-5;
                    
                    disp('Performing Collaborative Clustering ...');
                    tic;
                    [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea(1:NS, :)', Kc, Kr, option);
                    toc;
                    
                    %%  Normalize Dictionaries
                    sY = sqrt(sum(Y.^2,2));
                    sA = sqrt(sum(A.^2));
                    Y = Y ./ repmat(sY,1,size(Y,2));
                    A = A ./ repmat(sA,size(A,1),1);
                    S = repmat(sA',1, size(S,2)) .* S .* repmat(sY',size(S,1),1);
                    
                    
                    
                    %%  Shrink S
                    scale = numSbj * numFea/(Kc*Kr);
                    Var_S = max(S.^2- scale*lambda, eps);
                    shrkS = lambda*sqrt(2 * numSbj * numFea ./ Var_S);
                    S_bar = S ;
                    Sr =  max(abs(S_bar)-shrkS, 0).* sign(S_bar);%
                    
                    SbjFea_new = A*Sr*Y;
                    %             minSbjFea_new = repmat(min(SbjFea_new),size(SbjFea_new,1),1);
                    %             maxSbjFea_new = repmat(max(SbjFea_new),size(SbjFea_new,1),1);
                    %             nSbjFea_new = (SbjFea_new-minSbjFea_new) ./ max(eps,maxSbjFea_new-minSbjFea_new);
                    nSbjFea = SbjFea_new';
                    %%
                    [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea', Kc, Kr, option);
                    
                    %             %%  Normalize Dictionaries
                    %             sY = sqrt(sum(Y.^2,2));
                    %             sA = sqrt(sum(A.^2));
                    %             Y = Y ./ repmat(sY,1,size(Y,2));
                    %             A = A ./ repmat(sA,size(A,1),1);
                    %             S = repmat(sA',1, size(S,2)) .* S .* repmat(sY',size(S,1),1);
                    %             %%
                    
                    kLab = litekmeans(Y',Kc,'Replicates',30);
                    [~,aLab] = max(A,[],2);
                    clustermatr1(kk, :) = aLab;
                    clustermatr2(kk, randid(1:NS)) = kLab;
                    %                 if 67+kk<=numSbj
                    %                     clustermatr2(kk+1, 1+kk:67+kk) = aLab;
                    %                 else
                    %                     clustermatr2(kk+1, 1+kk:numSbj) = aLab(1:numSbj-kk);
                    %                     clustermatr2(kk+1, 1:67-numSbj+kk) = aLab(numSbj-kk+1:67);
                    %                 end
                    
                    
                end
                
                
                rind1 = [];
                rind2 = [];
                clustermatr2(clustermatr2==0)=Krv+1;
                for id1 = 1:NS-1
                    for id2 = id1+1:NS
                        rind1 = [ rind1, RandIndex( clustermatr1(id1, :), clustermatr1(id2, :) ) ];
                        rind2 = [ rind2, RandIndex( clustermatr2(id1, :), clustermatr2(id2, :) ) ];
                    end
                end
                rind_vec1 = [rind_vec1 mean(rind1)];
                rind_vec2 = [rind_vec2 mean(rind2)];
            end
            
            save(['rind_vec1_' num2str(rep) '.mat'],'rind_vec1');
            save(['rind_vec2_' num2str(rep) '.mat'],'rind_vec2');
        end
        
        %         %% plot clustering results
        %         figure; colormap('parula');hold on; box on;
        %         imagesc(nSbjFea(skInd,saInd)); axis image; colorbar; caxis([min(nSbjFea(:)),max(nSbjFea(:))]); xlabel('Feature ID'); ylabel('Sbj ID');
        %         cpt = find(diff(skLab)~=0) + 1;
        %         for cpi=1:length(cpt)
        %             plot([0,size(nSbjFea,2)+1],[cpt(cpi),cpt(cpi)],'r--','LineWidth',2);
        %         end
        %
        %         cpa = find(diff(saLab)~=0);
        %         for cpi=1:length(cpa)
        %             plot([cpa(cpi),cpa(cpi)],[1,size(nSbjFea,1)],'r--','LineWidth',2);
        %         end
        %         saveas(gcf, ['cluster_result_merged_Kc',num2str(Kc),'_Kr',num2str(Kr) '.png']);
        %         %%
        %         a(end, start_col:2:end_col) = array2table(aLab');
        %         a(1:end-1, end) = array2table(kLab);
        %         af = [a(:, 1:4), a(:, start_col:2:end_col), a(:,end)];
        %         writetable(af,['cc_result_','Kc',num2str(Kc),'_Kr',num2str(Kr),'.csv']);
        %
        %         outName = ['cc_result_merged','Kc',num2str(Kc),'_Kr',num2str(Kr),'.mat'];
        %         save(outName,'nSbjFea','ndSbjFea','aLab','kLab');
        
    end
end

rindvec = rind_vec2(1:21);
rindvec = rindvec + rind_vec2(22:42);
rindvec = rindvec + rind_vec2(43:63);
rindvec = rindvec./3;
[val, ind] = max(rindvec);
lambda_rng(ind)

disp('Finished.');
