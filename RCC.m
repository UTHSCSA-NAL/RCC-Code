
% clc;
clear;
close all;

addpath(genpath('nmfv1_4'))

a = dataset('XLSFile','merged_MCI.csv');
a = dataset2table(a);
start_col = 6; end_col = 63;%157;%

sbjFea = table2array(a(1:end-1, start_col:end_col)); %, a(:, 266:523)double(brain);

sumSbjFea = repmat(sum(sbjFea, 2), 1, size(sbjFea,2));
nSbjFea = sbjFea ./ sumSbjFea;
nSbjFea = nSbjFea-min(nSbjFea(:));
nSbjFea = nSbjFea./max(nSbjFea(:));
nSbjFea(:,std(nSbjFea)==0) = [];

expec = 0;
[numSbj, numFea] = size(nSbjFea);
% clustering
Kcv = 3;  % number of subject sub-cluster
Krv = 9;% number of feature sub-cluster
for i = 1:length(Kcv)
    
    cls = [];
    
    for j = 1:length(Krv)
        
        close all;
        Kc = Kcv(i); Kr = Krv(j);
        
        lambda = 2^-12;
        option.orthogonal = [1,1];
        option.iter = 300000;%40000
        option.dis = 0;
        option.residual = 1e-5;
        option.tof = 1e-5;
        for inneriter = 1:1
            disp('Performing Collaborative Clustering ...');
            tic;
            rand('twister',73);
            [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea', Kc, Kr, option);
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
            nSbjFea2 = SbjFea_new';
            %%
            %             rand('twister',17);
            [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea2', Kc, Kr, option);
            
        end
        %         rand('twister',37);
        kLab = litekmeans(Y',Kc,'Replicates',30);
        %kLab = litekmeans(nSbjFea,Kc,'Replicates',20);     % clustering based on
        %original features
        [weight,aLab] = max(A,[],2);
        
        [skLab,skInd] = sort(kLab);
        [saLab,saInd] = sort(aLab);
        
        %% plot clustering results
        figure; colormap('parula');hold on; box on;
        imagesc(nSbjFea2(skInd,saInd)); axis image; colorbar; caxis([min(nSbjFea2(:)),max(nSbjFea2(:))]); xlabel('Feature ID'); ylabel('Sbj ID');
        cpt = find(diff(skLab)~=0) + 1;
        for cpi=1:length(cpt)
            plot([0,size(nSbjFea2,2)+1],[cpt(cpi),cpt(cpi)],'r--','LineWidth',2);
        end
        
        cpa = find(diff(saLab)~=0);
        for cpi=1:length(cpa)
            plot([cpa(cpi),cpa(cpi)],[1,size(nSbjFea2,1)],'r--','LineWidth',2);
        end
        saveas(gcf, ['cluster_result_merged_Kc',num2str(Kc),'_Kr',num2str(Kr) '.png'])
        
        cl = kLab';
        
        a(end, start_col:end_col) = array2table(aLab');
        a(1:end-1, end) = array2table(kLab);
        af = [a(:, 1), a(:, start_col:end_col), a(:,end)];
        writetable(af,['rcc_result_','Kc',num2str(Kc),'_Kr',num2str(Kr),'.csv']);
        
    end
    
    %     save('fealab.mat','weight', 'aLab');
    
    outName = ['U_kLab',num2str(Kc) '_' num2str(Kr), '.mat'];
    save(outName,'cl');
    
end

disp('Finished.');
