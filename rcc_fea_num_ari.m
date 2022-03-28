
% clc;
clear;
close all;

addpath(genpath('nmfv1_4'))

fn_txt        =   'Results_randidx.txt';
fd_txt        =   fopen( fullfile(fn_txt), 'at');

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
NS = round( 1 * numSbj ) ;%0.8

Nrep = 30;

Kcv = 2:5;  % number of subject sub-cluster
Krv = 4:11;%11:19;% number of feature sub-cluster
% clustering
for i = 1:length(Kcv)
    
    Kc = Kcv(i); 
    fprintf(fd_txt, '\n%d\t', Kc);
    
    for j = 1:length(Krv)
        
        SMtr = [];
        clustermatr1 = zeros(Nrep, numSbj);
        
        for k = 1:Nrep
            
            randid = randperm(numSbj);
            nSbjFea2 = nSbjFea(randid(1:numSbj), :);
            
            close all;
            Kr = Krv(j);
            
            lambda = 2^-9.5;
            option.orthogonal = [1,1];
            option.iter = 300000;%40000
            option.dis = 0;
            option.residual = 1e-5;
            option.tof = 1e-5;
            
            for inneriter = 1:1
                disp('Performing Collaborative Clustering ...');
                tic;
                %             rand('twister',7);
                [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea2', Kc, Kr, option);
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
                
                [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(SbjFea_new, Kc, Kr, option);
                
            end
            
            kLab = litekmeans(Y',Kc,'Replicates',30);
            
            clustermatr1(k, randid(1:NS)) = kLab;
            
        end

        rind1 = [];
        clustermatr1(clustermatr1==0)=Kc+1;
        for id1 = 1:Nrep-1
            for id2 = id1+1:Nrep
                rind1 = [ rind1, RandIndex( clustermatr1(id1, :), clustermatr1(id2, :) ) ];
            end
        end
        
        fprintf(fd_txt, '%2.3f\t', mean(rind1));
 
    end
    
end

disp('Finished.');
fclose all;