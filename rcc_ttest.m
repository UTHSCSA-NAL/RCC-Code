a1 = dataset('XLSFile','merged.csv');
a1 = dataset2table(a1);
start_col = 2; end_col = 59;%
sbjFea1 = table2array(a1(1:end-1, start_col: end_col)); %, a(:, 266:523)double(brain);

% sumSbjFea = repmat(sum(sbjFea1, 2), 1, size(sbjFea1,2));
% nSbjFea1 = sbjFea1 ./ sumSbjFea;
% nSbjFea1 = nSbjFea1./max(nSbjFea1(:));
% deno1 = table2array(a1(1:end-1, 9));
% denom1 = repmat(deno1, 1, size(sbjFea1, 2));
% denom1(denom1==0) = eps;
% sbjFea1 = sbjFea1./ denom1;

a2 = dataset('XLSFile','rcc_result_Kc3_Kr4.csv');
a2 = dataset2table(a2);
sbjFea2 = table2array(a2(1:end-1, start_col: end_col)); %, a(:, 266:523)double(brain);

% sumSbjFea = repmat(sum(sbjFea2, 2), 1, size(sbjFea2,2));
% nSbjFea2 = sbjFea2 ./ sumSbjFea;
% sbjFea2 = nSbjFea2 ./ max(nSbjFea1(:));
% deno2 = table2array(a2(1:end-1, 9));
% denom2 = repmat(deno2, 1, size(sbjFea2, 2));
% denom2(denom2==0) = eps;
% sbjFea2 = sbjFea2./ denom2;

nc = 3;
klab4_g2 = load(['U_kLab' num2str(nc) '_4.mat']);


% % G1&G2
tmp = zeros(nc, end_col-start_col+1);
tmp_coe = tmp;
tmp1 = cell(nc, 1);
tmp_coe1 = tmp1;
k = 0;
for anchor = 1
    for test = 1:nc
        k = k+1;
        tmp1(k) = cellstr(['PVal_C' num2str(anchor) 'P' num2str(test)]);
        tmp_coe1(k) = cellstr(['TSta_C' num2str(anchor) 'P' num2str(test)]);
        for i = 1:1-start_col+end_col
            x = (sbjFea1(:, i));
            y = (sbjFea2(klab4_g2.cl==test, i));
            [h1, p, ci, coe] = ttest2(x,y);
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p);
            tmp(k, i) = adj_p;
            tmp_coe(k, i) = coe.tstat;
        end
    end
end

c42=nc;
tmp1 = cell2table(tmp1);
a1(1:c42, start_col:end_col) = array2table(tmp);
af_p1 = [tmp1, a1(1:c42, start_col:end_col)];

tmp_coe1 = cell2table(tmp_coe1);
a1(1:c42, start_col:end_col) = array2table(tmp_coe);
af_c1 = [tmp_coe1, a1(1:c42, start_col:end_col)];

%Within Group2
c42 = nchoosek(nc,2);
tmp = zeros(c42, end_col-start_col+1);
tmp_coe = tmp;
tmp1 = cell(c42, 1);
tmp_coe1 = tmp1;
k = 0;
for anchor = 1:nc
    for test = anchor+1:nc
        k = k+1;
        tmp1(k) = cellstr(['PVal_P' num2str(anchor) 'P' num2str(test)]);
        tmp_coe1(k) = cellstr(['TSta_P' num2str(anchor) 'P' num2str(test)]);
        for i = start_col:end_col
            x = table2array(a2(klab4_g2.cl==anchor, i));
            y = table2array(a2(klab4_g2.cl==test, i));
            [h1, p, ci, coe] = ttest2(x,y);
            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p);
            tmp(k, i-start_col+1) = adj_p;
            tmp_coe(k, i-start_col+1) = coe.tstat;
        end
    end
end

tmp1 = cell2table(tmp1);
a1(1:c42, start_col:end_col) = array2table(tmp);
af_p2 = [tmp1, a1(1:c42, start_col:end_col)];
% writetable([af_p1; af_p2],['PVal_' num2str(nc) '.csv']);

tmp_coe1 = cell2table(tmp_coe1);
a1(1:c42, start_col:end_col) = array2table(tmp_coe);
af_c2 = [tmp_coe1, a1(1:c42, start_col:end_col)];
% writetable([af_c1; af_c2],['TStats_' num2str(nc) '.csv']);
af_p1.Properties.VariableNames(1) = {'name'};
af_p2.Properties.VariableNames(1) = {'name'};
af_c1.Properties.VariableNames(1) = {'name'};
af_c2.Properties.VariableNames(1) = {'name'};
pt = [af_p1; af_p2; af_c1; af_c2];

pt = rows2vars(pt);

writetable(pt,['C:\Users\LiuHan\Downloads\ROI_Mapping_Combined\PVal_TStat_binomalized' num2str(nc) '.csv']);

