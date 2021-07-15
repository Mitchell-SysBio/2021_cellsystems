%% Temporal dynamics of synthetic consortia
% this script takes the correlation of assemblies' composition [grown in 
% nutrient poor media (M9) or nutrient rich media (LB)] at day
% 14 and plots a heatmap of the correlation coefficients and a bar plot of 
% the sorted correlation coefficients.
% Plots sup figure 9D-E
% 2021/04/27
%% Load data
load corrM9.mat
load corrLB.mat
%% Barplot of sorted correlation coefficients
pearson = [pearsonLB(:); pearsonM9(:)];
[y, pearsonInx] = sort(pearson);
colorBi = pearsonInx;
colorBi(pearsonInx>36) = 2 ;
colorBi(pearsonInx<=36) = 1 ;
figure;
b = bar(y);
set(gca,'ylim',[0.4 1.1])
grid on; box on;
ylabel('pearson correlation coefficient')
xlabel({'community pairs' '(sorted by pearson correlation coefficient)'})
set(gcf,'position',[440,526,744,272])
%% Heatmap of correlation coefficient separated by media 
figure;
imagesc(pearsonLB,[0.5 1]); colormap(gray);
xtickangle(90);
axis square;
title('pearson correlation coefficient - LB');
colorbar
figure;
imagesc(pearsonM9,[0.5 1]); colormap(gray);
xtickangle(90);
axis square;
title('pearson correlation coefficient - M9');
colorbar