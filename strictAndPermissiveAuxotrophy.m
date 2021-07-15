%% Strict and Permissive Auxotrophy
% analyzes data of deletion strains grown in nutrient poor media or
% nutrient poor media conditioned by a wild-type strain
% Plots figure 3B-C
% 2020/10/11
%%
load cmAUC.mat
map = [0.4940 0.1840 0.5560
    0.5 0.5 0.5
    0.2968 1.0000 0.5097
    0.6350 0.0780 0.1840
    0.2143 0.4156 1.0000];
%% scatter plot
% plots figure 3B
x=mmAUC;
y=cmAUC;
textCord = y>=0.5;
figure;
hold on
for i =1:length(x)
    plot(x(i),y(i),'o','MarkerFaceColor',map(hitAnnotations(i,1),:),...
    'MarkerEdgeColor',map(hitAnnotations(i,1),:),'MarkerSize',9)
end
plot([min([x; y]) 2.5],[min([x; y]) 2.5], '-k')
plot(0,y(67),'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',12)
text([x(textCord)+0.03],y(textCord),hitNames(textCord),'FontSize',12,'FontName','avenir')
title('Area under the curve at 11 h','FontSize',18,'FontName','avenir')
xlabel('AUC - minimal media','FontSize',16,'FontName','avenir')
ylabel('AUC - conditioned media','FontSize',16,'FontName','avenir')
axis square
grid on
%% collaboration score filtered first and normalized
ideal = [min(mmAUC) max(cmAUC)];
eucDist = zeros(length(cmAUC),1);
for iHit=1:length(cmAUC)
    eucDist(iHit,1) = sqrt((mmAUC(iHit,1)-ideal(1,1)).^2 + (cmAUC(iHit,1)-ideal(1,2)).^2);
end
eucDistFilltered = eucDist(cmAUC-mmAUC>0,1);
score = (max(eucDistFilltered)-eucDistFilltered)./max(eucDistFilltered);
hitNames2 = hitNames(cmAUC-mmAUC>0,1);
hitAnnotations2 = hitAnnotations(cmAUC-mmAUC>0,1);
[score inxEuc] = sort(score);
hitOrdNames = hitNames2(inxEuc);
hitOrdAnn = hitAnnotations2(inxEuc);
X = categorical(hitOrdNames);
X = reordercats(X,hitOrdNames);
%% Plot sorted collaboration score
% Plots figure 3C
figure
hold on
for k = 1:size(score,1)
    cord = hitOrdAnn(k,:);
    c = bar(k,score(k,1),'FaceColor',map(cord,:));
end
box on
grid on
ylabel('collaboration score filtered')
set(gca,'ylim',[0 1.05])
%save('validationEucDist4.mat','hitOrdNames','score','hitOrdAnn')
%% Load screen Results
[~,CMgene] = xlsread("OD5 - log2 and pval.xlsx", 'A2:A3545'); % this information is in supplementary table 1, sheet 'Cond. media 5X10^9 cells per mL'
CMfc = xlsread("OD5 - log2 and pval.xlsx", 'C2:C3545');
%find coordinates of hits in R output
CMcord = zeros(length(hitOrdNames),1);
for iGene = 1:length(hitOrdNames)
    CMcord(iGene) = find(strcmp(hitOrdNames(iGene),CMgene));
end
hitFC = CMfc(CMcord,:);
%% Scatter of collaboration score and screen fold-change
% Plots figure 3C inset
figure
for k = 1:size(score,1)
    cord = hitOrdAnn(k,:);
    plot(hitFC(k,1),score(k,1),'o','MarkerFaceColor',map(cord,:),...
        'MarkerEdgeColor',map(cord,:),'MarkerSize',14);
    hold on
end
%text(hitFC,score,hitOrdNames)
%set(gca,'ylim',[1 3])
xlabel('collaboration score')
grid on;box on; axis square
[r p] = corr(hitFC,score);
%%
score = (max(eucDistFilltered)-eucDist)./max(eucDistFilltered);
% save('eucDistAll.mat', 'score','hitNames')