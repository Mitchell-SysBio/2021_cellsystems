%% Cross-feeding in strain pairs
% this script analyzes absorbance data from co-cultures to identify pairs
% that improve their growth in co-culture relative to individual cultures.
% It uses the area under the growth curve to assess cooperation
% Plots figure 4A-E
%% Load 
load('complementation.mat');
mStrain=38;
stNameHM{11} = 'icd'; % change first letter to lowercase
stNameHM{20} = 'panZ'; % change annotated name 'yhhK' to its synonym 'panZ'
strainNames{21} = 'panZ';
categories = [1 2 3 9 10 12 13 14 15 16 17 18 28 35 36 37 38 11 27 29 30 ...
34 31 32 33 4 5 6 7 8 19 20 21 22 23 24 25 26];
%% Max OD alphabetically and by category
[stNameHMsort, inx] = sort(stNameHM);
maxODsort = maxODMatrix(inx,inx);
% figure %heatmap
% h = heatmap(stNameHMsort(categories),stNameHMsort(categories),maxODsort(categories,categories));
% h.Colormap = gray;
% h.ColorLimits = [0 1]
% h.Title = strcat('Max OD sorted + cat');
%% Check replicates correlation
% figure
% k = 1;
% hold on
% for i=1:length(maxODsort)
%     for j=(i+1):length(maxODsort)
%         plot(maxODsort(i,j),maxODsort(j,i),'o')
%         junk(k,1) = maxODsort(i,j);
%         junk(k,2) = maxODsort(j,i);
%         k = k+1;
%     end
% end
% figure
% plot(junk(:,1),junk(:,2),'o')
% [r p] = corr(junk(:,1),junk(:,2));

%% Area Matrix
% figure 4B
areaMat = zeros(mStrain,mStrain);
for i = 1:mStrain
    iStrain=strains(i);
    curCor1 = find(coordinates1 == iStrain);
    for j = 1:mStrain
        jStrain=strains(j);
        curCor2 = find(coordinates2 == jStrain);
        curInx = intersect(curCor1, curCor2);
        areaMat(i,j) = trapz(timeArray,mean(OD600(:,curInx),2));
    end
end
areaMatOrd = areaMat(inx,inx);
figure 
h = heatmap(stNameHMsort(categories),stNameHMsort(categories),areaMatOrd(categories,categories));
h.Colormap = gray;
% h.ColorLimits = [0 5];
h.Title = strcat('Area alph & cat');
%% Area Difference Matrix
areaMatrix = zeros(mStrain,mStrain);
for i = 1:mStrain
    iStrain=strains(i);
    curCor1 = find(coordinates1 == iStrain);
    for j = 1:mStrain
        jStrain=strains(j);
        curCor2 = find(coordinates2 == jStrain);
        curInx = intersect(curCor1, curCor2);
        areaMatrix(i,j) = trapz(timeArray,mean(OD600(:,curInx),2))-...
            max(min(trapz(timeArray,OD600(:,curCor1))),min(trapz(timeArray,OD600(:,curCor2))));
    end
end
%% histogram of AUC difference values
% figure 4C
figure
histogram(areaMatrix(:,:),25,'normalization','probability', 'Orientation', 'horizontal')
title('AUC difference values - normalized by row & col')
ylabel('proportion')
xlabel({'Improvement by collaboration' 'normalized AUC diference'})
set(gca,'ylim',[-0.21 10.1]);
set(gcf,'position',[440,219,310,579])
box on; grid on
%% histogram absolute area values
% sup figure 5b
areaTotMatrix = zeros(mStrain,mStrain);
for i = 1:mStrain
    iStrain=strains(i);
    curCor1 = find(coordinates1 == iStrain);
    for j = 1:mStrain
        jStrain=strains(j);
        curCor2 = find(coordinates2 == jStrain);
        curInx = intersect(curCor1, curCor2);
        areaTotMatrix(i,j) = trapz(timeArray,mean(OD600(:,curInx),2));
    end
end
areaSelfPairs = zeros(38,1);
for i = 1:38
    areaSelfPairs(i,1) = areaTotMatrix(i,i);
end
junk = areaTotMatrix;
for i = 1:38
    junk(i,i) = NaN;
end

figure;
hold on;
histogram(areaSelfPairs,'BinEdges',0:0.5:10,'Normalization','Probability', 'Orientation', 'horizontal','LineWidth',1,'FaceColor','k','EdgeAlpha',0.5)
histogram(junk,'BinEdges',0:0.5:10,'Normalization','Probability', 'Orientation', 'horizontal','LineWidth',1,'FaceColor','r','EdgeColor','r','EdgeAlpha',0.5)
box on
set(gcf,'position',[440,265,248,533])
legend('singleStrain','pairs')
xlabel('proportion')
ylabel('AUC')
grid on
%% make area difference matrix symmetric
% makes symmetric by taking the mean of the two replicates
areaIntMAtrix=zeros(mStrain,mStrain);
for iRow = 1:nStrain-2
    for iCol = 1:nStrain-2
        areaIntMatrix(iRow,iCol) = mean([areaMatrix(iRow,iCol) areaMatrix(iCol,iRow)]);
    end
end
%% heatmap alphabetically and by category
% plot sorted by function and alphabetically
areaOrd = areaIntMatrix(inx,inx);
figure 
h = heatmap(stNameHMsort(categories),stNameHMsort(categories),areaOrd(categories,categories));
h.Colormap = gray;
% h.ColorLimits = [0 5];
h.Title = strcat('Area difference - alph & cat');
%% Make binary heatmap with threshold
areaOrd = areaIntMatrix(inx,inx);
areaBiMat = areaOrd;
areaBiMat(areaBiMat < 2) = 0;
areaBiMat(areaBiMat >= 2) = 1;
% figure
% ha = heatmap(stNameHMsort(categories),stNameHMsort(categories),areaBiMat(categories,categories));
% ha.Colormap = pink;
% ha.ColorLimits = [0 1];
% ha.CellLabelColor = 'none'
% ha.Title = strcat('Area difference - Binary');
%% find number of interactions
degree = zeros(mStrain,1);
for i=1:mStrain
    degree(i,1) = sum(areaBiMat(i,:)>0);    
end
% stNameFilt = stNameHM(degree>0);
areaMatFilt = areaMatrix;
areaMatFilt(degree<=0,:) = [];
areaMatFilt(:,degree<=0) = [];
areaMatFilt(areaMatFilt<0) = 0;
% writetable(T, 'areaInteractionEdgesRCnormalizedT2.xlsx');
% save('areaThreshold2.mat','areaIntMatrix', 'inx', 'stNameHMsort', 'categories','degree','stNameHM')


%% Plot example growth curves
gcStrain = [20 21]; % individual strains to plot
for i = 1:length(gcStrain)
     figure
hold on;
    plot(timeArray, selfPairs(:,gcStrain(i)),'o-k','MarkerFaceColor','k');
    set(gca,'ylim',[0 max(OD600,[],'all')]);
    title(strcat('Δ',strainNames(gcStrain(i))))
end
%% growth curves replicates
% get mean absorbance for all pairs to plot
meanOD=cell(40,40);
sdOD=cell(40,40);
k=1;
for iStrain = 1:nStrain
    for jStrain = 1:nStrain
        curCor1 = find(coordinates1 == iStrain);
        curCor2 = find(coordinates1 == jStrain);
        curCor3 = find(coordinates2 == iStrain);
        curCor4 = find(coordinates2 == jStrain);
        curInx = unique([intersect(curCor1, curCor4) intersect(curCor2, curCor3)]);
        meanOD{iStrain,jStrain} = nanmean(OD600(:,curInx),2);
        sdOD{iStrain,jStrain} = nanstd(OD600(:,curInx),[],2);
        k=k+1;
    end
end
meanCordinates = [];
for i=1:nStrain
    meanCordinates = [meanCordinates linspace(i,40,nStrain-i+1)];
end
junk = [];
for i=1:nStrain
    junk = [junk linspace(i,i,nStrain-i+1)];
end
meanCordinates=[meanCordinates; junk];
clear junk

%% plot mean pair absorbance over time
s1 = 20;
s2 = 21;
cordStrain = intersect(find(meanCordinates(2,:)==s1),find(meanCordinates(1,:)==s2));
figure
hold on;
y = meanOD{meanCordinates(2,cordStrain),meanCordinates(1,cordStrain)};
plot(timeArray,y,'o-','MarkerFaceColor',[1.0000    0.9062         0],'MarkerEdgeColor',   [1.0000    0.9062         0],'Color',[1.0000    0.9062         0])
set(gca,'ylim',[0 max(OD600,[],'all')]);
title(strcat('Δ',strainNames(meanCordinates(2,cordStrain)),' - ','Δ', strainNames(meanCordinates(1,cordStrain))))

