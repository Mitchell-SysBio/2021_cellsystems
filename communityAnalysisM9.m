%% Temporal dynamics of synthetic consortia - Nutrient poor media
% this script takes the RPM of members in assemblies grown over time in 
% nutrient poor media (M9) and tracks and plots the composition over time.
% It also calculates the entropy (shannon diversity) of the assemblies over 
% time.
% Plots figure 5A-B, 5D; sup figure 9A, 9C
% 2020/09/15
%% user definitions
fileToSave = 'RPM_COUNTS_QS10.mat'; % name of output file
load dataset_QS202; 

nScreens = 25 ; % 12 communities * 4 time points
nConditions = 24;
nCommunities = 6;
junk = [3 ; 4 ; 5 ; 6 ; 7 ; 8 ; 9 ; 10 ; 11 ; 12 ; 13 ;14];
junk = [junk (junk + 24)];
junk = [junk (junk(:,1) + 36)];
junk = [junk (junk(:,1) + 12)];
comInx = junk(1:6,:);
timePoints = [1 2 3 4];
[inxBio] = xlsread('functionAnnotations.xlsx','F2:F39');
% clear junk
%% make sensible names
for j=1:nConditions % should match size of dataset
    i = comInx(j);
    curFileName = dataset(i).fileName;
    newStr = split(curFileName,["_",'.']); % split by _ or .
    datasetName{j} = newStr{1};
end
%% calculate RPM, COUNTS and compactRPM
RPM = []; % master matrix to hold all RPMs
COUNTS = [];
labels = {}; % labels for RPM matrix
for j = 1:nConditions
    iDataset = comInx(j);
    RPM = [RPM,dataset(iDataset).hits.RPM];
    COUNTS = [COUNTS,dataset(iDataset).hits.counts];
    labels = {labels{:},datasetName{j}};
end
%%
% comMembers = dataset(1).hits.names(RPM(:,1)>500);
[~, comMembers] = xlsread('functionAnnotations.xlsx','A2:A39');

%% filter by RPM
identifiedBarcodeRPM = cell(1,nConditions);
barcodesRPM = cell(1,nConditions);
  for i=1:nConditions
    barcWithCounts = RPM(:,i)>100;
    identifiedBarcodeRPM{1,i} = dataset(i).hits.names(barcWithCounts);
    barcodesRPM{1,i} = RPM(barcWithCounts,i);
end 
figure
for i=1:nConditions
    x = categorical(identifiedBarcodeRPM{1,i});
    subplot(4,6,i)
    bar(x,barcodesRPM{1,i})
    set(gca,'ylim',[0 10^6])
    title(labels(i))
end
set(gcf,'position',[1,200,1440,598])
%% Get RPM values for all community members
comRPM = zeros(size(comMembers,1),nCommunities*2);
k=1;
for iCom = 1:nCommunities
    for iTime = 1:size(timePoints,2)
        comCord = comInx(iCom,iTime);
        for iMemb = 1:size(comMembers,1)
            member = comMembers{iMemb,1};
            memCord = find(strcmp(member,dataset(comCord).hits.names));
            comRPM(iMemb,k) = RPM(memCord,(iTime-1)*6+iCom);    
        end
        k=k+1;
    end
end
%% Normalize values
normComRPM = comRPM./sum(comRPM,1);
normComRPM = normComRPM(inxBio,:);
%% plot all communities throughout the day
load colormap.mat
bioColormap = funcColormap(inxBio,:);
figure
for iCom = 1:nCommunities
    subplot(1,6,iCom)
    bB = bar(timePoints,normComRPM(:,(iCom-1)*4+1:(iCom-1)*4+4),'stacked','EdgeColor','none');
    hold on
    for k = 1:size(comMembers,1)
        bB(k).FaceColor = bioColormap(k,:);
    end
    set(gca,'ylim',[0 1])
end
%% community last timepoint - Sup figure 9C
comLabels = {'1','2','3','4','5','6'};
figure; cB = bar(normComRPM(:,[4:4:24])','stack','EdgeColor','none');
cB(1).FaceColor = 'flat';
for k = 1:size(comMembers,1)
cB(k).FaceColor = bioColormap(k,:);
end
set(gca,'ylim',[0 1],'xtick',[1:6],'xticklabel',comLabels)
axis square;
legend(comMembers(inxBio))
%% pie chart last timepoint
aaa = [4:4:24];
figure;
k=1;
for i=1:length(aaa)
    subplot(1,6,k)
    pie(normComRPM(:,aaa(i)));
    k=k+1;
end
%% entropy
entroControl = -38*1/38*log2(1/38);
entroM9=zeros(nCommunities,numel(timePoints));
iCom=1;
for jCom = 1:nCommunities
    for iTime = 1:numel(timePoints)
        entroCom = 0;
        for iMemb = 1:length(comRPM)
            if normComRPM(iMemb,iCom) == 0
                continue
            end
            int = normComRPM(iMemb,iCom)*log2(normComRPM(iMemb,iCom));
            entroCom = entroCom -int;
        end
        entroM9(jCom,iTime) = entroCom;
        iCom = iCom+1;
    end
end
%% plot entropy
normEntropy = entroM9./entroControl;
mod = 'b*(1*exp(-x*1))+(1-b)';
x = [0 1 3 7 14]';
y = [1 mean(normEntropy)]';
pseudo_x = linspace(min(x),max(x),100);
f1 = fit(x,y,mod);
marker = {'d','s','*','p','h','^'};
figure
hold on
for i=1:4
    plot([1 3 7 14],normEntropy(i,:),marker{i},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',8)
end
errorbar([1 3 7 14],mean(normEntropy),std(normEntropy),'ok','MarkerFaceColor','k','MarkerSize',8)
plot(0,1,'ok','MarkerFaceColor','k','MarkerSize',8)
plot(pseudo_x,f1(pseudo_x),'k');
set(gca,'ylim',[0 1.1],'xlim',[0 14])
grid on; box on;
xlabel('time (days)')
ylabel('normalized entropy')
save entroM9.mat normEntropy
%% plot all communities throughout time - strain area
% modCord = [randCord(1:4) randCord(end-1:end) randCord(5:end-2)];
% newColormap = junk(modCord,:);
figure
for iCom = 1:nCommunities
    subplot(1,6,iCom)
    aP = area([0 1 3 7 14],[linspace(1/38,1/38,38); normComRPM(:,(iCom-1)*4+1:(iCom-1)*4+4)'],'EdgeColor','none');
    hold on
    for k = 1:size(comMembers,1)
        aP(k).FaceColor = bioColormap(k,:);
    end
    set(gca,'ylim',[0 1],'xlim',[0 14])
    axis square
end
legend(comMembers(inxBio),'NumColumns',2)
set(gcf,'position',[440,310,928,456])

