%% Temporal dynamics of synthetic consortia - Nutrient rich media
% this script takes the RPM of members in assemblies grown over time in 
% nutrient rich media (LB) and tracks and plots the composition over time. 
% It also calculates the entropy (shannon diversity) of the assemblies over 
% time. Then, it plots them along the entropy of assemblies grown in
% nutrient poor media (M9).
% Plots sup figure 9A-C
% 2021/04/24
%% user definitions
load RPM_COUNTS_QS10.mat; 

nCommunities = 6;
nConditions = 24; % 6 communities * 4 time points

comInx = [1:6:24; 2:6:24; 3:6:24; 4:6:24; 5:6:24; 6:6:24]; % index of communities through time
timePoints = [1 2 3 4];

[inxBio] = xlsread('functionAnnotations.xlsx','F2:F37'); % index of auxotrophs based on their function

% create a colormap with random colors from the built in jet colormap
color = jet(38);
randCord = randperm(38);
myColormap = color(randCord,:);
%% Get RPM values for all community members
comRPM = zeros(size(members,1),nCommunities*2);
k=1;
for iCom = 1:nCommunities
    for iTime = 1:size(timePoints,2)
        comCord = comInx(iCom,iTime);
        for iMemb = 1:size(members,1)
            str = members{iMemb,1};
            memCord = find(strcmp(str,dataset(comCord).hits.names));
            comRPM(iMemb,k) = RPM(memCord,comCord);
        end
        k=k+1;
    end
end
clear dataset
%% Normalize values to get proportion of each auxotroph
normRPM = comRPM./sum(comRPM,1);
normRPM = normRPM(inxBio,:); % sorts members based on their function

%% plot all communities throughout time - strain area
figure
for iCom = 1:nCommunities
    subplot(1,6,iCom)
       aP = area([1 3 7 14],[normRPM(:,(iCom-1)*4+1:(iCom-1)*4+4)'],'EdgeColor','none');
       hold on
    for k = 1:size(members,1)
        aP(k).FaceColor = myColormap(k,:);
    end
    axis square;
    set(gca,'ylim',[0 1],'xlim',[1 14])
end
    legend(members(inxBio))
%% entropy
entroControl = -38*1/38*log2(1/38);
entroCom=zeros(nCommunities,numel(timePoints));
iCom=1;
for jCom = 1:nCommunities
    for iTime = 1:numel(timePoints)
        entro = 0;
        for iMemb = 1:length(comRPM)
            if normRPM(iMemb,iCom) == 0
                continue
            end
            int = normRPM(iMemb,iCom)*log2(normRPM(iMemb,iCom));
            entro = entro -int;
        end
        entroCom(jCom,iTime) = entro;
        iCom = iCom+1;
    end
end

%% plot normalized entropy
normEntropyLB = entroCom./entroControl;
figure; hold on
plot([0 0],[entroControl/entroControl entroControl/entroControl+0.01],'dk','MarkerSize',8)
for iCom = 1:nCommunities
    plot([1 3 7 14],normEntropyLB(iCom,:),'dk','MarkerSize',8)
end
set(gca,'ylim',[0 1.1],'xlim',[0 14])
xlabel('time (days)')
ylabel('normalized entropy')
grid on; box on;

%% plot entropy M9 and LB
load entroM9.mat
figure; hold on;
plot(0,1,'ok','MarkerFaceColor','k','MarkerSize',8)
plot(0,1,'d','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8)
plot([1 3 7 14],normEntropy,'ok','MarkerSize',8)
plot([1 3 7 14],normEntropyLB,'d','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',8)
grid on; box on;
set(gca,'ylim',[0 1.1],'xlim',[0 14])

%% fit entropy of consortia in M9 and LB
mod = 'b*(1*exp(-x*1))+(1-b)';
mod2 = 'a*x^2-b*x+c';
x = [0 1 3 7 14]';
y1 = [1 mean(normEntropy)]';
y2 = [1 mean(normEntropyLB)]';
pseudo_x = linspace(min(x),max(x),100);
f1 = fit(x,y1,mod);
f2 = fit(x,y2,mod2);

%% Plot individual consortia, mean of consortia, and fit
figure; hold on;
errorbar([0 1 3 7 14],[1 mean(normEntropy)],std([[1 1 1 1 1 1]',normEntropy]),'ok','MarkerFaceColor','k','MarkerSize',8)
errorbar([0 1 3 7 14],[1 mean(normEntropyLB)],std([[1 1 1 1 1 1]', normEntropyLB]),'d','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8)
plot([1 3 7 14],normEntropy,'ok','MarkerSize',8)
plot([1 3 7 14],normEntropyLB,'d','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',8)
plot(pseudo_x,f2(pseudo_x),'--','Color',[0.5 0.5 0.5]);
plot(pseudo_x,f1(pseudo_x),'k');
grid on; box on;
set(gca,'ylim',[0 1.1],'xlim',[0 14])
xlabel('time (days)')
ylabel('normalized entropy')
%% community last timepoint
comLabels = {'1','2','3','4','5','6'};
figure; cB = bar(normRPM(:,[4:4:24])','stack','EdgeColor','none');
for k = 1:size(members,1)
cB(k).FaceColor = myColormap(k,:);
end
set(gca,'ylim',[0 1],'xtick',[1:6],'xticklabel',comLabels)
axis square;
legend(members(inxBio))
set(gcf,'position',[440,310,800,456])