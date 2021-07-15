%% Random Simulation of Synthetic Consortia
% this script calculates an emprical p-value for getting an assembly
% supported by our pairwise interaction matrix given the 38 starting strains
% Plots sup figure 8
% 2021/03/10
%% load complementation data
load('areaThreshold2.mat'); % load matrix with interaction information from complementation experiment
% the next two lines sort by function and alphabetically
intMatrix = areaIntMatrix(inx,inx);
intMatrix = intMatrix(categories,categories);
% the next two lines make the matrix binary
intMatrix(intMatrix < 2) = 0;
intMatrix(intMatrix >= 2) = 1; 
stNameSort = stNameHMsort(categories); % sort the names of the auxotrophs
% remove diagonal
for i = 1:length(intMatrix)
        intMatrix(i,i) = 0;
end 
clear areaIntMatrix degree stNameHM stNameHMsort
%% Get p-value for network size 7
% calculates the p-value of obtaining a fully connected network composed of
% 7 auxotrophs
nStrains = 7;
nSamples = 100000;

sample = zeros(nSamples,nStrains);
for i=1:nSamples
    sample(i,:) = randperm(38,nStrains);
end

score = zeros(1,nSamples);
for iSample = 1:nSamples
    comComp = sample(iSample,:);
    int = zeros(nStrains,nStrains);
    for i = 1:numel(comComp)
        s1 = comComp(1,i);
        for j = 1:nStrains
            s2 = comComp(1,j);
            int(i,j) = intMatrix(s1,s2);
        end
        if sum(sum(int,2)>0) >= nStrains
            score(1,iSample) = 1;
        end
    end
end
pval7 = sum(score)/nSamples;
%% Get p-value for every network 
% calculates the p-value of obtaining a fully connected network composed of
% 1 to 38 auxotrophs

nSamples = 100000;
pval = zeros(1,length(intMatrix));
for nStrains = 1:length(intMatrix)
    score = zeros(1,nSamples);
    sample = zeros(nSamples,nStrains);
    for i=1:nSamples
        sample(i,:) = randperm(38,nStrains);
    end
    for iSample = 1:nSamples
        comComp = sample(iSample,:);
        int = zeros(nStrains,nStrains);
        for i = 1:nStrains
            s1 = comComp(1,i);
            for j = 1:numel(comComp)
                s2 = comComp(1,j);
                int(i,j) = intMatrix(s1,s2);
            end
            if sum(sum(int,2)>0) >= nStrains
                score(1,iSample) = 1;
            end
        end
    end
    pval(1,nStrains) = sum(score)/nSamples;
end
% plot p-value vs network size
figure;
plot([2:length(intMatrix)],pval(2:end),'or','MarkerFaceColor','r','MarkerSize',10)
xlabel('network size','FontSize',18)
ylabel('p-value','FontSize',18)
grid on; box on; axis square

%% weigh p-value by probability of getting N strain community
for i = 1:length(intMatrix)
    totNetworks(i) = nchoosek(38,i);
end
for i = 2:length(intMatrix)
    weighTot(i-1,1) = pval(i)*totNetworks(i)/sum(totNetworks(2:end));
end

%plot p-value weighted by all networks possible
figure;
bar([2:length(intMatrix)],weighTot')
xlabel('network size','FontSize',18)
ylabel('p-value weighted by total number of networks per network size','FontSize',18)
grid on; box on; axis square

%plot log10 number of networks by network size
figure;plot([2:length(intMatrix)],log10(totNetworks(2:end)),'o','MarkerEdgeColor',...
    'k','MarkerFaceColor','k','MarkerSize',10)
xlabel('network size','FontSize',18)
ylabel('log10(number of networks per network size)','FontSize',18)
grid on; box on; axis square

%%
empirical_pval = sum(weighTot);
