%% Interaction network by functional group
% takes the pairwise interactions, pools them together by their
% functional group and builds a network
% makes input to plot figure 4F using cytoscape
% 2020/11/13
%% load variables
load areaThreshold2.mat
areaBiMat = areaIntMatrix(inx,inx);
areaBiMat(areaBiMat < 1.5) = 0;
areaBiMat(areaBiMat >= 1.5) = 1;
areaBiMat = areaBiMat(categories,categories);
% catName(1:17,1) = {'amino acid'};
% catName(18:22,1) = {'carbon'};
% catName(23:25,1) = {'nucleotide'};
% catName(26:38,1) = {'vitamin'};
%% Pools interactions by group
aa = [1:17];
nuc = [23:24];
carb = [18:22];
bio = [26:30];
nad = [31:33];
pan = [34:37];
funCord = {aa, nuc, carb, bio, nad, pan};
for iFun = 1:length(funCord)
    curCord1 = funCord{iFun};
    for jFun = 1:length(funCord)
        curCord2 = funCord{jFun};
        funIntMat(iFun,jFun) = sum(areaBiMat(curCord1,curCord2),'all');
        funIntMat(jFun,iFun) = sum(areaBiMat(curCord1,curCord2),'all');
    end    
end
%% Universe of interactions
% calculates the number of possible interactions of strains within one
% group with strains within another group, for all of the groups

for iFun = 1:length(funCord)
    curCord1 = funCord{iFun};
    for jFun = 1:length(funCord)
        curCord2 = funCord{jFun};
        funUniverseMat(iFun,jFun) = length(curCord1)*length(curCord2);
        funUniverseMat(jFun,iFun) = length(curCord1)*length(curCord2);
    end    
end
%% Normalize number of interactions by total possible interactions
normFunIntMat = funIntMat./funUniverseMat;
%% Creates network
funNames = {'aa','nuc','carb','bio','nad','pan'};
fg = graph(normFunIntMat,funNames);
T = fg.Edges;
% plot(fg)
% writetable(T, 'areaInteractionGoupsT15RelToAll.xlsx'); % input for
% plotting using cytoscape