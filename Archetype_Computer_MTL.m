function H=Archetype_Computer_MTL(Data)

% Function for computing MTL archetype scores on new data

Data=Data./sum(Data,1); %ensures proper normalization
W=readtable('Archetype_Gene_Scores.csv'); % provides the archetype weights for each gene
W=W{2:end,2:end};
H=NMF_New_Weights(Data,W,12);