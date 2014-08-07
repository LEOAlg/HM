%% EM - Algorithm for implementation of heirarchical Bayes net
clear all;
format bank;
sep = filesep();
dataAddress = strcat('..',sep,'data',sep);
Z{1} = importdata(strcat(dataAddress,'AllPower.txt'));
Z{2} = importdata(strcat(dataAddress,'AllPerformance.txt'));
X = importdata(strcat(dataAddress,'X.txt'));


%% Parameters
[n,m] = size(Z{1});

accuracy = zeros(m,1);

%% Sample

numSamples = 50;
id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
id2 = randperm(n); id2 = id2(1:numSamples); %random points 


%% Create folder
Y_nameId = 1; %% Change this variable for performance = 1/system-power = 2;

parfor Y_nameId = 1:2,
    i = 1;
	[ accuracy, w_pred ] = splitEM( X,Z,Y_nameId,id1,i );
    wl{Y_nameId} = w_pred;
end
