%% EM - Algorithm for implementation of heirarchical Bayes net
clear all;
format bank;
loaddata;
clearvars -except Z filename Y_name sep normalization destAddress ...
                  saveDataAddress epsilon unionSet fontsize X;
Y_nameId = 1; %% Change this variable for performance = 1/system-power = 2;
y = Z{1,Y_nameId};  

%% Parameters
n = 1024;
m = length(filename);
pi = 1;
tau = 1;
epsilon = epsilon(:,Y_nameId);
iteration_limit = 50;
accuracy = zeros(m,1);
model = 'Quadratic';
whichstats = {'yhat','r','mse','rsquare','adjrsquare','beta','leverage','cookd','tstat'};

%% Random sampling
% W = sparse(ones(n,1)); % Your matrix
% K = ceil(noise*n); % Number of elements to pull out of the matrix
% I = randperm(numel(W));
% W(I(1:K)) = 0; % The random entries

W = zeros(n,1);
numSamples = 50;
id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
% points from convex hull
id2 = unionSet(1:ceil(length(unionSet)/numSamples):length(unionSet)); % uniform points % works for kmeans% not for nn

randid = randperm(length(unionSet));
id3 = sort(unionSet(randid(1:numSamples))); %random points 

imppoints = [226.00
        255.00
        256.00
        482.00
        511.00
        512.00
        994.00
       1023.00
       1024.00];
id4 = [sort(unionSet(randid(1:(numSamples - length(imppoints))))); imppoints]; %random points 

W(id1) = 1;

Y_name{Y_nameId}
numSamples = sum(W)
%% Create folder
folderName = strcat(destAddress,Y_name{Y_nameId},'_',num2str(numSamples),...
            '_',normalization{Y_nameId }); % string used for file name as well
%mkdir(folderName);
%% LOG File 
[idum,hostname]= system('hostname');
LogfileName = strcat('LOG_EMbayes',hostname,'.txt');
fileID = fopen(LogfileName,'at+');
str = strcat('------- \t',datestr(clock),'\t', hostname,'\t',Y_name{Y_nameId}...
            ,'-------- Normalization',normalization{Y_nameId },'\n');
fprintf(fileID, str);
fprintf(fileID,'i \t filename \t accuracy  \t elapsedTime \t n \t samples \t m \t epsilon \t iteration_limit \t error\n');
file_entry1 = [n, numSamples, m, iteration_limit];
fclose(fileID);
%%
variance = 0;
%%
%DifficultApps{1,:} = [3,6,9,14,21,24];
DifficultApps{1,:} = [14,21];
DifficultApps{2,:} = [14,15,21];
DifficultApps = DifficultApps(Y_nameId,:);
DifficultApps = DifficultApps{1};

for i = 7
%for j = 1:numel(DifficultApps)
 %   i = DifficultApps(j);
    %% initializing the vector with missing values
    Y_known = W.*y{i};
    y_em = y; 
    y_em{i} = Y_known;
    I = diag(W.*ones(n,1));
    %% Start EM algorithm 
    fprintf('i = %f, App = %s\n',i, filename{i});
    SupportData = cell2mat(y); 
    SupportData(:,i) = [];
    Old.mu = mean(SupportData,2);    
    Old.C = speye(n);
    Old.sigma = 1; 
    
    %
    Y_sample = y{i}(id1);        
    [a,b] = size(X);
    X_sample = X(id1,:);
    stats = regstats(Y_sample,X_sample,model,whichstats);
    Old.mu = x2fx(X, model)*stats.beta;   
    %
    
    
    startMain = tic;   
    [return_em, it] = runEM4(Old, W, cell2mat(y_em), i ,epsilon(i), iteration_limit, y);
    w_pred{i} = return_em.w;
    residual = w_pred{i} - y{i};
    rss = sum(residual.*residual); 
    tss = sum((y{i}-mean(y{i})).*(y{i}-mean(y{i})));
    residualsquare = 1-rss/tss;
    accuracy(i) = 1-(1-residualsquare)*(n-1)/(n-2);
    
    %% Plot

        
%         figure(24*i+1+1);
%         mat2 = return_em.Cl;    mat2 = mat2(:,:,i); 
%         mat2 = diag(diag(mat2))^(-0.5)*mat2*diag(diag(mat2))^(-0.5);
%         imagesc(mat2);colorbar; colormap(jet);
%         title('Cl(7)');
%         covmat2{1} = mat2;
%         figure(i*2);
%         mat1 = return_em.em_t.C;    
%         mat1 = diag(diag(mat1))^(-0.5)*mat1*diag(diag(mat1))^(-0.5);
%         imagesc(mat1);colorbar; colormap(jet);
%         title('Sigma');
%         covmat1{i} = mat1;


%     close all;
%     h = figure;
%     hold on;
%     plot(y{i},'r','LineWidth',2.5);
%     plot(w_pred{i},'b','LineWidth',2.5); 
%     scatter(1:length(Y_known),Y_known,'c','fill');
%     title(filename{i}, 'FontSize', fontsize);
%     hl = legend('Data','Estimate');
%     set(hl, 'FontSize',fontsize);
%     ylabel( Y_name{Y_nameId},'FontSize', fontsize);
%     xlabel('Configuration index','FontSize', fontsize);
%     annotation('textbox',[0.15 0.8 0.1 0.1], 'String', strcat('Accuracy = ',... 
%                 sprintf('%5.2f\n',accuracy(i))), 'FontSize', fontsize);
%      saveas(h, strcat(folderName,sep, filename{i}),'png');
%      saveas(h, strcat(folderName,sep, filename{i}),'fig');
%     hold off;
    
    %% Log file
    elapsedTime = toc(startMain);
    fileID = fopen(LogfileName,'at+');
    fprintf(fileID, strcat(num2str(i),'\t',filename{i},'\t'));
    file_entry = [ accuracy(i), it, elapsedTime];
    dlmwrite(LogfileName, full([file_entry,file_entry1]),'-append','delimiter','\t');  
    fclose(fileID);
    fprintf('accuracy = %f, time = %f\n\n', accuracy(i), elapsedTime);
    %pause;
end
wl = w_pred;
fileID = fopen(LogfileName,'at+');
str = strcat('                  ******END******               \n');
fprintf(fileID, str);
fclose(fileID);

%% Save variables
variableName = strcat(saveDataAddress,Y_name{Y_nameId},'_',...
                num2str(numSamples),'_',normalization{Y_nameId },'.mat');

%save(variableName, 'wl','accuracy');	
 
