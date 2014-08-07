%% Load data and variable names
clear all;
loaddata;
format bank;
clearvars -except filename Y_name sep Z normalization destAddress dataSize ...
          wlperformance wlpower fontsize saveDataAddress X Y_original energy ...
          LEO fontsize color* marker* legend* normValues;

numSamples = 20;
performanceId = 1;
systemPowerId = 2;
linewidth = 2;
idlePower = 75;
idlePerformance = 0;
%% Initialize other variables
n = dataSize;
m = length(filename);

truePerfWl = cell2mat(Z{1,1});
truePowerWl = cell2mat(Z{1,2});

Reg = '';
folderName = strcat(saveDataAddress, Reg,Y_name{performanceId},'_',...
                    num2str(numSamples),'_',normalization{performanceId},'.mat');
load(folderName); 
estPerfWl = cell2mat(wl);

folderName = strcat(saveDataAddress,Reg, Y_name{systemPowerId},'_',...
                    num2str(numSamples),'_',normalization{systemPowerId},'.mat');
load(folderName);
estPowerWl = cell2mat(wl);

Reg = 'Reg';
folderName = strcat(saveDataAddress, Reg,Y_name{performanceId},'_',...
                    num2str(numSamples),'_',normalization{performanceId},'.mat');
load(folderName); 
RegPerfWl = cell2mat(wl);
folderName = strcat(saveDataAddress,Reg, Y_name{systemPowerId},'_',...
                    num2str(numSamples),'_',normalization{systemPowerId},'.mat');
load(folderName);
RegPowerWl = cell2mat(wl);

Reg = 'QuadReg';
folderName = strcat(saveDataAddress, Reg,Y_name{performanceId},'_',...
                    num2str(numSamples),'_',normalization{performanceId},'.mat');
load(folderName); 
QuadRegPerfWl = cell2mat(wl);
folderName = strcat(saveDataAddress,Reg, Y_name{systemPowerId},'_',...
                    num2str(numSamples),'_',normalization{systemPowerId},'.mat');
load(folderName);
QuadRegPowerWl = cell2mat(wl);

%% Renormalize data

for i = 1:m
    t = normValues{1}(:,i);  
    truePerformance{i}= (truePerfWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);    
    estPerformance{i} = (estPerfWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);    
    RegPerformance{i} = (RegPerfWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);
    QuadRegPerformance{i} = (QuadRegPerfWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);
    
    t = normValues{2}(:,i); 
    truePower{i} = (truePowerWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);    
    estPower{i} = (estPowerWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);
    RegPower{i} = (RegPowerWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);
    QuadRegPower{i} = (QuadRegPowerWl(:,i)/100 -1)*(max(t) - min(t)) + min(t);
    
    truePerformance{i}(n+1) = idlePerformance;    
    estPerformance{i}(n+1) = idlePerformance;    
    RegPerformance{i}(n+1) = idlePerformance;
    QuadRegPerformance{i}(n+1) = idlePerformance;
    
    truePower{i}(n+1) = idlePower;    
    estPower{i}(n+1) = idlePower;    
    RegPower{i}(n+1) = idlePower;
    QuadRegPower{i}(n+1) = idlePower;
end
%% Create folder for plots
folderName = strcat(destAddress,'LP_',num2str(numSamples),'_',normalization{1}); 
mkdir(folderName);
n = n+1;
options=optimset('Display', 'off');
for i = 1:m
    %% Solve the linear program    
    W = max(truePerformance{i}); 
    W1 = min(truePerformance{i}); 
    w = W1:(W - W1)/100:W;
    T = 1.1; 
    Aeq = ones(1, n);
    beq = T;
    lb = zeros(n ,1);
    ub = T*ones(n ,1);
%     KY{i} = unique(convhull(truePerformance{ i } , truePower{ i }),'stable');
%     Ky{i} = unique(convhull(estPerformance{ i } , estPower{ i }),'stable');
%     KyReg{i} = unique(convhull(RegPerformance{ i } , RegPower{ i }),'stable');
    for j = 1:length(w)
           %% For true data
        f = truePower{i};
        A = -truePerformance{i}';
        b = -w(j);       
        [x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
        trueEnergy(i,j) = fval; 
        %% For estimated data
        f = estPower{i};
        A = [-truePerformance{i}';-estPerformance{i}'];
        b = [-w(j),-w(j)];
        [xest,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);    
        estEnergy(i,j) = truePower{i}'*xest; 
        if(xest >1)            
            pause;
        end
        
        %% For simple regression data
        f = RegPower{i};
        A = [-truePerformance{i}';-RegPerformance{i}'];
        b = [-w(j),-w(j)];
        [xest,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);    
        RegEnergy(i,j) = truePower{i}'*xest; 
        if(xest >1)
            pause;
        end
        
        %% For simple regression data
        f = QuadRegPower{i};
        A = [-truePerformance{i}';-QuadRegPerformance{i}'];
        b = [-w(j),-w(j)];
        [xest,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);    
        QuadRegEnergy(i,j) = truePower{i}'*xest; 
        if(xest >1)
            pause;
        end
        %% Race to idle Time
        Trun = T*min(w(j)/truePerformance{i}(1024),1); 
        Tidle = max(T - Trun,0);
        %idleEnergy(i,j) = Trun*truePower{i}(1024) + Tidle*75; % may be min power is better
        idleEnergy(i,j) = Trun*truePower{i}(1024) + Tidle*min(truePower{i}); % may be min power is better

        
    end
%     A = [X(KY{i},:),truePerformance{ i }(KY{i}) , truePower{ i }(KY{i})];
%     HullConfig{i} = A;
    %% Plot 
    close all;
    clear f j;
    h = figure;     
    hold on;      
    x = w/max(w);    
    M = energy(1,i)*[estEnergy(i,:); QuadRegEnergy(i,:); RegEnergy(i,:);   trueEnergy(i,:); idleEnergy(i,:)];    
    mean(M,2)
    [s1,s2] = size(M);
    for j = 1:s1
        a =  M(j,:);
        g = plot( x, a );   
        f(j) = plot( x(1:20:length(x)), a(1:20:length(x)),markerType{j}(2));  
        set([f(j),g],'Color',colorMat(j,:),'LineWidth',linewidth,'MarkerSize',markerSize);
    end
    
    legend(f,legendL5a,'Location','NorthWest','FontSize',fontsize);
    title(filename{i}, 'FontSize', fontsize);
    ylabel('Energy(in J)','FontSize', fontsize);
    xlabel('Utilization','FontSize', fontsize);
    saveas(h, strcat(folderName,sep, filename{i}) ,'png');
    hold off;  

    %pause;
end
variableName = strcat(saveDataAddress,'LP_',num2str(numSamples),'_',normalization{1},'.mat');
save(variableName, 'estEnergy', 'QuadRegEnergy', 'RegEnergy',   'trueEnergy', 'idleEnergy');
%%


