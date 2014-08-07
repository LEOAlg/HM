clear all;
sep = filesep();
dataAddress = strcat('..',sep,'datafromharper',sep);
saveDataAddress = strcat('SavedData',sep);
%destAddress = strcat('..',sep,'WriteUp',sep,'Figures',sep);

destAddress = strcat('..',sep,'..',sep,'..',sep,'Spring 2014',sep,...
                     'ASPLOS-15',sep,'figures',sep);

filename{1} = 'blackscholes' ;
filename{2} = 'bodytrack';
filename{3} = 'cfd';
filename{4} = 'fluidanimate';
filename{5} =  'HOP';
filename{6} = 'Jacobi' ;
filename{7} = 'Kmeans';
filename{8} = 'lud'; 
filename{9} =  'nn';
filename{10} = 'particlefilter';
filename{11} = 'PLSA';
filename{12} = 'Scalparc';
filename{13} =  'swaptions';
filename{14} =  'Swish';
filename{15} = 'vips';
filename{16} = 'x264'; 
filename{17} = 'apr'; 
filename{18} = 'btree'; 
filename{19} = 'semphy';
filename{20} = 'streamcluster';
filename{21} = 'svmrfe'; 
filename{22} = 'backprop'; 
filename{23} = 'bfs'; 
filename{24} = 'Kmeansnf'; 
filename{25} = 'filebound';

Y_name{1} = 'performance';
Y_name{2} = 'system-power';
Y_name{3} = 'system-energy';

X_name{1} = '#core, frequency';
X_name{2} = '#core, frequency, #memory controller';   

%suffix = '1socket';
suffix = '32threads';
prefix = '';
%prefix = 'Rapl';

numfiles = length( filename );
numY = length(Y_name);
numplots = numfiles*numY;
dataSize = 1024;
epsilon = 10000*ones(numfiles,numY);
id1 = [3,7,8,9,21,24];
epsilon(id1,1) = 40000;
id2 = [15];
epsilon(id2,2) = 40000;

fontsize = 25;
%% Extract data to variable Y.
n = 1024;
   energy = [];
   numSamples = 20;
id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
for i = 1:numfiles     
    Fileplotname = strcat( dataAddress, filename{i}, sep, prefix,...
                               filename{i}, suffix, '.results');
    A = importdata( Fileplotname );
    [a,b] = size(A);
        % removing the on-chip power and energy
    energy(:,i) = A(1:n,b);
    if (b == 7)            
        A(:,5) = [];
    elseif (b == 8)
        A(:,5) = [];
        
        A(:,5) = [];
    end
    power(:,i) = A(1:n,5);
    for j = 1:numY  
        Y = A(:,j+3);        
        a = Y(1:dataSize);
        Y_original{j}{i} = a;
        
        b = a(id1);
        %Z{j}{i} = (100/(max(a) - min(a)))*(a - min(a)) + 100;   	     
        %normalization{j} = '100timesmaxminY'; 
        %Z{j}{i} = (100/(a(n) - a(1)))*(a - a(1)) + 100;   	     
        %normalization{j} = '100timeslowhighY'; 
        Z{j}{i} = (100/(max(b) - min(b)))*(a - min(b)) + 100;   		     
        normalization{j} = '100timesmaxminYest';
        normValues{j}(:,i) = b;
    end 
   
end
time = energy./power;
X = A(1:1024,1:3); % original data
unionSet = importdata( strcat( saveDataAddress, 'Unionestconvexhull2.txt' ));

%% Color and legend
colorMat = [rgb('LawnGreen'); rgb('Crimson');rgb( 'Goldenrod');...
            rgb('Blue'); rgb('Black')]; %gymbk v>^so

markerSize = 6;
markers = 'v>^so';
for(i = 1:length(markers))
    markerType{i} = strcat('-',markers(i));
end

hatchSymbol = '\x-+/';
LEO = 'LEO';
legendL4 = {LEO,   'Online','Offline','Optimal'};
legendL4b ={LEO,  'Online','Offline','Data-points'};
legend4c = {LEO,   'Online','Offline','Race-to-idle'};
legendL5a = legendL4 ;  legendL5a{5} ='Race-to-idle';
legendL5b = legendL4 ;  legendL5b{5} ='Data-points';