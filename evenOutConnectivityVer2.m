function [ numToAdd, numToDel ] = evenOutConnectivityVer2( numPreList, locPost, boost, connString )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numPreList = double(numPreList);

%% get mean connectivity
globalMeanConnectivity = mean( numPreList( locPost(:,1) > 500 &...
                                           locPost(:,1) < 1000 &...
                                           locPost(:,2) > 500 &...
                                           locPost(:,2) < 1000 ) );

%%
a=0;
b=1;
targets = [[0 0 0 a a a a 0 0 0]  [0 0 0 a a a a 0 0 0]  [0 0 0 a a a a 0 0 0];  ...
           [0 0 0 0 a a 0 0 0 0]  [0 0 0 0 a a 0 0 0 0]  [0 0 0 0 a a 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [a 0 0 0 0 0 0 0 0 b]  [b 0 0 0 0 0 0 0 0 a]  [a 0 0 0 0 0 0 0 0 b];  ...
           [a a 0 0 0 0 0 0 b b]  [b b 0 0 0 0 0 0 a a]  [a a 0 0 0 0 0 0 b b];  ...
           [a a 0 0 0 0 0 0 b b]  [b b 0 0 0 0 0 0 a a]  [a a 0 0 0 0 0 0 b b];  ...
           [a 0 0 0 0 0 0 0 0 b]  [b 0 0 0 0 0 0 0 0 a]  [a 0 0 0 0 0 0 0 0 b];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [a 0 0 0 0 0 0 0 0 b]  [b 0 0 0 0 0 0 0 0 a]  [a 0 0 0 0 0 0 0 0 b];  ...
           [a a 0 0 0 0 0 0 b b]  [b b 0 0 0 0 0 0 a a]  [a a 0 0 0 0 0 0 b b];  ...
           [a a 0 0 0 0 0 0 b b]  [b b 0 0 0 0 0 0 a a]  [a a 0 0 0 0 0 0 b b];  ...
           [a 0 0 0 0 0 0 0 0 b]  [b 0 0 0 0 0 0 0 0 a]  [a 0 0 0 0 0 0 0 0 b];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [0 0 0 0 a a 0 0 0 0]  [0 0 0 0 a a 0 0 0 0]  [0 0 0 0 a a 0 0 0 0];  ...
           [0 0 0 a a a a 0 0 0]  [0 0 0 a a a a 0 0 0]  [0 0 0 a a a a 0 0 0];  ...
           ...
           [0 0 0 a a a a 0 0 0]  [0 0 0 a a a a 0 0 0]  [0 0 0 a a a a 0 0 0];  ...
           [0 0 0 0 a a 0 0 0 0]  [0 0 0 0 a a 0 0 0 0]  [0 0 0 0 a a 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [a 0 0 0 0 0 0 0 0 b]  [b 0 0 0 0 0 0 0 0 a]  [a 0 0 0 0 0 0 0 0 b];  ...
           [a a 0 0 0 0 0 0 b b]  [b b 0 0 0 0 0 0 a a]  [a a 0 0 0 0 0 0 b b];  ...
           [a a 0 0 0 0 0 0 b b]  [b b 0 0 0 0 0 0 a a]  [a a 0 0 0 0 0 0 b b];  ...
           [a 0 0 0 0 0 0 0 0 b]  [b 0 0 0 0 0 0 0 0 a]  [a 0 0 0 0 0 0 0 0 b];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0];  ...
           ];

targets = globalMeanConnectivity*ones(30,30) + boost * targets;

numToAdd = zeros(size(numPreList));
numToDel = zeros(size(numPreList));
           
nx = 30;
ny = 30;
for i = 1 : nx
    %report on progress
    disp(['Updating ' connString ' connectivity region row ' num2str(i)]);
    
    for j = 1 : ny
        %get region bounds
        lb = (i-1) / nx * 1500;
        rb = i / nx * 1500;
        bb = (j-1) / ny * 1500;
        tb = j / ny * 1500;
        
        target = targets(j,i);
        
        %pick out cells in region
        picks = find( locPost(:,1) >= lb &...
                      locPost(:,1) < rb &...
                      locPost(:,2) >= bb &...
                      locPost(:,2) < tb );               
        
        %calculate total and mean number connections
        totalNumConnections = sum(numPreList(picks));
        meanNumConnections = totalNumConnections / length(picks);
        
        if meanNumConnections < target        
            %add to make up difference if necessary     
            numToAdd(picks)=randRound(numPreList(picks) * target / meanNumConnections)  - numPreList(picks);
        end
        
        %delete to make up difference if necessary
        if meanNumConnections > target            
            numToDel(picks)=numPreList(picks) - randRound(numPreList(picks) * target / meanNumConnections);
        end                        
    end
end

numToAdd = max(0,numToAdd);
numToDel = max(0,numToDel);


end %function end

