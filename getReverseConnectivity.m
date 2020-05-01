function [ aPostBList, numAPostBList ] = getReverseConnectivity( bPostAList, numBPostAList, numBs, oppListMaxSize )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes heren
% make numFBPreInh and fBPreInh
numAs = size(bPostAList,1);

aPostBList = zeros(numBs,oppListMaxSize);
numAPostBList = zeros(numBs,1);
for i = 1 : numAs
    aCell = i;
    bCells = bPostAList(aCell,1:numBPostAList(aCell));
    
    numAPostBList(bCells) = numAPostBList(bCells)+1;
    for j = 1 : numBPostAList(aCell)
        bCell = bCells(j);
        aPostBList(bCell,numAPostBList(bCell)) = aCell;
    end
end

end

