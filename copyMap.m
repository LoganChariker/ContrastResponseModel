function [ newMap ] = copyMap( oldMap )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
newMap=containers.Map();
ks = keys(oldMap);
vals = values(oldMap);
for i = 1 : length(oldMap)
    key = ks{i};
    val = vals{i};
    
    newMap(key)=val;
end

end

