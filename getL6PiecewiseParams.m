function [ outxs, outys, outslopes ] = getL6PiecewiseParams( xs, ys )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if xs(end) > 900
    outxs = [xs];
    outys = [ys];
    outslopes = ( outys(2:end) - outys(1:end-1) ) ./ ( outxs(2:end) - outxs(1:end-1) );
else
    outxs = [xs, 1000];
    outys = [ys, ys(end) + (ys(end)-ys(end-1))/(xs(end)-xs(end-1))*(1000-xs(end))];
    outslopes = ( outys(2:end) - outys(1:end-1) ) ./ ( outxs(2:end) - outxs(1:end-1) );
end
end

