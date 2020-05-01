function [ output ] = randRound( m )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r = rand(size(m));
f = floor(m);
c = m-f;

output = round(double(f) + double(r < c));

end

