function [ output ] = rowmap( func, A )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if size(A,1) > 0
    firstResult = func(A(1,:));
    output = zeros(size(A,1),size(firstResult,2));
    if size(A,1)>0
        output(1,:) = firstResult;
        for i = 2 : size(A,1)
            output(i,:) = func(A(i,:));
        end
    end
else
    output = [];
end

end

