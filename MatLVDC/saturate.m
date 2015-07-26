function [ x ] = saturate(x, lowerBound, upperBound)
%SATURATE(x, lowerBound, upperBound) returns the saturated input value

x = max(x, lowerBound);
x = min(x, upperBound);

end

