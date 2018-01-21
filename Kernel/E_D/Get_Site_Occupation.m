function [ n ] = Get_Site_Occupation(W,B)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


n = zeros(length(B(1,:)),1);
for count = 1:length(B(1,:))
    n(count) = W'*diag(B(:,count))*W;   
end


end

