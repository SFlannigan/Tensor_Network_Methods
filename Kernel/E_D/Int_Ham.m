function [ H_Int ] = Int_Ham( B,U_A,U_B )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

coord=[];
value=[];
Int_A = zeros(length(B(:,1)),length(B(1,:)));
Int_B = zeros(length(B(:,1)),length(B(1,:)));

for v = 1:length(B(:,1))
    for u = 1:length(B(1,:))
        if (u/2 - floor(u/2))~=0
            Int(v,u) = U_A(u)/2 * B(v,u)*(B(v,u)-1);
        else
            Int(v,u) = U_B(u)/2 * B(v,u)*(B(v,u)-1);
        end
    end 
    coord = [coord,v];
    value = [value,sum(Int(v,:))];
end

H_Int = sparse(coord,coord,value,length(B(:,1)),length(B(:,1)));

end

