function [H_Diag] = Onsite_Ham(B,E)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

coord=[];
value=[];

for v = 1:length(B(:,1))
    onsite=0;
    for u = 1:length(B(1,:))
        onsite=onsite+E(u)*B(v,u);
    end
    coord = [coord,v];
    value = [value,onsite];
end

H_Diag = sparse(coord,coord,value,length(B(:,1)),length(B(:,1)));

end

