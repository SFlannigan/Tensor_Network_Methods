function [F] = Binary_Search(T,T_check)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

check=0;
counter = floor(length(T)/2);
remain = 0;
add = floor(counter/2);
while check==0
    if T_check>T(counter)
        counter = counter + add;
        if counter >length(T)
            counter=length(T);
        end
        add = floor((add+2*remain)/2);
        remain = add/2-floor(add/2);
    elseif T_check<T(counter)
        counter = counter - add;
        if counter <1
            counter=1;
        end
        add = floor((add+2*remain)/2);
        remain = add/2-floor(add/2);
    elseif T_check==T(counter)
        F=counter;
        check = 1;
    end
end

end

