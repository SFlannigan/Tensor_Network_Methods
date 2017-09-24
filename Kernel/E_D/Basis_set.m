function [B] = Basis_set( N,M )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nk=nchoosek(M+N-1,N);
B = zeros(nk,M);

B(1,1) = N;
count = 2;
while B(end,end) ~= N
    for count2 = 1:M-1
        if B(count-1,count2) ~=0
            if all(B(count-1,(count2+1):M-1) == 0)
                if count2 ~=1
                    B(count,1:(count2-1)) = B(count-1,1:(count2-1));
                end
                B(count,count2) = B(count-1,count2) - 1;
                
                B(count,count2+1) = N - sum(B(count,:));
                if count2+1 ~= M
                    B(count,(count2+2):end) = 0;
                end
            end
        end  
    end
    count = count+1;
end

end

