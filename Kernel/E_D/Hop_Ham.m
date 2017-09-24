function [ H_kin ] = Hop_Ham(B,t,tdash,Boundary_Conditions)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

p = zeros(1,length(B(1,:)));

r_coord=[];
c_coord=[];
value=[];
for iter = 1:length(B(1,:))
    p(iter) = 100*iter+3;
end

T = sum(p(ones(length(B(:,1)),1),:).^(0.5).*B,2);
[T,ind]=sort(T);
for v = 1:length(B(:,1))
    for count = 1:(length(B(1,:)))
        if B(v,count) ~=0
            A = transpose(B(v,:));
            C=0;
            if (count/2 - floor(count/2))~=0
                % Odd sites = A sites
                A = A(:,ones(1,1));
                A(count,1) = B(v,count)-1;
                
                if count ~= (length(B(1,:))-1)
                    C(1) = -t*sqrt((B(v,count+1)+1)*(B(v,count)));
                    A(count+1,1) = B(v,count+1)+1;
                else
                    C(1) = -t*sqrt((B(v,count+1)+1)*(B(v,count)));
                    A(count+1,1) = B(v,count+1)+1;
                end
                
            else
                %Even sites - B sites
                A = A(:,ones(1,1));
                A(count,1) = B(v,count)-1;
                if count ~= length(B(1,:))
                    C(1) = -tdash*sqrt((B(v,count+1)+1)*(B(v,count)));
                    A(count+1,1) = B(v,count+1)+1;
                else
                    if strcmp(Boundary_Conditions, 'Periodic')
                        C(1) = -tdash*sqrt((B(v,1)+1)*(B(v,count)));
                        A(1,1) = B(v,1)+1;
                    else
                        C(1) = 0;
                        A(1,1) = B(v,1)+1;
                    end
                end
            end
            
            T_r = sum(transpose(p(ones(length(A(1,:)),1),:)).^(0.5).*A);
            for counter = 1:length(A(1,:))
                L = find(T==T_r(counter));
                u = ind(L);
                r_coord = [r_coord,u];
                c_coord = [c_coord,v];
                value = [value,C(counter)];
            end
        end
    end
end

H_kin = sparse(r_coord,c_coord,value,length(B(:,1)),length(B(:,1)));
H_kin_conj = sparse(c_coord,r_coord,conj(value),length(B(:,1)),length(B(:,1)));

H_kin = H_kin+H_kin_conj;
end

