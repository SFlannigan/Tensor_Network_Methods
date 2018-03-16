function [ H_kin ] = Hop_Ham_SawTooth(B,t,tdash)
% % 
            %                       A   A   A
            %                    J'/ \ / \ /
            %  SAWTOOTH LATTICE - B---B---B
            %                       J
            %

p = zeros(1,length(B(1,:)));

r_coord=[];
c_coord=[];
value=[];
for iter = 1:length(B(1,:))
    p(iter) = 100*iter+3;
end

T = sum(p(ones(length(B(:,1)),1),:).^(0.5).*B,2);
[T,ind]=sort(T);

if (length(B(1,:))/2 - floor(length(B(1,:))/2)) == 0
    for v = 1:length(B(:,1))
        for count = 1:(length(B(1,:)))
            if B(v,count) ~=0
                A = transpose(B(v,:));
                C=0;
                if (count/2 - floor(count/2))~=0
                    % Odd sites = B sites
                    A = A(:,ones(1,2));
                    A(count,1) = B(v,count)-1;
                    A(count,2) = B(v,count)-1;

                    if count ~= (length(B(1,:))-1)
                        C(1) = -tdash*sqrt((B(v,count+1)+1)*(B(v,count)));
                        A(count+1,1) = B(v,count+1)+1;
                        C(2) = -t*sqrt((B(v,count+2)+1)*(B(v,count)));
                        A(count+2,2) = B(v,count+2)+1;  
                    else
                        C(1) = -tdash*sqrt((B(v,count+1)+1)*(B(v,count)));
                        A(count+1,1) = B(v,count+1)+1;
                        C(2) = 0;
                        A(1,2) = B(v,1)+1;
                    end
                    
                else
                    %Even sites - A sites
                    A = A(:,ones(1,1));
                    A(count,1) = B(v,count)-1;
                    if count ~= length(B(1,:))
                        C(1) = -tdash*sqrt((B(v,count+1)+1)*(B(v,count)));
                        A(count+1,1) = B(v,count+1)+1;
                    else
                        C(1) = 0;
                        A(1,1) = B(v,1)+1;
                    end
                end
                
                T_r = sum(transpose(p(ones(length(A(1,:)),1),:)).^(0.5).*A);
                for counter = 1:length(A(1,:))
%                     [F] = Binary_Search(T,T_r(counter));
                    F = find(T==T_r(counter));
                    u = ind(F);
                    r_coord = [r_coord,u];
                    c_coord = [c_coord,v];
                    value = [value,C(counter)];
                end
            end
        end    
    end
else
    for v = 1:length(B(:,1))
        for count = 1:(length(B(1,:)))
            if B(v,count) ~=0
                A = transpose(B(v,:));
                C=0;
                if (count/2 - floor(count/2))~=0
                    % Odd sites = A sites
                    A = A(:,ones(1,4));
                    A(count,1) = B(v,count)-1;
                    A(count,2) = B(v,count)-1;
                    if count ~= (length(B(1,:)))
                        C(1) = -tdash*sqrt((B(v,count+1)+1)*(B(v,count)));
                        A(count+1,1) = B(v,count+1)+1;
                        C(2) = -t*sqrt((B(v,count+2)+1)*(B(v,count)));
                        A(count+2,2) = B(v,count+2)+1;  
                    else
                        C(1) = 0;
                        A(1,1) = B(v,1)+1;
                        C(2) = 0;
                        A(1,2) = B(v,1)+1;
                    end
                    
                else
                    %Even sites - B sites
                    A = A(:,ones(1,2));
                    A(count,1) = B(v,count)-1;
                    C(1) = -tdash*sqrt((B(v,count+1)+1)*(B(v,count)));
                    A(count+1,1) = B(v,count+1)+1;
                end
                
                T_r = sum(transpose(p(ones(length(A(1,:)),1),:)).^(0.5).*A);
                for counter = 1:length(A(1,:))
%                     [F] = Binary_Search(T,T_r(counter));
                    F = find(T==T_r(counter));
                    u = ind(F);
                    r_coord = [r_coord,u];
                    c_coord = [c_coord,v];
                    value = [value,C(counter)];
                end
            end
        end    
    end
end

H_kin = sparse(r_coord,c_coord,value,length(B(:,1)),length(B(:,1)));
H_kin_conj = sparse(c_coord,r_coord,conj(value),length(B(:,1)),length(B(:,1)));

H_kin = H_kin+H_kin_conj;

end

