classdef mps_var
    % Matrix product state formulism that conserves particle number for the
    % Bose-Hubbard model.
    %
    % We use the cell functionality for each tensor at every site. The position
    % within each 2D cell dictates the local dmension and the number of particles
    % to the right (or left).
    %
    % Includes the algorothms, 2-site TDVP, 2-site DMRG, 1-site DMRG,
    % 1-site TDVP (but this latter algorithm does not work, but I am not
    % sure why...)
    %
    % To do:
    %       I think that I need to include more detailed comments on the
    %       time evolution procedures in the new formalism. I think I
    %       should produce a set of notes that goes into the details.
    %
    
    properties
        % Maximum bond dimension.
        Dmax
        
        % Matrix data
        %
        % Structure array of Cells. Each element of structure corresponds
        % to each lattice site and is a 2D Cell array - one dimension is
        % the local dimension and the other is the number of particles to
        % the right (or left).
        %
        % 2D matrices are held in each element the Cell array's
        % corresponding to the two bond dimenions
        data
        
        % Threshold of what is considered zero
        zero_thres
        
        % Number of particles
        N
        
        % Maximum number of particles allowed on a single site
        N_max
        
        % Local dimension
        d
        
        % Vector containing chosen positions of initial particles
        P_Positions
        
        % Local creation and Annihilation operators
        b_dag
        b
        
    end
    
    methods
        function mps_var=mps_var(M,Dmax,N,N_max)
            % Class constructor
            
            mps_var.Dmax=Dmax;
            mps_var.data=cell(1,M);
            mps_var.zero_thres=1e-15;
            mps_var.N=N;
            mps_var.N_max=N_max;
            mps_var.d=N_max+1;
            mps_var.b = sqrt(diag(1:N_max,1));
            mps_var.b_dag = sqrt(diag(1:N_max,-1));
            
        end
        
        function mps_var=set_Particle_Position(mps_var,Vec)
            % Vec == vector containing the positions of the N particles.
            
            if numel(Vec)==mps_var.N
                mps_var.P_Positions=Vec;
            else
                disp(['Particle Position Vector is the Wrong Length.' ...
                    ' Using Random Distribution.']);
            end
        end
        
        function mps_var = Canonicalise_1s(mps_var,Canon,Norm)
           % Single site Canonicalisation procedure
           % 
           % Canon == 'L-R' For Left Canonicalisation
           %       == 'R-L' For Right Canonicalisation
           %
           % Norm == 'true' to normalise the mps
           
           M=size(mps_var.data,2);
           if strcmp(Norm,'true')
               if strcmp(Canon,'L-R')
                   for m = 1:(M)
                       mps_var=mps_var.Do_SVD_Single_Tensor(m,Canon,Norm);
                   end
               elseif strcmp(Canon, 'R-L')
                   for m = M:-1:1
                       mps_var=mps_var.Do_SVD_Single_Tensor(m,Canon,Norm);
                   end
               else
                   error('Must choose canonicalisation');
               end
           else
               if strcmp(Canon,'L-R')
                   for m = 1:(M-1)
                       mps_var=mps_var.Do_SVD_Single_Tensor(m,Canon,Norm);
                   end
               elseif strcmp(Canon, 'R-L')
                   for m = M:-1:2
                       mps_var=mps_var.Do_SVD_Single_Tensor(m,Canon,Norm);
                   end
               else
                   error('Must choose canonicalisation');
               end
           end
        end
        
        function mps_var = Do_SVD_Single_Tensor(mps_var,m,Canon,Norm)
           
            if strcmp(Canon,'L-R')
                [u,v]=size(mps_var.data{m}{1});
                [~,NR]=size(mps_var.data{m});
                NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Do SVD
                Q=cell(mps_var.d,NR);
                S=cell(1,NR);
                for NR_iter=1:NR
                    d_vec=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                    if ~isempty(d_vec)
                        [U,Sv,V]=svd(cell2mat(mps_var.data{m}(d_vec,NR_iter)));
                        
                        Sv(abs(Sv)<mps_var.zero_thres)=0;
                        U(abs(U)<mps_var.zero_thres)=0;
                        V(abs(V)<mps_var.zero_thres)=0;
                        
                        dd_vec = ((d_vec(1)-1)*u+1):(d_vec(end)-1)*u+u;
                        
                        U_temp=zeros(mps_var.d*u,min(size(U,2)/length(d_vec)*mps_var.d,v));
                        U_temp(dd_vec,1:min(size(U,2),v))=U(:,1:min(size(U,2),v));
                        
                        cell_vec=u*ones(mps_var.d,1);
                        Q(:,NR_iter)=mat2cell(U_temp,cell_vec,min(size(U,2)/length(d_vec)*mps_var.d,v));
                        
                        Sv_temp=zeros(size(U_temp,2),v);
                        Sv_temp(1:min(size(U,2),v),:)=Sv(1:min(size(U,2),v),:);
                        
                        S{1,NR_iter}=Sv_temp*(V');
                    end
                end
                S(cellfun(@isempty,S)) = {zeros(v,v)};
                if strcmp(Norm,'true')
                    norm_val=0;
                    for NR_iter=1:NR
                        norm_val=norm_val+S{1,NR_iter}*S{1,NR_iter}';
                    end
                    
                    norm_val=norm_val(1);
                    S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
                end
                Q(cellfun(@isempty,Q)) = {zeros(u,v)};
                mps_var.data{m}=Q;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if m ~= size(mps_var.data,2)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Apply Singular values to next tensor
                    [~,NR_1]=size(mps_var.data{m+1});
                    
                    temp_E = cell(mps_var.d,NR_1);
                    for k = 1:mps_var.d
                        track=1;
                        for k_2 = (mps_var.N+1-k+1):-1:1
                            if k_2 <= NR_1
                                temp_E{k,k_2} = S{1,mps_var.N+1-track+1};
                            end
                            track=track+1;
                        end
                    end
                    temp_E(cellfun(@isempty,temp_E)) = {zeros(v,v)};
                    
                    mps_var.data{m+1} = cellfun(@(y,z) y*z, temp_E,mps_var.data{m+1}, 'UniformOutput', false);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            elseif strcmp(Canon, 'R-L')
                [u,v]=size(mps_var.data{m}{1});
                [~,NR]=size(mps_var.data{m});
                NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Convert mth Tensor to N_Left Labelling
                temp_T = cell(mps_var.d,NR);
                temp_IT = cell(mps_var.d,NL);
                expect = max(NR,NL):-1:1;
                for k = 1:mps_var.d
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2>mps_var.N+1-NR-k+1 && k_2<=NL
                            temp_T(k,k_2) = {expect(1,max(NR,NL)-track+1)};
                        end
                        track=track+1;
                    end
                    for k_2 = max(1,mps_var.N+1-NR-k+2):min(NL,(mps_var.N+2-k))
                        temp_IT(k,k_2) = {k_2};
                    end
                end
                
                temp = mps_var.data{m};
                
                R_mm = cell(mps_var.d,NL);
                for alpha = 1:mps_var.d
                    R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Do SVD
                Q=cell(mps_var.d,NL);
                S=cell(1,NL);
                for NL_iter=1:NL
                    d_vec=max(1,mps_var.N+3-NL_iter-NR):min(mps_var.d,mps_var.N+2-NL_iter);
                    
                    if ~isempty(d_vec)
                        Mat = cell2mat(R_mm(d_vec,NL_iter)); % (d u),v
                        Mat = reshape(Mat,u,length(d_vec),v); % u,d,v
                        Mat = permute(Mat,[1 3 2]); % u,v,d
                        Mat = reshape(Mat,u,length(d_vec)*v); % u,(d v)
                        
                        [U,Sv,V]=svd(Mat);
                        V=V';
                        
                        Sv(abs(Sv)<mps_var.zero_thres)=0;
                        U(abs(U)<mps_var.zero_thres)=0;
                        V(abs(V)<mps_var.zero_thres)=0;
                        
                        dd_vec = ((d_vec(1)-1)*v+1):(d_vec(end)-1)*v+v;
                        
                        V_temp=zeros(min(size(V,1)/length(d_vec)*mps_var.d,u),mps_var.d*v);
                        V_temp(1:min(size(V,1),u),dd_vec)=V(1:min(size(V,1),u),:); % u,(v d)
                        
                        Mat = reshape(V_temp,u,v,mps_var.d); % u,v,d
                        Mat = permute(Mat,[1 3 2]); %u,d,v
                        Mat = reshape(Mat,u*mps_var.d,v);% (d u),v
                        
                        cell_vec=u*ones(mps_var.d,1);
                        Q(:,NL_iter)=mat2cell(Mat,cell_vec,v);
                        
                        Sv_temp=zeros(u,size(V_temp,1));
                        Sv_temp(:,1:min(size(V,1),u))=Sv(:,1:min(size(V,1),u));
                        
                        S{1,NL_iter}=U*Sv_temp;
                    end
                end
                S(cellfun(@isempty,S)) = {zeros(u,u)};
                if strcmp(Norm,'true')
                    norm_val=0;
                    for NL_iter=1:NL
                        norm_val=norm_val+S{1,NL_iter}'*S{1,NL_iter};
                    end
                    
                    norm_val=norm_val(1);
                    S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
                end
                Q(cellfun(@isempty,Q)) = {zeros(u,v)};
                R_mm=Q;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Convert mth tensor back into N_Right labelling
                temp = cell(mps_var.d,NR);
                for alpha = 1:mps_var.d
                    temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
                end
                temp(cellfun(@isempty,temp)) = {zeros(u,v)};
                
                mps_var.data{m}=temp;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if m ~= 1
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Convert (m-1)th Tensor to N_Left Labelling
                    [~,NR_1]=size(mps_var.data{m-1});
                    NL_1=min(mps_var.N+1,(mps_var.d-1)*(m-2)+1);
                    
                    temp_T = cell(mps_var.d,NR_1);
                    temp_IT = cell(mps_var.d,NL_1);
                    expect = max(NR_1,NL_1):-1:1;
                    for k = 1:mps_var.d
                        track=1;
                        for k_2 = (mps_var.N+1-k+1):-1:1
                            if k_2>mps_var.N+1-NR_1-k+1 && k_2<=NL_1
                                temp_T(k,k_2) = {expect(1,max(NR_1,NL_1)-track+1)};
                            end
                            track=track+1;
                        end
                        for k_2 = max(1,mps_var.N+1-NR_1-k+2):min(NL_1,(mps_var.N+2-k))
                            temp_IT(k,k_2) = {k_2};
                        end
                    end
                    
                    temp = mps_var.data{m-1};
                    
                    R_mm_1 = cell(mps_var.d,NL_1);
                    for alpha = 1:mps_var.d
                        R_mm_1(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
                    end
                    R_mm_1(cellfun(@isempty,R_mm_1)) = {zeros(size(mps_var.data{m-1}{1}))};
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Apply Singular values to next tensor
                    temp_E = cell(mps_var.d,NL_1);
                    for k = 1:mps_var.d
                        track=1;
                        for k_2 = (mps_var.N+1-k+1):-1:1
                            if k_2 <= NL_1
                                temp_E{k,k_2} = S{1,mps_var.N+1-track+1};
                            end
                            track=track+1;
                        end
                    end
                    temp_E(cellfun(@isempty,temp_E)) = {zeros(u,u)};
                    
                    R_mm_1 = cellfun(@(y,z) y*z, R_mm_1,temp_E, 'UniformOutput', false);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Convert (m-1)th tensor back into N_Right labelling
                    temp = cell(mps_var.d,NR_1);
                    for alpha = 1:mps_var.d
                        temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm_1(alpha,cell2mat(temp_IT(alpha,:)));
                    end
                    temp(cellfun(@isempty,temp)) = {zeros(size(mps_var.data{m-1}{1}))};
                    
                    mps_var.data{m-1}=temp;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            else
                error('Must choose canonicalisation');
            end
                   
        end
        
        function N = Full_Norm(mps_var)
            % Find Overlap of MPS (Normalisation)
            
            M = size(mps_var.data,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Overlap
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            expect=cell(1,mps_var.N+1);
            expect(cellfun(@isempty,expect)) = {1};
            for n = 1:M
                [a,g] = size(mps_var.data{n});
                [alpha,gamma] = size(mps_var.data{n}{end,1});
                
                temp_E = cell(a,g);
                for k = 1:a
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2 <= g
                            temp_E{k,k_2} = expect{1,mps_var.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                
                temp = cellfun(@(x,y,z) x'*y*z, mps_var.data{n},temp_E,mps_var.data{n}, 'UniformOutput', false);
                
                expect=cell(1,g);
                expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                
                for i=1:a
                    expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                end
                if n<M
                    if size(expect,2)<mps_var.d
                        expect(:,size(expect,2)+1:mps_var.d)={zeros(gamma,gamma)};
                    end
                end
            end
            N=cell2mat(expect);
        end
        
        function mps_var = Increase_bond_dim(mps_var,Bond_Dim)
            % Increase Maximum bond dimension. Adds Zeros to the MPS entries.
            
             M = size(mps_var.data,2);
             
             for m = 1:M
                 [~,g] = size(mps_var.data{m});
                 Nl = min((mps_var.d-1)*(m-1)+1,mps_var.N+1);
                 [alpha,beta] = size(mps_var.data{m}{1,1});
                 
                 if m==1
                     BL=1;
                 else
                     BL = max(alpha,Bond_Dim);
                 end
                 if m==M
                     BR=1;
                 else
                     BR = max(beta,Bond_Dim);
                 end
                 
                 for count = 1:mps_var.d
                     for Num_r = 1:g
                         temp = mps_var.data{m}{count,Num_r};
                         temp(:,beta+1:BR) = zeros(alpha,length(beta+1:BR));
                         temp(alpha+1:BL,:) = zeros(length(alpha+1:BL),BR);
                         mps_var.data{m}{count,Num_r}=temp;
                     end
                 end
                 
             end
            
        end        
        
        function mps_var = set_rand_product_state_for_Variational_Algorithms(mps_var,Bond_Dimension)
            % Set initial random product state.
            %
            % N == Number of particles.
            % Nmax == Maximum number of particles allowed on a single site.
            
            M = size(mps_var.data,2);
            
            mps_var.Dmax=Bond_Dimension;
            
            if M*mps_var.N_max < mps_var.N
                error('Too many particles!');
            end
            
            % Generate N random integers between 1 and M.
            % These correspond to an intial random location for each particle.
            if isempty(mps_var.P_Positions)
                R = randi(M,1,mps_var.N);
                mps_var.P_Positions=R;
            else
                R = mps_var.P_Positions;
            end
            R = sort(R);
            
            A = zeros(M,1);
            C=0;
            check=0;
            for m = 1:M
                A(m) = C;
                C=0;
                for n = 1:mps_var.N
                    if R(n) == m
                        A(m) = A(m)+1; % Total number of particles on site m
                    end
                end
                % If the intial random locations give more particles than are
                % allowed (dictated by N_max) then we carry the extra particles
                % over to the next site.
                if A(m)>mps_var.N_max
                    C=A(m)-mps_var.N_max;
                    A(m)=mps_var.N_max;
                    if m==M
                        % We have reached the end of the lattice and we cannot
                        % carry over the extra particles.
                        % We have C extra particles that do not have a lattice site.
                        % We need to apply an extra algorithm (below) to sort them into
                        % sites that have spaces.
                        check=1;
                    end
                end
                
                u=min(mps_var.d^(m-1),min(mps_var.d^(M-m+1),Bond_Dimension));
                v=min(mps_var.d^m,min(mps_var.d^(M-m),Bond_Dimension));
                
                Mat = rand(Bond_Dimension)+1i*rand(Bond_Dimension);
                [Mat,~]=qr(Mat);
                
                Mat=Mat(1:u,1:v);
                
                if m == 1
                    B = cell(mps_var.d,mps_var.N+1);
                    B{A(m)+1,mps_var.N-A(m)+1}=Mat;
                    B(cellfun(@isempty,B)) = {zeros(u,v)};
                elseif m==M
                    B=cell(mps_var.d,1);
                    B{A(m)+1,1} = Mat;
                    B(cellfun(@isempty,B)) = {zeros(u,v)};
                else
                    xr = min(mps_var.N+1,mps_var.N_max*(M-m)+1);
                    B = cell(mps_var.d,xr);
                    B{A(m)+1,mps_var.N+1-sum(A(1:m))}=Mat;
                    B(cellfun(@isempty,B)) = {zeros(u,v)};
                end
                mps_var.data{m} = B;
                
            end
            
            % Apply extra algorithm to sort extra particles into free lattice sites
            if check==1
                m=1;
                % We alter the values of A(m) to ensure the all particles have a
                % lattice.
                while C ~=0
                    D=A(m)+C;
                    if D > mps_var.N_max
                        C = D-mps_var.N_max;
                        D=mps_var.N_max;
                    else
                        C=0;
                    end
                    A(m)=D;
                    m=m+1;
                end
                
                % With the new A(m) we construct the MPS
                for m = 1:M
                    u=min(mps_var.d^(m-1),min(mps_var.d^(M-m+1),Bond_Dimension));
                    v=min(mps_var.d^m,min(mps_var.d^(M-m),Bond_Dimension));
                    
                    Mat = rand(Bond_Dimension)+1i*rand(Bond_Dimension);
                    [Mat,~]=qr(Mat);
                    
                    Mat=Mat(1:u,1:v);
                    
                    if m == 1
                        B = cell(mps_var.d,mps_var.N+1);
                        B{A(m)+1,mps_var.N-A(m)+1}=Mat;
                        B(cellfun(@isempty,B)) = {zeros(u,v)};
                    elseif m==M
                        B=cell(mps_var.d,1);
                        B{A(m)+1,1} = Mat;
                        B(cellfun(@isempty,B)) = {zeros(u,v)};
                    else
                        xr = min(mps_var.N+1,mps_var.N_max*(M-m)+1);
                        B = cell(mps_var.d,xr);
                        B{A(m)+1,mps_var.N+1-sum(A(1:m))}=Mat;
                        B(cellfun(@isempty,B)) = {zeros(u,v)};
                    end
                    mps_var.data{m} = B;
                end
            end
        end
      
        function mps_var = DMRG(mps_var,H,Sweeps)
            % Apply DMRG procedure
            
            mps_var=mps_var.Canonicalise_1s('R-L','true');
                        
            M = size(mps_var.data,2);
            R = cell(1,M);
            L = cell(1,M);
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Initial Left Effective Hamiltonian
            %%%%%%%%%%%%%%%%%%%%%%%
            INIT_temp = cell(1,mps_var.N+1);
            INIT_temp(cellfun(@isempty,INIT_temp)) = {1};
            INIT{1} = INIT_temp;
            for m = M:-1:2
                if m==M
                    L{m} = mps_var.Add_Left_Eff_Ham(H,INIT,m);
                else
                    L{m} = mps_var.Add_Left_Eff_Ham(H,L{m+1},m);
                end
            end
            
            R_INIT_temp = cell(1,mps_var.N+1);
            R_INIT_temp(cellfun(@isempty,R_INIT_temp)) = {1};
            R_INIT{1}=R_INIT_temp;
                        
            %%%%%%%%%%%%%%%%%%%%%%%
            % Begin Algorithm
            %%%%%%%%%%%%%%%%%%%%%%%
            for s = 1:Sweeps
                tic
                %%%%%%%%%%%%%%%%%%%%%%%
                % Sweep Left to Right
                %%%%%%%%%%%%%%%%%%%%%%%
                for m = 1:(M-1)
                    if m == 1
                        H_Eff=mps_var.Create_Eff_Ham(H,R_INIT,L{m+1},m);
                    else
                        H_Eff=mps_var.Create_Eff_Ham(H,R{m-1},L{m+1},m);
                    end
                    
                    mps_var=mps_var.Variationally_Optimise_Tensor(H_Eff,m,'L-R');
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    % Update Right Effective Hamiltonian
                    %%%%%%%%%%%%%%%%%%%%%%%
                    if m ==1
                        R{m} = mps_var.Add_Right_Eff_Ham(H,R_INIT,m);
                    else
                        R{m} = mps_var.Add_Right_Eff_Ham(H,R{m-1},m);
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%
                % Sweep Right to Left
                %%%%%%%%%%%%%%%%%%%%%%%
                for m = M:-1:2
                    if m == M
                        H_Eff=mps_var.Create_Eff_Ham(H,R{m-1},INIT,m);
                    else
                        H_Eff=mps_var.Create_Eff_Ham(H,R{m-1},L{m+1},m);
                    end
                    
                    mps_var=mps_var.Variationally_Optimise_Tensor(H_Eff,m,'R-L');
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    % Update Left Effective Hamiltonian
                    %%%%%%%%%%%%%%%%%%%%%%%
                    if m ==M
                        L{m} = mps_var.Add_Left_Eff_Ham(H,INIT,m);
                    elseif m ~= 1
                        L{m} = mps_var.Add_Left_Eff_Ham(H,L{m+1},m);
                    end
                    
                end
                                                
                Check_Norm = abs(mps_var.Full_Norm);
                P_Num=0;
                for m =1:M
                    P_Num=P_Num+mps_var.Site_Site_Particle_Corr(m,m);
                end
                
                Energy=mps_var.Get_Energy(H);
                Var=mps_var.Get_Var(H);
                
                disp(['Sweep=',num2str(s),'/',num2str(Sweeps),' -- CPU Time=',num2str(toc)]);
                disp(['Norm=',num2str(Check_Norm),' -- P_Num=',num2str(P_Num)]);
                disp(['Energy=',num2str(Energy),' -- Var=',num2str(Var)]);
                                
            end
            
        end
        
        function mps_var = Variationally_Optimise_Tensor(mps_var,H_Eff,m,Sweep_Direction)
            
            % H_Eff == (u NR' d' v),(u NR d v)
            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find Smallest eigenvector
            d_cell=cell(1,NR);
            Vec=cell(1,NR);
            List=zeros(mps_var.d,NR);
            for NR_iter = 1:NR
                d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                Vec(d_cell{1,NR_iter},NR_iter)=mps_var.data{m}(d_cell{1,NR_iter},NR_iter);
                List(d_cell{1,NR_iter},NR_iter)=ones(length(d_cell{1,NR_iter}),1);
            end
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            % Vec == d,NR
            Vec = reshape(Vec,mps_var.d*NR,1); % (NR d),1
            Vec = cell2mat(Vec(Ind,1)); % (NR d u),v
                
            Vec = reshape(Vec,u,nnz(Len),v); %u,(NR d),v
            Vec = permute(Vec,[3,2,1]); %v,(NR d),u
            Vec = reshape(Vec,u*nnz(Len)*v,1); %(u NR d v),1
            
            opts.v=Vec;
            try
                [Vec,~]=eigs(H_Eff,1,'smallestreal',opts); %(u NR d v),1
            catch
                [Vec,~]=eigs(eye(size(H_Eff,1))+H_Eff,1,'smallestreal',opts); %(u NR d v),1
            end
            
            Vec = reshape(Vec,v,nnz(Len),u); % v,(NR d),u
            Vec = permute(Vec,[3,2,1]); % u,(NR d),v
            Vec = reshape(Vec,u*nnz(Len),v); % (NR d u),v
            
            row_vec=u*ones(nnz(Len),1);
            Vec = mat2cell(Vec,row_vec,v); % (NR d),1
            
            New_Vec=cell(mps_var.d*NR,1);
            New_Vec(Ind,1)=Vec;
            New_Vec=reshape(New_Vec,mps_var.d,NR);
            New_Vec(cellfun(@isempty,New_Vec)) = {zeros(u,v)};
            mps_var.data{m}=New_Vec;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do SVD
            mps_var=mps_var.Do_SVD_Single_Tensor(m,Sweep_Direction,'false');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function mps_var = DMRG_2s(mps_var,H,Sweeps)
            % Apply 2 Site DMRG procedure
            
            mps_var=mps_var.Canonicalise_1s('R-L','true');
                        
            M = size(mps_var.data,2);
            R = cell(1,M);
            L = cell(1,M);
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Initial Left Effective Hamiltonian
            %%%%%%%%%%%%%%%%%%%%%%%
            INIT_temp = cell(1,mps_var.N+1);
            INIT_temp(cellfun(@isempty,INIT_temp)) = {1};
            INIT{1} = INIT_temp;
            for m = M:-1:3
                if m==M
                    L{m} = mps_var.Add_Left_Eff_Ham(H,INIT,m);
                else
                    L{m} = mps_var.Add_Left_Eff_Ham(H,L{m+1},m);
                end
            end
            
            R_INIT_temp = cell(1,mps_var.N+1);
            R_INIT_temp(cellfun(@isempty,R_INIT_temp)) = {1};
            R_INIT{1}=R_INIT_temp;
                        
            %%%%%%%%%%%%%%%%%%%%%%%
            % Begin Algorithm
            %%%%%%%%%%%%%%%%%%%%%%%
            for s = 1:Sweeps
                tic
                %%%%%%%%%%%%%%%%%%%%%%%
                % Sweep Left to Right
                %%%%%%%%%%%%%%%%%%%%%%%
                Total_error=0;
                for m = 1:(M-2)
                    if m == 1
                        H_Eff=mps_var.Create_Eff_Ham_2s(H,R_INIT,L{m+2},m);
                    else
                        H_Eff=mps_var.Create_Eff_Ham_2s(H,R{m-1},L{m+2},m);
                    end
                    
                    [mps_var,Error]=mps_var.Variationally_Optimise_Tensor_2s(H_Eff,m,'L-R');
                    Total_error=Total_error+Error;
                    %%%%%%%%%%%%%%%%%%%%%%%
                    % Update Right Effective Hamiltonian
                    %%%%%%%%%%%%%%%%%%%%%%%
                    if m ==1
                        R{m} = mps_var.Add_Right_Eff_Ham(H,R_INIT,m);
                    else
                        R{m} = mps_var.Add_Right_Eff_Ham(H,R{m-1},m);
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%
                % Sweep Right to Left
                %%%%%%%%%%%%%%%%%%%%%%%
                for m = (M-1):-1:2
                    if m == (M-1)
                        H_Eff=mps_var.Create_Eff_Ham_2s(H,R{m-1},INIT,m);
                    elseif m==1
                        H_Eff=mps_var.Create_Eff_Ham_2s(H,R_INIT,L{m+2},m);
                    else
                        H_Eff=mps_var.Create_Eff_Ham_2s(H,R{m-1},L{m+2},m);
                    end
                    
                    [mps_var,Error]=mps_var.Variationally_Optimise_Tensor_2s(H_Eff,m,'R-L');
                    Total_error=Total_error+Error;
                    %%%%%%%%%%%%%%%%%%%%%%%
                    % Update Left Effective Hamiltonian
                    %%%%%%%%%%%%%%%%%%%%%%%
                    if m ==(M-1)
                        L{m+1} = mps_var.Add_Left_Eff_Ham(H,INIT,m+1);
                    elseif m ~= 1
                        L{m+1} = mps_var.Add_Left_Eff_Ham(H,L{m+2},m+1);
                    end
                    
                end
                                                
                Check_Norm = abs(mps_var.Full_Norm);
                P_Num=0;
                for m =1:M
                    P_Num=P_Num+mps_var.Site_Site_Particle_Corr(m,m);
                end
                
                Energy=mps_var.Get_Energy(H);
                Var=mps_var.Get_Var(H);
                
                disp(['Sweep=',num2str(s),'/',num2str(Sweeps),' -- CPU Time=',num2str(toc)]);
                disp(['Turncation Error=',num2str(Total_error)]);
                disp(['Norm=',num2str(Check_Norm),' -- P_Num=',num2str(P_Num)]);
                disp(['Energy=',num2str(Energy),' -- Var=',num2str(Var)]);
                                
            end
            
        end
        
        function [mps_var,Total_error] = Variationally_Optimise_Tensor_2s(mps_var,H_Eff,m,Sweep_Direction)
            
            [~,NR_2] = size(mps_var.data{m+1});
            
            Nl_2 = min(mps_var.N+1,(mps_var.d-1)*m+1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the tensor for the two-sites so we can apply the operator.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} to N_Right labelling
            %%%%%%%%%%%%%%%%%%
            temp_T = cell(mps_var.d,NR_2);
            temp_IT = cell(mps_var.d,Nl_2);
            expect = max(NR_2,Nl_2):-1:1;
            for k = 1:mps_var.d
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2>mps_var.N+1-NR_2-k+1 && k_2<=Nl_2
                        temp_T(k,k_2) = {expect(1,max(NR_2,Nl_2)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_var.N+1-NR_2-k+2):min(Nl_2,(mps_var.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_var.data{m+1};
            
            R_mm = cell(mps_var.d,Nl_2);
            for alpha = 1:mps_var.d
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            [nrows, ncols] = cellfun(@size, temp);
            beta = max(max(nrows));
            v = max(max(ncols));
            R_mm(cellfun(@isempty,R_mm)) = {zeros(beta,v)};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Collect entries of (L_m X R_mm)
            % Construct Two-Site Tensor
            NR1=size(mps_var.data{m},2);
            NL2=min(mps_var.N+1,(mps_var.d-1)*(m)+1);
            NR1_vec=NR1:-1:1;
            NL2_vec=1:NL2;
            endR=NR1+1-(mps_var.N+1-NL2+1);
            startL=mps_var.N+1-NR1+1;
            startR=endR-length(NL2_vec)+startL;
            
            Col_vec=cell(1,min(length(NR1_vec),length(NL2_vec)));
            track=1;
            for n=startR:min(length(NR1_vec),length(NL2_vec))
                Col_vec{track}=[NR1_vec(n),NL2_vec(track-1+startL)];
                track=track+1;
            end
            
            NR=size(mps_var.data{m},2);
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            NR_min=mps_var.N+1-NL+1;
            
            d_cell=cell(1,NR);
            for NR_iter = 1:NR
               d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
            end
            
            dd_vec=(1:mps_var.d)'-1;
            NL_vec=(mps_var.N:-1:0);
            NL_vec=NL_vec(ones(mps_var.d,1),1:NR);
            NR_vec_temp=(0:mps_var.N)';
            NR_vec=zeros(mps_var.d,NR);
            for NR_iter=1:NR
                NR_vec(1:size(NR_iter:min(length(NR_vec_temp),(NR_iter+mps_var.d-1)),2),NR_iter)=NR_vec_temp(NR_iter:min(end,(NR_iter+mps_var.d-1)));
            end

            ind=find(NR_vec<(NR_min-1));
            NR_vec(ind)=0;
            NL_vec(ind)=0;
            
            List=-dd_vec(:,ones(1,size(NL_vec,2))) + NR_vec + NL_vec;
            List(List~=mps_var.N)=0;
            List(List==mps_var.N)=1;
            
            NR_vec=(NR_vec+1).*List;
            
            NR2=size(mps_var.data{m+1},2);
            NL2=min(mps_var.N+1,(mps_var.d-1)*(m)+1);
            NR2_min=mps_var.N+1-NL2+1;
            
            d_cell2=cell(1,NR2);
            for NR2_iter = 1:NR2
               d_cell2{1,NR2_iter}=max(1,mps_var.N+3-NR2_iter-NL2):min(mps_var.d,mps_var.N+2-NR2_iter);
            end
            
            NL2_vec=(mps_var.N:-1:0);
            NL2_vec=NL2_vec(ones(mps_var.d,1),1:NR2);
            NR2_vec_temp=(0:mps_var.N)';
            NR2_vec=zeros(mps_var.d,NR2);
            for NR2_iter=1:NR2
                NR2_vec(1:size(NR2_iter:min(length(NR_vec_temp),(NR2_iter+mps_var.d-1)),2),NR2_iter)=NR2_vec_temp(NR2_iter:min(end,(NR2_iter+mps_var.d-1)));
            end

            ind=find(NR2_vec<(NR2_min-1));
            NR2_vec(ind)=0;
            NL2_vec(ind)=0;
            
            List2=-dd_vec(:,ones(1,size(NL2_vec,2))) + NR2_vec + NL2_vec;
            List2(List2~=mps_var.N)=0;
            List2(List2==mps_var.N)=1;
            
            NL2_vec=(NL2_vec+1).*List2;
                       
            shift=NL2_vec-dd_vec(:,ones(1,size(NL2_vec,2)));
            shift(shift<0)=0;
            
            List2L=zeros(mps_var.d,NL2);
            NL_Left_Vec=zeros(mps_var.d,NL2);
            for iter=1:mps_var.d
                shift_L=shift(iter,:);
                shift_L(shift_L==0)=[];
                start_shift=find(List2(iter,:)~=0);
                List2L(iter,shift_L)=List2(iter,start_shift);
                NL_Left_Vec(iter,shift_L)=NL2_vec(iter,start_shift);
            end
            
            dd_List=zeros(mps_var.d^2,size(Col_vec,2));
            List1=List(:,end:-1:1);
            dd_cell=cell(1,size(Col_vec,2));
            for n = 1:size(Col_vec,2)
                dd_List(:,n)=kron(List1(:,startR+n-1),List2L(:,startL+n-1));
                dd_cell{n}=find(dd_List(:,n)~=0);
            end
            dd_map=[kron(dd_vec,ones(mps_var.d,1)),kron(ones(mps_var.d,1),dd_vec)]+1;
            
            % dd_List == (d1 d2),(NR1 NL2)
            Len_dd=reshape(dd_List,mps_var.d^2*size(Col_vec,2),1); %(NR1 NL2 d1 d2)
            Ind=find(Len_dd~=0);
            
            Theta_Mat=cell(size(dd_List));
            
            for N_iter=1:size(Col_vec,2)
                for d_iter=dd_cell{N_iter}'
                    if ~isempty(Theta_Mat{d_iter,N_iter})
                        Theta_Mat{d_iter,N_iter}=Theta_Mat{d_iter,N_iter}+mps_var.data{m}{dd_map(d_iter,1),Col_vec{1,N_iter}(1)}*R_mm{dd_map(d_iter,2),Col_vec{1,N_iter}(2)};
                    else
                        Theta_Mat{d_iter,N_iter}=mps_var.data{m}{dd_map(d_iter,1),Col_vec{1,N_iter}(1)}*R_mm{dd_map(d_iter,2),Col_vec{1,N_iter}(2)};
                    end
                end
            end
            % Theta_Mat == (d1 d2),(NR1 NL2) -- {Theta_Mat} == u,v
            
            [u,~] = size(mps_var.data{m}{1});
            [~,v] = size(mps_var.data{m+1}{1,1});
            
            Theta_Mat = reshape(Theta_Mat,mps_var.d^2*size(Col_vec,2),1); %(NR1 NL2 d1 d2),1
            Theta_Mat = Theta_Mat(Ind,1);
            Theta_Mat = cell2mat(Theta_Mat); % (NR1 NL2 d1 d2 u),v
            Theta_Mat = reshape(Theta_Mat,u,nnz(Len_dd),v); % u,(NR1 NL2 d1 d2),v
            Theta_Mat = permute(Theta_Mat,[3 2 1]); % v,(NR1 NL2 d1 d2),u
            Theta_Mat = reshape(Theta_Mat,u*nnz(Len_dd)*v,1);% (u NR1 NL2 d1 d2 v),1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Theta_Mat'*Theta_Mat
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find smallest eigenvector
            % H_Eff == (u' NR1' NL2' d1' d2' v'),(u NR1 NL2 d1 d2 v)

            opts.v=Theta_Mat;
            try
                [Theta_Mat,~]=eigs(H_Eff,1,'smallestreal',opts); %(u NR d v),1
            catch
                [Theta_Mat,~]=eigs(eye(size(H_Eff,1))+H_Eff,1,'smallestreal',opts); %(u NR d v),1
            end
            
%             [V,E]=eig(H_Eff)

            Theta_Mat = reshape(Theta_Mat,v,nnz(Len_dd),u); %v,(NR1 NL2 d1 d2),u
            Theta_Mat = permute(Theta_Mat,[3 2 1]); %  u,(NR1 NL2 d1 d2),v
            Theta_Mat = reshape(Theta_Mat,u*nnz(Len_dd),v); % (NR1 NL2 d1 d2 u),v

            row_vec=u*ones(nnz(Len_dd),1);
            Theta_Mat = mat2cell(Theta_Mat,row_vec,v); % (NR1 NL2 d1 d2),1
            
            Theta_Mat_New=cell(mps_var.d^2*size(Col_vec,2),1);
            Theta_Mat_New(Ind,1)=Theta_Mat; % (NR1 NL2 d1 d2),1
            
            Theta_Mat_New=reshape(Theta_Mat_New,mps_var.d^2,size(Col_vec,2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do SVD's
            R_mm = cell(mps_var.d,NL2);
            L_m = cell(mps_var.d,NR);

            r_min=1;
            Norm=0;
            Total_error=0;
            for count2 = 1:size(Theta_Mat_New,2)
                Theta = cell2mat(Theta_Mat_New(:,count2));
                
                d_1=nnz(List(:,Col_vec{count2}(1)));
                d_2=nnz(List2L(:,Col_vec{count2}(2)));
                
                % Theta == (d_1 d_2 alpha),(gamma)
                Theta = reshape(Theta,u,d_2,d_1,v); %(a),(l+1),(l),(g)
                Theta = permute(Theta,[1 3 4 2]); %(a),(l),(g),(l+1)
                Theta = reshape(Theta,u*d_1,v*d_2); %(l a),(l+1 g)
                Theta = reshape(Theta,u*d_1,v*d_2); %(l a),(l+1 g)
                [A,S,V]=svd(Theta);
                
                S((abs(S))/S(1,1)<mps_var.zero_thres)=0;
                
                r = nnz(S);
                if r>r_min
                    r_min=r;
                    if r_min>mps_var.Dmax
                        r_min=mps_var.Dmax;
                    end
                end
                if r==0
                    r=1;
                    S = S(1:r,1:r);
                    Error=0;
                elseif r>mps_var.Dmax
                    r=mps_var.Dmax;
                    if any(size(S(r+1:end,r+1:end))==1)
                        Error=sum(S(r+1:end,r+1:end).^2);
                    else
                        Error=sum(diag(S(r+1:end,r+1:end)).^2);
                    end
                    S = S(1:r,1:r);
                else
                    Error=0;
                    S = S(1:r,1:r);
                end
%                 diag(S)'*diag(S)
                Norm = Norm+diag(S)'*diag(S);
                if strcmp(Sweep_Direction,'L-R')
                    A = A(:,1:r);
                    V = S*(V(:,1:r)');
                else
                    A = A(:,1:r)*S;
                    V = V(:,1:r)';
                end
                
                Total_error=Total_error+Error;
                
                % A == (l a),(r)
                A = reshape(A,u,d_1,r); %(a),(l),(r)
                A = permute(A,[1 3 2]); %(a),(r),(l)
                track=1;
                for count = unique(dd_map(dd_cell{count2},1))'
                    L_m(count,Col_vec{count2}(1)) = {A(:,:,track)};
                    track=track+1;
                end
                
                % V == (r),(l+1 g)
                V = reshape(V,r,v,d_2); %(r),(g),(l+1)
                track=1;
                for count = unique(dd_map(dd_cell{count2},2))'
                    R_mm(count,Col_vec{count2}(2)) = {V(:,:,track)};
                    track=track+1;
                end
            end
            R_mm(cellfun(@isempty,R_mm)) = {zeros(r_min,v)};
            L_m(cellfun(@isempty,L_m)) = {zeros(u,r_min)};
                        
            [nrows, ~] = cellfun(@size, L_m);
            u = max(max(nrows));
            [~, ncol] = cellfun(@size, R_mm);
            v = max(max(ncol));
            
            L_m = cellfun(@(x) cell_resize_func(x,[u,r_min]),L_m,'UniformOutput', false);
            R_mm = cellfun(@(x) cell_resize_func(x,[r_min,v]),R_mm,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update First Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mps_var.data{m}=L_m;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update Second Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} back to N_Left labelling
            %%%%%%%%%%%%%%%%%%
            
            temp = cell(mps_var.d,NR_2);
            for alpha = 1:mps_var.d
                temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(r_min,v)};
            
            mps_var.data{m+1}=temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            function Out = cell_resize_func(X,r)
                % Function used to resize MPS
                % tensor.
                
                [a,b2]= size(X);
                Out=X;
                if a < r(1)
                    Out(a+1:r(1),:) = zeros(r(1)-a,size(Out,2));
                end
                if b2 < r(2)
                    Out(:,b2+1:r(2)) = zeros(size(Out,1),r(2)-b2);
                end
                
            end
        end

        function mps_var = TDVP(mps_var,H,dt,T)
            % Apply TDVP procedure
            
            % Ensure begining tensor is right canonical and normalised.
            mps_var=mps_var.Canonicalise_1s('R-L','true');
            
            time_steps=floor(T/abs(dt));
            
            M = size(mps_var.data,2);
            R=cell(1,M);
            L = cell(1,M);
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Initial Left Effective Hamiltonian
            %%%%%%%%%%%%%%%%%%%%%%%
            INIT_temp = cell(1,mps_var.N+1);
            INIT_temp(cellfun(@isempty,INIT_temp)) = {1};
            INIT{1}=INIT_temp;
            for m = M:-1:2
                if m==M
                    L{m} = mps_var.Add_Left_Eff_Ham(H,INIT,m);
                else
                    L{m} = mps_var.Add_Left_Eff_Ham(H,L{m+1},m);
                end
            end
            
            R_INIT_temp = cell(1,mps_var.N+1);
            R_INIT_temp(cellfun(@isempty,R_INIT_temp)) = {1};
            R_INIT{1}=R_INIT_temp;
                        
            %%%%%%%%%%%%%%%%%%%%%%%
            % Begin Algorithm
            %%%%%%%%%%%%%%%%%%%%%%%
            for t = 1:time_steps
                
                %%%%%%%%%%%%%%%%%%%%%%%
                % Sweep Left to Right
                %%%%%%%%%%%%%%%%%%%%%%%
                for m = 1:(M-1)
                    
                    if m == 1
                        H_Eff=mps_var.Create_Eff_Ham(H,R_INIT,L{m+1},m);
                    else
                        H_Eff=mps_var.Create_Eff_Ham(H,R{m-1},L{m+1},m);
                    end
                    
                    [mps_var,S]=mps_var.Apply_H_Eff_Right(H_Eff,dt/2,m,'false');

                    %%%%%%%%%%%%%%%%%%%%%%%
                    % Update Right Effective Hamiltonian
                    %%%%%%%%%%%%%%%%%%%%%%%
                    if m ==1
                        R{m} = mps_var.Add_Right_Eff_Ham(H,R_INIT,m);
                    else
                        R{m} = mps_var.Add_Right_Eff_Ham(H,R{m-1},m);
                    end

                    K_Eff=mps_var.Create_Eff_K_Right(H_Eff,m);
%                     K_Eff=mps_var.Create_Eff_K_Right(R{m},L{m+1});
                    mps_var=mps_var.Apply_K_Eff_Right(S,K_Eff,dt/2,m,'false');
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%
                % Sweep Right to Left
                %%%%%%%%%%%%%%%%%%%%%%%
                for m = M:-1:1
                    
                    if m == 1
                        H_Eff=mps_var.Create_Eff_Ham(H,R_INIT,L{m+1},m);
                    elseif m == M
                        H_Eff=mps_var.Create_Eff_Ham(H,R{m-1},INIT,m);
                    else
                        H_Eff=mps_var.Create_Eff_Ham(H,R{m-1},L{m+1},m);
                    end
                    
                    if m == M
                        [mps_var,S]=mps_var.Apply_H_Eff_Left(H_Eff,dt,m,'false');
                    else
                        [mps_var,S]=mps_var.Apply_H_Eff_Left(H_Eff,dt/2,m,'false');
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    % Update Left Effective Hamiltonian
                    %%%%%%%%%%%%%%%%%%%%%%%
                    if m ==M
                        L{m} = mps_var.Add_Left_Eff_Ham(H,INIT,m);
                    elseif m ~= 1
                        L{m} = mps_var.Add_Left_Eff_Ham(H,L{m+1},m);
                    end
                    
                    if m ~=1              
                        K_Eff=mps_var.Create_Eff_K_Left(H_Eff,m);
%                         K_Eff=mps_var.Create_Eff_K_Left(R,{m-1},L{m});
                        mps_var=mps_var.Apply_K_Eff_Left(S,K_Eff,dt/2,m,'false');
                    end
                end
                        
                mps_var=mps_var.Canonicalise_1s('R-L','true');
                Check_Norm = mps_var.Full_Norm;
                P_Num=0;
                for m =1:M
                    P_Num=P_Num+mps_var.Site_Site_Particle_Corr(m,m);
                end
                
                Energy=mps_var.Get_Energy(H);
                Var=mps_var.Get_Var(H);
                                
                disp(['Time=',num2str(t*dt)]);
                disp(['Norm=',num2str(Check_Norm),' -- P_Num=',num2str(P_Num)]);
                disp(['Energy=',num2str(Energy),' -- Var=',num2str(Var)]);
            end
        end
                           
        function mps_var = TDVP_2s(mps_var,H,dt,T,Trotter_Splitting)
            % Apply 2-site TDVP procedure
            
            % Ensure begining tensor is right canonical and normalised.
            mps_var=mps_var.Canonicalise_1s('R-L','true');
            
            time_steps=floor(T/abs(dt));
            
            M = size(mps_var.data,2);
            R=cell(1,M);
            L = cell(1,M);
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Initial Left Effective Hamiltonian
            %%%%%%%%%%%%%%%%%%%%%%%
            INIT_temp = cell(1,mps_var.N+1);
            INIT_temp(cellfun(@isempty,INIT_temp)) = {1};
            INIT{1}=INIT_temp;
            for m = M:-1:3
                if m==M
                    L{m} = mps_var.Add_Left_Eff_Ham(H,INIT,m);
                else
                    L{m} = mps_var.Add_Left_Eff_Ham(H,L{m+1},m);
                end
            end
            
            R_INIT_temp = cell(1,mps_var.N+1);
            R_INIT_temp(cellfun(@isempty,R_INIT_temp)) = {1};
            R_INIT{1}=R_INIT_temp;
                        
            %%%%%%%%%%%%%%%%%%%%%%%
            % Begin Algorithm
            %%%%%%%%%%%%%%%%%%%%%%%
            tau=dt;
            if strcmp(Trotter_Splitting,'true')
                dt_vec=tau*[1/12,1/12,1/12,-1/6,1/12,0,1/12,0,1/12,0,1/12,1/12,1/12,1/12,0,1/12,0,1/12,0,1/12,-1/6,1/12,1/12,1/12];
            else
                dt_vec=tau*[1/2,1/2];
            end
            Sweep_Direction='L-R';
            tic
            for t = 1:time_steps
                for dt=dt_vec
                    if strcmp(Sweep_Direction,'L-R')
                        %%%%%%%%%%%%%%%%%%%%%%%
                        % Sweep Left to Right
                        %%%%%%%%%%%%%%%%%%%%%%%
                        Sweep_Direction='R-L';
                        Total_error=0;
                        for m = 1:(M-1)
                            
                            if m == 1
                                H_Eff=mps_var.Create_Eff_Ham_2s(H,R_INIT,L{m+2},m);
                            elseif m==(M-1)
                                H_Eff=mps_var.Create_Eff_Ham_2s(H,R{m-1},INIT,m);
                            else
                                H_Eff=mps_var.Create_Eff_Ham_2s(H,R{m-1},L{m+2},m);
                            end
                            
                            [mps_var,Error]=mps_var.Apply_H_Eff_2s(H_Eff,dt,m,'L-R');
                            Total_error=Total_error+Error;
                            %%%%%%%%%%%%%%%%%%%%%%%
                            % Update Right Effective Hamiltonian
                            %%%%%%%%%%%%%%%%%%%%%%%
                            if m ==1
                                R{m} = mps_var.Add_Right_Eff_Ham(H,R_INIT,m);
                            else
                                R{m} = mps_var.Add_Right_Eff_Ham(H,R{m-1},m);
                            end
                            
                            if m~=(M-1)
                                H_Eff=mps_var.Create_Eff_Ham(H,R{m},L{m+2},m+1);
                                [mps_var]=mps_var.Apply_1S_H_Right(H_Eff,-dt,m+1);
                            end
                        end
                        
                    else
                        Sweep_Direction='L-R';
                        %%%%%%%%%%%%%%%%%%%%%%%
                        % Sweep Right to Left
                        %%%%%%%%%%%%%%%%%%%%%%%
                        for m = (M-1):-1:1
                            
                            if m == 1
                                H_Eff=mps_var.Create_Eff_Ham_2s(H,R_INIT,L{m+2},m);
                            elseif m == (M-1)
                                H_Eff=mps_var.Create_Eff_Ham_2s(H,R{m-1},INIT,m);
                            else
                                H_Eff=mps_var.Create_Eff_Ham_2s(H,R{m-1},L{m+2},m);
                            end
                            
                            [mps_var,Error]=mps_var.Apply_H_Eff_2s(H_Eff,dt,m,'R-L');
                            Total_error=Total_error+Error;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%
                            % Update Left Effective Hamiltonian
                            %%%%%%%%%%%%%%%%%%%%%%%
                            if m ==(M-1)
                                L{m+1} = mps_var.Add_Left_Eff_Ham(H,INIT,m+1);
                            elseif m ~= 1
                                L{m+1} = mps_var.Add_Left_Eff_Ham(H,L{m+2},m+1);
                            end
                            
                            if m~=1
                                H_Eff=mps_var.Create_Eff_Ham(H,R{m-1},L{m+1},m);
                                [mps_var]=mps_var.Apply_1S_H_Left(H_Eff,-dt,m);
                            end
                        end
                    end
                    
                end
                        
%                 mps_var=mps_var.Canonicalise_1s('R-L','true');
                Check_Norm = mps_var.Full_Norm;
                P_Num=0;
                for m =1:M
                    P_Num=P_Num+mps_var.Site_Site_Particle_Corr(m,m);
                end
                
                Energy=mps_var.Get_Energy(H);
                Var=mps_var.Get_Var(H);
                                
                disp(['Time=',num2str(t*tau),' -- CPU Time=',num2str(toc)]);
                disp(['Truncation Error=',num2str(Total_error)]);
                disp(['Norm=',num2str(Check_Norm),' -- P_Num=',num2str(P_Num)]);
                disp(['Energy=',num2str(Energy),' -- Var=',num2str(Var)]);
            end
        end
      
        function H_Eff = Create_Eff_Ham_2s(mps_var,H,R,L,m)
            
            NR1=size(mps_var.data{m},2);
            NL2=min(mps_var.N+1,(mps_var.d-1)*(m)+1);
            
            NR1_vec=NR1:-1:1;
            NL2_vec=1:NL2;

            endR=NR1+1-(mps_var.N+1-NL2+1);
            startL=mps_var.N+1-NR1+1;
            
            startR=endR-length(NL2_vec)+startL;
            
            Col_vec=cell(1,min(length(NR1_vec),length(NL2_vec)));
            track=1;
            for n=startR:min(length(NR1_vec),length(NL2_vec))
                Col_vec{track}=[NR1_vec(n),NL2_vec(track-1+startL)];
                track=track+1;
            end
            
            NR=size(mps_var.data{m},2);
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            NR_min=mps_var.N+1-NL+1;
            [kappa1,kappa2]=size(H.data{m});
            
            d_cell=cell(1,NR);
            for NR_iter = 1:NR
               d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
            end
            
            dd_vec=(1:mps_var.d)'-1;
            NL_vec=(mps_var.N:-1:0);
            NL_vec=NL_vec(ones(mps_var.d,1),1:NR);
            NR_vec_temp=(0:mps_var.N)';
            NR_vec=zeros(mps_var.d,NR);
            for NR_iter=1:NR
                NR_vec(1:size(NR_iter:min(length(NR_vec_temp),(NR_iter+mps_var.d-1)),2),NR_iter)=NR_vec_temp(NR_iter:min(end,(NR_iter+mps_var.d-1)));
            end

            ind=find(NR_vec<(NR_min-1));
            NR_vec(ind)=0;
            NL_vec(ind)=0;
            
            List=-dd_vec(:,ones(1,size(NL_vec,2))) + NR_vec + NL_vec;
            List(List~=mps_var.N)=0;
            List(List==mps_var.N)=1;
            
            NR_vec=(NR_vec+1).*List;
            
            NR2=size(mps_var.data{m+1},2);
            NL2=min(mps_var.N+1,(mps_var.d-1)*(m)+1);
            NR2_min=mps_var.N+1-NL2+1;
            [~,kappa4]=size(H.data{m+1});
            
            d_cell2=cell(1,NR2);
            for NR2_iter = 1:NR2
               d_cell2{1,NR2_iter}=max(1,mps_var.N+3-NR2_iter-NL2):min(mps_var.d,mps_var.N+2-NR2_iter);
            end
            
            NL2_vec=(mps_var.N:-1:0);
            NL2_vec=NL2_vec(ones(mps_var.d,1),1:NR2);
            NR2_vec_temp=(0:mps_var.N)';
            NR2_vec=zeros(mps_var.d,NR2);
            for NR2_iter=1:NR2
                NR2_vec(1:size(NR2_iter:min(length(NR_vec_temp),(NR2_iter+mps_var.d-1)),2),NR2_iter)=NR2_vec_temp(NR2_iter:min(end,(NR2_iter+mps_var.d-1)));
            end

            ind=find(NR2_vec<(NR2_min-1));
            NR2_vec(ind)=0;
            NL2_vec(ind)=0;
            
            List2=-dd_vec(:,ones(1,size(NL2_vec,2))) + NR2_vec + NL2_vec;
            List2(List2~=mps_var.N)=0;
            List2(List2==mps_var.N)=1;
            
            NL2_vec=(NL2_vec+1).*List2;
            
            shift=NL2_vec-dd_vec(:,ones(1,size(NL2_vec,2)));
            shift(shift<0)=0;
            
            List2L=zeros(mps_var.d,NL2);
            NL_Left_Vec=zeros(mps_var.d,NL2);
            for iter=1:mps_var.d
                shift_L=shift(iter,:);
                shift_L(shift_L==0)=[];
                start_shift=find(List2(iter,:)~=0);
                List2L(iter,shift_L)=List2(iter,start_shift);
                NL_Left_Vec(iter,shift_L)=NL2_vec(iter,start_shift);
            end
            
            dd_List=zeros(mps_var.d^2,size(Col_vec,2));
            List1=List(:,end:-1:1);
            dd_cell=cell(1,size(Col_vec,2));
            for n = 1:size(Col_vec,2)
                dd_List(:,n)=kron(List1(:,startR+n-1),List2L(:,startL+n-1));
                dd_cell{n}=find(dd_List(:,n)~=0);
            end
            dd_map=[kron(dd_vec,ones(mps_var.d,1)),kron(ones(mps_var.d,1),dd_vec)]+1;
            
            % dd_List == (d1 d2),(NR1 NL2)
            Len_dd=reshape(dd_List,mps_var.d^2*size(Col_vec,2),1); %(NR1 NL2 d1 d2)
            Ind=find(Len_dd~=0);
            
            Theta = cell(nnz(Len_dd),1);
            for k1=1:kappa1
                for k4=1:kappa4
                    
                    for k2=1:kappa2
                        if ~all(all(H.data{m}{k1,k2}==0)) && ~all(all(H.data{m+1}{k2,k4}==0))
                            if all(all(H.data{m}{k1,k2}==1))
                                O1=eye(mps_var.d);
                            else
                                O1=H.data{m}{k1,k2};
                            end
                            if all(all(H.data{m+1}{k2,k4}==1))
                                O2=eye(mps_var.d);
                            else
                                O2=H.data{m+1}{k2,k4};
                            end
                            O=kron(O1,O2);
                            
                            H_temp=zeros(mps_var.d^2,mps_var.d^2,size(Col_vec,2),size(Col_vec,2));
                            
                            NR1_shift=H.N_track{m}{k1,k2}(2);
                            
                            for N_iter=1:size(Col_vec,2)
                                if N_iter-NR1_shift>0 && N_iter-NR1_shift<=size(Col_vec,2)
                                    H_temp(dd_cell{N_iter-NR1_shift},dd_cell{N_iter},N_iter-NR1_shift,N_iter)=O(dd_cell{N_iter-NR1_shift},dd_cell{N_iter});
                                end
                            end
                            
                            % H_temp == (d1' d2'),(d1 d2),(NR1' NL2'),(NR1 NL2)
                            H_temp=permute(H_temp,[1 3 2 4]); % (d1' d2'),(NR1' NL2'),(d1 d2),(NR1 NL2)
                            H_temp=reshape(H_temp,mps_var.d^2*size(Col_vec,2),size(Col_vec,2)*mps_var.d^2); %(NR1' NL2' d1' d2'),(NR1 NL2 d1 d2)
                            
                            trackL=1;
                            for N_iter=1:size(Col_vec,2)
                                for d_iter=dd_cell{N_iter}'
                                    pos=(N_iter-1)*(mps_var.d^2)+d_iter;
                                    % Theta == (u' NR1' NL2' d1' d2' v'),(u NR1 NL2 d1 d2 v)
                                    if ~isempty(Theta{trackL,1})
                                        Theta{trackL,1} = Theta{trackL,1}+kron(kron(R{k1}{NR_vec(dd_map(d_iter,1),Col_vec{1,N_iter}(1))},H_temp(pos,Ind)),transpose(L{k4}{NL_Left_Vec(dd_map(d_iter,2),Col_vec{1,N_iter}(2))}));
                                    else
                                        Theta{trackL,1} = kron(kron(R{k1}{NR_vec(dd_map(d_iter,1),Col_vec{1,N_iter}(1))},H_temp(pos,Ind)),transpose(L{k4}{NL_Left_Vec(dd_map(d_iter,2),Col_vec{1,N_iter}(2))}));
                                    end
                                    trackL=trackL+1;
                                end
                            end
                            
                        end
                    end
                end
            end  
            % Theta == (NR1' NL2' d1' d2'),1 -- {Theta} == (u' v'),(u NR1 NL2 d1 d2 v)
            [u,~]=size(mps_var.data{m}{1});
            [~,v]=size(mps_var.data{m+1}{1});
            Theta(cellfun(@isempty,Theta)) = {zeros(u*v,u*nnz(Len_dd)*v)};
            H_Eff=cell2mat(Theta); % (NR1' NL2' d1' d2' u' v'),(u NR1 NL2 d1 d2 v)
            H_Eff = reshape(H_Eff,v,u,nnz(Len_dd),u*nnz(Len_dd)*v); %v',u',(NR1' NL2' d1' d2'),(u NR1 NL2 d1 d2 v)
            H_Eff = permute(H_Eff,[1 3 2 4]); % v',(NR1' NL2' d1' d2'),u',(u NR1 NL2 d1 d2 v)
            H_Eff = reshape(H_Eff,u*nnz(Len_dd)*v,u*nnz(Len_dd)*v); %(u' NR1' NL2' d1' d2' v'),(u NR1 NL2 d1 d2 v)
            H_Eff(abs(H_Eff)<mps_var.zero_thres)=0;
            H_Eff=sparse(H_Eff);
        end    
        
        function [mps_var,Total_error] = Apply_H_Eff_2s(mps_var,H_Eff,dt,m,Sweep_Direction)

            [~,NR_2] = size(mps_var.data{m+1});
            
            Nl_2 = min(mps_var.N+1,(mps_var.d-1)*m+1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the tensor for the two-sites so we can apply the operator.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} to N_Right labelling
            %%%%%%%%%%%%%%%%%%
            temp_T = cell(mps_var.d,NR_2);
            temp_IT = cell(mps_var.d,Nl_2);
            expect = max(NR_2,Nl_2):-1:1;
            for k = 1:mps_var.d
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2>mps_var.N+1-NR_2-k+1 && k_2<=Nl_2
                        temp_T(k,k_2) = {expect(1,max(NR_2,Nl_2)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_var.N+1-NR_2-k+2):min(Nl_2,(mps_var.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_var.data{m+1};
            
            R_mm = cell(mps_var.d,Nl_2);
            for alpha = 1:mps_var.d
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            [nrows, ncols] = cellfun(@size, temp);
            beta = max(max(nrows));
            v = max(max(ncols));
            R_mm(cellfun(@isempty,R_mm)) = {zeros(beta,v)};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Collect entries of (L_m X R_mm)
            % Construct Two-Site Tensor
            NR1=size(mps_var.data{m},2);
            NL2=min(mps_var.N+1,(mps_var.d-1)*(m)+1);
            NR1_vec=NR1:-1:1;
            NL2_vec=1:NL2;
            endR=NR1+1-(mps_var.N+1-NL2+1);
            startL=mps_var.N+1-NR1+1;
            startR=endR-length(NL2_vec)+startL;
            
            Col_vec=cell(1,min(length(NR1_vec),length(NL2_vec)));
            track=1;
            for n=startR:min(length(NR1_vec),length(NL2_vec))
                Col_vec{track}=[NR1_vec(n),NL2_vec(track-1+startL)];
                track=track+1;
            end
            
            NR=size(mps_var.data{m},2);
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            NR_min=mps_var.N+1-NL+1;
            
            d_cell=cell(1,NR);
            for NR_iter = 1:NR
               d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
            end
            
            dd_vec=(1:mps_var.d)'-1;
            NL_vec=(mps_var.N:-1:0);
            NL_vec=NL_vec(ones(mps_var.d,1),1:NR);
            NR_vec_temp=(0:mps_var.N)';
            NR_vec=zeros(mps_var.d,NR);
            for NR_iter=1:NR
                NR_vec(1:size(NR_iter:min(length(NR_vec_temp),(NR_iter+mps_var.d-1)),2),NR_iter)=NR_vec_temp(NR_iter:min(end,(NR_iter+mps_var.d-1)));
            end

            ind=find(NR_vec<(NR_min-1));
            NR_vec(ind)=0;
            NL_vec(ind)=0;
            
            List=-dd_vec(:,ones(1,size(NL_vec,2))) + NR_vec + NL_vec;
            List(List~=mps_var.N)=0;
            List(List==mps_var.N)=1;
            
            NR_vec=(NR_vec+1).*List;
            
            NR2=size(mps_var.data{m+1},2);
            NL2=min(mps_var.N+1,(mps_var.d-1)*(m)+1);
            NR2_min=mps_var.N+1-NL2+1;
            
            d_cell2=cell(1,NR2);
            for NR2_iter = 1:NR2
               d_cell2{1,NR2_iter}=max(1,mps_var.N+3-NR2_iter-NL2):min(mps_var.d,mps_var.N+2-NR2_iter);
            end
            
            NL2_vec=(mps_var.N:-1:0);
            NL2_vec=NL2_vec(ones(mps_var.d,1),1:NR2);
            NR2_vec_temp=(0:mps_var.N)';
            NR2_vec=zeros(mps_var.d,NR2);
            for NR2_iter=1:NR2
                NR2_vec(1:size(NR2_iter:min(length(NR_vec_temp),(NR2_iter+mps_var.d-1)),2),NR2_iter)=NR2_vec_temp(NR2_iter:min(end,(NR2_iter+mps_var.d-1)));
            end

            ind=find(NR2_vec<(NR2_min-1));
            NR2_vec(ind)=0;
            NL2_vec(ind)=0;
            
            List2=-dd_vec(:,ones(1,size(NL2_vec,2))) + NR2_vec + NL2_vec;
            List2(List2~=mps_var.N)=0;
            List2(List2==mps_var.N)=1;
            
            NL2_vec=(NL2_vec+1).*List2;
                       
            shift=NL2_vec-dd_vec(:,ones(1,size(NL2_vec,2)));
            shift(shift<0)=0;
            
            List2L=zeros(mps_var.d,NL2);
            NL_Left_Vec=zeros(mps_var.d,NL2);
            for iter=1:mps_var.d
                shift_L=shift(iter,:);
                shift_L(shift_L==0)=[];
                start_shift=find(List2(iter,:)~=0);
                List2L(iter,shift_L)=List2(iter,start_shift);
                NL_Left_Vec(iter,shift_L)=NL2_vec(iter,start_shift);
            end
            
            dd_List=zeros(mps_var.d^2,size(Col_vec,2));
            List1=List(:,end:-1:1);
            dd_cell=cell(1,size(Col_vec,2));
            for n = 1:size(Col_vec,2)
                dd_List(:,n)=kron(List1(:,startR+n-1),List2L(:,startL+n-1));
                dd_cell{n}=find(dd_List(:,n)~=0);
            end
            dd_map=[kron(dd_vec,ones(mps_var.d,1)),kron(ones(mps_var.d,1),dd_vec)]+1;
            
            % dd_List == (d1 d2),(NR1 NL2)
            Len_dd=reshape(dd_List,mps_var.d^2*size(Col_vec,2),1); %(NR1 NL2 d1 d2)
            Ind=find(Len_dd~=0);
            
            Theta_Mat=cell(size(dd_List));
            
            for N_iter=1:size(Col_vec,2)
                for d_iter=dd_cell{N_iter}'
                    if ~isempty(Theta_Mat{d_iter,N_iter})
                        Theta_Mat{d_iter,N_iter}=Theta_Mat{d_iter,N_iter}+mps_var.data{m}{dd_map(d_iter,1),Col_vec{1,N_iter}(1)}*R_mm{dd_map(d_iter,2),Col_vec{1,N_iter}(2)};
                    else
                        Theta_Mat{d_iter,N_iter}=mps_var.data{m}{dd_map(d_iter,1),Col_vec{1,N_iter}(1)}*R_mm{dd_map(d_iter,2),Col_vec{1,N_iter}(2)};
                    end
                end
            end
            % Theta_Mat == (d1 d2),(NR1 NL2) -- {Theta_Mat} == u,v
            
            [u,~] = size(mps_var.data{m}{1});
            [~,v] = size(mps_var.data{m+1}{1,1});
            
            Theta_Mat = reshape(Theta_Mat,mps_var.d^2*size(Col_vec,2),1); %(NR1 NL2 d1 d2),1
            Theta_Mat = Theta_Mat(Ind,1);
            Theta_Mat = cell2mat(Theta_Mat); % (NR1 NL2 d1 d2 u),v
            Theta_Mat = reshape(Theta_Mat,u,nnz(Len_dd),v); % u,(NR1 NL2 d1 d2),v
            Theta_Mat = permute(Theta_Mat,[3 2 1]); % v,(NR1 NL2 d1 d2),u
            Theta_Mat = reshape(Theta_Mat,u*nnz(Len_dd)*v,1);% (u NR1 NL2 d1 d2 v),1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Apply Effective Hamiltonian
            % H_Eff == (u' NR1' NL2' d1' d2' v'),(u NR1 NL2 d1 d2 v)

            Theta_Mat=expv(-1i*dt,H_Eff,Theta_Mat); %(u NR' d' v),1
            
%             U_Eff=expm(-1i*H_Eff*dt); % (u NR' d' v),(u NR d v)
%             Theta_Mat=U_Eff*Theta_Mat; %(u NR' d' v),1

%             Theta_Mat = Theta_Mat - 1i*dt*H_Eff*Theta_Mat;

            Theta_Mat = reshape(Theta_Mat,v,nnz(Len_dd),u); %v,(NR1 NL2 d1 d2),u
            Theta_Mat = permute(Theta_Mat,[3 2 1]); %  u,(NR1 NL2 d1 d2),v
            Theta_Mat = reshape(Theta_Mat,u*nnz(Len_dd),v); % (NR1 NL2 d1 d2 u),v

            row_vec=u*ones(nnz(Len_dd),1);
            Theta_Mat = mat2cell(Theta_Mat,row_vec,v); % (NR1 NL2 d1 d2),1
            
            Theta_Mat_New=cell(mps_var.d^2*size(Col_vec,2),1);
            Theta_Mat_New(Ind,1)=Theta_Mat; % (NR1 NL2 d1 d2),1
            
            Theta_Mat_New=reshape(Theta_Mat_New,mps_var.d^2,size(Col_vec,2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do SVD's
            R_mm = cell(mps_var.d,NL2);
            L_m = cell(mps_var.d,NR);

            r_min=1;
            Norm=0;
            Total_error=0;
            for count2 = 1:size(Theta_Mat_New,2)
                Theta = cell2mat(Theta_Mat_New(:,count2));
                
                d_1=nnz(List(:,Col_vec{count2}(1)));
                d_2=nnz(List2L(:,Col_vec{count2}(2)));
                
                % Theta == (d_1 d_2 alpha),(gamma)
                Theta = reshape(Theta,u,d_2,d_1,v); %(a),(l+1),(l),(g)
                Theta = permute(Theta,[1 3 4 2]); %(a),(l),(g),(l+1)
                Theta = reshape(Theta,u*d_1,v*d_2); %(l a),(l+1 g)
                Theta = reshape(Theta,u*d_1,v*d_2); %(l a),(l+1 g)
                [A,S,V]=svd(Theta);
                
                S((abs(S))/S(1,1)<mps_var.zero_thres)=0;
                
                r = nnz(S);
                if r>r_min
                    r_min=r;
                    if r_min>mps_var.Dmax
                        r_min=mps_var.Dmax;
                    end
                end
                if r==0
                    r=1;
                    S = S(1:r,1:r);
                    Error=0;
                elseif r>mps_var.Dmax
                    r=mps_var.Dmax;
                    if any(size(S)==1)
                        Error=sum(S(r+1:end,r+1:end).^2);
                    else
                        Error=sum(diag(S(r+1:end,r+1:end)).^2);
                    end
                    S = S(1:r,1:r);
                else
                    Error=0;
                    S = S(1:r,1:r);
                end
%                 diag(S)'*diag(S)
                Norm = Norm+diag(S)'*diag(S);
                if strcmp(Sweep_Direction,'L-R')
                    A = A(:,1:r);
                    V = S*(V(:,1:r)');
                else
                    A = A(:,1:r)*S;
                    V = V(:,1:r)';
                end
                
                Total_error=Total_error+Error;
                
                % A == (l a),(r)
                A = reshape(A,u,d_1,r); %(a),(l),(r)
                A = permute(A,[1 3 2]); %(a),(r),(l)
                track=1;
                for count = unique(dd_map(dd_cell{count2},1))'
                    L_m(count,Col_vec{count2}(1)) = {A(:,:,track)};
                    track=track+1;
                end
                
                % V == (r),(l+1 g)
                V = reshape(V,r,v,d_2); %(r),(g),(l+1)
                track=1;
                for count = unique(dd_map(dd_cell{count2},2))'
                    R_mm(count,Col_vec{count2}(2)) = {V(:,:,track)};
                    track=track+1;
                end
            end
            R_mm(cellfun(@isempty,R_mm)) = {zeros(r_min,v)};
            L_m(cellfun(@isempty,L_m)) = {zeros(u,r_min)};
                        
            [nrows, ~] = cellfun(@size, L_m);
            u = max(max(nrows));
            [~, ncol] = cellfun(@size, R_mm);
            v = max(max(ncol));
            
            L_m = cellfun(@(x) cell_resize_func(x,[u,r_min]),L_m,'UniformOutput', false);
            R_mm = cellfun(@(x) cell_resize_func(x,[r_min,v]),R_mm,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update First Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mps_var.data{m}=L_m;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update Second Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} back to N_Left labelling
            %%%%%%%%%%%%%%%%%%
            
            temp = cell(mps_var.d,NR_2);
            for alpha = 1:mps_var.d
                temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(r_min,v)};
            
            mps_var.data{m+1}=temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            function Out = cell_resize_func(X,r)
                % Function used to resize MPS
                % tensor.
                
                [a,b2]= size(X);
                Out=X;
                if a < r(1)
                    Out(a+1:r(1),:) = zeros(r(1)-a,size(Out,2));
                end
                if b2 < r(2)
                    Out(:,b2+1:r(2)) = zeros(size(Out,1),r(2)-b2);
                end
                
            end
        
        end
        
        function [mps_var] = Apply_1S_H_Right(mps_var,H_Eff,dt,m)
            
            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            d_cell=cell(1,NR);
            Vec=cell(1,NR);
            List=zeros(mps_var.d,NR);
            for NR_iter = 1:NR
                d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                Vec(d_cell{1,NR_iter},NR_iter)=mps_var.data{m}(d_cell{1,NR_iter},NR_iter);
                List(d_cell{1,NR_iter},NR_iter)=ones(length(d_cell{1,NR_iter}),1);
            end
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            % Vec == d,NR
            Vec = reshape(Vec,mps_var.d*NR,1); % (NR d),1
            Vec = cell2mat(Vec(Ind,1)); % (NR d u),v
                
            Vec = reshape(Vec,u,nnz(Len),v); %u,(NR d),v
            Vec = permute(Vec,[3,2,1]); %v,(NR d),u
            Vec = reshape(Vec,u*nnz(Len)*v,1); %(u NR d v),1
            
            Vec=expv(-1i*dt,H_Eff,Vec); %(u NR' d' v),1
            
%             U_Eff=expm(-1i*H_Eff*dt); % (u NR' d' v),(u NR d v)
%             Vec=U_Eff*Vec; %(u NR' d' v),1

%             Vec = Vec - 1i*dt*H_Eff*Vec;
            
            Vec = reshape(Vec,v,nnz(Len),u); % v,(NR d),u
            Vec = permute(Vec,[3,2,1]); % u,(NR d),v
            Vec = reshape(Vec,u*nnz(Len),v); % (NR d u),v
            
            row_vec=u*ones(nnz(Len),1);
            Vec = mat2cell(Vec,row_vec,v); % (NR d),1
            
            New_Vec=cell(mps_var.d*NR,1);
            New_Vec(Ind,1)=Vec;
            New_Vec=reshape(New_Vec,mps_var.d,NR);
            New_Vec(cellfun(@isempty,New_Vec)) = {zeros(u,v)};
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do QR decomposition
            Q=cell(mps_var.d,NR);
            S=cell(1,NR);
            for NR_iter=1:NR
                d_vec=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                if ~isempty(d_vec)                    
                    [U,Sv]=qr(cell2mat(New_Vec(d_vec,NR_iter)));
                    
                    Sv(abs(Sv)<mps_var.zero_thres)=0;
                    U(abs(U)<mps_var.zero_thres)=0;
                    
                    dd_vec = ((d_vec(1)-1)*u+1):(d_vec(end)-1)*u+u;
                    
                    U_temp=zeros(mps_var.d*u,min(size(U,2)/length(d_vec)*mps_var.d,v));
                    U_temp(dd_vec,1:min(size(U,2),v))=U(:,1:min(size(U,2),v));
                    
                    cell_vec=u*ones(mps_var.d,1);
                    Q(:,NR_iter)=mat2cell(U_temp,cell_vec,min(size(U,2)/length(d_vec)*mps_var.d,v));
                    
                    Sv_temp=zeros(size(U_temp,2),v);
                    Sv_temp(1:min(size(U,2),v),:)=Sv(1:min(size(U,2),v),:);

                    S{1,NR_iter}=Sv_temp;
                end
            end
%             S(cellfun(@isempty,S)) = {zeros(v,v)};
%             if strcmp(Norm,'true')
%                 norm_val=0;
%                 for NR_iter=1:NR
%                     norm_val=norm_val+S{1,NR_iter}*S{1,NR_iter}';
%                 end
%                 
%                 norm_val=norm_val(1);
%                 S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
%             end
            Q(cellfun(@isempty,Q)) = {zeros(u,v)};
            mps_var.data{m}=Q;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if m~=size(mps_var.data,2)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply Singular values to next tensor
                [~,NR_1]=size(mps_var.data{m+1});
                
                temp_E = cell(mps_var.d,NR_1);
                for k = 1:mps_var.d
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2 <= NR_1
                            temp_E{k,k_2} = S{1,mps_var.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(v,v)};
                
                mps_var.data{m+1} = cellfun(@(y,z) y*z, temp_E,mps_var.data{m+1}, 'UniformOutput', false);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        end
        
        function [mps_var] = Apply_1S_H_Left(mps_var,H_Eff,dt,m)

            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            d_cell=cell(1,NR);
            Vec=cell(1,NR);
            List=zeros(mps_var.d,NR);
            for NR_iter = 1:NR
                d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                Vec(d_cell{1,NR_iter},NR_iter)=mps_var.data{m}(d_cell{1,NR_iter},NR_iter);
                List(d_cell{1,NR_iter},NR_iter)=ones(length(d_cell{1,NR_iter}),1);
            end
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            % Vec == d,NR
            Vec = reshape(Vec,mps_var.d*NR,1); % (NR d),1
            Vec = cell2mat(Vec(Ind,1)); % (NR d u),v
                
            Vec = reshape(Vec,u,nnz(Len),v); %u,(NR d),v
            Vec = permute(Vec,[3,2,1]); %v,(NR d),u
            Vec = reshape(Vec,u*nnz(Len)*v,1); %(u NR d v),1
            
            Vec=expv(-1i*dt,H_Eff,Vec); %(u NR' d' v),1
            
%             U_Eff=expm(-1i*H_Eff*dt); % (u NR' d' v),(u NR d v)
%             Vec=U_Eff*Vec; %(u NR' d' v),1

%             Vec = Vec - 1i*dt*H_Eff*Vec;

            Vec = reshape(Vec,v,nnz(Len),u); % v,(NR d),u
            Vec = permute(Vec,[3,2,1]); % u,(NR d),v
            Vec = reshape(Vec,u*nnz(Len),v); % (NR d u),v
            
            row_vec=u*ones(nnz(Len),1);
            Vec = mat2cell(Vec,row_vec,v); % (NR d),1
            
            New_Vec=cell(mps_var.d*NR,1);
            New_Vec(Ind,1)=Vec;
            New_Vec=reshape(New_Vec,mps_var.d,NR);
            New_Vec(cellfun(@isempty,New_Vec)) = {zeros(u,v)};
            mps_var.data{m}=New_Vec;
            
            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Convert mth Tensor to N_Left Labelling
            temp_T = cell(mps_var.d,NR);
            temp_IT = cell(mps_var.d,NL);
            expect = max(NR,NL):-1:1;
            for k = 1:mps_var.d
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2>mps_var.N+1-NR-k+1 && k_2<=NL
                        temp_T(k,k_2) = {expect(1,max(NR,NL)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_var.N+1-NR-k+2):min(NL,(mps_var.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_var.data{m};
            
            R_mm = cell(mps_var.d,NL);
            for alpha = 1:mps_var.d
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do QR decomposition
            Q=cell(mps_var.d,NL);
            S=cell(1,NL);
            if m~=1
                for NL_iter=1:NL
                    d_vec=max(1,mps_var.N+3-NL_iter-NR):min(mps_var.d,mps_var.N+2-NL_iter);
                    if ~isempty(d_vec)
                        Mat = cell2mat(R_mm(d_vec,NL_iter)); % (d u),v
                        Mat = reshape(Mat,u,length(d_vec),v); % u,d,v
                        Mat = permute(Mat,[1 3 2]); % u,v,d
                        Mat = reshape(Mat,u,length(d_vec)*v); % u,(d v)
                        Mat = transpose(Mat); %(d v),u
                        
                        [V,R]=qr(Mat); % V=(d v),(d v)
                        
                        Sv=transpose(R);
                        V=transpose(V);
                        
                        Sv(abs(Sv)<mps_var.zero_thres)=0;
                        V(abs(V)<mps_var.zero_thres)=0;
                        
                        dd_vec = ((d_vec(1)-1)*v+1):(d_vec(end)-1)*v+v;
                        
                        V_temp=zeros(min(size(V,1)/length(d_vec)*mps_var.d,u),mps_var.d*v);
                        V_temp(1:min(size(V,1),u),dd_vec)=V(1:min(size(V,1),u),:); % u,(v d)
                        
                        Mat = reshape(V_temp,u,v,mps_var.d); % u,v,d
                        Mat = permute(Mat,[1 3 2]); %u,d,v
                        Mat = reshape(Mat,u*mps_var.d,v);% (d u),v
                        
                        cell_vec=u*ones(mps_var.d,1);
                        Q(:,NL_iter)=mat2cell(Mat,cell_vec,v);
                        
                        Sv_temp=zeros(u,size(V_temp,1));
                        Sv_temp(:,1:min(size(V,1),u))=Sv(:,1:min(size(V,1),u));
                        
                        S{1,NL_iter}=Sv_temp;
                    end
                end
                S(cellfun(@isempty,S)) = {zeros(u,u)};
%                 if strcmp(Norm,'true')
%                     norm_val=0;
%                     for NL_iter=1:NL
%                         norm_val=norm_val+S{1,NL_iter}'*S{1,NL_iter};
%                     end
%                     
%                     norm_val=norm_val(1);
%                     S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
%                 end
                Q(cellfun(@isempty,Q)) = {zeros(u,v)};
                R_mm=Q;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Convert mth tensor back into N_Right labelling
            temp = cell(mps_var.d,NR);
            for alpha = 1:mps_var.d
                temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(u,v)};
            
            mps_var.data{m}=temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if m~=1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Convert (m-1)th Tensor to N_Left Labelling
                [~,NR_1]=size(mps_var.data{m-1});
                NL_1=min(mps_var.N+1,(mps_var.d-1)*(m-2)+1);
                
                temp_T = cell(mps_var.d,NR_1);
                temp_IT = cell(mps_var.d,NL_1);
                expect = max(NR_1,NL_1):-1:1;
                for k = 1:mps_var.d
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2>mps_var.N+1-NR_1-k+1 && k_2<=NL_1
                            temp_T(k,k_2) = {expect(1,max(NR_1,NL_1)-track+1)};
                        end
                        track=track+1;
                    end
                    for k_2 = max(1,mps_var.N+1-NR_1-k+2):min(NL_1,(mps_var.N+2-k))
                        temp_IT(k,k_2) = {k_2};
                    end
                end
                
                temp = mps_var.data{m-1};
                
                R_mm_1 = cell(mps_var.d,NL_1);
                for alpha = 1:mps_var.d
                    R_mm_1(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
                end
                R_mm_1(cellfun(@isempty,R_mm_1)) = {zeros(size(mps_var.data{m-1}{1}))};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply Singular values to next tensor
                temp_E = cell(mps_var.d,NL_1);
                for k = 1:mps_var.d
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2 <= NL_1
                            temp_E{k,k_2} = S{1,mps_var.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(u,u)};
                
                R_mm_1 = cellfun(@(y,z) y*z, R_mm_1,temp_E, 'UniformOutput', false);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Convert (m-1)th tensor back into N_Right labelling
                temp = cell(mps_var.d,NR_1);
                for alpha = 1:mps_var.d
                    temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm_1(alpha,cell2mat(temp_IT(alpha,:)));
                end
                temp(cellfun(@isempty,temp)) = {zeros(size(mps_var.data{m-1}{1}))};
                
                mps_var.data{m-1}=temp;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        end       
        
        function L = Add_Left_Eff_Ham(mps_var,H,L,m)
            
            [dd,NR] = size(mps_var.data{m});
            Nl = min(mps_var.N+1,(mps_var.d-1)*m+1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Change current Tensor to N_left labelling
            temp_T = cell(dd,NR);
            temp_IT = cell(dd,Nl);
            expect = max(NR,Nl):-1:1;
            for k = 1:dd
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2>mps_var.N+1-NR-k+1 && k_2<=Nl
                        temp_T(k,k_2) = {expect(1,max(NR,Nl)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_var.N+1-NR-k+2):min(Nl,(mps_var.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_var.data{m};
            
            R_mm = cell(dd,Nl);
            for alpha = 1:dd
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            [nrows, ncols] = cellfun(@size, temp);
            a = max(max(nrows));
            g = max(max(ncols));
            R_mm(cellfun(@isempty,R_mm)) = {zeros(a,g)};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Apply the MPO H{m} so site m
            Store=mps_var.Apply_MPO_Single_Site_Left_Form(R_mm,H,m);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Contract Store with exisiting L{m+1}
            temp=cell(size(L,1),1);
            for n = 1:size(L,1)
                temp_E = cell(size(R_mm));
                for k = 1:dd
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2 <= size(R_mm,2)
                            temp_E{k,k_2} = L{n}{1,mps_var.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(g,g)};
                
                temp{n,1} = temp_E;
            end
            
            % H{m} == a,b  -- {mps_var.data{m}} == u,v
            % Store == a,b -- {Store} == d,NR -- {{Store}} == u,v
            % temp == b,1  -- {temp} == d,NR  -- {{temp}} == v,v'
            
            L_mm=cell(size(Store,1),1);
            for alpha=1:size(Store,1)
                for beta=1:size(Store,2)
                    if ~isempty(Store{alpha,beta})
                        if ~isempty(L_mm{alpha,1})
                            Val=cellfun(@(x,y) x*y, Store{alpha,beta}, temp{beta,1}, 'UniformOutput', false);
                            L_mm{alpha,1}=cellfun(@plus, L_mm{alpha,1}, Val, 'UniformOutput', false);
                        else
                            L_mm{alpha,1}=cellfun(@(x,y) x*y, Store{alpha,beta}, temp{beta,1}, 'UniformOutput', false);
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Finally contract conjugate of mth tensor
            L=cell(size(L_mm));
            for n = 1:size(L_mm,1)
                temp=cellfun(@(x,y) x*y',L_mm{n,1},R_mm, 'UniformOutput', false);
                
                expect=cell(1,size(temp,2));
                expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                for i=1:dd
                    expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                end
                
                L{n} = expect;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function Store = Apply_MPO_Single_Site_Left_Form(mps_var,R_mm,H,m)
            % Apply the MPO H to the MPS
            
            M = size(mps_var.data,2);
            
            [k1,k2]=size(H.data{m});
            [alpha,beta]=size(R_mm{1,1});
            Nl = min((mps_var.d-1)*(M-m)+1,mps_var.N+1);
            
            [a,g] = size(R_mm);
            Store=cell(k1,k2);
            for a_count = 1:k1
                for b_count = 1:k2
                    if all(all(H.data{m}{a_count,b_count}==1))
                        Store{a_count,b_count}=R_mm;
                    elseif ~all(all(H.data{m}{a_count,b_count}==0))
                        O=H.data{m}{a_count,b_count};
                        
                        temp = cell(a,g);
                        temp(cellfun(@isempty,temp)) = {zeros(alpha,beta)};
                        for l = 1:mps_var.d
                            for N_R = max(1,(g-Nl-l+2)):max(1,min(g,mps_var.N+1-l+1))
                                for count2 = (-mps_var.d+l):(l-1)
                                    temp{l-count2,N_R} = temp{l-count2,N_R}+O(l-count2,l)*R_mm{l,N_R};
                                end
                            end
                        end
                        
                        
                        if H.N_track{m}{a_count,b_count}(1)>0
                            m_temp = cell(a,g);
                            for m1 = (1+H.N_track{m}{a_count,b_count}(1)):g
                                m_temp(:,m1) = temp(:,m1-H.N_track{m}{a_count,b_count}(1));
                            end
                            m_temp(cellfun(@isempty,m_temp)) = {zeros(alpha,beta)};
                            temp = m_temp;
                        elseif H.N_track{m}{a_count,b_count}(1)<0
                            m_temp = cell(a,g);
                            for m1 = 1:(g+H.N_track{m}{a_count,b_count}(1))
                                m_temp(:,m1) = temp(:,m1-H.N_track{m}{a_count,b_count}(1));
                            end
                            m_temp(cellfun(@isempty,m_temp)) = {zeros(alpha,beta)};
                            temp=m_temp;
                        else
                            temp_E = cell(a,g);
                            for k = 1:a
                                track=1;
                                for k_2 = (mps_var.N+1-k+1):-1:1
                                    if k_2 <= g
                                        temp_E{k,k_2} = eye(alpha);
                                    end
                                    track=track+1;
                                end
                            end
                            temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                            temp=cellfun(@(x,y) x*y,temp_E,temp,'UniformOutput', false);
                        end
                        Store{a_count,b_count}=temp;
                    end
                    
                end
            end  
            Zero_cell=cell(mps_var.d,g);
            Zero_cell(cellfun(@isempty,Zero_cell)) = {zeros(alpha,beta)};
            Store(cellfun(@isempty,Store)) = {Zero_cell};
        end
        
        function R = Add_Right_Eff_Ham(mps_var,H,R,m)

            L_mm = mps_var.data{m};
            [dd,NR] = size(L_mm);
            [nrows, ~] = cellfun(@size, L_mm);
            a = max(max(nrows));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Apply the MPO H{m} so site m
            Store=mps_var.Apply_MPO_Single_Site_Right_Form(L_mm,H,m);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Contract Store with exisiting R{m-1}
            temp=cell(1,size(R,2));
            for n = 1:size(R,2)
                temp_E = cell(size(L_mm));
                for k = 1:mps_var.d
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2 <= NR
                            temp_E{k,k_2} = R{n}{1,mps_var.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(a,a)};
                
                temp{1,n} = temp_E;
            end
            
            % H{m} == a,b  -- {mps_var.data{m}} == u,v
            % Store == a,b -- {Store} == d,NR -- {{Store}} == u,v
            % temp == 1,a  -- {temp} == d,NR  -- {{temp}} == u',u
            
            R_mm=cell(1,size(Store,2));
            for beta=1:size(Store,2)
                for alpha=1:size(Store,1)
                    if ~isempty(Store{alpha,beta})
                        if ~isempty(R_mm{1,beta})
                            Val=cellfun(@(x,y) x*y, temp{1,alpha},Store{alpha,beta}, 'UniformOutput', false);
                            R_mm{1,beta}=cellfun(@plus, R_mm{1,beta}, Val, 'UniformOutput', false);
                        else
                            R_mm{1,beta}=cellfun(@(x,y) x*y, temp{1,alpha},Store{alpha,beta}, 'UniformOutput', false);
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Finally contract conjugate of mth tensor
            Z=cell(dd,NR);
            Z(cellfun(@isempty,Z)) = {zeros(size(L_mm{1}))};
            R_mm(cellfun(@isempty,R_mm)) = {Z};
            R=cell(size(R_mm));
            for n = 1:size(R_mm,2)
                temp=cellfun(@(x,y) x'*y,L_mm,R_mm{n}, 'UniformOutput', false);
                
                expect=cell(1,size(temp,2));
                expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                for i=1:dd
                    expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                end
                
                R{n} = expect;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function Store = Apply_MPO_Single_Site_Right_Form(mps_var,L_mm,H,m)
            % Apply the MPO H to the MPS
                        
            [k1,k2]=size(H.data{m});
            [alpha,beta]=size(L_mm{1,1});
            Nl = min((mps_var.d-1)*(m-1)+1,mps_var.N+1);
            
            [a,g] = size(L_mm);
            Store=cell(k1,k2);
            for a_count = 1:k1
                for b_count = 1:k2
                    if all(all(H.data{m}{a_count,b_count}==1))
                        Store{a_count,b_count}=L_mm;
                    elseif ~all(all(H.data{m}{a_count,b_count}==0))
                        O=H.data{m}{a_count,b_count};
                        
                        temp = cell(a,g);
                        temp(cellfun(@isempty,temp)) = {zeros(alpha,beta)};
                        for l = 1:mps_var.d
                            for N_R = max(1,(g-Nl-l+2)):max(1,min(g,mps_var.N+1-l+1))
                                for count2 = (-mps_var.d+l):(l-1)
                                    temp{l-count2,N_R} = temp{l-count2,N_R}+O(l-count2,l)*L_mm{l,N_R};
                                end
                            end
                        end
                        
                        if H.N_track{m}{a_count,b_count}(2)>0
                            m_temp = cell(a,g);
                            for m1 = (1+H.N_track{m}{a_count,b_count}(2)):g
                                m_temp(:,m1) = temp(:,m1-H.N_track{m}{a_count,b_count}(2));
                            end
                            m_temp(cellfun(@isempty,m_temp)) = {zeros(alpha,beta)};
                            temp = m_temp;
                        elseif H.N_track{m}{a_count,b_count}(2)<0
                            m_temp = cell(a,g);
                            for m1 = 1:(g+H.N_track{m}{a_count,b_count}(2))
                                m_temp(:,m1) = temp(:,m1-H.N_track{m}{a_count,b_count}(2));
                            end
                            m_temp(cellfun(@isempty,m_temp)) = {zeros(alpha,beta)};
                            temp=m_temp;
                        else
                            temp_E = cell(a,g);
                            for k = 1:a
                                track=1;
                                for k_2 = (mps_var.N+1-k+1):-1:1
                                    if k_2 <= g
                                        temp_E{k,k_2} = eye(alpha);
                                    end
                                    track=track+1;
                                end
                            end
                            temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                            temp=cellfun(@(x,y) x*y,temp_E,temp,'UniformOutput', false);
                        end
                        Store{a_count,b_count}=temp;
                    end
                    
                end
            end 
            Zero_cell=cell(mps_var.d,g);
            Zero_cell(cellfun(@isempty,Zero_cell)) = {zeros(alpha,beta)};
            Store(cellfun(@isempty,Store)) = {Zero_cell};
        end        
        
        function H_Eff = Create_Eff_Ham(mps_var,H,R,L,m)
            
            NR=size(mps_var.data{m},2);
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            NR_min=mps_var.N+1-NL+1;
            [kappa1,kappa2]=size(H.data{m});
            
            d_cell=cell(1,NR);
            for NR_iter = 1:NR
               d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
            end
            
            dd_vec=(1:mps_var.d)'-1;
            NL_vec=(mps_var.N:-1:0);
            NL_vec=NL_vec(ones(mps_var.d,1),1:NR);
            NR_vec_temp=(0:mps_var.N)';
            NR_vec=zeros(mps_var.d,NR);
            for NR_iter=1:NR
                NR_vec(1:size(NR_iter:min(length(NR_vec_temp),(NR_iter+mps_var.d-1)),2),NR_iter)=NR_vec_temp(NR_iter:min(end,(NR_iter+mps_var.d-1)));
            end

            ind=find(NR_vec<(NR_min-1));
            NR_vec(ind)=0;
            NL_vec(ind)=0;
            
            List=-dd_vec(:,ones(1,size(NL_vec,2))) + NR_vec + NL_vec;
            List(List~=mps_var.N)=0;
            List(List==mps_var.N)=1;
            
            NR_vec=(NR_vec+1).*List;
            NL_vec=(NL_vec+1).*List;
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            Theta = cell(nnz(Len),1);
            for k1=1:kappa1
                for k2=1:kappa2
                    if ~all(all(H.data{m}{k1,k2}==0))
                        if all(all(H.data{m}{k1,k2}==1))
                            O=eye(mps_var.d);
                        else
                            O=H.data{m}{k1,k2};
                        end
                        
                        H_temp=zeros(mps_var.d,mps_var.d,NR,NR);
                        
                        NR_shift=H.N_track{m}{k1,k2}(2);
                        if NR_shift==0
                            for NR_iter=1:NR
                                H_temp(d_cell{NR_iter},d_cell{NR_iter},NR_iter,NR_iter)=(O(d_cell{NR_iter},d_cell{NR_iter}));
                            end
                        else
                            for NR_iter=1:NR
                                if NR_iter+NR_shift>0 && NR_iter+NR_shift<=NR
                                    H_temp(d_cell{NR_iter+NR_shift},d_cell{NR_iter},NR_iter+NR_shift,NR_iter)=(O(d_cell{NR_iter+NR_shift},d_cell{NR_iter}));
                                end
                            end
                        end
                        % H_temp == d',d,NR',NR
                        H_temp=permute(H_temp,[1 3 2 4]); % d',NR',d,NR
                        H_temp=reshape(H_temp,mps_var.d*NR,NR*mps_var.d); %(NR' d'),(NR d)
                       
                        track=1;
                        for NR_iter=1:NR
                            for d_iter=d_cell{NR_iter}
                                pos=(NR_iter-1)*mps_var.d+d_iter;
                                % Theta == (u NR' d' v),(u NR d v)
                                if ~isempty(Theta{track,1})
                                    Theta{track,1} = Theta{track,1}+kron(kron(R{k1}{NR_vec(d_iter,NR_iter)},H_temp(pos,Ind)),transpose(L{k2}{NL_vec(d_iter,NR_iter)}));
                                else
                                    Theta{track,1} = kron(kron(R{k1}{NR_vec(d_iter,NR_iter)},H_temp(pos,Ind)),transpose(L{k2}{NL_vec(d_iter,NR_iter)}));
                                end
                                track=track+1;
                            end
                        end
                    end
                end
            end  
            % Theta == (NR d),1 -- {Theta} == (u v),(u NR d v)
            [u,v]=size(mps_var.data{m}{1});
            Theta(cellfun(@isempty,Theta)) = {zeros(u*v,u*nnz(Len)*v)};
            H_Eff=cell2mat(Theta); % (NR d u v),(u NR d v)
            H_Eff = reshape(H_Eff,v,u,nnz(Len),u*nnz(Len)*v); %v,u,(NR d),(u NR d v)
            H_Eff = permute(H_Eff,[1 3 2 4]); % v,(NR d),u,(u NR d v)
            H_Eff = reshape(H_Eff,u*nnz(Len)*v,u*nnz(Len)*v); %(u' NR' d' v'),(u NR d v)
            H_Eff(abs(H_Eff)<mps_var.zero_thres)=0;
            H_Eff=sparse(H_Eff);
        end
        
        function [mps_var,S] = Apply_H_Eff_Right(mps_var,H_Eff,dt,m,Norm)

            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            d_cell=cell(1,NR);
            Vec=cell(1,NR);
            List=zeros(mps_var.d,NR);
            for NR_iter = 1:NR
                d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                Vec(d_cell{1,NR_iter},NR_iter)=mps_var.data{m}(d_cell{1,NR_iter},NR_iter);
                List(d_cell{1,NR_iter},NR_iter)=ones(length(d_cell{1,NR_iter}),1);
            end
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            % Vec == d,NR
            Vec = reshape(Vec,mps_var.d*NR,1); % (NR d),1
            Vec = cell2mat(Vec(Ind,1)); % (NR d u),v
                
            Vec = reshape(Vec,u,nnz(Len),v); %u,(NR d),v
            Vec = permute(Vec,[3,2,1]); %v,(NR d),u
            Vec = reshape(Vec,u*nnz(Len)*v,1); %(u NR d v),1
            
            Vec=expv(-1i*dt,H_Eff,Vec); %(u NR' d' v),1
            
%             U_Eff=expm(-1i*H_Eff*dt); % (u NR' d' v),(u NR d v)
%             Vec=U_Eff*Vec; %(u NR' d' v),1

%             Vec = Vec - 1i*dt*H_Eff*Vec;
            
            Vec = reshape(Vec,v,nnz(Len),u); % v,(NR d),u
            Vec = permute(Vec,[3,2,1]); % u,(NR d),v
            Vec = reshape(Vec,u*nnz(Len),v); % (NR d u),v
            
            row_vec=u*ones(nnz(Len),1);
            Vec = mat2cell(Vec,row_vec,v); % (NR d),1
            
            New_Vec=cell(mps_var.d*NR,1);
            New_Vec(Ind,1)=Vec;
            New_Vec=reshape(New_Vec,mps_var.d,NR);
            New_Vec(cellfun(@isempty,New_Vec)) = {zeros(u,v)};
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do QR decomposition
            Q=cell(mps_var.d,NR);
            S=cell(1,NR);
            for NR_iter=1:NR
                d_vec=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                if ~isempty(d_vec)                    
                    [U,Sv]=qr(cell2mat(New_Vec(d_vec,NR_iter)));
                    
                    Sv(abs(Sv)<mps_var.zero_thres)=0;
                    U(abs(U)<mps_var.zero_thres)=0;
                    
                    dd_vec = ((d_vec(1)-1)*u+1):(d_vec(end)-1)*u+u;
                    
                    U_temp=zeros(mps_var.d*u,min(size(U,2)/length(d_vec)*mps_var.d,v));
                    U_temp(dd_vec,1:min(size(U,2),v))=U(:,1:min(size(U,2),v));
                    
                    cell_vec=u*ones(mps_var.d,1);
                    Q(:,NR_iter)=mat2cell(U_temp,cell_vec,min(size(U,2)/length(d_vec)*mps_var.d,v));
                    
                    Sv_temp=zeros(size(U_temp,2),v);
                    Sv_temp(1:min(size(U,2),v),:)=Sv(1:min(size(U,2),v),:);

                    S{1,NR_iter}=Sv_temp;
                end
            end
            S(cellfun(@isempty,S)) = {zeros(v,v)};
            if strcmp(Norm,'true')
                norm_val=0;
                for NR_iter=1:NR
                    norm_val=norm_val+S{1,NR_iter}*S{1,NR_iter}';
                end
                
                norm_val=norm_val(1);
                S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
            end
            Q(cellfun(@isempty,Q)) = {zeros(u,v)};
            mps_var.data{m}=Q;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
        end
        
        function K_Eff = Create_Eff_K_Right(mps_var,H_Eff,m)
            
            % H_Eff == %(u' NR' d' v'),(u NR d v)
            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            d_cell=cell(1,NR);
            Vec=cell(1,NR);
            List=zeros(mps_var.d,NR);
            for NR_iter = 1:NR
                d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                Vec(d_cell{1,NR_iter},NR_iter)=mps_var.data{m}(d_cell{1,NR_iter},NR_iter);
                List(d_cell{1,NR_iter},NR_iter)=ones(length(d_cell{1,NR_iter}),1);
            end
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            % Vec == d,NR
            Vec = reshape(Vec,mps_var.d*NR,1); % (NR d),1
            Vec = cell2mat(Vec(Ind,1)); % (NR d u),v
            
            Vec = reshape(Vec,u,nnz(Len),v); %u,(NR d),v
            Vec = permute(Vec,[3,2,1]); %v,(NR d),u
            Vec = reshape(Vec,v,u*nnz(Len)); %v1,(u NR d)
            
            H_temp = reshape(H_Eff,v,nnz(Len)*u,v,nnz(Len)*u); %v',(u' NR' d'),v,(u NR d)
            
            K_temp=ncon({H_temp,Vec,conj(Vec)},{[-3 1 -1 2],[-2 2],[-4 1]}); % v,v1,v',v1'
            K_temp=permute(reshape(K_temp,v^2,v^2),[2 1]); %(v1' v'),(v1 v)
            K_temp(abs(K_temp)<mps_var.zero_thres)=0;
            
%             K_temp=eye(v^2,v^2);

            K_Eff=cell(1,NR);
            K_Eff(cellfun(@isempty,K_Eff)) = {K_temp};
            
            
%             kappa = size(R,2);
%             NR = size(R{1},2);
%             NL = size(L{1},2);
% 
%             K_Eff=cell(1,NR);
%             for num = 1:min(NL,NR)
%                 for k=1:kappa
%                     if ~isempty(K_Eff{1,num-NL+mps_var.N+1})
%                         K_Eff{1,num-NL+mps_var.N+1}=K_Eff{1,num-NL+mps_var.N+1}+kron(R{k}{num-NL+mps_var.N+1},transpose(L{k}{NL-num+1}));
%                     else
%                         K_Eff{1,num-NL+mps_var.N+1} = kron(R{k}{num-NL+mps_var.N+1},transpose(L{k}{NL-num+1}));
%                     end
%                 end
%             end
%             u=size(R{1}{1},1);
%             v=size(L{1}{1},1);
%             K_Eff(cellfun(@isempty,K_Eff)) = {zeros(u*v,u*v)};
        end
        
        function mps_var = Apply_K_Eff_Right(mps_var,S,K_Eff,dt,m,Norm)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Time evolve Singular Values, S           
            S = cellfun(@(y,z) Apply_Cell_Expv(y,z,dt), K_Eff,S, 'UniformOutput', false);
            
            if strcmp(Norm,'true')
                NR=size(S,2);
                norm_val=0;
                for NR_iter=1:NR
%                     norm_val=norm_val+S{1,NR_iter}*S{1,NR_iter}';
                    norm_val=norm_val + sum(diag(S{1,NR_iter}'*S{1,NR_iter}));
                end
                
                norm_val=norm_val(1);
                S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Multiply Singular Values to next tensor            
            [u,~]=size(mps_var.data{m+1}{1});
            [~,NR]=size(mps_var.data{m+1});
            
            temp_E = cell(mps_var.d,NR);
            for k = 1:mps_var.d
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2 <= NR
                        temp_E{k,k_2} = S{1,mps_var.N+1-track+1};
                    end
                    track=track+1;
                end
            end
            temp_E(cellfun(@isempty,temp_E)) = {zeros(u,u)};
            
            mps_var.data{m+1} = cellfun(@(y,z) y*z, temp_E,mps_var.data{m+1}, 'UniformOutput', false);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            function Out = Apply_Cell_Expv(y,z,dt)
                % Do y*z
                % y==(a,b),(a,b)
                % z==a,b
                
                [alpha,beta]=size(z);
                
                vec = permute(z, [2,1]); %b,a
                vec = reshape(vec,alpha*beta,1); %(a,b)
                
                if ~(all(vec)==0)
                    
                    vec = expv(1i*dt,y,vec);
                    
%                     vec = vec + 1i*dt*y*vec;
                    
                    
                    %                 U_T = expm(1i*dt*y);
                    %                 vec = U_T*vec;
                    %
                    %                 check=abs(U_T*U_T'-eye(size(U_T)))>1e-10;
                    %                 nnz(check)
                    
                    Out = permute(reshape(vec,beta,alpha),[2 1]);
                else
                    Out=z;
                end
                
            end
            
        end
        
        function [mps_var,S] = Apply_H_Eff_Left(mps_var,H_Eff,dt,m,Norm)

            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            d_cell=cell(1,NR);
            Vec=cell(1,NR);
            List=zeros(mps_var.d,NR);
            for NR_iter = 1:NR
                d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                Vec(d_cell{1,NR_iter},NR_iter)=mps_var.data{m}(d_cell{1,NR_iter},NR_iter);
                List(d_cell{1,NR_iter},NR_iter)=ones(length(d_cell{1,NR_iter}),1);
            end
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            % Vec == d,NR
            Vec = reshape(Vec,mps_var.d*NR,1); % (NR d),1
            Vec = cell2mat(Vec(Ind,1)); % (NR d u),v
                
            Vec = reshape(Vec,u,nnz(Len),v); %u,(NR d),v
            Vec = permute(Vec,[3,2,1]); %v,(NR d),u
            Vec = reshape(Vec,u*nnz(Len)*v,1); %(u NR d v),1
            
            Vec=expv(-1i*dt,H_Eff,Vec); %(u NR' d' v),1
            
%             U_Eff=expm(-1i*H_Eff*dt); % (u NR' d' v),(u NR d v)
%             Vec=U_Eff*Vec; %(u NR' d' v),1

%             Vec = Vec - 1i*dt*H_Eff*Vec;

            Vec = reshape(Vec,v,nnz(Len),u); % v,(NR d),u
            Vec = permute(Vec,[3,2,1]); % u,(NR d),v
            Vec = reshape(Vec,u*nnz(Len),v); % (NR d u),v
            
            row_vec=u*ones(nnz(Len),1);
            Vec = mat2cell(Vec,row_vec,v); % (NR d),1
            
            New_Vec=cell(mps_var.d*NR,1);
            New_Vec(Ind,1)=Vec;
            New_Vec=reshape(New_Vec,mps_var.d,NR);
            New_Vec(cellfun(@isempty,New_Vec)) = {zeros(u,v)};
            mps_var.data{m}=New_Vec;
            
            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Convert mth Tensor to N_Left Labelling
            temp_T = cell(mps_var.d,NR);
            temp_IT = cell(mps_var.d,NL);
            expect = max(NR,NL):-1:1;
            for k = 1:mps_var.d
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2>mps_var.N+1-NR-k+1 && k_2<=NL
                        temp_T(k,k_2) = {expect(1,max(NR,NL)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_var.N+1-NR-k+2):min(NL,(mps_var.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_var.data{m};
            
            R_mm = cell(mps_var.d,NL);
            for alpha = 1:mps_var.d
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do QR decomposition
            Q=cell(mps_var.d,NL);
            S=cell(1,NL);
            if m~=1
                for NL_iter=1:NL
                    d_vec=max(1,mps_var.N+3-NL_iter-NR):min(mps_var.d,mps_var.N+2-NL_iter);
                    if ~isempty(d_vec)
                        Mat = cell2mat(R_mm(d_vec,NL_iter)); % (d u),v
                        Mat = reshape(Mat,u,length(d_vec),v); % u,d,v
                        Mat = permute(Mat,[1 3 2]); % u,v,d
                        Mat = reshape(Mat,u,length(d_vec)*v); % u,(d v)
                        Mat = transpose(Mat); %(d v),u
                        
                        [V,R]=qr(Mat); % V=(d v),(d v)
                        
                        Sv=transpose(R);
                        V=transpose(V);
                        
                        Sv(abs(Sv)<mps_var.zero_thres)=0;
                        V(abs(V)<mps_var.zero_thres)=0;
                        
                        dd_vec = ((d_vec(1)-1)*v+1):(d_vec(end)-1)*v+v;
                        
                        V_temp=zeros(min(size(V,1)/length(d_vec)*mps_var.d,u),mps_var.d*v);
                        V_temp(1:min(size(V,1),u),dd_vec)=V(1:min(size(V,1),u),:); % u,(v d)
                        
                        Mat = reshape(V_temp,u,v,mps_var.d); % u,v,d
                        Mat = permute(Mat,[1 3 2]); %u,d,v
                        Mat = reshape(Mat,u*mps_var.d,v);% (d u),v
                        
                        cell_vec=u*ones(mps_var.d,1);
                        Q(:,NL_iter)=mat2cell(Mat,cell_vec,v);
                        
                        Sv_temp=zeros(u,size(V_temp,1));
                        Sv_temp(:,1:min(size(V,1),u))=Sv(:,1:min(size(V,1),u));
                        
                        S{1,NL_iter}=Sv_temp;
                    end
                end
                S(cellfun(@isempty,S)) = {zeros(u,u)};
                if strcmp(Norm,'true')
                    norm_val=0;
                    for NL_iter=1:NL
                        norm_val=norm_val+S{1,NL_iter}'*S{1,NL_iter};
                    end
                    
                    norm_val=norm_val(1);
                    S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
                end
                Q(cellfun(@isempty,Q)) = {zeros(u,v)};
                R_mm=Q;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Convert mth tensor back into N_Right labelling
            temp = cell(mps_var.d,NR);
            for alpha = 1:mps_var.d
                temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(u,v)};
            
            mps_var.data{m}=temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function K_Eff = Create_Eff_K_Left(mps_var,H_Eff,m)
            
            [u,v]=size(mps_var.data{m}{1});
            [~,NR]=size(mps_var.data{m});
            NL=min(mps_var.N+1,(mps_var.d-1)*(m-1)+1);
            
            d_cell=cell(1,NR);
            Vec=cell(1,NR);
            List=zeros(mps_var.d,NR);
            for NR_iter = 1:NR
                d_cell{1,NR_iter}=max(1,mps_var.N+3-NR_iter-NL):min(mps_var.d,mps_var.N+2-NR_iter);
                Vec(d_cell{1,NR_iter},NR_iter)=mps_var.data{m}(d_cell{1,NR_iter},NR_iter);
                List(d_cell{1,NR_iter},NR_iter)=ones(length(d_cell{1,NR_iter}),1);
            end
            
            % List == d,NR
            Len=reshape(List,mps_var.d*NR,1); %(NR d),1
            Ind=find(Len~=0);
            
            % Vec == d,NR
            Vec = reshape(Vec,mps_var.d*NR,1); % (NR d),1
            Vec = cell2mat(Vec(Ind,1)); % (NR d u),v
            
            Vec = reshape(Vec,u,nnz(Len),v); %u,(NR d),v
            Vec = permute(Vec,[3,2,1]); %v,(NR d),u
            Vec = reshape(Vec,v*nnz(Len),u); %(NR d v),u1
            
            % H_Eff == %(u' NR' d' v'),(u NR d v)
            H_temp = reshape(H_Eff,v*nnz(Len),u,v*nnz(Len),u); %(NR' d' v'),u',(NR d v),u
            
            K_temp=ncon({H_temp,Vec,conj(Vec)},{[1 -4 2 -2],[2 -1],[1 -3]}); % u1,u,u1',u'
            K_temp=permute(reshape(K_temp,u^2,u^2),[2 1]); %(u' u1'),(u u1)
            K_temp(abs(K_temp)<mps_var.zero_thres)=0;
            
%             K_temp=zeros(u^2,u^2);
            
            K_Eff=cell(1,NL);
            K_Eff(cellfun(@isempty,K_Eff)) = {K_temp};
            
%             kappa = size(R,2);
%             NR = size(R{1},2);
%             NL = size(L{1},2);
%             
%             K_Eff=cell(1,NL);
%             for k=1:kappa
%                 for num = 1:min(NL,NR)
%                     if ~isempty(K_Eff{1,NL-num+1})
%                         K_Eff{1,NL-num+1}=K_Eff{1,NL-num+1}+kron(R{k}{num-NL+mps_var.N+1},transpose(L{k}{NL-num+1}));
%                     else
%                         K_Eff{1,NL-num+1} = kron(R{k}{num-NL+mps_var.N+1},transpose(L{k}{NL-num+1}));
%                     end
%                 end
%             end
%             u=size(R{1}{1},1);
%             v=size(L{1}{1},1);
%             K_Eff(cellfun(@isempty,K_Eff)) = {zeros(u*v,u*v)};
            
        end
        
        function mps_var = Apply_K_Eff_Left(mps_var,S,K_Eff,dt,m,Norm)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Time evolve Singular Values, S           
            S = cellfun(@(y,z) Apply_Cell_Expv(y,z,dt), K_Eff,S, 'UniformOutput', false);
            
            if strcmp(Norm,'true')
                NL=size(S,2);
                norm_val=0;
                for NL_iter=1:NL
%                     norm_val=norm_val+S{1,NL_iter}'*S{1,NL_iter};
                    norm_val=norm_val + sum(diag(S{1,NL_iter}'*S{1,NL_iter}));
                end
                
                norm_val=norm_val(1);
                S=cellfun(@(x) x/sqrt(diag(norm_val)),S, 'UniformOutput', false);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Convert (m-1)th tensor into N_Left labelling
            [dd,NR] = size(mps_var.data{m-1});
            Nl = min(mps_var.N+1,(mps_var.d-1)*(m-2)+1);

            temp_T = cell(dd,NR);
            temp_IT = cell(dd,Nl);
            expect_temp = max(NR,Nl):-1:1;
            for k = 1:dd
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2>mps_var.N+1-NR-k+1 && k_2<=Nl
                        temp_T(k,k_2) = {expect_temp(1,max(NR,Nl)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_var.N+1-NR-k+2):min(Nl,(mps_var.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_var.data{m-1};
            
            R_mm = cell(dd,Nl);
            for a = 1:dd
                R_mm(a,cell2mat(temp_IT(a,:))) = temp(a,cell2mat(temp_T(a,:)));
            end
            [nrows, ncols] = cellfun(@size, temp);
            u = max(max(nrows));
            v = max(max(ncols));
            R_mm(cellfun(@isempty,R_mm)) = {zeros(u,v)};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Multiply Singular Values to next tensor
            temp_E = cell(mps_var.d,Nl);
            for k = 1:mps_var.d
                track=1;
                for k_2 = (mps_var.N+1-k+1):-1:1
                    if k_2 <= Nl
                        temp_E{k,k_2} = S{1,mps_var.N+1-track+1};
                    end
                    track=track+1;
                end
            end
            temp_E(cellfun(@isempty,temp_E)) = {zeros(v,v)};
            
            R_mm = cellfun(@(y,z) y*z, R_mm,temp_E, 'UniformOutput', false);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Convert (m-1)th tensor back into N_Right labelling
            temp = cell(mps_var.d,NR);
            for a = 1:mps_var.d
                temp(a,cell2mat(temp_T(a,:))) = R_mm(a,cell2mat(temp_IT(a,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(size(mps_var.data{m-1}{1}))};
            
            mps_var.data{m-1}=temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            
            function Out = Apply_Cell_Expv(y,z,dt)
                % Do y*z
                % y==(a,b),(a,b)
                % z==a,b
                
                [alpha,beta]=size(z);
                
                vec = permute(z, [2,1]); %b,a
                vec = reshape(vec,alpha*beta,1); %(a,b)
                
                if ~(all(vec)==0)
                    vec = expv(1i*dt,y,vec);
                    
%                     vec = vec + 1i*dt*y*vec;

                    
                    %                 U_T = expm(1i*dt*y);
                    %                 vec = U_T*vec;
                    
                    Out = permute(reshape(vec,beta,alpha),[2,1]);
                else
                    Out=z;
                end
                
            end
            
        end
        
        function E = Get_Energy(mps_var,H)
            
            M = size(mps_var.data,2);
            
            R_INIT_temp = cell(1,mps_var.N+1);
            R_INIT_temp(cellfun(@isempty,R_INIT_temp)) = {1};
            R_INIT{1}=R_INIT_temp;
            
            R=cell(1,M);
            for m = 1:M
                if m ==1
                    R{m} = mps_var.Add_Right_Eff_Ham(H,R_INIT,m);
                else
                    R{m} = mps_var.Add_Right_Eff_Ham(H,R{m-1},m);
                end
            end
            E=cell2mat(R{M}{1});
        end
        
        function Var = Get_Var(mps_var,H)
            
            E=mps_var.Get_Energy(H);
            
            M = size(mps_var.data,2);
            
            R_INIT_temp = cell(1,mps_var.N+1);
            R_INIT_temp(cellfun(@isempty,R_INIT_temp)) = {1};
            R_INIT{1}=R_INIT_temp;
            
            R=cell(1,M);
            H_sq.data=cell(1,M);
            H_sq.N_track=cell(1,M);
            H_sq.d=H.d;
            for m = 1:M
                [alpha,beta]=size(H.data{m});
                H_sq.data{m}=cell(alpha^2,beta^2);
                
                for a1 = 1:alpha
                    for a2 = 1:alpha
                        for b1 = 1:beta
                            for b2 = 1:beta
                                if ~all(all(H.data{m}{a1,b1}==0)) && ~all(all(H.data{m}{a2,b2}==0))
                                    H_sq.N_track{m}{(a1-1)*alpha+a2,(b1-1)*beta+b2}=H.N_track{m}{a1,b1}+H.N_track{m}{a2,b2};
                                    H_sq.data{m}{(a1-1)*alpha+a2,(b1-1)*beta+b2}=H.data{m}{a1,b1}*H.data{m}{a2,b2};
                                else
                                    H_sq.data{m}{(a1-1)*alpha+a2,(b1-1)*beta+b2}=0;
                                end
                            end
                        end
                    end
                end
                
                if m ==1
                    R{m} = mps_var.Add_Right_Eff_Ham(H_sq,R_INIT,m);
                else
                    R{m} = mps_var.Add_Right_Eff_Ham(H_sq,R{m-1},m);
                end
            end
            Var=cell2mat(R{M}{1});
            
            Var=(Var-E^2)/Var;
        end
        
        function E = Site_Site_Particle_Corr(mps_var,site1,site2)
            % Find two site correlation function for N bosons in M site lattice
            %
            % b applied to site1
            % b_dag applied to site2
            % Then expectation value
            
            M = size(mps_var.data,2);
            O = mps_var.data;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Apply b to site1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,g] = size(mps_var.data{site1});
            [alpha,gamma] = size(mps_var.data{site1}{end,1});
            Nl = min((mps_var.d-1)*(site1-1)+1,mps_var.N+1);
            
            temp = cell(a,g);
            temp(cellfun(@isempty,temp)) = {zeros(alpha,gamma)};
            for l = 1:mps_var.d
                for N_R = max(1,(g-Nl-l+2)):max(1,min(g,mps_var.N+1-l+1))
                    for count2 = (-mps_var.d+l):(l-1)
                        temp{l-count2,N_R} = temp{l-count2,N_R}+mps_var.b(l-count2,l)*O{site1}{l,N_R};
                    end
                end
            end
            O{site1} = temp;
            
            for m_count = 1:(site1-1)
                [m_a,m_g]=size(O{m_count});
                m_temp = cell(m_a,m_g);
                for m1 = 1:m_g-1
                    m_temp(:,m1) = O{m_count}(:,m1+1);
                end
                m_temp(cellfun(@isempty,m_temp)) = {zeros(size(O{m_count}{2,2}))};
                O{m_count} = m_temp;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Apply b_dag to site2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,g] = size(mps_var.data{site2});
            [alpha,gamma] = size(mps_var.data{site2}{end,1});
            
            temp = cell(a,g);
            temp(cellfun(@isempty,temp)) = {zeros(alpha,gamma)};
            for l = 1:mps_var.d
                for N_R = 1:max(1,min(g,mps_var.N+1-l+1))
                    for count2 = (-mps_var.d+l):(l-1)
                        temp{l-count2,N_R} = temp{l-count2,N_R}+mps_var.b_dag(l-count2,l)*O{site2}{l,N_R};
                    end
                end
            end
            O{site2} = temp;
            
            for m_count = 1:(site2-1)
                [m_a,m_g]=size(O{m_count});
                m_temp = cell(m_a,m_g);
                for m1 = 2:m_g
                    m_temp(:,m1) = O{m_count}(:,m1-1);
                end
                m_temp(cellfun(@isempty,m_temp)) = {zeros(size(O{m_count}{2,2}))};
                O{m_count} = m_temp;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Expectation Value
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            expect=cell(1,mps_var.N+1);
            expect(cellfun(@isempty,expect)) = {1};
            for n = 1:M
                [a,g] = size(O{n});
                [alpha,gamma] = size(O{n}{end,1});
                
                temp_E = cell(a,g);
                for k = 1:a
                    track=1;
                    for k_2 = (mps_var.N+1-k+1):-1:1
                        if k_2 <= g
                            temp_E{k,k_2} = expect{1,mps_var.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                
                temp = cellfun(@(x,y,z) x'*y*z, mps_var.data{n},temp_E,O{n}, 'UniformOutput', false);
                
                expect=cell(1,g);
                expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                
                for i=1:a
                    expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                end
                if n<M
                    if size(expect,2)<mps_var.d
                        expect(:,size(expect,2)+1:mps_var.d)={zeros(gamma,gamma)};
                    end
                end
                
            end
            E=abs(cell2mat(expect));
        end
        
        function Corr = plot_Corr(mps_var,Corr)
            % Plot the results of the correlation function caluclation
            %
            % Either plots a 1D graph of the particle number distribution,
            % or a 2D surf plot of the site-site correlation function
            % depending on the dimension of the input matrix.
            
            M = max(size(Corr));
            dim = sum(size(Corr)>1);
            
            if dim == 2
                figure;
                surf(1:M,1:M,Corr)
                xlabel('Site, i','fontsize',20);
                ylabel('Site, J','fontsize',20);
                zlabel('Occupation Number','fontsize',20);
                set(gca,'fontsize',15)
                axis([1 M 1 M 0 max(max(Corr))]);
            elseif dim ==1
                figure;
                plot(1:M,Corr)
                xlabel('Site','fontsize',20);
                ylabel('Occupation Number','fontsize',20);
                set(gca,'fontsize',15)
                axis([1 M 0 max(max(Corr))]);
            else
                error('Unsuported Input!');
            end
        end

        function W_TEBD = Calc_State_Vector(mps_var,B)
            % Produce exact diagonalisation state vector from MPS
            
            W_TEBD=zeros(size(B,1),1);
            for n = 1:size(B,1)
                for m = 1:size(B,2)
                    d_iter = B(n,m)+1;
                    if d_iter<=mps_var.d
                        if m ~= size(B,2)
                            N_R = sum(B(n,(m+1):end))+1;
                        else
                            N_R=1;
                        end
                        if N_R <= max(mps_var.d*(size(B,2)-m),1)
                            if m==1
                                W=mps_var.data{m}{d_iter,N_R};
                            else
                                W=W*mps_var.data{m}{d_iter,N_R};
                            end
                        else
                            if m==1
                                W=zeros(size(mps_var.data{m}{1,1}));
                            else
                                W=W*zeros(size(mps_var.data{m}{1,1}));
                            end
                        end
                    else
                        if m==1
                            W=zeros(size(mps_var.data{m}{1,1}));
                        else
                            W=W*zeros(size(mps_var.data{m}{1,1}));
                        end
                    end
                end
                W_TEBD(n)=W;
            end
        end
    end
end

