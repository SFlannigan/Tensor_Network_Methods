classdef mps_cpn
    % Matrix product state formulism that conserves particle number for the
    % Bose-Hubbard model.
    %
    % We use the cell functionality for each tensor at every site. The position
    % within each 2D cell dictates the local dmension and the number of particles
    % to the right (or left).
    %
    % To do:
    %       I think that I need to include more detailed comments on the
    %       time evolution procedures in the new formalism. I think I
    %       should produce a set of notes that goes into the details.
    %

    
    
    properties (SetAccess=public)
    end
    
    
    
    properties (SetAccess=protected)
        
        % Maximum bond dimension.
        Dmax
        
        % Matrix data
        %
        % Structure array of Cells. Each element of structure corresponds
        % to each lattice site and is a 2D Cell array - one dimension is
        % the local dimension and the other is the number of particles to
        % the right (or left).
        %
        % !!! but you do not say anything below on how you change this
        % second dimension, i.e. when it count particle to the left and
        % when to the right. I assume you always count them to the right.
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
        
        % Local creation and Annihilation operators
        b_dag
        b
        
        % Vector containing Suzuki-Trotter time steps
        ST_Order
        
        % Vector containing exponential roots. Used with Runga-Kutta
        RK_Roots
        
        % Matrix Containing [lattice position,Number of applications] to 
        % implement SWAP gates
        Swap_Mat
        
        % Vector containing chosen positions of initial particles
        P_Positions
        
    end
    
    
    
    
    methods
        
        function mps_cpn=mps_cpn(M,Dmax,N,N_max)
            % Class constructor
            
            mps_cpn.Dmax=Dmax;
            mps_cpn.data=cell(1,M);
            mps_cpn.zero_thres=1e-15;
            mps_cpn.N=N;
            mps_cpn.N_max=N_max;
            mps_cpn.d=N_max+1;  % !!! in furute it should be variable as well
            mps_cpn.b = sqrt(diag(1:N_max,1));
            mps_cpn.b_dag = sqrt(diag(1:N_max,-1));
            
        end
        
        function mps_cpn=set_Particle_Position(mps_cpn,Vec)
            % Vec == vector containing the positions of the N particles.
            
            if numel(Vec)==mps_cpn.N
                mps_cpn.P_Positions=Vec;
            else
                display(['Particle Position Vector is the Wrong Length.' ...
                    ' Using Random Distribution.']);
            end
        end
        
        function mps_cpn=set_rand_product_state(mps_cpn)
            % Set initial random product state.
            %
            % N == Number of particles.
            % Nmax == Maximum number of particles allowed on a single site.
            
            M = size(mps_cpn.data,2);
            
            if M*mps_cpn.N_max < mps_cpn.N
                error('Too many particles!');
            end
            
            % Generate N random integers between 1 and M.
            % These correspond to an intial random location for each particle.
            if isempty(mps_cpn.P_Positions)
                R = randi(M,1,mps_cpn.N);
            else
                R = mps_cpn.P_Positions;
            end
            R = sort(R);
            
            % !!! this block of 40ish lines is bad because it needs either:
            % 1. names of variables that are selfexplanatory
            % 2. or just describe here what is being done until the
            % next empty line
            A = zeros(M,1);
            C=0;
            check=0;
            for m = 1:M
                A(m) = C;
                C=0;
                for n = 1:mps_cpn.N
                    if R(n) == m
                        A(m) = A(m)+1; % Total number of particles on site m
                    end
                end
                % If the intial random locations give more particles than are
                % allowed (dictated by N_max) then we carry the extra particles
                % over to the next site.
                if A(m)>mps_cpn.N_max
                    C=A(m)-mps_cpn.N_max;
                    A(m)=mps_cpn.N_max;
                    if m==M
                        % We have reached the end of the lattice and we cannot
                        % carry over the extra particles.
                        % We have C extra particles that do not have a lattice site.
                        % We need to apply an extra algorithm (below) to sort them into
                        % sites that have spaces.
                        check=1;
                    end
                end
                
                if m == 1
                    B = cell(mps_cpn.d,mps_cpn.N+1);
                    B{A(m)+1,mps_cpn.N-A(m)+1}=1;
                elseif m==M
                    B=cell(mps_cpn.d,1);
                    B{A(m)+1,1} = 1;
                else
                    xr = min(mps_cpn.N+1,mps_cpn.N_max*(M-m)+1);
                    B = cell(mps_cpn.d,xr);
                    B{A(m)+1,mps_cpn.N+1-sum(A(1:m))}=1;
                end
                B(cellfun(@isempty,B)) = {0};
                mps_cpn.data{m} = B;
                
            end
            
            % Apply extra algorithm to sort extra particles into free lattice sites
            if check==1
                m=1;
                % We alter the values of A(m) to ensure the all particles have a
                % lattice.
                while C ~=0
                    D=A(m)+C;
                    if D > mps_cpn.N_max
                        C = D-mps_cpn.N_max;
                        D=mps_cpn.N_max;
                    else
                        C=0;
                    end
                    A(m)=D;
                    m=m+1;
                end
                
                % With the new A(m) we construct the MPS
                for m = 1:M
                    if m == 1
                        B = cell(mps_cpn.d,mps_cpn.N+1);
                        B{A(m)+1,mps_cpn.N-A(m)+1}=1;
                    elseif m==M
                        B=cell(mps_cpn.d,1);
                        B{A(m)+1,1} = 1;
                    else
                        xr = min(mps_cpn.N+1,mps_cpn.N_max*(M-m)+1);
                        B = cell(mps_cpn.d,xr);
                        B{A(m)+1,mps_cpn.N+1-sum(A(1:m))}=1;
                    end
                    B(cellfun(@isempty,B)) = {0};
                    mps_cpn.data{m} = B;
                end
            end
            
            
%             % !!! how about this way? let's discuss
%             % #####################
%             if isempty(mps_cpn.P_Positions)
%                 R=randi(M,1,mps_cpn.N);
%             else
%                 R=mps_cpn.P_Positions;
%             end
%             if max(R>M) 
%                 error('wrong location of particles'); 
%             end
%             
%             A=zeros(M,1);
%             for iM=1:length(R)
%                A(R(iM))=A(R(iM))+1;
%             end
%             
%             particles_in_lattice=0; % to the right
%             for iM=M:-1:M
%                 mps_cpn.data{iM}=cell(mps_cpn.d,mps_cpn.N+1);
%                 mps_cpn.data{iM}{A(iM)+1,particles_in_lattice+1}=1;
%                 particles_in_lattice=particles_in_lattice+A(iM);
%             end
%             
%             if particles_in_lattice~=mps_cpn.N
%                 error('sanity check fail')
%             end
%             % #####################
        end
        
        function mps_cpn=set_bond_dim(mps_cpn,trunc)
            % Change the class bond dimension
            
            mps_cpn.Dmax=trunc;
            
        end
        
        function [mps_cpn,Total_error] = Canonicalisation_2s(mps_cpn,Sweep_Direction)
            % Canonicalise the MPS
            %
            % Sweep_Direction == 'L-R' - Right Canonicalised
            % Sweep_Direction == 'R-L' - Left Canonicalised
            
            M = size(mps_cpn.data,2);
            
            if strcmp(Sweep_Direction,'L-R')
                start=1;
                End_S = M-1;
                sign = 1;
            elseif strcmp(Sweep_Direction,'R-L')
                start = (M-1);
                End_S = 1;
                sign = -1;
            end
            
            Total_error=0;
            for m=start:sign:End_S
                [a_1,g_1] = size(mps_cpn.data{m});
                [a_2,g_2] = size(mps_cpn.data{m+1});
                
                Nl_2 = min(mps_cpn.N+1,(mps_cpn.d-1)*m+1);
                %%%%%%%%%%%%%%%%%%
                %Convert mps_cpn.data{m+1} to N_Right labelling (R_mm)
                %%%%%%%%%%%%%%%%%%
                temp_T = cell(a_2,g_2);
                temp_IT = cell(a_2,Nl_2);
                expect = max(g_2,Nl_2):-1:1;
                for k = 1:a_2
                    track=1;
                    for k_2 = (mps_cpn.N+1-k+1):-1:1
                        if k_2>mps_cpn.N+1-g_2-k+1 && k_2<=Nl_2
                            temp_T(k,k_2) = {expect(1,max(g_2,Nl_2)-track+1)};
                        end
                        track=track+1;
                    end
                    for k_2 = max(1,mps_cpn.N+1-g_2-k+2):min(Nl_2,(mps_cpn.N+2-k))
                        temp_IT(k,k_2) = {k_2};
                    end
                end
                
                temp = mps_cpn.data{m+1};
                
                R_mm = cell(a_2,Nl_2);
                for alpha = 1:a_2
                    R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
                end
                [nrows, ncols] = cellfun(@size, temp);
                beta = max(max(nrows));
                g = max(max(ncols));
                R_mm(cellfun(@isempty,R_mm)) = {zeros(beta,g)};
                
                %%%%%%%%%%%%%%%%%%
                %Collect entries of (L_m X R_mm) to apply the SVD's
                %%%%%%%%%%%%%%%%%%
                [nrows, ~] = cellfun(@size, mps_cpn.data{m});
                a = max(max(nrows));
                L_m = mps_cpn.data{m}(:,end:-1:1);
                
                d_1_s = (1:mps_cpn.N+1);
                d_1_s(d_1_s>mps_cpn.d)=mps_cpn.d;
                
                d_2_s = (mps_cpn.N+1):-1:1;
                d_2_s(d_2_s>mps_cpn.d)=mps_cpn.d;
                
                temp = cell(a_1*a_2,length((1+mps_cpn.N+1-g_1):Nl_2));
                track=1;
                for gamma = (1+mps_cpn.N+1-g_1):Nl_2
                    for counter = 1:d_1_s(gamma)
                        temp(((counter-1)*mps_cpn.d+1):(((counter-1)*mps_cpn.d)+d_2_s(gamma)),track) = {zeros(a,g)};
                    end
                    track=track+1;
                end
                               
                track=1;
                for gamma = (mps_cpn.N+1-g_1+1):Nl_2
                    for d_1 = 1:d_1_s(gamma)
                        for d_2 = 1:d_2_s(gamma)
                            if d_1+d_2-2<=mps_cpn.N
                                if ((d_1-1)*a_1+d_2)<=a_1*a_2 && track<=size(temp,2) && ~isempty(temp{(d_1-1)*a_1+d_2,track})
                                    temp{(d_1-1)*a_1+d_2,track} = temp{(d_1-1)*a_1+d_2,track} + L_m{d_1,track}*R_mm{d_2,gamma};
                                end
                            end
                        end
                    end
                    track=track+1;
                end
                
                temp_2 = cellfun(@(x) mps_cpn.func(x),temp,'UniformOutput', false);
                
                d_1_s=d_1_s((mps_cpn.N+1-g_1+1):Nl_2);
                d_2_s=d_2_s((mps_cpn.N+1-g_1+1):Nl_2);
                d_len = d_1_s.*d_2_s;
                
                Norm = zeros(max(d_len),min(Nl_2,g_1));
                for count = 1:min(Nl_2,g_1)
                    Norm(1:d_len(count),count) = cell2mat(temp_2(:,count));
                end
                
                Norm = sum(sum(Norm));
                temp = cellfun(@(x) x/sqrt(Norm),temp,'UniformOutput', false);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Do SVD's
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                R_mm = cell(mps_cpn.d,(mps_cpn.N+1));
                L_m = cell(mps_cpn.d,g_1);
                [a,~] = size(mps_cpn.data{m}{1,1});
                [~,g] = size(mps_cpn.data{m+1}{1,1});
                r_min=1;
                Norm=0;
                for count2 = 1:min(Nl_2,g_1)
                    Theta = cell2mat(temp(:,count2));
                    d_1 = d_1_s(count2);
                    d_2 = d_2_s(count2);
                    
                    % Theta == (d_1 d_2 alpha),(gamma)
                    Theta = reshape(Theta,a,d_2,d_1,g); %(a),(l+1),(l),(g)
                    Theta = permute(Theta,[1 3 4 2]); %(a),(l),(g),(l+1)
                    Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g)
                    Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g)
                    [A,S,V]=svd(Theta);
                    
                    S((abs(S))/S(1,1)<mps_cpn.zero_thres)=0;
                    
                    r = nnz(S);
                    if r>r_min
                        r_min=r;
                    end
                    if r==0
                        r=1;
                        S = S(1:r,1:r);
                        Error=0;
                    elseif r>mps_cpn.Dmax
                        r=mps_cpn.Dmax;
                        if any(size(S)==r+1)
                            Error=sum(S(r+1:end,r+1:end).^2);
                        else
                            Error=sum(diag(S(r+1:end,r+1:end)).^2);
                        end
                        S = S(1:r,1:r);
                    else
                        Error=0;
                        S = S(1:r,1:r);
                    end
                    Norm = Norm+diag(S)'*diag(S);
                    if diag(S)'*diag(S) ~=0
                        if strcmp(Sweep_Direction,'L-R')
                            A = A(:,1:r);
                            V = S*V(:,1:r)';
                        else
                            A = A(:,1:r)*S;
                            V = V(:,1:r)';
                        end
                    else
                        if strcmp(Sweep_Direction,'L-R')
                            A = A(:,1:r);
                            V = S*V(:,1:r)';
                        else
                            A = A(:,1:r)*S;
                            V = V(:,1:r)';
                        end
                    end
                    Total_error=Total_error+Error;
                    
                    % A == (l a),(r)
                    A = reshape(A,a,d_1,r); %(a),(l),(r)
                    A = permute(A,[1 3 2]); %(a),(r),(l)
                    for count = 1:d_1
                        L_m(count,count2) = {A(:,:,count)};
                    end
                    
                    % V == (r),(l+1 g)
                    V = reshape(V,r,g,d_2); %(r),(g),(l+1)
                    for count = 1:d_2
                        R_mm(count,count2+(mps_cpn.N+1-g_1)) = {V(:,:,count)};
                    end
                end
                R_mm(cellfun(@isempty,R_mm)) = {0};
                L_m(cellfun(@isempty,L_m)) = {0};
                
                
                
                [nrows, ~] = cellfun(@size, L_m);
                a = max(max(nrows));
                [~, ncol] = cellfun(@size, R_mm);
                g = max(max(ncol));
                
                L_m = cellfun(@(x) mps_cpn.cell_resize_func(x,[a,r_min]),L_m,'UniformOutput', false);
                R_mm = cellfun(@(x) mps_cpn.cell_resize_func(x,[r_min,g]),R_mm,'UniformOutput', false);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Update First Site
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                mps_cpn.data{m}=L_m(:,end:-1:1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Update Second Site
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%
                %Convert mps_cpn.data{m+1} back to N_Left labelling
                %%%%%%%%%%%%%%%%%%
                
                temp = cell(a_2,g_2);
                for alpha = 1:a_2
                    temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
                end
                temp(cellfun(@isempty,temp)) = {zeros(r_min,g)};
                
                mps_cpn.data{m+1}=temp;
            end
            
        end
        
        function x = func(mps_cpn,x)
            % Function used to calculate normalisation of two site MPS
            % tensor.
            % Used by 'Canonicalisation_2s' and 'Time_Evolve'.
            %
            % !!! Create a directory with useful extra functions and put
            % this one there. This function does not even use mps_cpn.
            
            if ~isempty(x)
                x = sum(diag(x'*x));
            end
            
        end
        
        function Out = cell_resize_func(mps_cpn,X,r)
            % Function used to resize MPS
            % tensor.
            % Used by 'Canonicalisation_2s' and 'Time_Evolve'.
            %
            % !!! Create a directory with useful extra functions and put
            % this one there. This function does not even use mps_cpn.
            
            
            [a,b2]= size(X);
            Out=X;
            if a < r(1)
                Out(a+1:r(1),:) = zeros(r(1)-a,size(Out,2));
            end
            if b2 < r(2)
                Out(:,b2+1:r(2)) = zeros(size(Out,1),r(2)-b2);
            end
            
        end
        
        function E = Site_Site_Particle_Corr(mps_cpn,site1,site2)
            % Find two site correlation function for N bosons in M site lattice
            %
            % b applied to site1
            % b_dag applied to site2
            % Then expectation value
            
            M = size(mps_cpn.data,2);
            O = mps_cpn.data;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Apply b to site1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [a,g] = size(mps_cpn.data{site1});
            [alpha,gamma] = size(mps_cpn.data{site1}{end,1});
            Nl = min((mps_cpn.d-1)*(site1-1)+1,mps_cpn.N+1);
            
            temp = cell(a,g);
            temp(cellfun(@isempty,temp)) = {zeros(alpha,gamma)};
            for l = 1:mps_cpn.d
                for N_R = max(1,(g-Nl-l+2)):max(1,min(g,mps_cpn.N+1-l+1))
                    for count2 = (-mps_cpn.d+l):(l-1)
                        temp{l-count2,N_R} = temp{l-count2,N_R}+mps_cpn.b(l-count2,l)*O{site1}{l,N_R};
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
            [a,g] = size(mps_cpn.data{site2});
            [alpha,gamma] = size(mps_cpn.data{site2}{end,1});
            
            temp = cell(a,g);
            temp(cellfun(@isempty,temp)) = {zeros(alpha,gamma)};
            for l = 1:mps_cpn.d
                for N_R = 1:max(1,min(g,mps_cpn.N+1-l+1))
                    for count2 = (-mps_cpn.d+l):(l-1)
                        temp{l-count2,N_R} = temp{l-count2,N_R}+mps_cpn.b_dag(l-count2,l)*O{site2}{l,N_R};
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
            expect=cell(1,mps_cpn.N+1);
            expect(cellfun(@isempty,expect)) = {1};
            for n = 1:M
                [a,g] = size(O{n});
                [alpha,gamma] = size(O{n}{end,1});
                
                temp_E = cell(a,g);
                for k = 1:a
                    track=1;
                    for k_2 = (mps_cpn.N+1-k+1):-1:1
                        if k_2 <= g
                            temp_E{k,k_2} = expect{1,mps_cpn.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                
                temp = cellfun(@(x,y,z) x'*y*z, mps_cpn.data{n},temp_E,O{n}, 'UniformOutput', false);
                
                expect=cell(1,g);
                expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                
                for i=1:a
                    expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                end
                if n<M
                    if size(expect,2)<mps_cpn.d
                        expect(:,size(expect,2)+1:mps_cpn.d)={zeros(gamma,gamma)};
                    end
                end
                
            end
            E=abs(cell2mat(expect));
        end
        
        function N = Full_Norm(mps_cpn)
            % Find Overlap of MPS (Normalisation)
            
            M = size(mps_cpn.data,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Overlap
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            expect=cell(1,mps_cpn.N+1);
            expect(cellfun(@isempty,expect)) = {1};
            for n = 1:M
                [a,g] = size(mps_cpn.data{n});
                [alpha,gamma] = size(mps_cpn.data{n}{end,1});
                
                temp_E = cell(a,g);
                for k = 1:a
                    track=1;
                    for k_2 = (mps_cpn.N+1-k+1):-1:1
                        if k_2 <= g
                            temp_E{k,k_2} = expect{1,mps_cpn.N+1-track+1};
                        end
                        track=track+1;
                    end
                end
                temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                
                temp = cellfun(@(x,y,z) x'*y*z, mps_cpn.data{n},temp_E,mps_cpn.data{n}, 'UniformOutput', false);
                
                expect=cell(1,g);
                expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                
                for i=1:a
                    expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                end
                if n<M
                    if size(expect,2)<mps_cpn.d
                        expect(:,size(expect,2)+1:mps_cpn.d)={zeros(gamma,gamma)};
                    end
                end
            end
            N=cell2mat(expect);
        end
        
        function mps_cpn = set_Suzuki_Trotter_order(mps_cpn,order)
            % Set vecotr containing Suzuki-Trotter time steps for TEBD time
            % evolution
            if order ==  4
                % !!! You should also mark in which direction you go. Look
                % at Eqs. (6.43)-(6.45) in Johannes's PhD
                mps_cpn.ST_Order = [1/12,1/12,1/12,-1/6,1/12,0,1/12,0,1/12,0,1/12,1/12,1/12,1/12,0,1/12,0,1/12,0,1/12,-1/6,1/12,1/12,1/12];
            end
        end
        
        function H_2s = Loc_2s_Ham(mps_cpn,M,J,U,E)
            % Generate 2-site local Hamiltonian for TEBD algorithm.
 
            H_2s = cell(1,M-1);
            
            num = mps_cpn.b_dag*mps_cpn.b;
            id = eye(mps_cpn.d);
            for m = 1:(M-1)
                if m ==1
                    H_1 = U/2*num*(num - id) + E(m)*num;
                    H_2 = U/4*num*(num - id) + E(m+1)/2*num;
                elseif m==(M-1)
                    H_1 = U/4*num*(num - id) + E(m)/2*num;
                    H_2 = U/2*num*(num - id) + E(m+1)*num;
                else
                    H_1 = U/4*num*(num - id) + E(m)/2*num;
                    H_2 = U/4*num*(num - id) + E(m+1)/2*num;
                end
                H_2s{m} = kron(H_1,id)+kron(id,H_2)-J*kron(mps_cpn.b_dag,mps_cpn.b)-J*kron(mps_cpn.b,mps_cpn.b_dag);
            end
            
        end
        
        function [mps_cpn,Total_error] = TEBD_Local_2s_Gates(mps_cpn,dt,J,U,E)
            % TEBD time evolution. Advance MPS by one time step with a
            % local 2 site Hamiltonian.
            %
            % Simple 1D Lattice
            
            M = size(mps_cpn.data,2);
            
            H_2s = mps_cpn.Loc_2s_Ham(M,J,U,E);
            
            Total_error=0;
            Sweep_Direction = 'L-R';
            for t = 1:size(mps_cpn.ST_Order,2);
                
                tau = dt*mps_cpn.ST_Order(t);
                U_T = cell(size(H_2s));
                for m = 1:(M-1)
                    U_T{m} = expm(-1i*tau*H_2s{m});
                end
                
                [mps_cpn,Error] = mps_cpn.Time_Evolve(U_T,Sweep_Direction);
                Total_error = Total_error+Error;
                
                if strcmp(Sweep_Direction,'L-R')
                    Sweep_Direction = 'R-L';
                else
                    Sweep_Direction = 'L-R';
                end
            end
            
        end
        
        function [mps_cpn,Total_error] = Time_Evolve(mps_cpn,U,Sweep_Direction)
            % Apply two-site local operators to Entire MPS and truncate if
            % necessary
            
            M = size(mps_cpn.data,2);
            if strcmp(Sweep_Direction,'L-R')
                start=1;
                End_S = M-1;
                sign = 1;
            elseif strcmp(Sweep_Direction,'R-L')
                start = (M-1);
                End_S = 1;
                sign = -1;
            end
            
            Total_error=0;
            for m=start:sign:End_S
                [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U{m},m,Sweep_Direction);
                Total_error = Total_error+Error;
            end
            
        end
        
        function [mps_cpn,Total_error] = Loc_Op_2s_Apply(mps_cpn,U,m,Sweep_Direction)
            % Apply the 2-site local operator to sites m and m+1. Apply 
            % SVD and truncate if necessary.  
            
            [a_1,g_1] = size(mps_cpn.data{m});  % maybe better names
            [a_2,g_2] = size(mps_cpn.data{m+1});
            
            Nl_2 = min(mps_cpn.N+1,(mps_cpn.d-1)*m+1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the tensor for the two-sites so we can apply the operator.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} to N_Right labelling
            %%%%%%%%%%%%%%%%%%
            temp_T = cell(a_2,g_2);
            temp_IT = cell(a_2,Nl_2);
            expect = max(g_2,Nl_2):-1:1;
            for k = 1:a_2
                track=1;
                for k_2 = (mps_cpn.N+1-k+1):-1:1
                    if k_2>mps_cpn.N+1-g_2-k+1 && k_2<=Nl_2
                        temp_T(k,k_2) = {expect(1,max(g_2,Nl_2)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_cpn.N+1-g_2-k+2):min(Nl_2,(mps_cpn.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_cpn.data{m+1};
            
            R_mm = cell(a_2,Nl_2);
            for alpha = 1:a_2
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            [nrows, ncols] = cellfun(@size, temp);
            beta = max(max(nrows));
            g = max(max(ncols));
            R_mm(cellfun(@isempty,R_mm)) = {zeros(beta,g)};
            
            %%%%%%%%%%%%%%%%%%
            %Collect entries of (L_m X R_mm) to apply the SVD's
            %%%%%%%%%%%%%%%%%%
            [nrows, ~] = cellfun(@size, mps_cpn.data{m});
            a = max(max(nrows));
            L_m = mps_cpn.data{m}(:,end:-1:1);
            
            d_1_s = (1:mps_cpn.N+1);
            d_1_s(d_1_s>mps_cpn.d)=mps_cpn.d;
            
            d_2_s = (mps_cpn.N+1):-1:1;
            d_2_s(d_2_s>mps_cpn.d)=mps_cpn.d;
            
            % Make temp as a matrix!!!
            % Use cells only when you do not know the exact sizes of the
            % object. Everything below can be rewritten with matrix
            % multiplications, i.e. it will be a gazilion time faster.
            temp = cell(a_1*a_2,length((1+mps_cpn.N+1-g_1):Nl_2));
            track=1;
            for gamma = (1+mps_cpn.N+1-g_1):Nl_2
                for counter = 1:d_1_s(gamma)
                    temp(((counter-1)*mps_cpn.d+1):(((counter-1)*mps_cpn.d)+d_2_s(gamma)),track) = {zeros(a,g)};
                end
                track=track+1;
            end
            
            % think about d1_vec, d2_vec, NL_vec,...
            
            % !!! This should be vectorized, i.e. written with as few for-loops
            % as possible. Currently, it is the bottleneck of your code.
            % In order for me to understand it we need to discuss.
            track=1;
            for gamma = (mps_cpn.N+1-g_1+1):Nl_2
                for d_1 = 1:d_1_s(gamma)
                    for d_2 = 1:d_2_s(gamma)
                        if d_1+d_2-2<=mps_cpn.N
                            
                            % !!! double-check for indices
                            table=[];
                            for count = (-min(d_1,track)+1):(d_2-1)
                                if ((d_1-1)*a_1+d_2 + count*(mps_cpn.d-1))<=...
                                        a_1*a_2 && track+count<=size(temp,2) &&...
                                        ~isempty(temp{(d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),track+count})
                                    table=[table;(d_1-1)*a_1+d_2,count,track+count];
                                end
                            end                            
                            for iii = 1:size(table,1)
                                    temp{table(iii,1)+table(iii,2),table(iii,3)} =...
                                    temp{table(iii,1)+table(iii,2),table(iii,3)} ...
                                        +U(table(iii,1)+table(iii,2),table(iii,1))*L_m{d_1,track}*R_mm{d_2,gamma}; 
                            end                            
                        end
                    end
                end
                track=track+1;
            end
            
            temp_2 = cellfun(@(x) mps_cpn.func(x),temp,'UniformOutput', false);
            
            d_1_s=d_1_s((mps_cpn.N+1-g_1+1):Nl_2);
            d_2_s=d_2_s((mps_cpn.N+1-g_1+1):Nl_2);
            d_len = d_1_s.*d_2_s;
            
            Norm = zeros(max(d_len),min(Nl_2,g_1));
            for count = 1:min(Nl_2,g_1)
                Norm(1:d_len(count),count) = cell2mat(temp_2(:,count));
            end
            
            Norm = sum(sum(Norm));
            temp = cellfun(@(x) x/sqrt(Norm),temp,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do SVD's
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R_mm = cell(mps_cpn.d,(mps_cpn.N+1));
            L_m = cell(mps_cpn.d,g_1);
            [a,~] = size(mps_cpn.data{m}{1,1});
            [~,g] = size(mps_cpn.data{m+1}{1,1});
            r_min=1;
            Norm=0;
            Total_error=0;
            for count2 = 1:min(Nl_2,g_1)
                Theta = cell2mat(temp(:,count2));
                d_1 = d_1_s(count2);
                d_2 = d_2_s(count2);
                
                % Theta == (d_1 d_2 alpha),(gamma)
                Theta = reshape(Theta,a,d_2,d_1,g); %(a),(l+1),(l),(g)
                Theta = permute(Theta,[1 3 4 2]); %(a),(l),(g),(l+1)
                Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g)
                Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g) % !!! remove this line
                [A,S,V]=svd(Theta);
                
                S((abs(S))/S(1,1)<mps_cpn.zero_thres)=0;
                
                r = nnz(S);
                if r>r_min
                    r_min=r;
                    if r_min>mps_cpn.Dmax
                        r_min=mps_cpn.Dmax;
                    end
                end
                if r==0
                    r=1;
                    S = S(1:r,1:r);
                    Error=0;
                elseif r>mps_cpn.Dmax
                    r=mps_cpn.Dmax;
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
                Norm = Norm+diag(S)'*diag(S);
                if strcmp(Sweep_Direction,'L-R')
                    A = A(:,1:r);
                    V = S*V(:,1:r)';
                else
                    A = A(:,1:r)*S;
                    V = V(:,1:r)';
                end
                
                Total_error=Total_error+Error;
                
                % A == (l a),(r)
                A = reshape(A,a,d_1,r); %(a),(l),(r)
                A = permute(A,[1 3 2]); %(a),(r),(l)
                for count = 1:d_1
                    L_m(count,count2) = {A(:,:,count)};
                end
                
                % V == (r),(l+1 g)
                V = reshape(V,r,g,d_2); %(r),(g),(l+1)
                for count = 1:d_2
                    R_mm(count,count2+(mps_cpn.N+1-g_1)) = {V(:,:,count)};
                end
            end
            R_mm(cellfun(@isempty,R_mm)) = {0};
            L_m(cellfun(@isempty,L_m)) = {0};
            
            
            
            [nrows, ~] = cellfun(@size, L_m);
            a = max(max(nrows));
            [~, ncol] = cellfun(@size, R_mm);
            g = max(max(ncol));
            
            L_m = cellfun(@(x) mps_cpn.cell_resize_func(x,[a,r_min]),L_m,'UniformOutput', false);
            R_mm = cellfun(@(x) mps_cpn.cell_resize_func(x,[r_min,g]),R_mm,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update First Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mps_cpn.data{m}=L_m(:,end:-1:1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update Second Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} back to N_Left labelling
            %%%%%%%%%%%%%%%%%%
            
            temp = cell(a_2,g_2);
            for alpha = 1:a_2
                temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(r_min,g)};
            
            mps_cpn.data{m+1}=temp;
        end
        
        function [mps_cpn,Total_error] = SWAP_Loc_Op_2s_Apply(mps_cpn,U,m,Sweep_Direction)
            % Apply the 2-site local operator to sites m and m+1 and then SWAP the
            % position of sites m and m+1. Apply SVD and truncate if necessary.  

            [a_1,g_1] = size(mps_cpn.data{m});
            [a_2,g_2] = size(mps_cpn.data{m+1});
            
            Nl_2 = min(mps_cpn.N+1,(mps_cpn.d-1)*m+1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the tensor for the two-sites so we can apply the operator.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} to N_Right labelling
            %%%%%%%%%%%%%%%%%%
            temp_T = cell(a_2,g_2);
            temp_IT = cell(a_2,Nl_2);
            expect = max(g_2,Nl_2):-1:1;
            for k = 1:a_2
                track=1;
                for k_2 = (mps_cpn.N+1-k+1):-1:1
                    if k_2>mps_cpn.N+1-g_2-k+1 && k_2<=Nl_2
                        temp_T(k,k_2) = {expect(1,max(g_2,Nl_2)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_cpn.N+1-g_2-k+2):min(Nl_2,(mps_cpn.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_cpn.data{m+1};
            
            R_mm = cell(a_2,Nl_2);
            for alpha = 1:a_2
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            [nrows, ncols] = cellfun(@size, temp);
            beta = max(max(nrows));
            g = max(max(ncols));
            R_mm(cellfun(@isempty,R_mm)) = {zeros(beta,g)};
            
            %%%%%%%%%%%%%%%%%%
            %Collect entries of (L_m X R_mm) to apply the SVD's
            %%%%%%%%%%%%%%%%%%
            [nrows, ~] = cellfun(@size, mps_cpn.data{m});
            a = max(max(nrows));
            L_m = mps_cpn.data{m}(:,end:-1:1);
            
            d_1_s = (1:mps_cpn.N+1);
            d_1_s(d_1_s>mps_cpn.d)=mps_cpn.d;
            
            d_2_s = (mps_cpn.N+1):-1:1;
            d_2_s(d_2_s>mps_cpn.d)=mps_cpn.d;
            
            temp = cell(a_1*a_2,length((1+mps_cpn.N+1-g_1):Nl_2));
            track=1;
            for gamma = (1+mps_cpn.N+1-g_1):Nl_2
                for counter = 1:d_1_s(gamma)
                    temp(((counter-1)*mps_cpn.d+1):(((counter-1)*mps_cpn.d)+d_2_s(gamma)),track) = {zeros(a,g)};
                end
                track=track+1;
            end
            
            track=1;
            for gamma = (mps_cpn.N+1-g_1+1):Nl_2
                for d_1 = 1:d_1_s(gamma)
                    for d_2 = 1:d_2_s(gamma)
                        if d_1+d_2-2<=mps_cpn.N
                            for count = (-min(d_1,track)+1):(d_2-1)
                                if ((d_1-1)*a_1+d_2 + count*(mps_cpn.d-1))<=a_1*a_2 && track+count<=size(temp,2) && ~isempty(temp{(d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),track+count})
                                    temp{(d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),track+count} = temp{(d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),track+count} +U((d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),(d_1-1)*a_1+d_2)*L_m{d_1,track}*R_mm{d_2,gamma};
                                end
                            end
                        end
                    end
                end
                track=track+1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SWAP
            s_d = 0:mps_cpn.d-1;
            s_id = ones(1,mps_cpn.d);
            s_N_1 = kron(s_d,s_id)';
            s_N_2 = kron(s_id,s_d)';
            Sw = s_N_2-s_N_1;
            
            temp_temp = cell(a_1*a_2,length((1+mps_cpn.N+1-g_1):Nl_2));
            for ii = 1:size(temp,1)
                temp_temp(ii,:) = circshift(temp(ii,:),[0,Sw(ii)]);
            end
            
            temp_temp = reshape(temp_temp,a_1,a_2,size(temp_temp,2));
            temp_temp = permute(temp_temp,[2 1 3]);
            temp_temp = reshape(temp_temp,size(temp));
            
            ind = find(1==cellfun(@isempty,temp_temp)~=cellfun(@isempty,temp));
            if ~isempty(ind)
                temp_temp(ind) = {zeros(a,g)};
            end
            
            ind_2 = find(1==mod(cellfun(@isempty,temp_temp)+1,2)~=mod(cellfun(@isempty,temp)+1,2));
            if ~isempty(ind_2)
                temp_temp(ind_2) = {[]};
            end
            temp = temp_temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            temp_2 = cellfun(@(x) mps_cpn.func(x),temp,'UniformOutput', false);
            
            d_1_s=d_1_s((mps_cpn.N+1-g_1+1):Nl_2);
            d_2_s=d_2_s((mps_cpn.N+1-g_1+1):Nl_2);
            d_len = d_1_s.*d_2_s;
            
            Norm = zeros(max(d_len),min(Nl_2,g_1));
            for count = 1:min(Nl_2,g_1)
                Norm(1:d_len(count),count) = cell2mat(temp_2(:,count));
            end
            
            Norm = sum(sum(Norm));
            temp = cellfun(@(x) x/sqrt(Norm),temp,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do SVD's
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R_mm = cell(mps_cpn.d,(mps_cpn.N+1));
            L_m = cell(mps_cpn.d,g_1);
            [a,~] = size(mps_cpn.data{m}{1,1});
            [~,g] = size(mps_cpn.data{m+1}{1,1});
            r_min=1;
            Norm=0;
            Total_error=0;
            for count2 = 1:min(Nl_2,g_1)
                Theta = cell2mat(temp(:,count2));
                d_1 = d_1_s(count2);
                d_2 = d_2_s(count2);
                
                % Theta == (d_1 d_2 alpha),(gamma)
                Theta = reshape(Theta,a,d_2,d_1,g); %(a),(l+1),(l),(g)
                Theta = permute(Theta,[1 3 4 2]); %(a),(l),(g),(l+1)
                Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g)
                Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g)
                [A,S,V]=svd(Theta);
                
                S((abs(S))/S(1,1)<mps_cpn.zero_thres)=0;
                
                r = nnz(S);
                if r>r_min
                    r_min=r;
                    if r_min>mps_cpn.Dmax
                        r_min=mps_cpn.Dmax;
                    end
                end
                if r==0
                    r=1;
                    S = S(1:r,1:r);
                    Error=0;
                elseif r>mps_cpn.Dmax
                    r=mps_cpn.Dmax;
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
                Norm = Norm+diag(S)'*diag(S);
                if strcmp(Sweep_Direction,'L-R')
                    A = A(:,1:r);
                    V = S*V(:,1:r)';
                else
                    A = A(:,1:r)*S;
                    V = V(:,1:r)';
                end
                
                Total_error=Total_error+Error;
                
                % A == (l a),(r)
                A = reshape(A,a,d_1,r); %(a),(l),(r)
                A = permute(A,[1 3 2]); %(a),(r),(l)
                for count = 1:d_1
                    L_m(count,count2) = {A(:,:,count)};
                end
                
                % V == (r),(l+1 g)
                V = reshape(V,r,g,d_2); %(r),(g),(l+1)
                for count = 1:d_2
                    R_mm(count,count2+(mps_cpn.N+1-g_1)) = {V(:,:,count)};
                end
            end
            R_mm(cellfun(@isempty,R_mm)) = {0};
            L_m(cellfun(@isempty,L_m)) = {0};
            
            
            
            [nrows, ~] = cellfun(@size, L_m);
            a = max(max(nrows));
            [~, ncol] = cellfun(@size, R_mm);
            g = max(max(ncol));
            
            L_m = cellfun(@(x) mps_cpn.cell_resize_func(x,[a,r_min]),L_m,'UniformOutput', false);
            R_mm = cellfun(@(x) mps_cpn.cell_resize_func(x,[r_min,g]),R_mm,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update First Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mps_cpn.data{m}=L_m(:,end:-1:1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update Second Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} back to N_Left labelling
            %%%%%%%%%%%%%%%%%%
            
            temp = cell(a_2,g_2);
            for alpha = 1:a_2
                temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(r_min,g)};
            
            mps_cpn.data{m+1}=temp;
        end
        
        function [mps_cpn,Total_error] = SWAP_2s_Only(mps_cpn,m,Sweep_Direction)
            % SWAP the position of sites m and m+1. Apply SVD and truncate 
            % if necessary.
            
            [a_1,g_1] = size(mps_cpn.data{m});
            [a_2,g_2] = size(mps_cpn.data{m+1});
            
            Nl_2 = min(mps_cpn.N+1,(mps_cpn.d-1)*m+1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the tensor for the two-sites so we can apply the operator.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} to N_Right labelling
            %%%%%%%%%%%%%%%%%%
            temp_T = cell(a_2,g_2);
            temp_IT = cell(a_2,Nl_2);
            expect = max(g_2,Nl_2):-1:1;
            for k = 1:a_2
                track=1;
                for k_2 = (mps_cpn.N+1-k+1):-1:1
                    if k_2>mps_cpn.N+1-g_2-k+1 && k_2<=Nl_2
                        temp_T(k,k_2) = {expect(1,max(g_2,Nl_2)-track+1)};
                    end
                    track=track+1;
                end
                for k_2 = max(1,mps_cpn.N+1-g_2-k+2):min(Nl_2,(mps_cpn.N+2-k))
                    temp_IT(k,k_2) = {k_2};
                end
            end
            
            temp = mps_cpn.data{m+1};
            
            R_mm = cell(a_2,Nl_2);
            for alpha = 1:a_2
                R_mm(alpha,cell2mat(temp_IT(alpha,:))) = temp(alpha,cell2mat(temp_T(alpha,:)));
            end
            [nrows, ncols] = cellfun(@size, temp);
            beta = max(max(nrows));
            g = max(max(ncols));
            R_mm(cellfun(@isempty,R_mm)) = {zeros(beta,g)};
            
            %%%%%%%%%%%%%%%%%%
            %Collect entries of (L_m X R_mm) to apply the SVD's
            %%%%%%%%%%%%%%%%%%
            [nrows, ~] = cellfun(@size, mps_cpn.data{m});
            a = max(max(nrows));
            L_m = mps_cpn.data{m}(:,end:-1:1);
            
            d_1_s = (1:mps_cpn.N+1);
            d_1_s(d_1_s>mps_cpn.d)=mps_cpn.d;
            
            d_2_s = (mps_cpn.N+1):-1:1;
            d_2_s(d_2_s>mps_cpn.d)=mps_cpn.d;
            
            temp = cell(a_1*a_2,length((1+mps_cpn.N+1-g_1):Nl_2));
            track=1;
            for gamma = (1+mps_cpn.N+1-g_1):Nl_2
                for counter = 1:d_1_s(gamma)
                    temp(((counter-1)*mps_cpn.d+1):(((counter-1)*mps_cpn.d)+d_2_s(gamma)),track) = {zeros(a,g)};
                end
                track=track+1;
            end
            
            id = eye(mps_cpn.d^2);
            
            track=1;
            for gamma = (mps_cpn.N+1-g_1+1):Nl_2
                for d_1 = 1:d_1_s(gamma)
                    for d_2 = 1:d_2_s(gamma)
                        if d_1+d_2-2<=mps_cpn.N
                            for count = (-min(d_1,track)+1):(d_2-1)
                                if ((d_1-1)*a_1+d_2 + count*(mps_cpn.d-1))<=a_1*a_2 && track+count<=size(temp,2) && ~isempty(temp{(d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),track+count})
                                    temp{(d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),track+count} = temp{(d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),track+count} +id((d_1-1)*a_1+d_2 + count*(mps_cpn.d-1),(d_1-1)*a_1+d_2)*L_m{d_1,track}*R_mm{d_2,gamma};
                                end
                            end
                        end
                    end
                end
                track=track+1;
            end
            
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SWAP
            s_d = 0:mps_cpn.d-1;
            s_id = ones(1,mps_cpn.d);
            s_N_1 = kron(s_d,s_id)';
            s_N_2 = kron(s_id,s_d)';
            Sw = s_N_2-s_N_1;
            
            temp_temp = cell(a_1*a_2,length((1+mps_cpn.N+1-g_1):Nl_2));
            for ii = 1:size(temp,1)
                temp_temp(ii,:) = circshift(temp(ii,:),[0,Sw(ii)]);
            end
            
            temp_temp = reshape(temp_temp,a_1,a_2,size(temp_temp,2));
            temp_temp = permute(temp_temp,[2 1 3]);
            temp_temp = reshape(temp_temp,size(temp));
            
            ind = find(1==cellfun(@isempty,temp_temp)~=cellfun(@isempty,temp));
            if ~isempty(ind)
                temp_temp(ind) = {zeros(a,g)};
            end
            
            ind_2 = find(1==mod(cellfun(@isempty,temp_temp)+1,2)~=mod(cellfun(@isempty,temp)+1,2));
            if ~isempty(ind_2)
                temp_temp(ind_2) = {[]};
            end
            temp = temp_temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            temp_2 = cellfun(@(x) mps_cpn.func(x),temp,'UniformOutput', false);
            
            d_1_s=d_1_s((mps_cpn.N+1-g_1+1):Nl_2);
            d_2_s=d_2_s((mps_cpn.N+1-g_1+1):Nl_2);
            d_len = d_1_s.*d_2_s;
            
            Norm = zeros(max(d_len),min(Nl_2,g_1));
            for count = 1:min(Nl_2,g_1)
                Norm(1:d_len(count),count) = cell2mat(temp_2(:,count));
            end
            
            Norm = sum(sum(Norm));
            temp = cellfun(@(x) x/sqrt(Norm),temp,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do SVD's
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R_mm = cell(mps_cpn.d,(mps_cpn.N+1));
            L_m = cell(mps_cpn.d,g_1);
            [a,~] = size(mps_cpn.data{m}{1,1});
            [~,g] = size(mps_cpn.data{m+1}{1,1});
            r_min=1;
            Norm=0;
            Total_error=0;
            for count2 = 1:min(Nl_2,g_1)
                Theta = cell2mat(temp(:,count2));
                d_1 = d_1_s(count2);
                d_2 = d_2_s(count2);
                
                % Theta == (d_1 d_2 alpha),(gamma)
                Theta = reshape(Theta,a,d_2,d_1,g); %(a),(l+1),(l),(g)
                Theta = permute(Theta,[1 3 4 2]); %(a),(l),(g),(l+1)
                Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g)
                Theta = reshape(Theta,a*d_1,g*d_2); %(l a),(l+1 g)
                [A,S,V]=svd(Theta);
                
                S((abs(S))/S(1,1)<mps_cpn.zero_thres)=0;
                
                r = nnz(S);
                if r>r_min
                    r_min=r;
                    if r_min>mps_cpn.Dmax
                        r_min=mps_cpn.Dmax;
                    end
                end
                if r==0
                    r=1;
                    S = S(1:r,1:r);
                    Error=0;
                elseif r>mps_cpn.Dmax
                    r=mps_cpn.Dmax;
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
                Norm = Norm+diag(S)'*diag(S);
                if strcmp(Sweep_Direction,'L-R')
                    A = A(:,1:r);
                    V = S*V(:,1:r)';
                else
                    A = A(:,1:r)*S;
                    V = V(:,1:r)';
                end
                
                Total_error=Total_error+Error;
                
                % A == (l a),(r)
                A = reshape(A,a,d_1,r); %(a),(l),(r)
                A = permute(A,[1 3 2]); %(a),(r),(l)
                for count = 1:d_1
                    L_m(count,count2) = {A(:,:,count)};
                end
                
                % V == (r),(l+1 g)
                V = reshape(V,r,g,d_2); %(r),(g),(l+1)
                for count = 1:d_2
                    R_mm(count,count2+(mps_cpn.N+1-g_1)) = {V(:,:,count)};
                end
            end
            R_mm(cellfun(@isempty,R_mm)) = {0};
            L_m(cellfun(@isempty,L_m)) = {0};
            
            
            
            [nrows, ~] = cellfun(@size, L_m);
            a = max(max(nrows));
            [~, ncol] = cellfun(@size, R_mm);
            g = max(max(ncol));
            
            L_m = cellfun(@(x) mps_cpn.cell_resize_func(x,[a,r_min]),L_m,'UniformOutput', false);
            R_mm = cellfun(@(x) mps_cpn.cell_resize_func(x,[r_min,g]),R_mm,'UniformOutput', false);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update First Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mps_cpn.data{m}=L_m(:,end:-1:1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update Second Site
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%
            %Convert mps_cpn.data{m+1} back to N_Left labelling
            %%%%%%%%%%%%%%%%%%
            
            temp = cell(a_2,g_2);
            for alpha = 1:a_2
                temp(alpha,cell2mat(temp_T(alpha,:))) = R_mm(alpha,cell2mat(temp_IT(alpha,:)));
            end
            temp(cellfun(@isempty,temp)) = {zeros(r_min,g)};
            
            mps_cpn.data{m+1}=temp;
        end
        
        function Corr = plot_Corr(mps_cpn,Corr)
            % Plot the results of the correlation function caluclation
            %
            % Either plots a 1D graph of the particle number distribution,
            % or a 2D surf plot of the site-site correlation function
            % depending on the dimension of the input matrix.
            %
            % !!! Create a directory with useful extra functions and put
            % this one there. This function does not even use mps_cpn.
            
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

        function Corr = plot_Corr_Y_Junc(mps_cpn,Corr,M_L,M_R1,M_R2)
            % Plot the results of the correlation function caluclation
            %
            % Either plots a 1D graph of the particle number distribution,
            % or a 2D surf plot of the site-site correlation function
            % depending on the dimension of the input matrix.
            %
            % !!! Create a directory with useful extra functions and put
            % this one there. This function does not even use mps_cpn.
            
            M = max(size(Corr));
            dim = sum(size(Corr)>1);
            
            if dim ==1
                figure
                ax1 = subplot(3,2,3);
                plot(1:M_L+1,Corr(1:M_L+1),'b')
                xlabel('Lattice Site','fontsize',16);
                ylabel('Population','fontsize',16);
                axis([1 M_L+1 0 1])
                ax2 = subplot(3,2,2);
                plot(M_L+2:(M_L+1+M_R1),Corr(M_L+2:(M_L+1+M_R1)),'r')
                xlabel('Lattice Site','fontsize',16);
                ylabel('Population','fontsize',16);
                axis([M_L+2 (M_L+1+M_R1) 0 1])
                ax3 = subplot(3,2,6);
                plot(M_L+2:(M_L+1+M_R2),Corr(M_L+2+M_R1:(M_L+1+M_R1+M_R2)),'g')
                xlabel('Lattice Site','fontsize',16);
                ylabel('Population','fontsize',16);
                axis([M_L+2 (M_L+1+M_R2) 0 1])
            else
                error('Unsuported Input!');
            end
        end
        
        function [mps_cpn,Total_error] = TEBD_Local_2s_Gates_SawTooth(mps_cpn,dt,J,Jdash,U,E)
            % TEBD time evolution. Advance MPS by one time step with a
            % local 2 site Hamiltonian.
            %
            % For the SawTooth lattice:
            %                       A   A   A
            %                    J'/ \ / \ /
            %  SAWTOOTH LATTICE - B---B---B
            %                       J
            %
            % Applies hopping between B siters using SWAP gates
            
            M = size(mps_cpn.data,2);
            
            H_1 = mps_cpn.Loc_2s_Ham(M,Jdash,U,E);
            H_2 = mps_cpn.Loc_2s_Ham(M/2,J,U,E);

            Total_error=0;
            Sweep_Direction = 'L-R';
            for t = 1:size(mps_cpn.ST_Order,2);
                
                tau = dt*mps_cpn.ST_Order(t);
                U_T_1 = cell(size(H_1));
                U_T_2 = cell(size(H_2));
                for m = 1:(M-1)
                    U_T_1{m} = expm(-1i*tau*H_1{m});
                    if m <= size(H_2,2);
                        U_T_2{m} = expm(-1i*tau*H_2{m});
                    end
                end
                
                [mps_cpn,Error] = mps_cpn.Time_Evolve_SawTooth(U_T_1,U_T_2,Sweep_Direction);
                Total_error = Total_error+Error;
                
                if strcmp(Sweep_Direction,'L-R')
                    Sweep_Direction = 'R-L';
                else
                    Sweep_Direction = 'L-R';
                end
            end
            
        end
        
        function [mps_cpn,Total_error] = Time_Evolve_SawTooth(mps_cpn,U_T_1,U_T_2,Sweep_Direction)
            % Apply two-site local operators to Entire MPS and truncate if
            % necessary
            %
            % For the SawTooth lattice:
            %                       A   A   A
            %                    J'/ \ / \ /
            %  SAWTOOTH LATTICE - B---B---B
            %                       J
            %
            % Applies hopping between B siters using SWAP gates
            
            M = size(mps_cpn.data,2);
            if strcmp(Sweep_Direction,'L-R')
                Total_error=0;
                for m=1:2:(M-1)
                    if m ~= M-1
                        [mps_cpn,Error] = mps_cpn.SWAP_Loc_Op_2s_Apply(U_T_1{m},m,Sweep_Direction);
                        Total_error = Total_error+Error;
                        [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_T_2{(m+1)/2},m+1,Sweep_Direction);
                        Total_error = Total_error+Error;
                        [mps_cpn,Error] = mps_cpn.SWAP_2s_Only(m,Sweep_Direction);
                        Total_error = Total_error+Error;
                        [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_T_1{m+1},m+1,Sweep_Direction);
                        Total_error = Total_error+Error;
                    else
                        [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_T_1{m},m,Sweep_Direction);
                        Total_error = Total_error+Error;
                    end
                end
            elseif strcmp(Sweep_Direction,'R-L')
                Total_error=0;
                for m=M:-2:2
                    if m==M
                        [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_T_1{m-1},m-1,Sweep_Direction);
                        Total_error = Total_error+Error;
                    else
                        [mps_cpn,Error] = mps_cpn.SWAP_Loc_Op_2s_Apply(U_T_1{m},m,Sweep_Direction);
                        Total_error = Total_error+Error;
                        [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_T_2{m/2},m-1,Sweep_Direction);
                        Total_error = Total_error+Error;
                        [mps_cpn,Error] = mps_cpn.SWAP_2s_Only(m,Sweep_Direction);
                        Total_error = Total_error+Error;
                        [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_T_1{m-1},m-1,Sweep_Direction);
                        Total_error = Total_error+Error;
                    end
                end
            end
        end
        
        function [mps_cpn,Total_error] = TEBD_Local_2s_Gates_Y_Junction(mps_cpn,dt,J,U,E,M_L,M_R1,M_R2)
            % TEBD time evolution. Advance MPS by one time step with a
            % local 2 site Hamiltonian.
            %
            %  Two Right Branches are connected by a long range interaction term
            %  implemented with SWAP gates
            %
            %  MPS Layout:
            %  LEFT_BRANCH - CENTRAL_SITE - RIGHT_BRANCH_1 - RIGHT_BRANCH_2
            %
            % Many SWAP gates are implemented in the same operation to connect the
            % Central site with the second right branch.
            
            M = size(mps_cpn.data,2);
            
            H_LR1 = mps_cpn.Loc_2s_Ham(M_L+M_R1+1,J,U,E(1:M_L+M_R1+1));
            H_R2 = mps_cpn.Loc_2s_Ham(M_R2+1,J,U,E(M_L+M_R1+1:end));

            num = mps_cpn.b_dag*mps_cpn.b;
            id = eye(mps_cpn.d);
            H_1 = U/4*num*(num - id) + E(M_L)/2*num;
            H_2 = U/6*num*(num - id) + E(M_L+1)/3*num;
            H_Bond{1} = kron(H_1,id)+kron(id,H_2)-J*kron(mps_cpn.b_dag,mps_cpn.b)-J*kron(mps_cpn.b,mps_cpn.b_dag);

            H_1 = U/6*num*(num - id) + E(M_L+1)/3*num;
            H_2 = U/4*num*(num - id) + E(M_L+2)/2*num;
            H_Bond{2} = kron(H_1,id)+kron(id,H_2)-J*kron(mps_cpn.b_dag,mps_cpn.b)-J*kron(mps_cpn.b,mps_cpn.b_dag);

            
            Total_error=0;
            Sweep_Direction = 'L-R';
            for t = 1:size(mps_cpn.ST_Order,2)
                
                tau = dt*mps_cpn.ST_Order(t);
                U_LR1 = cell(size(H_LR1));
                U_R2 = cell(1,M-1);
                track=2;
                for m = 1:(M-1)
                    if m > M_L+M_R1+1
                        U_R2{m} = expm(-1i*tau*H_R2{track});
                        track=track+1;
                    elseif m < M_L+M_R1+1
                        U_LR1{m} = expm(-1i*tau*H_LR1{m});
                    end
                end
                U_Bond=cell(1,2);
                for m = 1:2
                    U_Bond{m}=expm(-1i*tau*H_Bond{m});
                end
                
                [mps_cpn,Error] = mps_cpn.Time_Evolve_Y_Junction(U_LR1,U_Bond,U_R2,Sweep_Direction,M_L,M_R1,M_R2);
                Total_error = Total_error+Error;
                
                if strcmp(Sweep_Direction,'L-R')
                    Sweep_Direction = 'R-L';
                else
                    Sweep_Direction = 'L-R';
                end
            end
            
        end
        
        function [mps_cpn,Total_error] = Time_Evolve_Y_Junction(mps_cpn,U_LR1,U_Bond,U_R2,Sweep_Direction,M_L,M_R1,M_R2)
            % Apply two-site local operators to Entire MPS and truncate if
            % necessary
            %
            %  Two Right Branches are connected by a long range interaction term
            %  implemented with SWAP gates
            %
            %  MPS Layout:
            %  LEFT_BRANCH - CENTRAL_SITE - RIGHT_BRANCH_1 - RIGHT_BRANCH_2
            %
            % Many SWAP gates are implemented in the same operation to connect the
            % Central site with the second right branch.
            
            
            if strcmp(Sweep_Direction,'L-R')
                Total_error=0;
                for m=1:(M_L-1)
                    [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_LR1{m},m,Sweep_Direction);
                    Total_error = Total_error+Error;
                end
                [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_Bond{1},M_L,Sweep_Direction);
                Total_error = Total_error+Error;
                [mps_cpn,Error] = mps_cpn.SWAP_Loc_Op_2s_Apply(U_Bond{2},M_L+1,Sweep_Direction);
                Total_error = Total_error+Error;
                m=M_L+2;
                for sw = 1:M_R1-1
                    [mps_cpn,Error] = mps_cpn.SWAP_2s_Only(m,Sweep_Direction);
                    Total_error = Total_error+Error;
                    m=m+1;
                end
                [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_Bond{2},m,Sweep_Direction);
                Total_error = Total_error+Error;
                m=m-1;
                for sw = 1:(M_R1-1)
                    [mps_cpn,Error] = mps_cpn.SWAP_2s_Only(m,'R-L');
                    Total_error = Total_error+Error;
                    m=m-1;
                end
                [mps_cpn,Error] = mps_cpn.SWAP_2s_Only(m,Sweep_Direction);
                Total_error = Total_error+Error;
                
                for m = (M_L+2):(M_L+M_R1)
                    [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_LR1{m},m,Sweep_Direction);
                    Total_error = Total_error+Error;
                end
                [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(eye(mps_cpn.d^2),M_L+M_R1+1,Sweep_Direction);
                Total_error = Total_error+Error;
                
                for m = (M_L+M_R1+2):(M_L+M_R1+M_R2)
                    [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_R2{m},m,Sweep_Direction);
                    Total_error = Total_error+Error;
                end
            elseif strcmp(Sweep_Direction,'R-L')
                Total_error=0;
                for m = (M_L+M_R1+M_R2):-1:(M_L+M_R1+2)
                    [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_R2{m},m,Sweep_Direction);
                    Total_error = Total_error+Error;
                end
                [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(eye(mps_cpn.d^2),M_L+M_R1+1,Sweep_Direction);
                Total_error = Total_error+Error;
                
                for m = (M_L+M_R1):-1:(M_L+2)
                    [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_LR1{m},m,Sweep_Direction);
                    Total_error = Total_error+Error;
                end

                [mps_cpn,Error] = mps_cpn.SWAP_Loc_Op_2s_Apply(U_Bond{2},M_L+1,'L-R');
                Total_error = Total_error+Error;
                m=M_L+2;
                for sw = 1:M_R1-1
                    [mps_cpn,Error] = mps_cpn.SWAP_2s_Only(m,'L-R');
                    Total_error = Total_error+Error;
                    m=m+1;
                end
                [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_Bond{2},m,Sweep_Direction);
                Total_error = Total_error+Error;
                m=m-1;
                for sw = 1:M_R1
                    [mps_cpn,Error] = mps_cpn.SWAP_2s_Only(m,Sweep_Direction);
                    Total_error = Total_error+Error;
                    m=m-1;
                end
                [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_Bond{1},M_L,Sweep_Direction);
                Total_error = Total_error+Error;
                for m=(M_L-1):-1:1
                    [mps_cpn,Error] = mps_cpn.Loc_Op_2s_Apply(U_LR1{m},m,Sweep_Direction);
                    Total_error = Total_error+Error;
                end
            end
        end
        
        function Energy = Get_Energy_2s(mps_cpn,J,U,E)
            
            
            M = size(mps_cpn.data,2);
            
            H_2s = mps_cpn.Loc_2s_Ham(M,J,U,E);
            Energy=0;
            
            for m=1:(M-1)
                mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
                
                [temp_mps_cpn,~] = mps_cpn.Loc_Op_2s_Apply(H_2s{m},m,'L-R');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculate Overlap
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                expect=cell(1,temp_mps_cpn.N+1);
                expect(cellfun(@isempty,expect)) = {1};
                for n = 1:M
                    [a,g] = size(temp_mps_cpn.data{n});
                    [alpha,gamma] = size(temp_mps_cpn.data{n}{end,1});
                    
                    temp_E = cell(a,g);
                    for k = 1:a
                        track=1;
                        for k_2 = (temp_mps_cpn.N+1-k+1):-1:1
                            if k_2 <= g
                                temp_E{k,k_2} = expect{1,temp_mps_cpn.N+1-track+1};
                            end
                            track=track+1;
                        end
                    end
                    temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                    
                    temp = cellfun(@(x,y,z) x'*y*z, mps_cpn.data{n},temp_E,temp_mps_cpn.data{n}, 'UniformOutput', false);
                    
                    expect=cell(1,g);
                    expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                    
                    for i=1:a
                        expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                    end
                    if n<M
                        if size(expect,2)<temp_mps_cpn.d
                            expect(:,size(expect,2)+1:temp_mps_cpn.d)={zeros(gamma,gamma)};
                        end
                    end
                end
                Energy=Energy+cell2mat(expect);
                
            end
        end
        
        function Var = Get_Var_2s(mps_cpn,J,U,E)
            
            Energy = mps_cpn.Get_Energy_2s(J,U,E);
            
            M = size(mps_cpn.data,2);
            
            H_2s = mps_cpn.Loc_2s_Ham(M,J,U,E);
            Var=0;
            
            for m=1:(M-1)
                mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
                
                [temp_mps_cpn_2,~] = mps_cpn.Loc_Op_2s_Apply(H_2s{m},m,'L-R');
                
                for mm = 1:(M-1)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{mm}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_2s{mm},mm,'L-R');
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Calculate Overlap
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    expect=cell(1,temp_mps_cpn.N+1);
                    expect(cellfun(@isempty,expect)) = {1};
                    for n = 1:M
                        [a,g] = size(temp_mps_cpn.data{n});
                        [alpha,gamma] = size(temp_mps_cpn.data{n}{end,1});
                        
                        temp_E = cell(a,g);
                        for k = 1:a
                            track=1;
                            for k_2 = (temp_mps_cpn.N+1-k+1):-1:1
                                if k_2 <= g
                                    temp_E{k,k_2} = expect{1,temp_mps_cpn.N+1-track+1};
                                end
                                track=track+1;
                            end
                        end
                        temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};
                        
                        temp = cellfun(@(x,y,z) x'*y*z, mps_cpn.data{n},temp_E,temp_mps_cpn.data{n}, 'UniformOutput', false);
                        
                        expect=cell(1,g);
                        expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                        
                        for i=1:a
                            expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                        end
                        if n<M
                            if size(expect,2)<temp_mps_cpn.d
                                expect(:,size(expect,2)+1:temp_mps_cpn.d)={zeros(gamma,gamma)};
                            end
                        end
                    end
                    Var=Var+cell2mat(expect);
                end
            end
            
            Var=(Var-Energy^2)/Var;
            
        end
        
        function Energy = Get_Energy_2s_Y_Junction(mps_cpn,J,U,E,M_L,M_R1,M_R2)

            M = size(mps_cpn.data,2);
            
            temp_H_LR1 = mps_cpn.Loc_2s_Ham(M_L+M_R1+1,J,U,E(1:M_L+M_R1+1));
            temp_H_R2 = mps_cpn.Loc_2s_Ham(M_R2+1,J,U,E(M_L+M_R1+1:end));
            
            H_LR1 = cell(size(temp_H_LR1));
            H_R2 = cell(1,M-1);
            track=2;
            for m = 1:(M-1)
                if m > M_L+M_R1+1
                    H_R2{m} = temp_H_R2{track};
                    track=track+1;
                elseif m < M_L+M_R1+1
                    H_LR1{m} = temp_H_LR1{m};
                end
            end
            num = mps_cpn.b_dag*mps_cpn.b;
            id = eye(mps_cpn.d);
            H_1 = U/4*num*(num - id) + E(M_L)/2*num;
            H_2 = U/6*num*(num - id) + E(M_L+1)/3*num;
            H_Bond{1} = kron(H_1,id)+kron(id,H_2)-J*kron(mps_cpn.b_dag,mps_cpn.b)-J*kron(mps_cpn.b,mps_cpn.b_dag);

            H_1 = U/6*num*(num - id) + E(M_L+1)/3*num;
            H_2 = U/4*num*(num - id) + E(M_L+2)/2*num;
            H_Bond{2} = kron(H_1,id)+kron(id,H_2)-J*kron(mps_cpn.b_dag,mps_cpn.b)-J*kron(mps_cpn.b,mps_cpn.b_dag);

            
            Energy=0;
            
            for m=1:(M_L-1)
                mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
                [temp_mps_cpn,~] = mps_cpn.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                Energy=Energy+calc_overlap(M,temp_mps_cpn,mps_cpn);
            end
            m=M_L;
            mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
            [temp_mps_cpn,~] = mps_cpn.Loc_Op_2s_Apply(H_Bond{1},m,'L-R');
            Energy=Energy+calc_overlap(M,temp_mps_cpn,mps_cpn);
            m=M_L+1;
            mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
            [temp_mps_cpn,~] = mps_cpn.Loc_Op_2s_Apply(H_Bond{2},m,'L-R');
            Energy=Energy+calc_overlap(M,temp_mps_cpn,mps_cpn);
            
            m=M_L+1;
            for sw = 1:M_R1
                [mps_cpn,~] = mps_cpn.SWAP_2s_Only(m,'L_R');
                m=m+1;
            end
            mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
            [temp_mps_cpn,~] = mps_cpn.Loc_Op_2s_Apply(H_Bond{2},m,'L_R');
            
            for sw = 1:M_R1
                [mps_cpn,~] = mps_cpn.SWAP_2s_Only(m,'R-L');
                [temp_mps_cpn,~] = temp_mps_cpn.SWAP_2s_Only(m,'R-L');
                m=m-1;
            end
            Energy=Energy+calc_overlap(M,temp_mps_cpn,mps_cpn);
            
            for m = (M_L+2):(M_L+M_R1)
                mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
                [temp_mps_cpn,~] = mps_cpn.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                Energy=Energy+calc_overlap(M,temp_mps_cpn,mps_cpn);
            end
            
            for m = (M_L+M_R1+2):(M_L+M_R1+M_R2)
                mps_cpn.Dmax=size(mps_cpn.data{m}{1,1},2);
                [temp_mps_cpn,~] = mps_cpn.Loc_Op_2s_Apply(H_R2{m},m,'L-R');
                Energy=Energy+calc_overlap(M,temp_mps_cpn,mps_cpn);
            end
            
            
            function Y=calc_overlap(M,temp_mps_cpn,mps_cpn)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculate Overlap
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                expect=cell(1,temp_mps_cpn.N+1);
                expect(cellfun(@isempty,expect)) = {1};
                for n = 1:M
                    if any(size(temp_mps_cpn.data{n}{1,1})<size(mps_cpn.data{n}{1,1}))
                        [r1,r2]=size(mps_cpn.data{n}{1,1});
                        temp_mps_cpn.data{n}=cellfun(@(x) temp_mps_cpn.cell_resize_func(x,[r1,r2]),temp_mps_cpn.data{n},'UniformOutput', false);
                    end
                    
                    [a,g] = size(temp_mps_cpn.data{n});
                    [alpha,gamma] = size(temp_mps_cpn.data{n}{end,1});
                    
                    temp_E = cell(a,g);
                    for k = 1:a
                        track=1;
                        for k_2 = (temp_mps_cpn.N+1-k+1):-1:1
                            if k_2 <= g
                                temp_E{k,k_2} = expect{1,temp_mps_cpn.N+1-track+1};
                            end
                            track=track+1;
                        end
                    end
                    temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};

                    temp = cellfun(@(x,y,z) x'*y*z, mps_cpn.data{n},temp_E,temp_mps_cpn.data{n}, 'UniformOutput', false);
                    
                    expect=cell(1,g);
                    expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                    
                    for i=1:a
                        expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                    end
                    if n<M
                        if size(expect,2)<temp_mps_cpn.d
                            expect(:,size(expect,2)+1:temp_mps_cpn.d)={zeros(gamma,gamma)};
                        end
                    end
                end
                Y=cell2mat(expect);
            end
            
        end
        
        function [Var,Energy] = Get_Var_2s_Y_Junction(mps_cpn,J,U,E,M_L,M_R1,M_R2)

            Energy = mps_cpn.Get_Energy_2s_Y_Junction(J,U,E,M_L,M_R1,M_R2);
            
            M = size(mps_cpn.data,2);
            
            temp_H_LR1 = mps_cpn.Loc_2s_Ham(M_L+M_R1+1,J,U,E(1:M_L+M_R1+1));
            temp_H_R2 = mps_cpn.Loc_2s_Ham(M_R2+1,J,U,E(M_L+M_R1+1:end));
            
            H_LR1 = cell(size(temp_H_LR1));
            H_R2 = cell(1,M-1);
            track=2;
            for m = 1:(M-1)
                if m > M_L+M_R1+1
                    H_R2{m} = temp_H_R2{track};
                    track=track+1;
                elseif m < M_L+M_R1+1
                    H_LR1{m} = temp_H_LR1{m};
                end
            end
            
            Var=0;
            
            for mm=1:(M_L+1)
                mps_cpn.Dmax=size(mps_cpn.data{mm}{1,1},2);
                [temp_mps_cpn_2,~] = mps_cpn.Loc_Op_2s_Apply(H_LR1{mm},mm,'L-R');
                for m=1:(M_L+1)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
                
                m=M_L+1;
                for sw = 1:M_R1
                    [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'L_R');
                    m=m+1;
                end
                temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{M_L+1},m,'L_R');
                
                for sw = 1:M_R1
                    [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'R-L');
                    [temp_mps_cpn,~] = temp_mps_cpn.SWAP_2s_Only(m,'R-L');
                    m=m-1;
                end
                Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                
                for m = (M_L+2):(M_L+M_R1)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
                
                for m = (M_L+M_R1+2):(M_L+M_R1+M_R2)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_R2{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
            end

            mm=M_L+1;
            for sw = 1:M_R1
                [mps_cpn,~] = mps_cpn.SWAP_2s_Only(mm,'L_R');
                mm=mm+1;
            end
            mps_cpn.Dmax=size(mps_cpn.data{mm}{1,1},2);
            [temp_mps_cpn_2,~] = mps_cpn.Loc_Op_2s_Apply(H_LR1{M_L+1},mm,'L_R');
            
            for sw = 1:M_R1
                [mps_cpn,~] = mps_cpn.SWAP_2s_Only(mm,'R-L');
                [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(mm,'R-L');
                mm=mm-1;
            end
            for m=1:(M_L+1)
                temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
            end
            
            m=M_L+1;
            for sw = 1:M_R1
                [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'L_R');
                m=m+1;
            end
            temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
            [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{M_L+1},m,'L_R');
            
            for sw = 1:M_R1
                [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'R-L');
                [temp_mps_cpn,~] = temp_mps_cpn.SWAP_2s_Only(m,'R-L');
                m=m-1;
            end
            Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
            
            for m = (M_L+2):(M_L+M_R1)
                temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
            end
            
            for m = (M_L+M_R1+2):(M_L+M_R1+M_R2)
                temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_R2{m},m,'L-R');
                Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
            end

            for mm = (M_L+2):(M_L+M_R1)
                mps_cpn.Dmax=size(mps_cpn.data{mm}{1,1},2);
                [temp_mps_cpn_2,~] = mps_cpn.Loc_Op_2s_Apply(H_LR1{mm},mm,'L-R');
                for m=1:(M_L+1)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
                
                m=M_L+1;
                for sw = 1:M_R1
                    [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'L_R');
                    m=m+1;
                end
                temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{M_L+1},m,'L_R');
                
                for sw = 1:M_R1
                    [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'R-L');
                    [temp_mps_cpn,~] = temp_mps_cpn.SWAP_2s_Only(m,'R-L');
                    m=m-1;
                end
                Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                
                for m = (M_L+2):(M_L+M_R1)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
                
                for m = (M_L+M_R1+2):(M_L+M_R1+M_R2)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_R2{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
            end
            
            for mm = (M_L+M_R1+2):(M_L+M_R1+M_R2)
                mps_cpn.Dmax=size(mps_cpn.data{mm}{1,1},2);
                [temp_mps_cpn_2,~] = mps_cpn.Loc_Op_2s_Apply(H_R2{mm},mm,'L-R');
                for m=1:(M_L+1)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
                
                m=M_L+1;
                for sw = 1:M_R1
                    [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'L_R');
                    m=m+1;
                end
                temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{M_L+1},m,'L_R');
                
                for sw = 1:M_R1
                    [temp_mps_cpn_2,~] = temp_mps_cpn_2.SWAP_2s_Only(m,'R-L');
                    [temp_mps_cpn,~] = temp_mps_cpn.SWAP_2s_Only(m,'R-L');
                    m=m-1;
                end
                Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                
                for m = (M_L+2):(M_L+M_R1)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_LR1{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
                
                for m = (M_L+M_R1+2):(M_L+M_R1+M_R2)
                    temp_mps_cpn_2.Dmax=size(temp_mps_cpn_2.data{m}{1,1},2);
                    [temp_mps_cpn,~] = temp_mps_cpn_2.Loc_Op_2s_Apply(H_R2{m},m,'L-R');
                    Var=Var+calc_overlap(M,temp_mps_cpn,mps_cpn);
                end
                
                Var=(Var-Energy^2)/Var;
            end
            
            
            function Y=calc_overlap(M,temp_mps_cpn,mps_cpn)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculate Overlap
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                expect=cell(1,temp_mps_cpn.N+1);
                expect(cellfun(@isempty,expect)) = {1};
                for n = 1:M
                    if any(size(temp_mps_cpn.data{n}{1,1})<size(mps_cpn.data{n}{1,1}))
                        [r1,r2]=size(mps_cpn.data{n}{1,1});
                        temp_mps_cpn.data{n}=cellfun(@(x) temp_mps_cpn.cell_resize_func(x,[r1,r2]),temp_mps_cpn.data{n},'UniformOutput', false);
                    end
                    
                    [a,g] = size(temp_mps_cpn.data{n});
                    [alpha,gamma] = size(temp_mps_cpn.data{n}{end,1});
                    
                    temp_E = cell(a,g);
                    for k = 1:a
                        track=1;
                        for k_2 = (temp_mps_cpn.N+1-k+1):-1:1
                            if k_2 <= g
                                temp_E{k,k_2} = expect{1,temp_mps_cpn.N+1-track+1};
                            end
                            track=track+1;
                        end
                    end
                    temp_E(cellfun(@isempty,temp_E)) = {zeros(alpha,alpha)};

                    temp = cellfun(@(x,y,z) x'*y*z, mps_cpn.data{n},temp_E,temp_mps_cpn.data{n}, 'UniformOutput', false);
                    
                    expect=cell(1,g);
                    expect(cellfun(@isempty,expect)) = {zeros(size(temp{1,1}))};
                    
                    for i=1:a
                        expect(1,:) = cellfun(@plus, expect, temp(i,:), 'UniformOutput', false);
                    end
                    if n<M
                        if size(expect,2)<temp_mps_cpn.d
                            expect(:,size(expect,2)+1:temp_mps_cpn.d)={zeros(gamma,gamma)};
                        end
                    end
                end
                Y=cell2mat(expect);
            end
            
        end
        
        function W_TEBD = Calc_State_Vector(mps_cpn,B)
            % Produce exact diagonalisation state vector from MPS
            
            W_TEBD=zeros(size(B,1),1);
            for n = 1:size(B,1)
                
                for m = 1:size(B,2)
                    mps_cpn.d = B(n,m)+1;
                    if m ~= size(B,2)
                        N_R = sum(B(n,(m+1):end))+1;
                    else
                        N_R=1;
                    end
                    
                    if m ==1
                        W=mps_cpn.data{m}{mps_cpn.d,N_R};
                    else
                        W=W*mps_cpn.data{m}{mps_cpn.d,N_R};
                    end
                end
                W_TEBD(n)=W;
            end
        end
        
        function mps_cpn = Set_Runga_Kutta_Order(mps_cpn,order)
            % Set the roots of exp(-i*dt*H). Used in the Runga-Kutta MPO
            % time evolution
            
            % coefficients of exp(x) in the reversed order
            p=1./factorial(order-1:-1:0);
            % roots
            mps_cpn.RK_Roots=roots(p);
            
        end
        
        function [mps_cpn,Total_error] = Runga_Kutta_Time_Evolve(mps_cpn,H,dt)
            % Runga-Kutta time evolution. Advance the MPS by one time step
            % with the operator exp(-1*dt*H).
            %
            % [1] - Juan Jose Garcia-Ripoll, 'Time evolution of Matrix 
            % Product States' IOP, (2006)
            
            order = size(mps_cpn.RK_Roots,1);
            
            Total_error=0;
            for c = 1:order
                H=H.Eff_Time_Evolve(dt,mps_cpn.RK_Roots(c));
                
                mps_cpn = mps_cpn.Apply_MPO(H);
                
                % Truncate and Canonicalise
                [mps_cpn,error] = Canonicalisation_2s(mps_cpn,'L-R');
                Total_error=Total_error+error;
            end
        end
        
        function mps_cpn = Apply_MPO(mps_cpn,H)
            % Apply the MPO H to the MPS
            
            M = size(mps_cpn.data,2);
            
            for m = 1:M
                [k1,k2]=size(H.data{m});
                Store = cell(k1,k2);
                [alpha,beta]=size(mps_cpn.data{m}{1,1});
                [~,NR]=size(mps_cpn.data{m});
                Nl = min((mps_cpn.d-1)*(m-1)+1,mps_cpn.N+1);
                
                [a,g] = size(mps_cpn.data{m});
                temp_Zeros = cell(a,g);
                temp_Zeros(cellfun(@isempty,temp_Zeros)) = {zeros(alpha,beta)};
                for a_count = 1:k1
                    for b_count = 1:k2
                        if all(all(H.data{m}{a_count,b_count}==0))
                            Store{a_count,b_count}=cell2mat(cellfun(@(x,y) x.*y, mps_cpn.data{m},temp_Zeros, 'UniformOutput', false));
                        elseif all(all(H.data{m}{a_count,b_count}==1))
                            Store{a_count,b_count}=cell2mat(mps_cpn.data{m});
                        else
                            O=H.data{m}{a_count,b_count};
                            
                            temp = cell(a,g);
                            temp(cellfun(@isempty,temp)) = {zeros(alpha,beta)};
                            for l = 1:mps_cpn.d
                                for N_R = max(1,(g-Nl-l+2)):max(1,min(g,mps_cpn.N+1-l+1))
                                    for count2 = (-mps_cpn.d+l):(l-1)
                                        temp{l-count2,N_R} = temp{l-count2,N_R}+O(l-count2,l)*mps_cpn.data{m}{l,N_R};
                                    end
                                end
                            end
                            
                            if H.N_track{m}{b_count}==1
                                m_temp = cell(a,g);
                                for m1 = 2:g
                                    m_temp(:,m1) = temp(:,m1-1);
                                end
                                m_temp(cellfun(@isempty,m_temp)) = {zeros(alpha,beta)};
                                temp = m_temp;
                            elseif H.N_track{m}{b_count}==-1
                                m_temp = cell(a,g);
                                for m1 = 1:g-1
                                    m_temp(:,m1) = temp(:,m1+1);
                                end
                                m_temp(cellfun(@isempty,m_temp)) = {zeros(alpha,beta)};
                                temp=m_temp;
                            end
                            
                            Store{a_count,b_count}=cell2mat(temp);
                        end
                        
                    end
                end
                
                Store = cell2mat(Store);
                % Store == (k1,d,alpha),(k2,NR,beta)
                Mat = reshape(Store,alpha,mps_cpn.d,k1,beta,NR,k2); %alpha,d,k1,beta,NR,k2
                Mat = permute(Mat,[1 3 4 6 2 5]); %alpha,k1,beta,k2,d,NR
                Mat = reshape(Mat,alpha*k1,beta*k2,mps_cpn.d,NR); %(k1 alpha),(k2 beta),d,NR
                
                temp = cell(mps_cpn.d,NR);
                for count = 1:mps_cpn.d
                    for Num_r = max(1,(g-Nl-count+2)):max(1,min(g,mps_cpn.N+1-count+1))
                        temp(count,Num_r) = {Mat(:,:,count,Num_r)};
                    end
                end
                temp(cellfun(@isempty,temp)) = {zeros(alpha*k1,beta*k2)};
                
                mps_cpn.data{m} = temp;
                
            end
            
        end
        
        function mps_cpn = Increase_bond_dim(mps_cpn,Bond_Dim)
            % Increase Maximum bond dimension. Adds Zeros to the MPS entries.
            
             M = size(mps_cpn.data,2);
             
             for m = 1:M
                 [~,g] = size(mps_cpn.data{m});
                 Nl = min((mps_cpn.d-1)*(m-1)+1,mps_cpn.N+1);
                 [alpha,beta] = size(mps_cpn.data{m}{1,1});
                 
                 BL = max(alpha,Bond_Dim);
                 BR = max(beta,Bond_Dim);
                 
                 for count = 1:mps_cpn.d
                     for Num_r = max(1,(g-Nl-count+2)):max(1,min(g,mps_cpn.N+1-count+1))
                         temp = mps_cpn.data{m}{count,Num_r};
                         temp(:,beta+1:BR) = zeros(alpha,length(beta+1:BR));
                         temp(alpha+1:BL,:) = zeros(length(alpha+1:BL),BR);
                         mps_cpn.data{m}{count,Num_r}=temp;
                     end
                 end
                 
             end
            
        end
        
        function mps_cpn = DMRG(mps_cpn,H,Sweeps)
            % Apply DMRG procedure
            
            M = size(mps_cpn.data,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Calculate Initial Left Effective Hamiltonian
            %%%%%%%%%%%%%%%%%%%%%%%
            L = cell(1,M);
            L{1} = 1;
            for m = 1:(M-1)
                L{m+1} = mps_cpn.Add_Left_Eff_Ham(H,L{m},m+1);
            end
            
        end
        
        function mps_cpn = Add_Left_Eff_Ham(mps_cpn,H,L,m)
            
%             [k1,k2] = size(H.data,
            
        end

    end
    
end

