classdef mpo_cpn
    % Creates a variety of different MPO's. The main purpose of this class
    % is to create an MPO that can be used in conjunction with the class 
    % "mpo_cpn" to produce algorithms (time-evolution and ground state 
    % searches) that conserve particle number in the Bose-Hubbard model.
    %
    % The class has two main properties: 
    %           1) The conventional MPO in the usual form
    %           2) A similar structure that defines whether each operator
    %           in the MPO increases, decreases or leaves unchanged the
    %           particle number
    
    properties (SetAccess=public)
    end

    properties (SetAccess=protected)
        
        % 1) The conventional MPO in the usual form
        data
        
        % 2) A similar structure that defines whether each operator
        % in the MPO increases, decreases or leaves unchanged the
        % particle number.
        %
        % If the MPO has bond dimension (alpha,beta), the bond dimension
        % of the N_track structure is (1,beta). Each beta defines whether
        % that entry multiplies, at some point to the right, an operator 
        % that increases or decreases the particle number - or alternatively
        % if it never does. 
        % This is important for the "mps_cpn" storage - as it has
        % place holders for different number of particles to the right and
        % left.
        N_track
        
        % Local Bond dimension
        d
        
    end
    
    methods
        
        function mpo_cpn=mpo_cpn(M,N_max)
            % Class constructor
            
            mpo_cpn.data=cell(1,M);
            mpo_cpn.N_track=cell(1,M);
            mpo_cpn.d=N_max+1;
            
        end
        
        function mpo_cpn=Simple_1D_Nearest(mpo_cpn,J,U,E)
            % Create MPO for a simple 1D lattice with only nearest
            % neighbour hopping contributions.
            
            M = size(mpo_cpn.data,2);
            N_max = mpo_cpn.d-1;
            
            b = sqrt(diag(1:N_max,1)); b_dag = sqrt(diag(1:N_max,-1));
            Num_0 = b_dag*b; id = eye(mpo_cpn.d); Z = 0;%zeros(mpo_cpn.d);
            
            mpo_cpn.data{1} = [{(E(1)*Num_0 + U/2*Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{1}];
            
            % N_track:
            % In this case "mpo_cpn.N_track" is the same for every site -
            % except the final site.
            % The first column of the mth MPO multiplies the identity of the
            % (m+1)th MPO. And this entry continues to multiply operators
            % that do not change the particle number. So nothing special is
            % required for this column.
            % The second column adds a particle, so must multiply an
            % annhilation operator at some point in the future.
            % Similar for the third column.
            % The final column is an identity. But it will multiply the
            % final row of the next MPO. But again, we do not need to take
            % into account an increase or decrease in particle number.
            %       For instance, take the multiplication of this entry 
            %       (identity) to the (4,2) entry of the (m+1)th MPO 
            %       (Creation Operator). This will increase the number to 
            %       the right. But at some point in the future this entry
            %       will multiply an annhilation operator, thus decreasing
            %       the number to the right again, leaving the MPS of this
            %       entry the same.
            mpo_cpn.N_track{1} = [{0},{-1},{1},{0}];
            
            mpo_cpn.data{M} = [{1};{b};{b_dag};{(E(M)*Num_0 + U/2 * Num_0*(Num_0-id))}];
            mpo_cpn.N_track{M} = [{0}];
            
            for m = 2:(M-1)
                mpo_cpn.data{m} = [{1},{Z},{Z},{Z};{b},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{1}];
                mpo_cpn.N_track{m} = [{0},{-1},{1},{0}];
            end
            
        end
        
        function mpo_cpn=Simple_1D_Next_Nearest(mpo_cpn,J,J_2,U,E)
            % Create MPO for a simple 1D lattice with nearest neighbour
            % and next nearest neighbour hopping contributions.
            
            M = size(mpo_cpn.data,2);
            N_max = mpo_cpn.d-1;
            
            b = sqrt(diag(1:N_max,1)); b_dag = sqrt(diag(1:N_max,-1));
            Num_0 = b_dag*b; id = eye(mpo_cpn.d); Z = zeros(mpo_cpn.d);
            
            mpo_cpn.data{1} = [{(E(1)*Num_0 + U/2*Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
            mpo_cpn.N_track{1} = [{0},{-1},{1},{0}];
            
            mpo_cpn.data{2} = [{id},{Z},{Z},{Z},{Z},{Z};{b},{Z},{Z},{Z},{J_2/J*id},{Z};{b_dag},{Z},{Z},{Z},{Z},{J_2/J*id};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{Z},{Z},{id}];
            mpo_cpn.N_track{2} = [{0},{-1},{1},{-1},{1},{0}];
            
            mpo_cpn.data{M} = [{id};{b};{b_dag};{b};{b_dag};{(E(M)*Num_0 + U/2 * Num_0*(Num_0-id))}];
            mpo_cpn.N_track{M} = [{0}];
            
            for m = 3:(M-1)
                mpo_cpn.data{m} = [{id},{Z},{Z},{Z};{b},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
                mpo_cpn.N_track{m} = [{0},{-1},{1},{-1},{1},{0}];
            end
            
        end
        
        function mpo_cpn=Simple_1D_Long_Range(mpo_cpn,J,Jdash,U,E)
            % Create MPO for a simple 1D lattice with exponentially
            % decreasing hopping terms.
            
            M = size(mpo_cpn.data,2);
            N_max = mpo_cpn.d-1;
            
            b = sqrt(diag(1:N_max,1)); b_dag = sqrt(diag(1:N_max,-1));
            Num_0 = b_dag*b; id = eye(mpo_cpn.d); Z = zeros(mpo_cpn.d);
            
            mpo_cpn.data{1} = [{(E(1)*Num_0 + U/2*Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
            mpo_cpn.N_track{1} = [{0},{-1},{1},{0}];
            
            mpo_cpn.data{M} = [{id};{b};{b_dag};{(E(M)*Num_0 + U/2 * Num_0*(Num_0-id))}];
            mpo_cpn.N_track{M} = [{0}];
            
            for m = 2:(M-1)
                mpo_cpn.data{m} = [{id},{Z},{Z},{Z};{b},{Jdash/J*id},{Z},{Z};{b_dag},{Z},{Jdash/J*id},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
                mpo_cpn.N_track{m} = [{0},{-1},{1},{0}];
            end
            
        end
        
        function mpo_cpn=SawTooth_1D(mpo_cpn,J,U,E)
            % Create MPO for a 1D SawTooth lattice.
            % 
            %                       A   A   A
            %                    J'/ \ / \ /
            %  SAWTOOTH LATTICE - B---B---B
            %                       J
            %
            
            M = size(mpo_cpn.data,2);
            N_max = mpo_cpn.d-1;
            
            b = sqrt(diag(1:N_max,1)); b_dag = sqrt(diag(1:N_max,-1));
            Num_0 = b_dag*b; id = eye(mpo_cpn.d); Z = zeros(mpo_cpn.d);
            
            mpo_cpn.data{1} = [{(E(1)*Num_0 + U/2*Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
            mpo_cpn.N_track{1} = [{0},{-1},{1},{0}];
            
            mpo_cpn.data{M} = [{id};{b};{b_dag};{(E(M)*Num_0 + U/2 * Num_0*(Num_0-id))}];
            mpo_cpn.N_track{M} = [{0}];
            
            for m = 2:(M-1)
                if (floor(m/2)-m/2)==0
                    mpo_cpn.data{m} = [{id},{Z},{Z},{Z};{b},{Jdash/J*id},{Z},{Z};{b_dag},{Z},{Jdash/J*id},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
                else
                    mpo_cpn.data{m} = [{id},{Z},{Z},{Z};{b},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
                end
                mpo_cpn.N_track{m} = [{0},{-1},{1},{0}];
            end
            
        end
        
        function mpo_cpn=Y_Junction(mpo_cpn,J,U,E,M_L,M_R1,M_R2)
            % Create MPO for a Y-Junction with only nearest
            % neighbour hopping contributions.
            %
            %  MPS/MPO Layout:
            %  LEFT_BRANCH - CENTRAL_SITE - RIGHT_BRANCH_1 - RIGHT_BRANCH_2
            
            M = M_L+M_R1+M_R2+1;
            N_max = mpo_cpn.d-1;
            
            b = sqrt(diag(1:N_max,1)); b_dag = sqrt(diag(1:N_max,-1));
            Num_0 = b_dag*b; id = eye(mpo_cpn.d); Z = zeros(mpo_cpn.d);
            
            mpo_cpn.data{1} = [{(E(1)*Num_0 + U/2*Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
            mpo_cpn.N_track{1} = [{0},{-1},{1},{0}];

            for m = 2:M_L
                mpo_cpn.data{m} = [{id},{Z},{Z},{Z};{b},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
                mpo_cpn.N_track{m} = [{0},{-1},{1},{0}];
            end
            
            mpo_cpn.data{M_L+1} = [{id},{Z},{Z},{Z},{Z},{Z};{b},{Z},{Z},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z},{Z},{Z};{(E(M_L+1)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{-J*b_dag},{-J*b},{id}];
            mpo_cpn.N_track{M_L+1} = [{0},{-1},{1},{-1},{1},{0}];
            
            for m = (M_L+2):M_L+M_R1
                mpo_cpn.data{m} = [{id},{Z},{Z},{Z},{Z},{Z};{b},{Z},{Z},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z},{Z},{Z};{Z},{Z},{Z},{id},{Z},{Z};{Z},{Z},{Z},{Z},{id},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{Z},{Z},{id}];
                mpo_cpn.N_track{m} = [{0},{-1},{1},{-1},{1},{0}];
            end
            
            mpo_cpn.data{M_L+M_R1+1} = [{id},{Z},{Z},{Z};{b},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z};{Z},{id},{Z},{Z};{Z},{Z},{id},{Z};{(E(M_L+M_R1)*Num_0 + U/2 * Num_0*(Num_0-id))},{Z},{Z},{id}];
            mpo_cpn.N_track{M_L+M_R1+1} = [{0},{-1},{1},{0}];
            
            for m = (M_L+M_R1+1+2):(M-1)
                mpo_cpn.data{m} = [{id},{Z},{Z},{Z};{b},{Z},{Z},{Z};{b_dag},{Z},{Z},{Z};{(E(m)*Num_0 + U/2 * Num_0*(Num_0-id))},{-J*b_dag},{-J*b},{id}];
                mpo_cpn.N_track{m} = [{0},{-1},{1},{0}];
            end
            
            mpo_cpn.data{M} = [{id};{b};{bdag};{(E(M)*Num_0 + U/2 * Num_0*(Num_0-id))}];
            mpo_cpn.N_track{M} = [{0}];
            
        end
        
        function mpo_cpn=Eff_Time_Evolve(mpo_cpn,dt,r)
            % Calculate Time evolution operator for Runga-Kutta Method
            % U = exp(-i*dt*H) ~ PROD_r (-i*dt*H-r)
            
            M = size(mpo_cpn.data,2);
            id=eye(mpo_cpn.d);
            
            for m = 1:M
                [alpha,beta]=size(mpo_cpn.data{m});
                
                mpo_cpn.data{m}{alpha,1} = mpo_cpn.data{m}{alpha,1}*(-1i*dt) - id*r/M;
                for count = 2:(beta-1)
                    mpo_cpn.data{m}{alpha,count}=-1i*dt*mpo_cpn.data{m}{alpha,count};
                end
            end

        end
        
    end
    
end

