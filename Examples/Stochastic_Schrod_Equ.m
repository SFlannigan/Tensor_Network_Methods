%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Two-Level Atom - Solution from Stochastic Schrodinger Equ. and 
%%  Continuous Measurement
%%
%%  Units: a = hbar = 1
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all


clear; clc; close all

save_data='true';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Parameters
count=1;
% for V_0 = 5:30%10%[10,16,20]
    V_0=10; % Lattice Potential Depth
    N=500; % Number of Real space points in Wannier function calculation
    a_s=350; % Characteristic Scattering Length in units of Bohr radius (a_0)
    a = 256e-9; % Lattice Spacing
    a_0 = 5.29e-11; % Bohr Radius
%     for a_s_loop = 350%-200:10:200
    a_s_loop=200;
        tic
        a_s = a_s_loop*a_0/a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plancks constant
        hbar = 1.0545718e-34;
        %Mass of Ceasium Atom
        M = 132.9*1.6726e-27;
        %Time Unit
        time_unit = M*a^2/hbar;
                
        s_0 = 1.00624;
        n_minus = 0.06;
        n_plus = 0.06;
        a_minus = -850; %a_0
        a_plus = 1060; %a_0
        
        if a_s_loop<0
            C_a = 4590*sinh(2*n_minus)/(sin(s_0*log(a_s_loop/a_minus))^2+sinh(n_minus)^2);
            L_3 = 3*C_a*a_s^4;
        elseif a_s_loop==0
            L_3 = 0;
        else
            C_a = 67.1*exp(-2*n_plus)*(cos(s_0*log(a_s_loop/a_plus))^2+sinh(n_plus)^2) + 16.8*(1-exp(-4*n_plus));
            L_3 = 3*C_a*a_s^4;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find Hubbard-Model Coefficients
        [J(count),U(count),gamma_3(count)] = Get_Coefficients(V_0,a_s,L_3,N);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_dag = diag(sqrt(1:4),-1);
b = diag(sqrt(1:4),1);
num = b_dag*b;
id = eye(5);
H_BH = -J(count)*kron(b_dag,b) -J(count)*kron(b,b_dag) + U(count)/2*kron(num*(num-id),id) + U(count)/2*kron(id,num*(num-id));
H = H_BH - 1i*gamma_3(count)/12*(kron(b_dag*b_dag*b_dag*b*b*b,id)+kron(id,b_dag*b_dag*b_dag*b*b*b));

Jump_1 = kron(b*b*b,id);
Jump_2 = kron(id,b*b*b);


time_step=1;
Tf=5000; %gamma*T

Total_Ex = 10000;

O = zeros(Tf/time_step+1,Total_Ex);
O_2 = zeros(Tf/time_step+1,Total_Ex);
O_3 = zeros(Tf/time_step+1,Total_Ex);
O_4 = zeros(Tf/time_step+1,Total_Ex);
time_plot=zeros(Tf/time_step+1,1);
for ex = 1:Total_Ex
tic
d = 5^2;
Psi = zeros(d,1); %(d'),(d)
Psi(13,1) = 1; % Begin in State of two particles on each site

time=0;
count=1;
while time<(Tf-time_step)

    if time < Tf-time_step
        Psi_N=Psi;
        n = rand;
        U = expm(-1i*H*time_step);
        for t = time:time_step:Tf
            Psi_temp = U*Psi;
            Psi = Psi_temp/sqrt(Psi_temp'*Psi_temp);
            Psi_N = U*Psi_N;
            
            N_1 = kron(num,id);
            N_2 = kron(id,num);
            
            O(count,ex) = abs(Psi'*N_1*Psi) + abs(Psi'*N_2*Psi);
            
            time_plot(count)=t;
            count=count+1;
            if Psi_N'*Psi_N<=n
                count=count-1;
                break
            end
            
        end
        time=t;
        
        p1 = Jump_1*Psi;
        P_1 = p1'*p1;
        
        p2 = Jump_2*Psi;
        P_2 = p2'*p2;
        
        if P_1>P_2
            Psi = p1/sqrt(P_1);
        else
            Psi = p2/sqrt(P_2);
        end
    end
end

if (floor(ex/20)-ex/20)==0
   display(['Iteration ',num2str(ex)]); 
end

end

O_av = mean(O,2);

figure; 
plot(time_plot,O_av);
xlabel('\gammat','fontsize',20);
% title(['Continuous Measurement - Total Runs=',num2str(Total_Ex),':    \Omega=',num2str(omega),'; \gamma=',num2str(gamma),'; \Delta=',num2str(delta)],'fontsize',20);
% legend('\rho_{22}','\rho_{11}','\rho_{12}')
toc

