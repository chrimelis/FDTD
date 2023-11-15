clear all;
%{
    Used for PLOTTING

    This File is used for ERROR analysis of the FDTD method
    We implement the FDTD method in 2D
    for the TM problem. 
    We compute each field component 4 times:
        1) The real field (using an extended space without reflections)
        2) The field after implementing P.E.C. boundary conditions
        3) The field  >>     >>         Mur absorbing b.c.(only the LEFT boundary the rest with P.E.C.)
        4) The field  >>     >>         PML (only the LEFT boundary the rest with P.E.C.)
%}
% Constants definition
e0 = 8.85418781 * 10^(-12);
m0 = 4 * pi * 10^(-7);
c = 1/sqrt(e0*m0);

f = 10^(10);
lambda = c/f;
Xmax = 10 * lambda; % Xmax = Ymax

dx = lambda / 10;   % dx = dy
% By CFL we obtain for free space
dtmax = dx/(c * sqrt(2));
p = 0.8; % percentage of allowable time step
dt = p * dtmax;

% The following commented out code, is used in case we want to observe what happens
% when the CFL stability condition is violated

% dt = 1.2 * dx/(c);

%{
    Xmax = 10lambda = 100 dx = Ymax
    The surviving field coordinates are (Ez, Hx, Hy)
    - Ez(x,y) x belongs to [0 dx, 100 dx] and y belongs to [0 dx, 100 dx]
    - Hx(x,y) x belongs to [0 dx, 100 dx] and y belongs to [0.5 dx, 99.5 dx]
    - Hy(x,y) x belongs to [0.5 dx, 99.5 dx] and y belongs to [0 dx, 100 dx]
    
    N indicates the # of nodes that exist in a row/column of the grid
    above the central node (Xmax/2,Ymax/2) <-> (N+1, N+1) 
%}
N = round(Xmax/2/dx);
if mod(Xmax, 2*dx) ~= 0
    fprintf('N = Xmax/(2*dx) MUST be an integer\n')
    return
end 

% field coordinates initialization

% We denote r_A the entity A for the extended computational space 
% >>   >>   m_A  >>    >>    >>  >>  field caused with Mur conditions(left wall only)
% >>   >>   p_A  >>    >>    >>  >>  field caused with PML conditions(left wall only)

ext_f = 2; % extension factor  

Ez = zeros(2*N+1, 2*N+1); r_Ez = zeros((2+2*ext_f)*N+1,(2+2*ext_f)*N+1); 
m_Ez = zeros(2*N+1, 2*N+1); p_Ez = zeros(2*N+1, 2*N+1);

Hx = zeros(2*N+1, 2*N);   r_Hx = zeros((2+2*ext_f)*N+1, (2+2*ext_f)*N);
m_Hx = zeros(2*N+1, 2*N); p_Hx = zeros(2*N+1, 2*N);

Hy = zeros(2*N, 2*N+1);   r_Hy = zeros((2+2*ext_f)*N, (2+2*ext_f)*N+1);
m_Hy = zeros(2*N, 2*N+1); p_Hy = zeros(2*N, 2*N+1);

% Material description
e(1:2*N+1, 1:2*N+1) = e0; r_e(1:(2+2*ext_f)*N+1, 1:(2+2*ext_f)*N+1) = e0;
sigma = zeros(2*N+1, 2*N+1); r_sigma = zeros((2+2*ext_f)*N+1, (2+2*ext_f)*N+1);
m = m0;

% Scatterer
x0 = Xmax/2 + 3*lambda; % x coordinate of the center of the cylinder
y0 = Xmax/2;    % y coordinate of the center of the cylinder
R = lambda;     % the radius of the cylinder
sig = 1.2;  % sigma: the conductivity
er = 3.4;   % the relative dielectric constant

for i = 1:length(e)
    for j = 1:length(e)
        if ( ( (i-1)*dx - x0 )^2 + ( (j-1)*dx - y0 )^2 ) <= R^2
            sigma(i,j) = sig;
            e(i,j) = e(i,j) * er;
        end 
    end
end
% extended space material definition
for i = 1:length(r_e)
    for j = 1:length(r_e)
        if (i > ext_f*N) && (i < (2+ext_f)*N + 1) && (j > ext_f * N) && (j < (2+ext_f)*N+1)
            if ( ( (i- ext_f * N - 1)*dx - x0 )^2 + ( (j - ext_f * N - 1)*dx - y0 )^2 ) <= R^2
                r_sigma(i,j) = sig;
                r_e(i,j) = r_e(i,j) * er;
            end
        end 
    end
end

T = 1/f;
sp_axis = dx*(0:2*N);   % space axis
r_sp_axis = (-ext_f*N:(2+ext_f)*N)*dx;  % extended space axis

n1 = round(3*T/dt);  % corresponds to time n1*dt  
n2 = round(10*T/dt); %  >>    >>     >>    n2*dt
n3 = round(12*T/dt); %  >>    >>     >>    n3*dt

% Build the coefficient matrices of the FDTD equations
% Again r_A denotes the entity A in the extended computational space
Ca = (e - 0.5*dt * sigma) ./ (e + 0.5*dt *sigma); 
r_Ca = (r_e - 0.5*dt * r_sigma) ./ (r_e + 0.5*dt * r_sigma); 

Cb = dt/dx ./ (e + 0.5*dt * sigma); 
r_Cb = dt/dx ./ (r_e + 0.5*dt * r_sigma);

Da = -dt/m/dx; % m {i, j+0.5}
Db = dt/m/dx; % m {i+0.5,j}

%{
    Boundary conditions preparation    
%}

% Mur Absorbing Boundary Conditions (ONLY LEFT BOUNDARY)
Ez1_t1 = zeros(1,2*N+1);    % m_Ez(1,j) 1 time step  ago
Ez1_t2 = zeros(1,2*N+1);    % m_Ez(1,j) 2 time steps ago
Ez2_t1 = zeros(1,2*N+1);    % m_Ez(2,j) 1 time step  ago
Ez2_t2 = zeros(1,2*N+1);    % m_Ez(2,j) 2 time steps ago

Hx1_t1 = zeros(1,2*N);    % m_Hx(1,j) 1 time step  ago
Hx1_t2 = zeros(1,2*N);    % m_Hx(1,j) 2 time steps ago
Hx2_t1 = zeros(1,2*N);    % m_Hx(2,j) 1 time step  ago
Hx2_t2 = zeros(1,2*N);    % m_Hx(2,j) 2 time steps ago
% Hy is not defined in the LEFT boundary (Yee cell)

% PML Absorbing Boundary Condition Setup(ONLY LEFT WALL)
Npml = 8;   % PML width
O = 2;  % order
Ref = 10^(-6);    % reflection coeff. for normal incidence

Hx_pml = zeros(Npml, 2*N);
Hy_pml = zeros(Npml, 2*N+1);
Ezx_pml = zeros(Npml, 2*N+1);
Ezy_pml = zeros(Npml, 2*N+1);
% PML conductivity matrices
se = -e0*c*log(Ref)/(2^(O+2)*dx*Npml^(O+1));
sh = se*m0/e0;  % by PML matching condition

sigmaE = zeros(1,Npml);
sigmaHy = zeros(1,Npml);
for i = 1:Npml
    sigmaE(i)  = se * ( (2*i+1)^(O+1) - (2*i-1)^(O+1) ); 
    sigmaHy(i) = sh * ( (2*(i-0.5)+1)^(O+1)-(2*(i-0.5)-1)^(O+1) );
end
sigmaE = fliplr(sigmaE);
sigmaHy = fliplr(sigmaHy);

% Build the coefficient matrices for PML equations
% for Ezx
Ca_pml = exp(1).^(-sigmaE * dt/e0);
Cb_pml = (1-Ca_pml)./(sigmaE * dx);
% for Hy
Day_pml=exp(1).^(-sigmaHy * dt/m0);
Dby_pml=(1-Day_pml)./(sigmaHy * dx);


e_pec = zeros(1,n3+2);  % time series of absolute error at node (1,N+1) PEC case
e_m = zeros(1,n3+2);    % time series of absolute error at node (1,N+1) MUR case
e_p = zeros(1,n3+2);    % time series of absolute error at node (1,N+1) PML case

time_steps = n3+1;
for n = 0:time_steps 
    % prepare left wall for Mur A.B.C.

    % m_Ez
    Ez1_t2 = Ez1_t1; 
    Ez1_t1 = m_Ez(1, 1:2*N+1);
    Ez2_t2 = Ez2_t1;
    Ez2_t1 = m_Ez(2, 1:2*N+1);
    % m_Hx
    Hx1_t2 = Hx1_t1; 
    Hx1_t1 = m_Hx(1, 1:2*N);
    Hx2_t2 = Hx2_t1;
    Hx2_t1 = m_Hx(2, 1:2*N);

    % update Ez
    for i = 2:2*N
        for j = 2:2*N
            % P.E.C. case
            Ez(i,j) = Ca(i,j) * Ez(i,j) + ...
            Cb(i,j) * (Hy(i,j) - Hy(i-1,j) + Hx(i,j-1) - Hx(i,j));
            % Mur & P.E.C. case
            m_Ez(i,j) = Ca(i,j) * m_Ez(i,j) + ...
            Cb(i,j) * (m_Hy(i,j) - m_Hy(i-1,j) + m_Hx(i,j-1) - m_Hx(i,j));
            % PML & P.E.C. case
            p_Ez(i,j) = Ca(i,j) * p_Ez(i,j) + ...
            Cb(i,j) * (p_Hy(i,j) - p_Hy(i-1,j) + p_Hx(i,j-1) - p_Hx(i,j)); 
        end
    end
    % add source
    Ez(N+1,N+1) = source(n*dt);
    m_Ez(N+1,N+1) = source(n*dt);
    p_Ez(N+1,N+1) = source(n*dt);

    % update r_Ez
    for i = 2:(2+2*ext_f)*N
        for j = 2:(2+2*ext_f)*N
            r_Ez(i,j) = r_Ca(i,j) * r_Ez(i,j) + ...
            r_Cb(i,j) * (r_Hy(i,j) - r_Hy(i-1,j) + r_Hx(i,j-1) - r_Hx(i,j)); 
        end
    end
    % add source
    r_Ez((ext_f+1)*N+1,(ext_f+1)*N+1) = source(n*dt);

    % Mur A.B.C.
    for j = 1:2*N+1
        if (j == 1) || (j == 2*N+1)
            % Use 1st order Mur ABC
            m_Ez(1,j) = Ez2_t1(j) - (dx - c*dt)/(dx + c*dt) * (m_Ez(2,j) - Ez1_t1(j));
        else
            % Use 2nd order Mur ABC
            m_Ez(1,j) = -Ez2_t2(j) - (dx - c*dt)/(dx + c*dt) * (m_Ez(2,j) + Ez1_t2(j)) + ...
            2*dx/(dx + c*dt) * (Ez1_t1(j) + Ez2_t1(j)) + (c*dt)^2 / (2*dx*(dx+ c*dt)) * ...
            ( Ez1_t1(j+1) - 2*Ez1_t1(j) + Ez1_t1(j-1) + Ez2_t1(j+1) - 2*Ez2_t1(j) + Ez2_t1(j-1) );
        end
    end

    % PML A.B.C.
    for j = 2:2*N % we exclude nodes (1,1) and (1,2*N+1) due to Hx indices
        p_Ez(1,j) = Ca(1,j) * p_Ez(1,j) + ...
        Cb(1,j) * (p_Hy(1,j) - Hy_pml(Npml,j) + p_Hx(1,j-1) - p_Hx(1,j));
    end
    for i = 2:Npml
        for j = 2:2*N
            % Ezx_pml
            Ezx_pml(i,j) = Ca_pml(i) * Ezx_pml(i,j) + ...
            Cb_pml(i) * (Hy_pml(i,j) - Hy_pml(i-1,j));
            % Ezy_pml
            Ezy_pml(i,j) = Ezy_pml(i,j) - dt/e0/dx * (p_Hx(i,j) - p_Hx(i,j-1));
        end
    end

    % update Hx
    for i = 2:2*N
        for j = 1:2*N
            % Only P.E.C. case
            Hx(i,j) = Hx(i,j) + Da * (Ez(i,j+1) - Ez(i,j));
            % Mur and P.E.C. case
            m_Hx(i,j) = m_Hx(i,j) + Da * (m_Ez(i,j+1) - m_Ez(i,j));
            % PML and P.E.C. case
            p_Hx(i,j) = p_Hx(i,j) + Da * (p_Ez(i,j+1) - p_Ez(i,j));
        end
    end

    % update r_Hx
    for i = 2:(2+2*ext_f)*N
        for j = 1:(2+2*ext_f)*N
            r_Hx(i,j) = r_Hx(i,j) + Da * (r_Ez(i,j+1) - r_Ez(i,j));
        end
    end
    % Mur A.B.C.
    for j = 1:2*N
        if (j == 1) || (j == 2*N)
            % Use 1st order Mur ABC
            m_Hx(1,j) = Hx2_t1(j) - (dx - c*dt)/(dx + c*dt) * (m_Hx(2,j) - Hx1_t1(j));
        else
            % Use 2nd order Mur ABC
            m_Hx(1,j) = -Hx2_t2(j) - (dx - c*dt)/(dx + c*dt) * (m_Hx(2,j) + Hx1_t2(j)) + ...
            2*dx/(dx + c*dt) * (Hx1_t1(j) + Hx2_t1(j)) + (c*dt)^2 / (2*dx*(dx+ c*dt)) * ...
            ( Hx1_t1(j+1) - 2*Hx1_t1(j) + Hx1_t1(j-1) + Hx2_t1(j+1) - 2*Hx2_t1(j) + Hx2_t1(j-1) );
        end
    end
    % PML A.B.C.
    for i = 2:Npml
        for j = 1:2*N
            Hx_pml(i,j) = Hx_pml(i,j) - dt/m/dx * (Ezx_pml(i,j+1) + Ezy_pml(i,j+1) - Ezx_pml(i,j) - Ezy_pml(i,j));
        end
    end


    % update Hy
    for i = 1:2*N
        for j = 2:2*N
            % Only P.E.C. case
            Hy(i,j) = Hy(i,j) + Db * (Ez(i+1,j) - Ez(i,j));
            % Mur and P.E.C. case
            m_Hy(i,j) = m_Hy(i,j) + Db * (m_Ez(i+1,j) - m_Ez(i,j));
            % PML and P.E.C. case
            p_Hy(i,j) = p_Hy(i,j) + Db * (p_Ez(i+1,j) - p_Ez(i,j));
        end
    end
    % update r_Hy
    for i = 1:(2+2*ext_f)*N
        for j = 2:(2+2*ext_f)*N
            r_Hy(i,j) = r_Hy(i,j) + Db * (r_Ez(i+1,j) - r_Ez(i,j));
        end
    end
    % PML A.B.C.
    for i = 1:(Npml-1)
        for j = 2:2*N
            Hy_pml(i,j) = Day_pml(i) * Hy_pml(i,j) + ...
            Dby_pml(i) * (Ezx_pml(i+1,j) + Ezy_pml(i+1,j) - Ezx_pml(i,j) - Ezy_pml(i,j));
        end
    end
    for j = 2:2*N
        Hy_pml(Npml,j) = Day_pml(Npml) * Hy_pml(Npml,j) + ...
        Dby_pml(Npml) * (p_Ez(1,j)  - Ezx_pml(Npml,j) - Ezy_pml(Npml,j));
    end

    % 1) true Ez @ time n*dt
    true_Ez = r_Ez( (ext_f*N+1):((2+ext_f)*N+1), (ext_f*N+1):((2+ext_f)*N+1));
    % 2) P.E.C. Ez is stored in Ez
    % 3) Mur (and P.E.C.) Ez is stored in m_Ez
    % 4) PML (and P.E.C.) Ez is stored in p_Ez

    quant = 10;
    a1 = abs(true_Ez(1,N+1-quant:N+1+quant) - Ez(1,N+1-quant:N+1+quant));
    a2 = abs(true_Ez(1,N+1-quant:N+1+quant) - m_Ez(1,N+1-quant:N+1+quant));
    a3 = abs(true_Ez(1,N+1-quant:N+1+quant) - p_Ez(1,N+1-quant:N+1+quant));

    if n == n1 || n == n2 || n == n3
        figure();clf;
        h3 = imagesc(sp_axis, sp_axis, Ez');
        set(gca, 'ydir', 'normal');
        hold on
        circle2(x0,y0,R);
        hold off
        title(sprintf('(P.E.C.) $E_z$ (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
        xlabel('x')
        ylabel('y')
        grid on
        colorbar;

        figure();clf;
        h4 = imagesc(sp_axis, sp_axis, m_Ez');
        set(gca, 'ydir', 'normal');
        hold on
        circle2(x0,y0,R);
        hold off
        title(sprintf('(MUR) $E_z$ (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
        xlabel('x')
        ylabel('y')
        grid on
        colorbar;

        figure();clf;
        h5 = imagesc(sp_axis, sp_axis, p_Ez');
        set(gca, 'ydir', 'normal');
        hold on
        circle2(x0,y0,R);
        hold off
        title(sprintf('(PML) $E_z$ (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
        xlabel('x')
        ylabel('y')
        grid on
        colorbar;

        figure();clf;
        h6 = imagesc(sp_axis, sp_axis, true_Ez');
        set(gca, 'ydir', 'normal');
        hold on
        circle2(x0,y0,R);
        hold off
        title(sprintf('(True Values) $E_z$ (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
        xlabel('x')
        ylabel('y')
        grid on
        colorbar;

        figure();clf; 
        stem(N+1-quant:N+1+quant, a1);
        axis tight;
        xlabel(sprintf('Index $j$ of Node $(1,j)$'), 'Interpreter', 'latex');
        ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
        title(sprintf('(P.E.C.) Absolute Error (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    

        figure();clf; 
        stem(N+1-quant:N+1+quant, a2);
        axis tight;
        title(sprintf('(MUR) Absolute Error (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
        xlabel(sprintf('Index $j$ of Node $(1,j)$'), 'Interpreter', 'latex');
        ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');

        figure();clf; 
        stem(N+1-quant:N+1+quant, a3);
        axis tight;
        title(sprintf('(PML) Absolute Error (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
        xlabel(sprintf('Index $j$ of Node $(1,j)$'), 'Interpreter', 'latex');
        ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
    end
    e_pec(n+1) = a1(quant + 1);
    e_m(n+1) =   a2(quant + 1);
    e_p(n+1) =   a3(quant + 1);
end

t_axis = (0:n3+1)*dt/T; % normalized time axis 

% P.E.C.
figure(); clf; plot(t_axis, e_pec);
ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
xlabel(sprintf('Normalized Time $t/T_0$'), 'Interpreter', 'latex');
title(sprintf('(P.E.C.) Absolute Error @ node $(1, %d)$',N+1), 'Interpreter', 'latex');
grid on; 
% Log P.E.C.
figure(); clf; plot(t_axis, e_pec);
ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
xlabel(sprintf('Normalized Time $t/T_0$'), 'Interpreter', 'latex');
title(sprintf('(P.E.C.) Absolute Error @ node $(1, %d)$',N+1), 'Interpreter', 'latex');
grid on; 
set(gca, 'YScale', 'log')

% 2nd order Mur
figure(); clf; plot(t_axis, e_m);
ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
xlabel(sprintf('Normalized Time $t/T_0$'), 'Interpreter', 'latex');
title(sprintf('(Mur) Absolute Error @ node $(1, %d)$',N+1), 'Interpreter', 'latex');
grid on; 
% Log MUR
figure(); clf; plot(t_axis, e_m);
ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
xlabel(sprintf('Normalized Time $t/T_0$'), 'Interpreter', 'latex');
title(sprintf('(Mur) Absolute Error @ node $(1, %d)$',N+1), 'Interpreter', 'latex');
grid on; 
set(gca, 'YScale', 'log')

% PML
figure(); clf; plot(t_axis, e_p);
ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
xlabel(sprintf('Normalized Time $t/T_0$'), 'Interpreter', 'latex');
title(sprintf('(PML) Absolute Error @ node $(1, %d)$',N+1), 'Interpreter', 'latex');
grid on; 
% Log PML
figure(); clf; plot(t_axis, e_p);
ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
xlabel(sprintf('Normalized Time $t/T_0$'), 'Interpreter', 'latex');
title(sprintf('(PML) Absolute Error @ node $(1, %d)$',N+1), 'Interpreter', 'latex');
grid on; 
set(gca, 'YScale', 'log')

% combined error graph
figure(); clf; plot(t_axis, e_m);
ylabel(sprintf('Absolute Error'), 'Interpreter', 'latex');
xlabel(sprintf('Normalized Time $t/T_0$'), 'Interpreter', 'latex');
title(sprintf('Absolute Error @ node $(1, %d)$',N+1), 'Interpreter', 'latex');
grid on; 
set(gca, 'YScale', 'log')
hold on
plot(t_axis, e_p);
hold off
legend('MUR','PML')