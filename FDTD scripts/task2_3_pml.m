clear all;
%{
    Used for VIDEO CAPTURING

    Here we implement the FDTD method in 2D
    for the TM problem. In this task we use
    the PML boundary conditions(ONLY LEFT WALL).
%}
e0 = 8.85418781 * 10^(-12);
m0 = 4 * pi * 10^(-7);
c = 1/sqrt(e0*m0);

f = 10^(10);
lambda = c/f;
Xmax = 10 * lambda; % Xmax = Ymax

dx = lambda / 10;   % dx = dy
% By CFL we obtain for free space
dtmax = dx/(c * sqrt(2));
p = 1; % percentage of allowable time step
dt = p * dtmax;

% The following is used in case we want to observe what happens
% when the CFL stability condition is violated

% dt = 1.2 * dx/(c);

%{
    Xmax = 10lambda = 100 dx = Ymax
    The surviving field coordinates are (Ez, Hx, Hy)
    - Ez(x,y) x belongs to [0 dx, 100 dx] and y belongs to [0 dx, 100 dx]
    - Hx(x,y) x belongs to [0 dx, 100 dx] and y belongs to [0.5 dx, 99.5 dx]
    - Hy(x,y) x belongs to [0.5 dx, 99.5 dx] and y belongs to [0 dx, 100 dx]
 
%}
N = round(Xmax/2/dx);
if mod(Xmax, 2*dx) ~= 0
    fprintf('N = Xmax/(2*dx) MUST be an integer\n')
    return
end 

% field coordinates initialization
Ez = zeros(2*N+1, 2*N+1);
Hx = zeros(2*N+1, 2*N);
Hy = zeros(2*N, 2*N+1);

% Material description
e(1:2*N+1, 1:2*N+1) = e0;
sigma = zeros(2*N+1, 2*N+1);
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


%{
    PML and PEC boundary conditions    
%}
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


T = 1/f;
sp_axis = dx*(0:2*N);   % space axis

n1 = round(3*T/dt);  % corresponds to time n1*dt  
n2 = round(10*T/dt); %  >>    >>     >>    n2*dt
n3 = round(12*T/dt); %  >>    >>     >>    n3*dt

% Build the coefficient matrices of the FDTD equations
Ca = (e - 0.5*dt * sigma) ./ (e + 0.5*dt *sigma);
Cb = dt/dx ./ (e + 0.5*dt * sigma);
Da = -dt/m/dx; % m {i, j+0.5}
Db = dt/m/dx; % m {i+0.5,j}


% Prepare figure setup
figure();clf;
h = imagesc(sp_axis,sp_axis,Ez');
set(gca, 'ydir', 'normal');
circle(x0,y0,R);
xlabel('x')
ylabel('y')
grid on
colorbar;


% Create a video writer object
writerObj = VideoWriter('animation.mp4', 'MPEG-4');
writerObj.FrameRate = 10;   % Set the frame rate of the video
% Open the VideoWriter object and create the animation
open(writerObj);

time_steps = n3+1;
for n = 0:time_steps 

    % update Ez
    for i = 2:2*N
        for j = 2:2*N
            Ez(i,j) = Ca(i,j) * Ez(i,j) + ...
            Cb(i,j) * (Hy(i,j) - Hy(i-1,j) + Hx(i,j-1) - Hx(i,j)); 
        end
    end
    % add source
    Ez(N+1,N+1) = source(n*dt);
    % PML A.B.C.
    for j = 2:2*N % we exclude nodes (1,1) and (1,2*N+1) due to Hx indices
        Ez(1,j) = Ca(1,j) * Ez(1,j) + ...
        Cb(1,j) * (Hy(1,j) - Hy_pml(Npml,j) + Hx(1,j-1) - Hx(1,j));
    end
    for i = 2:Npml
        for j = 2:2*N
            % Ezx_pml
            Ezx_pml(i,j) = Ca_pml(i) * Ezx_pml(i,j) + ...
            Cb_pml(i) * (Hy_pml(i,j) - Hy_pml(i-1,j));
            % Ezy_pml
            Ezy_pml(i,j) = Ezy_pml(i,j) - dt/e0/dx * (Hx(i,j) - Hx(i,j-1));
        end
    end

    % update Hx
    for i = 2:2*N
        for j = 1:2*N
            Hx(i,j) = Hx(i,j) + Da * (Ez(i,j+1) - Ez(i,j));
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
            Hy(i,j) = Hy(i,j) + Db * (Ez(i+1,j) - Ez(i,j));
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
        Dby_pml(Npml) * (Ez(1,j)  - Ezx_pml(Npml,j) - Ezy_pml(Npml,j));
    end
    
    set(h, 'CData', Ez');
    title(sprintf('Electric Field $E_z$ (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
    pause(0.1);
    frame = getframe(gcf);   % Capture the current frame of the figure
    writeVideo(writerObj, frame);   % Write the frame to the video file
end

% Close the VideoWriter object
close(writerObj);