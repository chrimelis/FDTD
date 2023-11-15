clear all;
%{
    This File is used for VIDEO CREATION
    Shows ONLY desired region

    Here we implement the FDTD method in 2D
    for the TM problem. 
        We compute the field in an EXTENDED SPACE
    but plot it only in the desired region.
        We run the program for a specific time interval such
    that no reflections occur from the P.E.C. boundaries
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

% We denote r_A the entity A for the extended computational space 
ext_f = 2; % extension factor  

Ez = zeros(2*N+1, 2*N+1); r_Ez = zeros((2+2*ext_f)*N+1,(2+2*ext_f)*N+1); 
Hx = zeros(2*N+1, 2*N);   r_Hx = zeros((2+2*ext_f)*N+1, (2+2*ext_f)*N);
Hy = zeros(2*N, 2*N+1);   r_Hy = zeros((2+2*ext_f)*N, (2+2*ext_f)*N+1);

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
r_sp_axis = (-ext_f*N:(2+ext_f)*N)*dx;  % extended space axis

%{
    PEC boundary conditions    
%}
T = 1/f;
sp_axis = dx*(0:2*N);   % space axis

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



% Prepare figure setup
figure();clf;
h = imagesc(sp_axis,sp_axis,r_Ez( (ext_f*N+1):((2+ext_f)*N+1), (ext_f*N+1):((2+ext_f)*N+1))');
set(gca, 'ydir', 'normal');
hold on
circle2(x0,y0,R);
rectangle('Position', [0, 0, Xmax, Xmax], 'EdgeColor', 'k', 'LineWidth',2);
hold off
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
    % update r_Ez
    for i = 2:(2+2*ext_f)*N
        for j = 2:(2+2*ext_f)*N
            r_Ez(i,j) = r_Ca(i,j) * r_Ez(i,j) + ...
            r_Cb(i,j) * (r_Hy(i,j) - r_Hy(i-1,j) + r_Hx(i,j-1) - r_Hx(i,j)); 
        end
    end
    % add source
    r_Ez((ext_f+1)*N+1,(ext_f+1)*N+1) = source(n*dt);

    % update Hx
    for i = 2:N
        for j = 1:2*N
            Hx(i,j) = Hx(i,j) + Da * (Ez(i,j+1) - Ez(i,j));
        end
    end
    % update r_Hx
    for i = 2:(2+2*ext_f)*N
        for j = 1:(2+2*ext_f)*N
            r_Hx(i,j) = r_Hx(i,j) + Da * (r_Ez(i,j+1) - r_Ez(i,j));
        end
    end
    

    % update Hy
    for i = 1:2*N
        for j = 2:2*N
            Hy(i,j) = Hy(i,j) + Db * (Ez(i+1,j) - Ez(i,j));
        end
    end
    % update r_Hy
    for i = 1:(2+2*ext_f)*N
        for j = 2:(2+2*ext_f)*N
            r_Hy(i,j) = r_Hy(i,j) + Db * (r_Ez(i+1,j) - r_Ez(i,j));
        end
    end
    
    set(h, 'CData', r_Ez( (ext_f*N+1):((2+ext_f)*N+1), (ext_f*N+1):((2+ext_f)*N+1))');
    title(sprintf('Electric Field $E_z$ (t = %.2f $T_0$)', n * dt/T), 'Interpreter', 'latex');    
    pause(0.001);
    frame = getframe(gcf);   % Capture the current frame of the figure
    writeVideo(writerObj, frame);   % Write the frame to the video file
end

% Close the VideoWriter object
close(writerObj);
