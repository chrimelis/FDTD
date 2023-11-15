clear all;

%{ 
    PART I
    Used for validating theoretical results for the FDTD method
    
    Definition of the wavenumber, frequency
    wavelength in free space(or air), and
    space and time discretization steps 
%} 
k = 200 * pi / 3;
lambda = 0.03;
f = 10^(10);
dx = lambda / 10;
dt = 7.07 * 10^(-12);
c = 3 * 10^8;

% The angles in degrees
theta = linspace(0,90, 9000);
% up_normalized = (arithmetic u_phase) / c
up_normalized = 1 / (pi * f * dt ) * asin( c*dt/dx * ...
    sqrt( (sin( k * cos(pi*theta/180) * dx/2 )).^2 + ...
          (sin( k * sin(pi*theta/180) * dx/2 )).^2 ) );

% Plot commands
figure();clf;
plot(theta,up_normalized);
title('Normalized Arithmetic Phase Velocity')
xlabel('\theta')
ylabel('u_{p,ar}/c')
grid on

