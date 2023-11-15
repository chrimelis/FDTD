function y = source(t)
% Hard (sinusoidal) source for the FDTD method
% returns the value of Ez of the central node of the mesh  
    f0 = 10^(10);
    y = sin(2*pi*f0*t);    
end