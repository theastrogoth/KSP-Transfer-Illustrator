function [E] = solveKepler(M,e)
% Solves Kepler's Equation for eccentric or hyperbolic anomaly via the
% Newton-Raphson method
% method
% Reference: https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
% Inputs:   M, true anomaly (rad)
%           e, eccentricity (-/-)
% Outputs:  E, eccentric anomaly (rad)
%       or  E, hyperbolic anomaly (rad)

    tolerance = 10E-6;
    iteration = 0;
    maxIteration = 1000;

if e==1     % parabolic case
    error('Parabolic case (e=1) not valid')
elseif e<1  % elliptical case
    if e < 0.8
        E = M;     %first guess before iterating
    else
        E = pi;
    end
    Eprev = (M+pi)*4;
    while abs((E-Eprev)/E)>tolerance
        Eprev = E;
        E = E - (E - e*sin(E) - M)/(1-e*cos(E));
        iteration = iteration+1;
        if iteration>maxIteration
            disp('Iteration limit exceeded')
            break
        end
    end
elseif e>1  % hyperbolic case
    if abs(M)>4*pi
        E = 4*pi*sign(M);
    else
        E = M;
    end
    Eprev = (M+pi)*4;
    while abs((E-Eprev)/E)>tolerance
        Eprev = E;
        E = E + (M - e*sinh(E) + E)/(e*cosh(E)-1);
        iteration = iteration+1;
        if iteration>maxIteration
            disp('Iteration limit exceeded')
            break
        end
    end    
end
end