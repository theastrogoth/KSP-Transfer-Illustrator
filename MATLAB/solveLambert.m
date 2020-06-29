function [transferOrbit,transferOrbitPC, planeChangeDV, planeChangeDT] = ...
    solveLambert(orbit1,orbit2,tStart,tFlight,planeChange)
% SolveLambert uses an iterative technique to solve Lambert's problem, 
% determining an orbit from two positions and a time interval
% Reference: http://www.braeunig.us/space/index.htm
% Inputs:   orbit1, Keplerian orbit at start (m)
%           orbit2, Keplerian orbit at finish (m)
%           tStart, flight start time (s)
%           tFlight, flight duration (s)
%           planeChange, true if a plane change occurs during transfer
% Outputs:  transferOrbit, Keplerian orbit from start to finish (or plane
%               change)
%           transferOrbitPC, Keplerian orbit from plane change to finish
if nargin<5
   planeChange = false; 
end
mu = orbit1.primary.mu;
    % gravity parameter of the orbited body
r1 = orbit1.getStateVector(tStart);
r2 = orbit2.getStateVector(tStart+tFlight);
    % positions of start and end bodies at departure and arrival times
if planeChange
   r2 = orbit1.toOrbitBasis(r2);
   r2OutOfPlane = r2(3);
   r2 = [r2(1);r2(2);0];
   r2 = orbit1.toPrimaryBasis(r2);
end
dnu = atan2(norm(cross(r1,r2)),dot(r1,r2));
    % change in true anomaly (rad)
theta1 = atan2(r1(2),r1(1));
theta2 = atan2(r2(2),r2(1));
    % angles in the ecliptic (xy) plane
dtheta = theta2-theta1;
while dtheta<0
    dtheta = dtheta+2*pi;
end
if dtheta > pi
    dnu = 2*pi-dnu;
end
% constants
k = norm(r1)*norm(r2)*(1-cos(dnu));
l = norm(r1) + norm(r2);
m = norm(r1)*norm(r2)*(1+cos(dnu));
% limiting p values (two parabolas)
pj = k/(l+(2*m)^0.5);
pjj = k/(l-(2*m)^0.5);
if dnu > pi
    pMin = 0;
    pMax = pjj;
else
    pMin = pj;
    pMax = Inf;
end
tolerance = 1E-6;
maxIteration = 200;
error = 1;
iteration = 0;
p = (pj+pjj)/2;
pNext = p;
while error>tolerance
    if iteration>maxIteration
        disp('Failure to converge')
        break
    end
    p = pNext;
    iteration = iteration+1;
    a = m*k*p/((2*m-(l^2))*(p^2)+2*k*l*p-k^2);
    f = 1-norm(r2)/p*(1-cos(dnu));
    g = norm(r1)*norm(r2)*sin(dnu)/(mu*p)^0.5;
    df = ((mu/p).^0.5)*tan(dnu/2)*((1-cos(dnu))/p - 1/norm(r1) - 1/norm(r2));
    if a>0
        sindE = -norm(r1)*norm(r2)*df/(mu*a)^0.5;
        cosdE = 1-norm(r1)/a*(1-f);
        dE = atan2(sindE,cosdE);
        while dE < 0
            dE = dE + 2*pi;
        end
        t = g+(((a^3)/mu)^0.5)*(dE-sindE);
        dtdp = -g/(2*p)-1.5*a*(t-g)*(k^2+(2*m-l^2)*p^2)/(m*k*p^2) + ...
            ((a^3/mu)^0.5)*2*k*sin(dE)/(p*(k-l*p));
    else
        dF = acosh(1-norm(r1)/a*(1-f));
        t = g+(((((-a)^3)/mu)^0.5)*(sinh(dF)-dF));
        dtdp = -g/(2*p)-1.5*a*(t-g)*(k^2+(2*m-l^2)*p^2)/(m*k*p^2) - ...
            (((-a)^3/mu)^0.5)*2*k*sinh(dF)/(p*(k-l*p));
    end
    error = abs(tFlight-t)/tFlight;
    pNext = p + (tFlight-t)/dtdp;
    if pNext < pMin
        pNext = pMin*1.000001;
    elseif pNext > pMax
        pNext = pMax*0.999999;
    end
end

v1 = (r2-f*r1)/g;
transferOrbit = Orbit.fromStateVector(r1,v1,tStart,orbit1.primary);

if planeChange
    transferAngle = transferOrbit.getAngleInEcliptic(tStart,r2);
    if transferAngle < 0
       transferAngle = transferAngle +2*pi; 
    end
    if transferAngle < pi/2
        thetaPC = 0;
    else
        thetaPC = transferAngle-pi/2;
    end
        % angle from ejection of plane change burn
    nuPC =  transferOrbit.getTrueAnomaly(tStart) + thetaPC;
    tPC = transferOrbit.getTime(nuPC);
    [rPC,vPCi] = transferOrbit.getStateVector(tPC);
        % state vector at plane change burn
    iPC = atan2(r2OutOfPlane,norm(r2));
        % inclination change needed
    vPCiEcliptic = transferOrbit.toOrbitBasis(vPCi);
    vPCfEcliptic = norm(vPCiEcliptic) * [cos(iPC)*vPCiEcliptic(1)/norm(vPCiEcliptic);
                                         cos(iPC)*vPCiEcliptic(2)/norm(vPCiEcliptic); 
                                         sin(iPC)];
%     vPCfEcliptic = vPCiEcliptic + [0;0;norm(vPCiEcliptic)*sin(iPC)];
    vPCf = transferOrbit.toPrimaryBasis(vPCfEcliptic);
        % velocity after plane change burn
    transferOrbitPC = Orbit.fromStateVector(rPC,vPCf,tPC,transferOrbit.primary);
    planeChangeDV = vPCf - vPCi;
    planeChangeDT = tPC - tStart;
else
    transferOrbitPC = transferOrbit;
    planeChangeDV = 0;
    planeChangeDT = 0;
end

end
