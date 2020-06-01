classdef Orbit
    properties (Access = public)
        a
        e
        inc
        Omega
        w
        Mo
        primary
    end
    methods (Access = public, Static)
        function orb = Orbit(a,e,inc,Omega,w,Mo,primary)
            if nargin > 0
                orb.a = a;          % semimajor axis (m)
                orb.e = e;          % eccentricity
                orb.inc = inc;      % inclination (rad)
                orb.Omega = Omega;  % longitude of ascending node (rad)
                orb.w = w;          % argument of periapsis (rad)
                orb.Mo = Mo;        % mean anomaly at epoch (rad)
                orb.primary = primary;
            end
        end
        function orb = fromStateVector(r,drdt,t,primary)
            orb = Orbit();
            if norm(r)==0
                error('Invalid position')
                    % mostly here in case this is used for the sun
            end
            orb.primary = primary;
            % constructs orbit from state vector, time and grav parameter
            orb.a = 1/(2/norm(r) - ((norm(drdt))^2)/primary.mu);
                % semimajor axis (m)
            h = cross(r,drdt);
                % orbital momentum vector (m^2/s)
            eVec = cross(drdt,h)/orb.primary.mu - r/norm(r);
                % eccentricity vector
            orb.e = norm(eVec);
                % eccentricity
            n = cross([0;0;1],h);
            if norm(n)==0
               n = eVec; 
            end
                % vector pointing toward ascending node
            if dot(r,drdt)>=0
                nu = real(acos(dot(eVec,r)/(norm(eVec)*norm(r))));
            else
                nu = real(2*pi - acos(dot(eVec,r)/(norm(eVec)*norm(r))));
            end
            if isnan(nu)
               nu = 0; 
            end
                % true anomaly (rad)
            orb.inc = acos(h(3)/norm(h));
                % inclination (rad)
            if n(2)>=0
                orb.Omega = acos(n(1)/norm(n));
            else
                orb.Omega = 2*pi - acos(n(1)/norm(n));
            end
            if isnan(orb.Omega)
                orb.Omega = 0;
            end
                % longitude of ascending node (rad)
            if eVec(3)>=0
                orb.w = real(acos(dot(n,eVec)/(norm(n)*norm(eVec))));
            else
                orb.w = 2*pi - real(acos(dot(n,eVec)/(norm(n)*norm(eVec))));
            end
            if isnan(orb.w)
                orb.w = 0;
            end
                % argument of periapsis (rad)
            if orb.e == 1   % parabolic case
                error('Parabolic case (e=1) not implemented')
            elseif orb.e<1  % elliptical case
                E = 2*atan(tan(nu/2)/sqrt((1+orb.e)/(1-orb.e)));
                    % eccentric anomaly (rad)
                M = E - orb.e*sin(E);
                    % mean anomaly (rad)
                orb.Mo = M - (sqrt(primary.mu/orb.a^3))*t;
                if orb.Mo<0 || orb.Mo>=2*pi
                    twopiRatio = orb.Mo/(2*pi);
                    orb.Mo = orb.Mo-2*pi*floor(twopiRatio);
                end
                    % mean anomaly at epoch (rad)
            elseif orb.e>1  % hyperbolic case
                F = acosh((orb.e+cos(nu))/(1+orb.e*cos(nu)));
                if nu <0
                    F=-F;
                end
                    % hyperbolic anomaly (rad)
                M = orb.e*sinh(F) - F;
                    % mean anomaly (rad)
                orb.Mo = M - (sqrt(primary.mu/-orb.a^3))*t;
            end
        end
    end
    methods (Access = public)
        function [r, drdt] = getStateVector(orb,t)
            if isempty(orb.primary)
                r = [0;0;0];
                drdt = [0;0;0];
                return
                    % mostly here for use with the sun
            end
            nu = orb.getTrueAnomaly(t);
            rMag = orb.a*(1-orb.e^2)/(1+orb.e*cos(nu));
                % distance to central body (m)
            o = rMag*[cos(nu);sin(nu);0];
                % position vector in orbital frame (periapsis on +x axis)
            phi = atan(orb.e*sin(nu)/(1+orb.e*cos(nu)));
                % flight path angle (rad) 
            drdtMag = sqrt(orb.mu*(2/rMag-1/orb.a));
                % speed (m/s)            
            dodt = drdtMag*[cos(nu+pi/2-phi);sin(nu+pi/2-phi);0];
                % velocity vector in orbital frame
            R1 =   [cos(-orb.Omega) -sin(-orb.Omega)	0;
                    sin(-orb.Omega) cos(-orb.Omega)     0;
                    0               0                   1];
                % rotate (z axis) to match longitude of ascending node
            R2 =   [1               0                   0;
                    0               cos(-orb.inc)       -sin(-orb.inc);
                    0               sin(-orb.inc)       cos(-orb.inc)];
                % rotate (x axis) to match inclination 
            R3 =   [cos(-orb.w)     -sin(-orb.w)        0;
                    sin(-orb.w)     cos(-orb.w)         0;
                    0               0                   1];
                % rotate (x axis) to match argument of periapsis
            r = (R3*R2*R1)'*o;
            drdt = (R3*R2*R1)'*dodt;
                % position and velocity vectors in inertial frame
        end
        function nu = getTrueAnomaly(orb, t)
            if orb.e==1     % parabolic case
                error('Parabolic case (e=1) not implemented')
            elseif orb.e<1  % elliptical case
                % Calculate Cartesian state vector at time t
                M = orb.Mo + t*(orb.mu/orb.a^3)^0.5; 
                    % mean anomaly (rad)
                if M<0 || M>=2*pi
                    twopiRatio = M/(2*pi);
                    M = M-2*pi*floor(twopiRatio);
                end
                E = solveKepler(M,orb.e); 
                    % eccentric anomaly (rad)
                nu = 2*atan2(sqrt(1+orb.e)*sin(E/2),sqrt(1-orb.e)*cos(E/2)); 
                    % true anomaly (rad)
            elseif orb.e>1  % hyperbolic case
                % Calculate Cartesian state vector at time t
                M = orb.Mo + t*(orb.mu/-orb.a^3)^0.5; 
                    % mean anomaly (rad)
                F = solveKepler(M,orb.e); 
                    % hyperbolic anomaly (rad)
                o = -[orb.a*(orb.e-cosh(F));orb.a*sqrt(orb.e^2-1)*sinh(F);0];
                    % position vector in orbital frame
                nu = atan2(o(2),o(1));
            end
        end
        function t = getTime(orb, nu)
            if orb.e == 1   % parabolic case
                error('Parabolic case (e=1) not implemented')
            elseif orb.e<1  % elliptical case
                E = 2*atan(tan(nu/2)/sqrt((1+orb.e)/(1-orb.e)));
                    % eccentric anomaly (rad)
                M = E - orb.e*sin(E);
                    % mean anomaly (rad)
                if M<0 || M>=2*pi
                    twopiRatio = M/(2*pi);
                    M = M-2*pi*floor(twopiRatio);
                end
                t = (M-orb.Mo)/sqrt(orb.mu/orb.a^3);
                while t<0
                   t = t+orb.period; 
                end
            elseif orb.e>1  % hyperbolic case
                F = acosh((orb.e+cos(nu))/(1+orb.e*cos(nu)));
                if nu <0
                    F=-F;
                end
                    % hyperbolic anomaly (rad)
                M = orb.e*sinh(F) - F;
                    % mean anomaly (rad)
                t = (M-orb.Mo)/sqrt(orb.mu/-orb.a^3);
            end
        end
        function gravParam = mu(orb)
            gravParam = orb.primary.mu;
        end
        function T = period(orb)
            T = 2*pi*sqrt(orb.a^3/orb.mu);
        end
        function theta = getAngleInEcliptic(orb,t,r2)
            r1 = orb.getStateVector(t);
            r1Ecliptic = orb.toOrbitBasis(r1);
            r2Ecliptic = orb.toOrbitBasis(r2);
            thetaR1 = atan2(r1Ecliptic(2),r1Ecliptic(1));
            thetaR2 = atan2(r2Ecliptic(2),r2Ecliptic(1));
            theta = thetaR2-thetaR1;
            if theta<-pi
                theta = theta + 2*pi;
            elseif theta>pi
                theta = theta - 2*pi;
            end
        end
        function [X,Y,Z] = getOrbitBasis(orb)
            if isempty(orb.primary)
               X = [1;0;0];
               Y = [0;1;0];
               Z = [0;0;1];
               return
                    % returns the celestial reference basis for the sun
            end
            [ro,vo] = orb.getStateVector(0);
                % position and velocity vectors at epoch
            Z = cross(ro,vo);
            Z = Z/norm(Z);
                % normal direction to orbital plane
            X = [1;0;0];        
                % Assumed celestial longitude
            X = X - dot(X,Z)*Z;
            X = X/norm(X);
                % projection of celestial longitude into orbital plane
            Y = cross(Z,X);
            Y = Y/norm(Y); 
                % 2nd basis vector (velocity direction if circular orbit)
        end
        function rOrb = toOrbitBasis(orb,rPrim)
            [X,Y] = orb.getOrbitBasis;
            rOrb = Transform3D(rPrim,X,Y,false);
        end
        function rPrim = toPrimaryBasis(orb,rOrb)
            [X,Y] = orb.getOrbitBasis;
            rPrim = Transform3D(rOrb,X,Y,true);
        end
        function [] = plotTrajectory(orb,t1,t2,color,startMarker,endMarker,lineType)
            if nargin <4
                color = 'w';
            end
            if nargin <5
               startMarker = 'none';
               endMarker = 'none';
            end
            if nargin <7
                lineType = '-';
            end
            t = linspace(t1,t2,1000);
            r = zeros(3,length(t));
            for ii = 1:length(t)
                if isempty(orb.primary.primary)
                    r(:,ii) = orb.getStateVector(t(ii));
                else
                    r(:,ii) = orb.primary.toPrimaryBasis(orb.getStateVector(t(ii)));
                end                    
            end
            hold on
            plot3(r(1,:),r(2,:),r(3,:),'LineStyle',lineType,'Color',color)
            plot3(r(1,1),r(2,1),r(3,1),'Marker',startMarker,'Color',color)
            plot3(r(1,end),r(2,end),r(3,end),'Marker',endMarker,'Color',color)
            hold off
            set(gca,'Color',[0.2,0.2,0.2])
            set(gca,'xticklabel',[])
            set(gca,'yticklabel',[])
            set(gca,'zticklabel',[])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'ztick',[])
            axis equal
        end
    end
end
