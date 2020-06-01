classdef Transfer < handle
   properties (Access = public)
      startOrbit            % starting parking orbit(assume equatorial, circular)
      endOrbit              % ending parking orbit (assume circular, not necessarily equatorial)
      startTime             % time of ejection from starting body's SOI (s)
      flightTime            % time interval between ejection and insertion (s)
      planeChange           % false if ballistic, true if plane change
      transferOrbit         % orbit between ejection and insertion
      transferOrbitPC       % transfer orbit after plane change burn
      ejectionTrajectory    % hyperbolic trajectory for leaving starting body
      insertionTrajectory   % hyperbolic trajectory for arriving at target
      ejectionDV            % delta v for ejection burn (m/s)
      insertionDV           % delta v for insertion burn (m/s)
      planeChangeDV         % delta v for plane change burn (m/s)
      ejectionDT            % delta t for ejection trajectory (s)
      insertionDT           % delta t for insertion trajectory (s)
      planeChangeDT         % time interval before plane change burn (s)
   end
   methods (Access = public, Static)
       function trs = Transfer(startOrbit,endOrbit,startTime,flightTime,planeChange)
           if nargin > 0
              if ~isequal(startOrbit.primary.primary,endOrbit.primary.primary)
                 error('start and end bodies do not orbit the same primary body') 
              end
              if nargin<5
                 trs.planeChange = false; 
              else
                  trs.planeChange = planeChange;
              end
              trs.startOrbit = startOrbit;
              trs.endOrbit = endOrbit;
              trs.startTime = startTime;
              trs.flightTime = flightTime;
              trs.getTransferParams; 
              trs.getEjectionParams;
              trs.getInsertionParams;
           end
       end
   end
   methods (Access = public)
       function sp = startPrimary(trs)
           sp = trs.startOrbit.primary;
       end
       function ep = endPrimary(trs)
           ep = trs.endOrbit.primary;           
       end
       function getTransferParams(trs)
           [trs.transferOrbit, trs.transferOrbitPC, trs.planeChangeDV, trs.planeChangeDT] = ...
               solveLambert(trs.startPrimary.orbit,trs.endPrimary.orbit,trs.startTime,trs.flightTime,trs.planeChange);
       end
       function getEjectionParams(trs)
           % assumes a circular parking orbit with no inclination
           % get parameters from parking orbit
           ro = trs.startOrbit.a;
                % parking orbit distance from gravity center
           mu = trs.startOrbit.mu;
           SOI = trs.startPrimary.SOI;
           [~,vp] = trs.startPrimary.getStateVector(trs.startTime);
           [~,vs] = trs.transferOrbit.getStateVector(trs.startTime);
           vsp = vs - vp;
           vsp = trs.startPrimary.toBodyBasis(vsp);
                % excess heliocentric velocity (rotated to body's orbital
                % reference frame)
           vo = sqrt(norm(vsp)^2 + 2*(mu/ro-mu/SOI));
                % velocity after ejection burn
           e = sqrt(1+2*(vo^2/2-mu/ro)*ro^2*vo^2/mu^2);
           a = 1/(2/ro-vo^2/mu);
                % hyperbolic orbit parameters  
           % describe positions at periapsis and SOI in the hyperbola's plane
           thetaSOI = acos(1/e*(a*(1-e^2)/SOI-1));
           phiSOI = atan(e*sin(thetaSOI)/(1+e*cos(thetaSOI)));
           thetaV = thetaSOI + pi/2 - phiSOI;
           vSOI = sqrt(mu*(2/SOI-1/a))*[cos(thetaV);sin(thetaV);0];
                % angles relative to periapsis
           roVector = ro*[1;0;0];
           voVector = vo*[0;1;0];
           no = [0;0;1];               % normal vector at periapsis BEFORE burn
                % vectors for periapsis (orbital reference frame)           
           % rotate coordinates to acheive appropriate ejection direction
           phi1 = atan2(vsp(3),sqrt(abs(norm(vSOI)^2-vSOI(1)^2-vsp(3)^2)));
           R1 =    [1  0          0         ;
                    0  cos(phi1) -sin(phi1) ;
                    0  sin(phi1)  cos(phi1)];
           R1vSOI = R1*vSOI;
                % rotate about x axis to match ejection z component
           theta2 = atan2(vsp(2),vsp(1)) - atan2(R1vSOI(2),R1vSOI(1));
           R2 =    [cos(theta2)  -sin(theta2) 0 ;
                    sin(theta2)   cos(theta2) 0 ;
                    0             0           1];
                % rotate about z axis to match ejection direction
           % apply transformations
           roVector = R2*R1*roVector;
           voVector = R2*R1*voVector;          
           % get trajectory parameters
           ET = Orbit.fromStateVector(roVector,voVector,0,trs.startPrimary);
           FSOI = sign(thetaSOI)*acosh((cos(thetaSOI)+e)/(1+e*cos(thetaSOI)));
           deltaM = e*sinh(FSOI) - FSOI - ET.Mo;
           trs.ejectionDT = deltaM * sqrt(-ET.a^3/ET.mu);
           trs.ejectionTrajectory = Orbit.fromStateVector(roVector,voVector,trs.getEjectionBurnTime,trs.startPrimary);
           % get prograde and normal components to burn
           vPark = cross(no,roVector);
           vPark = sqrt(mu/ro)*vPark/norm(vPark);
           trs.ejectionDV = voVector - vPark;
           % get parking orbit (should be equatorial and circular)
           trs.startOrbit = Orbit.fromStateVector(roVector,vPark,trs.getEjectionBurnTime,trs.startPrimary);
       end
       function getInsertionParams(trs)
           % assumes a circular parking orbit with the same inclination as
           % hyperbolic approach vector
           % get parameters from parking orbit
           ro = trs.endOrbit.a;
                % parking orbit distance from gravity center
           mu = trs.endOrbit.mu;
           SOI = trs.endPrimary.SOI;
           [~,vp] = trs.endPrimary.getStateVector(trs.startTime+trs.flightTime);
           if trs.planeChange
               [~,vs] = trs.transferOrbitPC.getStateVector(trs.startTime+trs.flightTime);
           else
               [~,vs] = trs.transferOrbit.getStateVector(trs.startTime+trs.flightTime);
           end
           vps = vp - vs;
           vps = trs.endPrimary.toBodyBasis(vps);
                % excess heliocentric velocity (rotated to body's
                % reference)
           vo = sqrt(norm(vps)^2 + 2*(mu/ro-mu/SOI));
                % velocity before insertion burn
           e = sqrt(1+2*(vo^2/2-mu/ro)*ro^2*vo^2/mu^2);
           a = 1/(2/ro-vo^2/mu);
                % hyperbolic orbit parameters  
           % describe positions at periapsis and SOI in the hyperbola's plane
           thetaSOI = -acos(1/e*(a*(1-e^2)/SOI-1));
           phiSOI = atan(e*sin(thetaSOI)/(1+e*cos(thetaSOI)));
           thetaV = thetaSOI + pi/2 - phiSOI;
           vSOI = sqrt(mu*(2/SOI-1/a))*[cos(thetaV);sin(thetaV);0];
                % angles relative to periapsis
           roVector = ro*[1;0;0];
           voVector = vo*[0;1;0];
           no = [0;0;1];               % normal vector at periapsis BEFORE burn
                % vectors for periapsis (orbital reference frame)           
           % rotate coordinates to acheive appropriate insertion direction
           phi1 = atan2(vps(3),sqrt(abs(norm(vSOI)^2-vSOI(1)^2-vps(3)^2)));
           R1 =    [1  0          0         ;
                    0  cos(phi1) -sin(phi1) ;
                    0  sin(phi1)  cos(phi1)];
           R1vSOI = R1*vSOI;
                % rotate about x axis to match insertion z component
           theta2 = atan2(vps(2),vps(1)) - atan2(R1vSOI(2),R1vSOI(1));
           R2 =    [cos(theta2)  -sin(theta2) 0 ;
                    sin(theta2)   cos(theta2) 0 ;
                    0             0           1];
                % rotate about z axis to match insertion direction
           % apply transformations
           roVector = R2*R1*roVector;
           voVector = R2*R1*voVector;  
           no = R2*R1*no;
           % get trajectory parameters
           IT = Orbit.fromStateVector(roVector,voVector,0,trs.endPrimary);
           FSOI = sign(thetaSOI)*acosh((cos(thetaSOI)+e)/(1+e*cos(thetaSOI)));
           deltaM = e*sinh(FSOI) - FSOI - IT.Mo;
           trs.insertionDT = deltaM * sqrt(-IT.a^3/IT.mu);
           trs.insertionTrajectory = Orbit.fromStateVector(roVector,voVector,trs.getInsertionBurnTime,trs.endPrimary);
           % get prograde and normal components to burn
           vPark = cross(no,roVector);
           vPark = sqrt(mu/ro)*vPark/norm(vPark);
           trs.insertionDV = voVector - vPark;
           % get parking orbit (should be equatorial and circular)
           trs.endOrbit = Orbit.fromStateVector(roVector,vPark,trs.getInsertionBurnTime,trs.endPrimary);
       end
       function plotTransferOrbit(trs,color)
           % plot start body, end body, and transfer orbits
           trs.endPrimary.plotTrajectory(trs.startTime,trs.startTime+trs.flightTime,'o','none')
           trs.startPrimary.plotTrajectory(trs.startTime,trs.startTime+trs.flightTime)
           if trs.planeChange
               trs.transferOrbit.plotTrajectory(trs.startTime,trs.startTime+trs.planeChangeDT,color,'o','^')
               trs.transferOrbitPC.plotTrajectory(trs.startTime+trs.planeChangeDT,trs.startTime+trs.flightTime,color,'^','x','--')
           else
               trs.transferOrbit.plotTrajectory(trs.startTime,trs.startTime+trs.flightTime,color,'o','x')
           end
           % set axis limits
           rMax = max([trs.startPrimary.orbit.a*(1+trs.startPrimary.orbit.e),...
               trs.endPrimary.orbit.a*(1+trs.endPrimary.orbit.e),...
               trs.transferOrbit.a*(1+trs.transferOrbit.e)]);     
           xlim([-1.25*rMax,1.25*rMax])
           ylim([-1.25*rMax,1.25*rMax])
           zlim([-1.25*rMax,1.25*rMax])
           % show angle between bodies at ejection in plane of start body
           hold on
           r1 = trs.startPrimary.getStateVector(trs.startTime);
           [X,Y] = trs.startPrimary.getBodyBasis;
           r2 = trs.endPrimary.getStateVector(trs.startTime);
           r1Ecliptic = Transform3D(r1,X,Y);
           r2Ecliptic = Transform3D(r2,X,Y);
           r2Plane = Transform3D([r2Ecliptic(1);r2Ecliptic(2);0],X,Y,true);
           thetaR1 = atan2(r1Ecliptic(2),r1Ecliptic(1));
           thetaR2 = atan2(r2Ecliptic(2),r2Ecliptic(1));
           phaseAngle = thetaR2-thetaR1;
           if phaseAngle<-pi
              phaseAngle = phaseAngle + 2*pi;  
              thetaR2 = thetaR2+2*pi;
           elseif phaseAngle>pi
              phaseAngle = phaseAngle - 2*pi; 
              thetaR1 = thetaR1+2*pi;
           end
           thetaArc = linspace(thetaR1,thetaR2,101);
           rArc = [cos(thetaArc); sin(thetaArc); 0*thetaArc]*rMax*1.1;
           rArc = Transform3D(rArc,X,Y,true);
           plot3(rArc(1,:),rArc(2,:),rArc(3,:),':','Color',color);
           plot3([r1(1,1),rArc(1,1)],[r1(2,1),rArc(2,1)],[r1(3,1),rArc(3,1)],':','Color',color)
           plot3([r2Plane(1),rArc(1,end)],[r2Plane(2),rArc(2,end)],[r2Plane(3),rArc(3,end)],':','Color',color)
           text(rArc(1,ceil(end/2))*1.1,rArc(2,ceil(end/2))*1.1,rArc(3,ceil(end/2))*1.1,[num2str(phaseAngle/pi*180,'%4.2f'),char(176)],'Color',color)
           % show primary body
           [xSphere,ySphere,zSphere] = sphere(20);
           xSphere = trs.transferOrbit.primary.eqr*xSphere;
           ySphere = trs.transferOrbit.primary.eqr*ySphere;
           zSphere = trs.transferOrbit.primary.eqr*zSphere;
           surface(xSphere,ySphere,zSphere,'FaceColor', 'none',...
               'EdgeColor',trs.transferOrbit.primary.color,'LineWidth',0.5)
           hold off
       end
       function plotEjection(trs,color)
           % show body's direction from primary and prograde direction
           hold on
           [rp,vp] = trs.startPrimary.getStateVector(trs.getEjectionBurnTime);
           rpDir = rp/norm(rp) * 2.5*trs.startOrbit.a;
           vpDir = vp/norm(vp) * 2.5*trs.startOrbit.a;
           quiver3(0,0,0,-rpDir(1),-rpDir(2),-rpDir(3),'Color',trs.startPrimary.color,'LineWidth',1.5)
           quiver3(0,0,0,vpDir(1),vpDir(2),vpDir(3),'Color',trs.startPrimary.color,'LineWidth',1.5)
           % show body
           [X,Y] = trs.startPrimary.getBodyBasis;
           [xSphere,ySphere,zSphere] = sphere(20);
           sz = size(xSphere);
           xSphere = trs.startPrimary.eqr*xSphere(:)';
           ySphere = trs.startPrimary.eqr*ySphere(:)';
           zSphere = trs.startPrimary.eqr*zSphere(:)';
           rSphere = Transform3D([xSphere;ySphere;zSphere],X,Y,true);
           surface(reshape(rSphere(1,:),sz),reshape(rSphere(2,:),sz),reshape(rSphere(3,:),sz),'FaceColor', 'none',...
               'EdgeColor',trs.startPrimary.color,'LineWidth',0.5)
           % show angle of burn location from planet's prograde direction
           ro = trs.ejectionTrajectory.getStateVector(trs.getEjectionBurnTime);
           roDirEcliptic = Transform3D(ro,X,Y);
           vpDirEcliptic = Transform3D(vp,X,Y);
           thetaRO = atan2(roDirEcliptic(2),roDirEcliptic(1));
           thetaVP = atan2(vpDirEcliptic(2),vpDirEcliptic(1));
           burnAngle = thetaRO - thetaVP; 
           if burnAngle<-pi
               burnAngle = burnAngle + 2*pi;
               thetaRO = thetaRO+2*pi;
           elseif burnAngle>pi
               burnAngle = burnAngle - 2*pi;
               thetaVP = thetaVP+2*pi;
           end
           thetaArc = linspace(thetaVP,thetaRO,101);
           rArc = [cos(thetaArc); sin(thetaArc); 0*thetaArc]*trs.startOrbit.a*1.5;
           rArc = Transform3D(rArc,X,Y,true);
           plot3(rArc(1,:),rArc(2,:),rArc(3,:),':','Color',color);
           plot3([0,rArc(1,1)],[0,rArc(2,1)],[0,rArc(3,1)],':','Color',color)
           plot3([0,rArc(1,end)],[0,rArc(2,end)],[0,rArc(3,end)],':','Color',color)
           text(rArc(1,ceil(end/2))*1.1,rArc(2,ceil(end/2))*1.1,rArc(3,ceil(end/2))*1.1,[num2str(burnAngle/pi*180,'%4.2f'),char(176)],'Color',color)
           % show direction of burn vector
           roPrimary = trs.startPrimary.toPrimaryBasis(ro);
           burnDir = trs.startPrimary.toPrimaryBasis(trs.ejectionDV/norm(trs.ejectionDV)) * 1.5*trs.startOrbit.a;
           quiver3(roPrimary(1),roPrimary(2),roPrimary(3),burnDir(1),burnDir(2),burnDir(3),'LineWidth',1.5,'Color',color)
           hold off
           % plot parking and ejection trajectories
           trs.ejectionTrajectory.plotTrajectory(trs.getEjectionBurnTime,trs.getEjectionBurnTime+trs.ejectionDT,color,'o','x')
           trs.startOrbit.plotTrajectory(trs.getEjectionBurnTime-trs.startOrbit.period,trs.getEjectionBurnTime,color,'none','none','--')  
           % set axis limits
           xlim([-5*trs.startOrbit.a,5*trs.startOrbit.a])
           ylim([-5*trs.startOrbit.a,5*trs.startOrbit.a])
           zlim([-5*trs.startOrbit.a,5*trs.startOrbit.a])
       end
       function plotInsertion(trs,color)
           % show body's direction from primary and direction of velocity
           hold on
           [rp,vp] = trs.endPrimary.getStateVector(trs.getInsertionBurnTime);
           rpDir = rp/norm(rp) * 2.5*trs.endOrbit.a;
           vpDir = vp/norm(vp) * 2.5*trs.endOrbit.a;
           quiver3(0,0,0,-rpDir(1),-rpDir(2),-rpDir(3),'Color',trs.endPrimary.color,'LineWidth',1.5)
           quiver3(0,0,0,vpDir(1),vpDir(2),vpDir(3),'Color',trs.endPrimary.color,'LineWidth',1.5)
           % show body
           [xSphere,ySphere,zSphere] = sphere(20);
           sz = size(xSphere);
           xSphere = trs.endPrimary.eqr*xSphere(:)';
           ySphere = trs.endPrimary.eqr*ySphere(:)';
           zSphere = trs.endPrimary.eqr*zSphere(:)';
           rSphere = Transform3D([xSphere;ySphere;zSphere],rpDir,vpDir);
           surface(reshape(rSphere(1,:),sz),reshape(rSphere(2,:),sz),reshape(rSphere(3,:),sz),'FaceColor', 'none',...
               'EdgeColor',trs.endPrimary.color,'LineWidth',0.5)
           % show direction of burn vector
           ro = trs.insertionTrajectory.getStateVector(trs.getInsertionBurnTime);
           burnDir = trs.insertionDV/norm(trs.insertionDV) * 1.5*trs.endOrbit.a;
           quiver3(ro(1),ro(2),ro(3),burnDir(1),burnDir(2),burnDir(3),'LineWidth',1.5,'Color',color)
           hold off
           % plot parking and ejection trajectories
           trs.insertionTrajectory.plotTrajectory(trs.getInsertionBurnTime-trs.insertionDT,trs.getInsertionBurnTime,color,'o','x')
           trs.endOrbit.plotTrajectory(trs.getInsertionBurnTime,trs.getInsertionBurnTime+trs.endOrbit.period,color,'none','none','--')
           % set axis limits
           xlim([-5*trs.endOrbit.a,5*trs.endOrbit.a])
           ylim([-5*trs.endOrbit.a,5*trs.endOrbit.a])
           zlim([-5*trs.endOrbit.a,5*trs.endOrbit.a])
       end
       function ejectBurnTime = getEjectionBurnTime(trs)
           ejectBurnTime = trs.startTime - trs.ejectionDT;
       end
       function insertBurnTime = getInsertionBurnTime(trs)
           insertBurnTime = trs.startTime + trs.flightTime + trs.insertionDT;
       end
   end
end