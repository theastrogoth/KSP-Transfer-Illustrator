classdef Body < handle
    properties (Access = public)
        name     % body's name
        eqr      % equatorial radius (m)
        mu       % gravitational parameter (m^3/s^2)
        SOI      % sphere of influence for patched conic (m)
        orbit    % Keplerian orbital elements
        color    % color chosen to represent body in plots
    end
    methods (Access = public, Static)
        function bd = Body(name,eqr,mu,SOI,orbit,color)
            if nargin > 0
                bd.name = name;
                bd.eqr = eqr;
                bd.mu = mu;
                bd.SOI = SOI;
                bd.orbit = orbit;
                bd.color = color;
            end
        end
    end
    methods (Access = public)
        function [] = plotTrajectory(bd,t1,t2,startMarker,endMarker,lineType)
            if nargin <4
               startMarker = 'none';
               endMarker = 'none';
            end
            if nargin <6
                lineType = '-';
            end
            bd.orbit.plotTrajectory(t1,t2,bd.color,startMarker,endMarker,lineType)
        end
        function prim = primary(bd)
            prim = bd.orbit.primary;
        end
        function [X,Y,Z] = getBodyBasis(bd)
           [X,Y,Z] = bd.orbit.getOrbitBasis; 
        end
        function rBod = toBodyBasis(bd,rPrim)
            rBod = bd.orbit.toOrbitBasis(rPrim);
        end
        function rPrim = toPrimaryBasis(bd,rBod)
            rPrim = bd.orbit.toPrimaryBasis(rBod);
        end
        function [r,drdt] = getStateVector(bd,t)
           [r,drdt] = bd.orbit.getStateVector(t); 
        end
        function T = period(bd)
           T = bd.orbit.period; 
        end
    end
end