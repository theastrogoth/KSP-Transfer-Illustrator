# -*- coding: utf-8 -*-

import math
import numpy as np
from numpy.linalg import norm
from orbit import Orbit
from body import Body

class Transfer:
    """Orbital transfer from a starting orbit to an ending orbit.
    
    Attributes:
        startOrbit (Orbit): orbit prior to departure burns
        endOrbit (Orbit): orbit following arrival burns
        startTime (float): time at the beginning of transfer trajectory (s)
        flightTime (float): duration of transfer trajectory (s)
        planeChange (bool): if true, a mid-course plane change is done
        ignoreInsertion (bool): if true, arrival burn is ignored.
        fixedEndTime (float): if given as input, the target position 
            will be specified by the position of the end orbit at this 
            time (s). Used for setting up successive transfers.
        transferOrbit (Orbit): orbital trajectory between start and end.
            If there is a plane change maneuver, this is the portion of the
            trajectory prior to the maneuver.
        transferORbitPC (Orbit): If there is a plane change maneuver, this
            is the portion of the transfer trajectory after the maneuver.
        ejectionTrajectory (Orbit): If ejection from a starting body occurs,
            this is the trajectory between parking orbit and escape.
        ejectionDV (array): 3D departure burn vector (m/s)
        insertionDV (float): magnitude of arrival burn (m/s)
        planeChangeDV (array): 3D plane change burn vector (m/s)
        ejectionDT (float): If ejection from a starting body occurs, this is
            the time interval between ejection burn and escape.
        insertionDT (float): If insertion at an ending body occurs, this is
            the time interval between encounter and insertion burn.
        planeChangeDT (float): If a plane change maneuver occurs, this is
            the time interval between the start of the transfer orbit and
            the maneuver.
            
    """
    
    def __init__(self, startOrbit, endOrbit, startTime, flightTime, 
                 planeChange = False, ignoreInsertion = False, 
                 fixedEndTime = None):
        
        # Assign input attributes
        self.startOrbit = startOrbit
        self.endOrbit = endOrbit
        self.startTime = startTime
        self.flightTime = flightTime
        self.planeChange = planeChange
        self.ignoreInsertion = ignoreInsertion
        self.fixedEndTime = fixedEndTime
        
        # These attributes get defined here but are filled in with methods
        self.transferOrbit = None
        self.transferOrbitPC = None
        self.ejectionTrajectory = None
        self.ejectionDV = 0
        self.insertionDV = 0
        self.planeChangeDV = 0
        self.ejectionDT = 0
        self.insertionDT = 0
        self.planeChangeDT = 0
        self.phaseAngle = 0
        self.ejectionBurnAngle = None
        
        # Calculate transfer, ejection, and insertion parameters
        self.get_transfer_details()
    
    
    @staticmethod
    def solve_lambert(startOrbit, endOrbit, startTime, flightTime,
                      planeChange = False, fixedEndTime = None,
                      tol = 1E-6, maxIt = 200):
        """Solves the Lambert problem to obtain a trajectory to the target.
        
        Args:
            startOrbit (Orbit): orbit prior to departure
            endOrbit (Orbit): orbit following arrival
            startTime (float): time at the beginning of transfer (s)
            flightTime (float): duration of transfer trajectory (s)
            planeChange (bool): if true, a mid-course plane change occurs
            fixedEndTime (float): if provided, fixes target location to end
                orbit's position at this time (s)
            tol (float): the maximum tolerance for iteration termination
            maxIt (int): the maximum number of iterations before breaking
        
        Returns:
            transferOrbit (Orbit): trajectory before plane change
            transferOrbitPC (Orbit): trajectory after plane change
            planeChangeDV (array): plane change burn vector (m/s)
            planeChangeDT (float): time between start and plane change (s)
        """
        
        # Set gravitational parameter for transfer orbit
        mu  = startOrbit.prim.mu
        
        # Get start position
        rStart = startOrbit.get_state_vector(startTime)[0]
        
        # Adjust end position based on inputs
        if fixedEndTime is None:
            rEnd = endOrbit.get_state_vector(startTime + flightTime)[0]
        else:
            # Choose the end position at a particular point in the orbit
            rEnd = endOrbit.get_state_vector(fixedEndTime)[0]
        if planeChange:
            # Rotate the target orbit to be coplanar with the starting one
            rEnd = startOrbit.from_primary_to_orbit_bases(rEnd)
            rEnd = np.array([rEnd[0],rEnd[1],0]) /                          \
                norm(np.array([rEnd[0],rEnd[1]])) * norm(rEnd);
            rEnd = startOrbit.from_orbit_to_primary_bases(rEnd)
        
        # Store magnitudes of position vectors for later use
        rStartMag = norm(rStart)
        rEndMag = norm(rEnd)
        
        # Get true anomaly change and angles in the ecliptic (x-y plane)
        dNu = math.atan2(norm(np.cross(rStart,rEnd)),np.dot(rStart,rEnd))
        thetaStart = math.atan2(rStart[1], rStart[0])
        thetaEnd = math.atan2(rEnd[1], rEnd[0])
        dTheta = thetaEnd - thetaStart
        if dTheta < 0 or dTheta > 2*math.pi:
            dTheta = dTheta - math.floor(dTheta/(2*math.pi))*2*math.pi
        if dTheta > math.pi:
            dNu = 2*math.pi - dNu
        
        # Set constants for p iteration
        k = rStartMag * rEndMag * (1-math.cos(dNu))
        L = rStartMag + rEndMag
        m = rStartMag * rEndMag * (1+math.cos(dNu))
        
        # Set bounds for p values
        pj = k / (L+math.sqrt(2*m))
        pjj = k / (L-math.sqrt(2*m))
        if dNu > math.pi:
            pMin = 0
            pMax = pjj
        else:
            pMin = pj
            pMax = math.inf
        
        # Initialize values prior to iteration
        it = 0
        err = tol+1
        p = (pj+pjj)/2
        pNext = p
        
        # Use Newton-p-iteration to minimize error for time of flight
        while err > tol:
            it = it+1
            if it > maxIt:
                break
            p = pNext
            a = m*k*p / ((2*m-L**2)*(p**2) + 2*k*L*p - k**2)
            f = 1 - rEndMag/p * (1 - math.cos(dNu))
            g = rStartMag * rEndMag * math.sin(dNu) / math.sqrt(mu*p)
            df = math.sqrt(mu/p)*math.tan(dNu/2)*((1-math.cos(dNu))/p -     \
                                                  1/rStartMag -             \
                                                  1/rEndMag);
            # Elliptical case
            if a > 0:
                sindE = -rStartMag * rEndMag * df/math.sqrt(mu*a)
                cosdE = 1 - rStartMag/a * (1-f)
                # Change in elliptical anomaly
                dE = math.atan2(sindE,cosdE)
                while dE < 0:
                    dE = dE + 2*math.pi
            
                # Time of flight and slope with respect to p
                t = g + math.sqrt((a**3)/mu) * (dE - sindE)
                dtdp = -g/(2*p) -                                           \
                        1.5*a*(t-g)*(k**2 + (2*m-L**2)*p**2) / (m*k*p**2) + \
                        math.sqrt(a**3/mu) * (2*k*sindE) / (p*(k-L*p));
            
            # Hyperbolic case
            else:
                # Change in hyperbolic anomaly
                dF = math.acosh(1 - rStartMag/a * (1-f))
                
                # Time of flight and slope with respect to p
                t = g + math.sqrt((-a)**3/mu) * (math.sinh(dF)  - dF)
                dtdp = -g/(2*p) -                                           \
                        1.5*a*(t-g)*(k**2 + (2*m-L**2)*p**2) / (m*k*p**2) - \
                        math.sqrt((-a)**3/mu) * (2*k*math.sinh(dF)) /       \
                                                (p*(k-L*p));
            
            # Compute error and next guess for p
            err = abs(flightTime-t)/flightTime
            pNext = p + (flightTime - t) / dtdp
            
            # If the next guess is outside of allowed bounds, use bisection
            if pNext < pMin:
                pNext = (p + pMin)/2
            elif pNext > pMax:
                pNext = (p + pMax)/2
        
        # From final p-iteration parameters, calculate velocity at the start
        # of the transfer orbit, and then define the transfer orbit
        vStart = (rEnd - f * rStart)/g
        transferOrbit = Orbit.from_state_vector(rStart,vStart, startTime,
                                                startOrbit.prim)
        
        # If a mid-course plane change maneuver will be used, a second
        # transfer orbit must be determined
        if planeChange:
            # Obtain the end position in its original plane
            if fixedEndTime is None:
                rEnd = endOrbit.get_state_vector(startTime+flightTime)[0]
            else:
                rEnd = endOrbit.get_state_vector(fixedEndTime)[0]
            
            # Get angle in the orbital plane between start and end
            transferAngle =                                                 \
                transferOrbit.get_angle_in_orbital_plane(startTime,rEnd)
            
            # The optimal angle for the plane change position is 90 degrees
            # before the end position. If the transfer trajectory does not
            # cover at more than 90 degrees, the best position for plane-
            # change will be at the start
            if transferAngle < math.pi/2:
                thetaPC = 0
            else:
                thetaPC = transferAngle - math.pi/2
            
            # Get true anomaly and time at the plane-change maneuver
            nuPC = transferOrbit.get_true_anomaly(startTime) + thetaPC
            tPC = transferOrbit.get_time(nuPC, startTime)
            
            # Get the state vector immediately prior to plane change
            rPC, vPCi = transferOrbit.get_state_vector(tPC)
            vPCiPlane = transferOrbit.from_primary_to_orbit_bases(vPCi)
            
            # Calculate the inclination change needed for the maneuver
            rPCPlane = transferOrbit.from_primary_to_orbit_bases(rPC)
            rEndPlane = transferOrbit.from_primary_to_orbit_bases(rEnd)
            nf = np.cross(rPCPlane,rEndPlane)
                # normal vector to plane after burn
            incPC = math.acos(nf[2] / norm(nf))
            
            # Rotate velocity vector prior to burn to get vector after burn
            vPCfPlane = np.array([math.cos(incPC) * vPCiPlane[0],           \
                                  math.cos(incPC) * vPCiPlane[1],           \
                                  math.sin(incPC) * norm(vPCiPlane)]);
            
            # Get the velocity after plane change in the primary bases
            vPCf = transferOrbit.from_orbit_to_primary_bases(vPCfPlane)
            
            # Define second part of transfer with position and velocity
            # vectors after the plane change maneuver
            transferOrbitPC = Orbit.from_state_vector(rPC,vPCf,tPC,
                                                      transferOrbit.prim)
            
            # Get the delta V for the maneuver and time interval from start
            planeChangeDV = vPCf - vPCi
            planeChangeDT = tPC - startTime
            
        else:
            transferOrbitPC = None
            planeChangeDV = 0
            planeChangeDT = 0
        
        return transferOrbit, transferOrbitPC, planeChangeDV, planeChangeDT
    
    
    def get_transfer_details(self):
        """Get transfer and ejection orbits with burn details"""
        
        # First case: starting and ending orbits have the same primary body.
        # No change of sphere of influence takes place.
        if self.startOrbit.prim == self.endOrbit.prim:
            self.transferOrbit, self.transferOrbitPC,                       \
                self.planeChangeDV, self.planeChangeDT =                    \
                    self.solve_lambert(self.startOrbit,                     \
                                       self.endOrbit,                       \
                                       self.startTime, self.flightTime,     \
                                       self.planeChange, self.fixedEndTime);
            
            # Get departure burn delta v
            vStart = self.startOrbit.get_state_vector(self.startTime)[1]
            vTrStart = self.transferOrbit.get_state_vector(self.startTime)[1]
            self.ejectionDV = vTrStart - vStart
            
            # Get arrival burn delta v
            if not self.ignoreInsertion:
                if self.planeChange:
                    vTrEnd = self.transferOrbitPC.get_state_vector(         \
                                    self.startTime + self.flightTime)[1];
                else:
                    vTrEnd = self.transferOrbit.get_state_vector(           \
                                    self.startTime + self.flightTime)[1]
                vEnd = self.endOrbit.get_state_vector(self.startTime +      \
                                                      self.flightTime)[1];
                self.insertionDV = vEnd - vTrEnd
            
            # Get phase angle
            self.phaseAngle = self.startOrbit.get_angle_in_orbital_plane(   \
                self.get_departure_burn_time(),                             \
                self.endOrbit.get_state_vector(                             \
                    self.get_departure_burn_time())[0]);
        
        # Second case: starting in orbit around a body and transferring to a
        # parking orbit around its primary body
        elif self.startOrbit.prim in self.endOrbit.prim.satellites:
            self.transferOrbit, self.transferOrbitPC,                       \
                self.planeChangeDV, self.planeChangeDT =                    \
                    self.solve_lambert(self.startOrbit.prim.orb,            \
                                       self.endOrbit,                       \
                                       self.startTime, self.flightTime,     \
                                       self.planeChange, self.fixedEndTime);
            
            self.get_ejection_details()
            
            # Get arrival burn delta v
            if not self.ignoreInsertion:
                if self.planeChange:
                    vTrEnd = self.transferOrbitPC.get_state_vector(         \
                                    self.startTime + self.flightTime)[1];
                else:
                    vTrEnd = self.transferOrbit.get_state_vector(           \
                                    self.startTime + self.flightTime)[1]
                vEnd = self.endOrbit.get_state_vector(self.startTime +      \
                                                      self.flightTime)[1];
                self.insertionDV = vEnd - vTrEnd
            
            # Get phase angle
            self.phaseAngle = self.startOrbit.prim.orb                      \
                .get_angle_in_orbital_plane(                                \
                    self.get_departure_burn_time(),                         \
                    self.endOrbit.get_state_vector(                         \
                        self.get_departure_burn_time())[0]);
        
        # Third case: starting in orbit around a body and transferring to a
        # parking orbit around one of its satellites
        elif self.endOrbit.prim in self.startOrbit.prim.satellites:
            self.transferOrbit, self.transferOrbitPC,                       \
                self.planeChangeDV, self.planeChangeDT =                    \
                    self.solve_lambert(self.startOrbit,                     \
                                       self.endOrbit.prim.orb,              \
                                       self.startTime, self.flightTime,     \
                                       self.planeChange, self.fixedEndTime);
            
            if not self.ignoreInsertion: self.get_insertion_details()
            
            # Get departure burn delta v
            vStart = self.startOrbit.get_state_vector(self.startTime)[1]
            vTrStart = self.transferOrbit.get_state_vector(self.startTime)[1]
            self.ejectionDV = vTrStart - vStart
            
            # Get phase angle
            self.phaseAngle = self.startOrbit.get_angle_in_orbital_plane(   \
                    self.get_departure_burn_time(),                         \
                    self.endOrbit.prim.orb.get_state_vector(                \
                        self.get_departure_burn_time())[0]);
        
        # Fourth case: starting in a parking orbit around a body,  and then
        # transfering to another body and parking there. Both bodies orbit
        # the same primary.
        elif self.startOrbit.prim.orb.prim == self.endOrbit.prim.orb.prim:
            self.transferOrbit, self.transferOrbitPC,                       \
                self.planeChangeDV, self.planeChangeDT =                    \
                    self.solve_lambert(self.startOrbit.prim.orb,            \
                                       self.endOrbit.prim.orb,              \
                                       self.startTime, self.flightTime,     \
                                       self.planeChange, self.fixedEndTime);
            
            self.get_ejection_details()
            if not self.ignoreInsertion: self.get_insertion_details()
        
            # Get phase angle
            self.phaseAngle = self.startOrbit.prim.orb                      \
                .get_angle_in_orbital_plane(                                \
                    self.get_departure_burn_time(),                         \
                    self.endOrbit.prim.orb.get_state_vector(                \
                        self.get_departure_burn_time())[0]);
    
    
    def get_ejection_details(self):
        """Get ejection trajectory with burn details."""
        
        # Assume that parking orbit is circular
        ro = self.startOrbit.a              # distance from body at burn
        mu = self.startOrbit.prim.mu
        rEscape = self.startOrbit.prim.soi  # distance from body at escape
        
        # Get velocity of primary body and velocity needed after excape
        vPrim = self.startOrbit.prim.orb.get_state_vector(self.startTime)[1]
        vTrans = self.transferOrbit.get_state_vector(self.startTime)[1]
        
        # Excess velocity needed at escape from primary's sphere of influence
        vRel = vTrans - vPrim
        vRel = self.startOrbit.prim.orb.from_primary_to_orbit_bases(vRel)
        
        # speed after ejection burn
        vo = math.sqrt(norm(vRel)**2 + 2*(mu/ro - mu/rEscape))
        
        # escape trajectory elements
        e = math.sqrt(1+2*(vo**2/2 - mu/ro) * ro**2 * vo**2 / mu**2)
        a = 1 / (2/ro - vo**2/mu)
        
        # Describe positions at the SOI escape in the hyperbolic
        # escape trajectory's orbital plane
        # true anomaly at escape
        try:
            thetaEscape = math.acos(1/e * (a*(1-e**2)/rEscape - 1))
        except ValueError:
            thetaEscape = math.acos(
                math.copysign(1, 1/e * (a*(1-e**2)/rEscape - 1)))
        # flight path angle at escape
        phiEscape = math.atan(e*math.sin(thetaEscape) /                     \
                              (1+e*math.cos(thetaEscape)));
        # velocity vector at escape in orbital reference bases
        vEscape = math.sqrt(mu * (2/rEscape - 1/a)) *                       \
            np.array([math.cos(thetaEscape + math.pi/2 - phiEscape),        \
                      math.sin(thetaEscape + math.pi/2 - phiEscape),        \
                      0]);
        
        # start orbit basis vectors
        Xo, Yo, Zo = self.startOrbit.get_basis_vectors()
        
        # post-burn position and velocity vectors at periapsis, in the orbit's
        # reference bases
        roVec = self.startOrbit.from_primary_to_orbit_bases(ro * Xo)
        voVec = self.startOrbit.from_primary_to_orbit_bases(vo * Yo)
        
        # Reperesent the escape velocity in the orbital reference bases
        vRel = self.startOrbit.from_primary_to_orbit_bases(vRel)
        vEscape = self.startOrbit.from_primary_to_orbit_bases(vEscape)
        
        # Rotate the ejection trajector to match the desired escape velocity
        # An assumption is made that the periapsis lies in the primary body's
        # x-y plane.
        # First rotate around x-axis to match z-component
        phi = math.atan2(vRel[2], math.sqrt(abs(norm(vEscape)**2 -          \
                                                vEscape[0]**2 - vRel[2]**2)));
        R1 = np.array([[1,              0,              0],                 \
                       [0,      math.cos(phi),      -math.sin(phi)],        \
                       [0,      math.sin(phi),       math.cos(phi)]]);
        R1vEscape = np.matmul(R1, vEscape)
        
        # Then rotate around z-axis to match ejection direction
        theta = math.atan2(vRel[1], vRel[0]) - math.atan2(R1vEscape[1],     \
                                                          R1vEscape[0]);
        R2 = np.array([[math.cos(theta),    -math.sin(theta),   0],         \
                       [math.sin(theta),     math.cos(theta),   0],         \
                       [0,              0,              1]]);
        
        # Apply rotations
        roVec = np.matmul(R2, np.matmul(R1, roVec))
        voVec = np.matmul(R2, np.matmul(R1, voVec))
        
        # Represent periapsis state vector in primary bases
        roVec = self.startOrbit.from_orbit_to_primary_bases(roVec)
        voVec = self.startOrbit.from_orbit_to_primary_bases(voVec)
        
        # Get burn vector
        vPark = np.cross(Zo,roVec)
        vPark = math.sqrt(mu/ro) * vPark/norm(vPark)
        self.ejectionDV = voVec - vPark;
        
        # adjust so that the mean anomaly at epoch is compatible with the 
        # transfer start time
        if e < 1:
            # elliptical case
            eccAnomEscape = 2*math.atan(math.tan(thetaEscape/2) /           \
                                  math.sqrt((1+e) / (1-e)))
            dMeanAnom = eccAnomEscape - e*math.sin(eccAnomEscape)
        else:
            # hyperbolic case
            hypAnomEscape = math.copysign(                                  \
                                math.acosh((math.cos(thetaEscape)+e) /      \
                                       (1 + e*math.cos(thetaEscape))),      \
                                           thetaEscape);
            dMeanAnom = e*math.sinh(hypAnomEscape)
        
        # Add the correct ejection trajectory and duration to the transfer
        self.ejectionDT = abs(dMeanAnom * math.sqrt((abs(a))**3/mu))
        self.ejectionTrajectory =                                           \
            Orbit.from_state_vector(roVec, voVec,                           \
                                    self.get_departure_burn_time(),         \
                                    self.startOrbit.prim);
        
        # Get angle from prograde of ejection burn
        progradeAngle =                                                     \
            self.startOrbit.prim.orb.get_angle_in_orbital_plane(            \
                0,                                                          \
                self.startOrbit.prim.orb.from_primary_to_orbit_bases(       \
                    self.startOrbit.prim.orb.get_state_vector(              \
                        self.get_departure_burn_time())[1]));
        burnAngle =                                                         \
            self.startOrbit.prim.orb.get_angle_in_orbital_plane(            \
                0,                                                          \
                self.ejectionTrajectory.get_state_vector(                   \
                    self.get_departure_burn_time())[0]);
        
        self.ejectionBurnAngle = Orbit.map_angle(burnAngle-progradeAngle)
        if self.ejectionBurnAngle > math.pi:
            self.ejectionBurnAngle = self.ejectionBurnAngle - 2*math.pi
    
    
    def get_insertion_details(self):
        """Get insertion burn delta V and time of burn."""
        
        # Assume that parking orbit is circular
        ro = self.endOrbit.a              # distance from body at burn
        mu = self.endOrbit.prim.mu
        rEnc = self.endOrbit.prim.soi     # distance from body at encounter
        
        # Get velocity of primary body and velocity at encounter
        vPrim = self.endOrbit.prim.orb.get_state_vector(self.startTime +    \
                                                        self.flightTime)[1];
        vTrans = self.transferOrbit.get_state_vector(self.startTime +       \
                                                     self.flightTime)[1];
        
        # Excess velocity at encounter at primary's sphere of influence
        vRel = vTrans - vPrim
        vRel = self.endOrbit.prim.orb.from_primary_to_orbit_bases(vRel)
        
        # speed after immediately before insertion burn
        vo = math.sqrt(norm(vRel)**2 + 2*(mu/ro - mu/rEnc))
        
        # We won't be concerned with direction when it comes to the insertion,
        # so we can subtract parking speed from speed before burn to get the
        # delta v needed for insertion (scalar)
        vPark = math.sqrt(mu/ro)
        self.insertionDV = vo - vPark
        
        # Get time interval between encounter and burn 
        e = math.sqrt(1+2*(vo**2/2 - mu/ro) * ro**2 * vo**2 / mu**2)
        a = 1 / (2/ro - vo**2/mu)
        try:
            thetaEnc = -math.acos(1/e * (a*(1-e**2)/rEnc - 1))
        except ValueError:
            thetaEnc = -math.acos(
                math.copysign(1, 1/e * (a*(1-e**2)/rEnc - 1)))
        if e < 1:
            # elliptical case
            eccAnomEnc = 2*math.atan(math.tan(thetaEnc/2) /                 \
                                  math.sqrt((1+e) / (1-e)))
            dMeanAnom = eccAnomEnc - e*math.sin(eccAnomEnc)
        else:
            # hyperbolic case
            hypAnomEnc = math.copysign(                                     \
                                math.acosh((math.cos(thetaEnc)+e) /         \
                                       (1 + e*math.cos(thetaEnc))),         \
                                           thetaEnc);
            dMeanAnom = e*math.sinh(hypAnomEnc) - hypAnomEnc
        
        self.insertionDT = abs(dMeanAnom * math.sqrt((abs(a))**3/mu))
    
    
    def get_departure_burn_time(self):
        """Get the time since epoch of departure burn.
        
        Returns:
            the time in seconds at departure burn
        """
        
        return self.startTime - self.ejectionDT
    
    
    def get_arrival_burn_time(self):
        """Get the time since epoch of arrival burn.
        
        Returns:
            the time in seconds at arrival burn
        """
        
        return self.startTime + self.flightTime + self.insertionDT
    
    
    def get_plane_change_time(self):
        """Get the time since epoch at the plane change maneuver.
        
        Returns:
            the time in seconds at plane change maneuver
        """
        
        return self.startTime + self.planeChangeDT
    
    
    def get_total_delta_V(self):
        """Get total delta V required for all parts of transfer.
        
        Returns:
            total delta v across all maneuvers (m/s)
        """
        
        return norm(self.ejectionDV) + norm(self.planeChangeDV) +           \
            norm(self.insertionDV)