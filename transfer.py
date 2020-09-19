import math
from random import gauss
import numpy as np
from numpy.linalg import norm
from orbit import Orbit
from body import Body
from copy import copy

class Transfer:
    """Orbital transfer from a starting orbit to an ending orbit.
    
    Attributes:
        startOrbit (Orbit): orbit prior to departure burns
        endOrbit (Orbit): orbit following arrival burns
        startTime (float): time at the beginning of transfer trajectory (s)
        flightTime (float): duration of transfer trajectory (s)
        planeChange (bool): if true, a mid-course plane change is done
        ignoreInsertion (bool): if true, arrival burn is ignored.
        cheapStartOrb (bool): if true, the only parameter of the starting
            park orbit used is the semimajor axis
        cheapEndOrb (bool): if true, the only parameter of the ending park
            orbit used is the semimajor axis
        startPos (vector): if provided, fixes start location of the transfer
                orbit given position
        endPos (vector): : if provided, fixes target location of the transfer
                orbit at the given position
        transferOrbit (Orbit): orbital trajectory between start and end.
            If there is a plane change maneuver, this is the portion of the
            trajectory prior to the maneuver.
        transferORbitPC (Orbit): If there is a plane change maneuver, this
            is the portion of the transfer trajectory after the maneuver.
        ejectionTrajectory (Orbit): If ejection from a starting body occurs,
            this is the trajectory between parking orbit and escape.
        insertionTrajectory (Orbit): If capture at a target body occurs,
            this is the trajectory between encounter and parking.
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
        convergenceFail (bool): If true, start and end positions for the
            Lambert problem did not converge via the genetic algorithm
            
    """
    
    def __init__(self, startOrbit, endOrbit, startTime, flightTime, 
                 planeChange = False, ignoreInsertion = False,
                 cheapStartOrb = False, cheapEndOrb = True):
        
        # Assign input attributes
        self.startOrbit = startOrbit
        self.endOrbit = endOrbit
        self.startTime = startTime
        self.flightTime = flightTime
        self.planeChange = planeChange
        self.ignoreInsertion = ignoreInsertion
        self.cheapStartOrb = cheapStartOrb
        self.cheapEndOrb = cheapEndOrb
        
        self.originalStartOrbit = copy(self.startOrbit)
        self.originalEndOrbit = copy(self.endOrbit)
        
        # These attributes get defined here but are filled in with methods
        self.startPos = None
        self.endPos = None
        self.transferOrbit = None
        self.transferOrbitPC = None
        self.ejectionTrajectory = None
        self.insertionTrajectory = None
        self.ejectionDV = 0
        self.insertionDV = 0
        self.planeChangeDV = 0
        self.ejectionDT = 0
        self.insertionDT = 0
        self.planeChangeDT = 0
        self.phaseAngle = 0
        self.ejectionBurnAngle = None
        self.convergenceFail = True
        
        # Calculate transfer, ejection, and insertion parameters
        self.get_transfer_details()
        # self.genetic_refine()
    
    
    # @staticmethod
    # def solve_lambert_izzo(startOrbit, endOrbit, startTime, flightTime,
    #                   planeChange = False, startPos = None, endPos = None,
    #                   tol = 1E-14, maxIt = 20):
    #     """Solves the Lambert problem to obtain a trajectory to the target.
        
    #     Args:
    #         startOrbit (Orbit): orbit prior to departure
    #         endOrbit (Orbit): orbit following arrival
    #         startTime (float): time at the beginning of transfer (s)
    #         flightTime (float): duration of transfer trajectory (s)
    #         planeChange (bool): if true, a mid-course plane change occurs
    #         startPos (vector): if provided, fixes start location at the
    #             given position
    #         endPos (vector): : if provided, fixes target location at the
    #             given position
    #         tol (float): the maximum tolerance for iteration termination
    #         maxIt (int): the maximum number of iterations before breaking
        
    #     Returns:
    #         transferOrbit (Orbit): trajectory before plane change
    #         transferOrbitPC (Orbit): trajectory after plane change
    #         planeChangeDV (array): plane change burn vector (m/s)
    #         planeChangeDT (float): time between start and plane change (s)
    #     """
        
    #     # Set gravitational parameter for transfer orbit
    #     mu = startOrbit.prim.mu
        
    #     # Get start position
    #     if startPos is None:
    #         rStart = startOrbit.get_state_vector(startTime)[0]
    #     else:
    #         rStart = startPos
        
    #     # Adjust end position based on inputs
    #     if endPos is None:
    #         rEnd = endOrbit.get_state_vector(startTime + flightTime)[0]
    #     else:
    #         rEnd = endPos
    #     if planeChange:
    #         # Rotate the target orbit to be coplanar with the starting one
    #         rEnd = startOrbit.from_primary_to_orbit_bases(rEnd)
    #         rEnd = np.array([rEnd[0],rEnd[1],0]) /                          \
    #             norm(np.array([rEnd[0],rEnd[1]])) * norm(rEnd);
    #         rEnd = startOrbit.from_orbit_to_primary_bases(rEnd)
        
    #     # Code adapted from Oldenhuis' Lambert Solver, Izzo's method
    #     # non-dimensional units
    #     r1 = norm(rStart)
    #     r1vec = rStart/r1;
    #     V = math.sqrt(mu/r1);
    #     r2vec = rEnd/r1;
    #     T = r1/V;                  
    
    
    @staticmethod
    def solve_lambert(startOrbit, endOrbit, startTime, flightTime,
                      planeChange = False, startPos = None, endPos = None,
                      tol = 1E-6, maxIt = 200):
        """Solves the Lambert problem to obtain a trajectory to the target.
        
        Args:
            startOrbit (Orbit): orbit prior to departure
            endOrbit (Orbit): orbit following arrival
            startTime (float): time at the beginning of transfer (s)
            flightTime (float): duration of transfer trajectory (s)
            planeChange (bool): if true, a mid-course plane change occurs
            startPos (vector): if provided, fixes start location at the
                given position
            endPos (vector): : if provided, fixes target location at the
                given position
            tol (float): the maximum tolerance for iteration termination
            maxIt (int): the maximum number of iterations before breaking
        
        Returns:
            transferOrbit (Orbit): trajectory before plane change
            transferOrbitPC (Orbit): trajectory after plane change
            planeChangeDV (array): plane change burn vector (m/s)
            planeChangeDT (float): time between start and plane change (s)
        """
        
        # Set gravitational parameter for transfer orbit
        mu = startOrbit.prim.mu
        
        # Get start position
        if startPos is None:
            rStart = startOrbit.get_state_vector(startTime)[0]
        else:
            rStart = startPos
        
        # Adjust end position based on inputs
        if endPos is None:
            rEnd = endOrbit.get_state_vector(startTime + flightTime)[0]
        else:
            rEnd = endPos
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
                print('Lambert solver failed to converge')
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
            if endPos is None:
                rEnd = endOrbit.get_state_vector(startTime + flightTime)[0]
            else:
                rEnd = endPos
            
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
            
            ## Calculate the inclination change needed for the maneuver
            # nTr = np.cross(rPC, vPCi)
            # nTr = nTr/norm(nTr)             # normal vector to pre-burn orbit
            
            # Get basis vectors for orbit after plane change
            zPC = np.cross(rPC, rEnd)
            zPC = zPC/norm(zPC)             # normal vector to post-burn orbit
            xPC= np.array([1, 0, 0])        # assumed celestial longitude
            xPC = xPC - np.dot(xPC,zPC) * zPC/norm(zPC)**2
            if norm(xPC) < 1E-15:
                xPC = np.array([0, math.copysign(1,zPC[0]), 0])
                zPC = np.array([math.copysign(1,zPC[0]), 0, 0])
            else:
                xPC = xPC / norm(xPC)
            yPC = np.cross(zPC,xPC)
            yPC = yPC/norm(yPC)
            
            # rotate velocity vector to new plane
            vPCf = transferOrbit.rotate_to_bases(vPCiPlane, xPC, yPC, True)
            
            # # normal vector to plane after burn
            # incPC = math.acos(np.dot(nTr,nTrPC))
            # if np.dot(nTrPC, vPCi) > 0:
            #     incPC = -incPC
            
            # # Rotate velocity vector prior to burn to get vector after burn
            # vPCfPlane = np.array([math.cos(incPC) * vPCiPlane[0],           \
            #                       math.cos(incPC) * vPCiPlane[1],           \
            #                       math.sin(incPC) * norm(vPCiPlane)]);
            
            # # Get the velocity after plane change in the primary bases
            # vPCf = transferOrbit.from_orbit_to_primary_bases(vPCfPlane)
            
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
                                       self.planeChange,                    \
                                       self.startPos, self.endPos);
            
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
            
            # Set start and end position for refinining
            self.startPos =                                                 \
                self.startOrbit.get_state_vector(self.startTime)[0];
            self.endPos =                                                   \
                self.endOrbit.get_state_vector(self.startTime +             \
                                               self.flightTime)[0];
        
        # Second case: starting in orbit around a body and transferring to a
        # parking orbit around its primary body
        elif self.startOrbit.prim in self.endOrbit.prim.satellites:
            self.transferOrbit, self.transferOrbitPC,                       \
                self.planeChangeDV, self.planeChangeDT =                    \
                    self.solve_lambert(self.startOrbit.prim.orb,            \
                                       self.endOrbit,                       \
                                       self.startTime, self.flightTime,     \
                                       self.planeChange,                    \
                                       self.startPos, self.endPos);
            
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
            
            # Set start and end position for refinining
            self.startPos =                                                 \
              self.startOrbit.prim.orb.get_state_vector(self.startTime)[0] +\
                self.ejectionTrajectory.get_state_vector(self.startTime)[0];
            self.endPos =                                                   \
                self.endOrbit.get_state_vector(self.startTime +             \
                                               self.flightTime)[0];
        
        # Third case: starting in orbit around a body and transferring to a
        # parking orbit around one of its satellites
        elif self.endOrbit.prim in self.startOrbit.prim.satellites:
            self.transferOrbit, self.transferOrbitPC,                       \
                self.planeChangeDV, self.planeChangeDT =                    \
                    self.solve_lambert(self.startOrbit,                     \
                                       self.endOrbit.prim.orb,              \
                                       self.startTime, self.flightTime,     \
                                       self.planeChange,                    \
                                       self.startPos, self.endPos);
            
            self.get_insertion_details()
            
            # Get departure burn delta v
            vStart = self.startOrbit.get_state_vector(self.startTime)[1]
            vTrStart = self.transferOrbit.get_state_vector(self.startTime)[1]
            self.ejectionDV = vTrStart - vStart
            
            # Get phase angle
            self.phaseAngle = self.startOrbit.get_angle_in_orbital_plane(   \
                    self.get_departure_burn_time(),                         \
                    self.endOrbit.prim.orb.get_state_vector(                \
                        self.get_departure_burn_time())[0]);
        
            # Set start and end position for refinining
            self.startPos =                                                 \
                self.startOrbit.get_state_vector(self.startTime)[0];
            self.endPos =                                                   \
              self.endOrbit.prim.orb.get_state_vector(                      \
                  self.startTime + self.flightTime)[0] +                    \
                self.insertionTrajectory.get_state_vector(                  \
                    self.startTime + self.flightTime)[0];
        
        # Fourth case: starting in a parking orbit around a body, and then
        # transfering to another body and parking there. Both bodies orbit
        # the same primary.
        elif self.startOrbit.prim.orb.prim == self.endOrbit.prim.orb.prim:
            self.transferOrbit, self.transferOrbitPC,                       \
                self.planeChangeDV, self.planeChangeDT =                    \
                    self.solve_lambert(self.startOrbit.prim.orb,            \
                                       self.endOrbit.prim.orb,              \
                                       self.startTime, self.flightTime,     \
                                       self.planeChange,                    \
                                       self.startPos, self.endPos);
            
            self.get_ejection_details()
            self.get_insertion_details()
        
            # Get phase angle
            self.phaseAngle = self.startOrbit.prim.orb                      \
                .get_angle_in_orbital_plane(                                \
                    self.get_departure_burn_time(),                         \
                    self.endOrbit.prim.orb.get_state_vector(                \
                        self.get_departure_burn_time())[0]);
            
            # Set start and end position for refinining
            self.startPos =                                                 \
              self.startOrbit.prim.orb.get_state_vector(self.startTime)[0] +\
                self.ejectionTrajectory.get_state_vector(self.startTime)[0];
            self.endPos =                                                   \
              self.endOrbit.prim.orb.get_state_vector(                      \
                  self.startTime + self.flightTime)[0] +                    \
                self.insertionTrajectory.get_state_vector(                  \
                    self.startTime + self.flightTime)[0];
        
        # Adjust phase angle to be within the range [-pi, pi]
        if self.phaseAngle < -math.pi:
            self.phaseAngle = self.phaseAngle + 2*math.pi
        elif self.phaseAngle > math.pi:
            self.phaseAngle = self.phaseAngle - 2*math.pi
    
    def get_ejection_details(self, tol = 0.1, maxIt = 50):
        """Get ejection trajectory with burn details."""
        
        # Assume that parking orbit is circular
        mu = self.startOrbit.prim.mu
        rEscape = self.startOrbit.prim.soi  # distance from body at escape
        
        # Get velocity of primary body and velocity needed after escape
        vPrim = self.startOrbit.prim.orb.get_state_vector(self.startTime)[1]
        vTrans = self.transferOrbit.get_state_vector(self.startTime)[1]

        err = tol + 1
        it = 0
        roNext = self.startOrbit.a
        while abs(err) > tol:
            it = it + 1
            if it > maxIt:
                # print('Ejection burn position failed to converge')
                # print(err)
                break
            
            # Periapsis radius of ejection trajectory (also burn position)
            ro = roNext
            
            # Excess velocity needed at escape from primary's sphere of influence
            vRel = vTrans - vPrim
            
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
            phiEscape = math.atan(e*math.sin(thetaEscape) /                 \
                                  (1+e*math.cos(thetaEscape)));
            # velocity vector at escape in orbital reference bases
            vEscape = math.sqrt(mu * (2/rEscape - 1/a)) *                   \
                np.array([math.cos(thetaEscape + math.pi/2 - phiEscape),    \
                          math.sin(thetaEscape + math.pi/2 - phiEscape),    \
                          0]);
            
            # start orbit basis vectors
            if not self.cheapStartOrb:
                Xo, Yo, Zo = self.startOrbit.get_basis_vectors()
            else:
                Xo = np.array([1,0,0])
                Yo = np.array([0,1,0])
                Zo = np.array([0,0,1])
            
            # post-burn position and velocity vectors at periapsis, in the orbit's
            # reference bases
            if not self.cheapStartOrb:
                roVec = self.startOrbit.from_primary_to_orbit_bases(ro * Xo)
                voVec = self.startOrbit.from_primary_to_orbit_bases(vo * Yo)
            else:
                roVec = ro * Xo
                voVec = vo * Yo
            
            # Reperesent the escape velocity in the orbital reference bases
            if not self.cheapStartOrb:
                vRel = self.startOrbit.from_primary_to_orbit_bases(vRel)
            
            # Rotate the ejection trajectory to match the desired escape velocity
            # An assumption is made that the periapsis lies in the primary body's
            # x-y plane.
            # First rotate around x-axis to match z-component
            phi = math.atan2(vRel[2], math.sqrt(abs(norm(vEscape)**2 -      \
                                             vEscape[0]**2 - vRel[2]**2)));
            R1 = np.array([[1,              0,              0],             \
                           [0,      math.cos(phi),      -math.sin(phi)],    \
                           [0,      math.sin(phi),       math.cos(phi)]]);
            R1vEscape = np.matmul(R1, vEscape)
            
            # Then rotate around z-axis to match ejection direction
            theta = math.atan2(vRel[1], vRel[0]) - math.atan2(R1vEscape[1], \
                                                              R1vEscape[0]);
            R2 = np.array([[math.cos(theta),    -math.sin(theta),   0],     \
                           [math.sin(theta),     math.cos(theta),   0],     \
                           [0,              0,              1]]);
            
            # Apply rotations
            roVec = np.matmul(R2, np.matmul(R1, roVec))
            voVec = np.matmul(R2, np.matmul(R1, voVec))
            
            # Represent periapsis state vector in primary bases
            if not self.cheapStartOrb:
                roVec = self.startOrbit.from_orbit_to_primary_bases(roVec)
                voVec =self.startOrbit.from_orbit_to_primary_bases(voVec)
            
            if self.cheapStartOrb:
                err = 0
            else:
                angleDiff = self.startOrbit.get_angle_in_orbital_plane(     \
                    self.startOrbit.get_time(0),roVec);
                roVecActual, vParkActual =                                  \
                    self.startOrbit.get_state_vector(                       \
                        self.startOrbit.get_time(angleDiff));
                prevErr = err
                err = norm(roVec) - norm(roVecActual)
                if abs(err)/abs(prevErr) > 0.9:
                    roNext = (norm(roVecActual) + ro)/2
                else:
                    roNext = norm(roVecActual)
        
        # Get burn vector
        if self.cheapStartOrb:
            Zo = np.matmul(R2, np.matmul(R1, Zo))
            vPark = np.cross(Zo,roVec)
            vPark = math.sqrt(mu/ro) * vPark/norm(vPark)
        else:
            vPark = vParkActual
        
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
            dMeanAnom = e*math.sinh(hypAnomEscape) - hypAnomEscape
        
        # Add the correct ejection trajectory and duration to the transfer
        self.ejectionDT = abs(dMeanAnom * math.sqrt((abs(a))**3/mu))
        self.ejectionTrajectory =                                           \
            Orbit.from_state_vector(roVec, voVec,                           \
                                    self.get_departure_burn_time(),         \
                                    self.startOrbit.prim);
        
        # Reset start parking orbit to match departure burn timing
        if self.cheapStartOrb:
            self.startOrbit =                                               \
                Orbit.from_state_vector(roVec,vPark,                        \
                                        self.get_departure_burn_time(),     \
                                        self.startOrbit.prim);
        
        # Get angle from prograde of ejection burn
        progradeAngle =                                                     \
            self.startOrbit.get_angle_in_orbital_plane(                     \
                0,                                                          \
                self.startOrbit.prim.orb.from_primary_to_orbit_bases(       \
                    self.startOrbit.prim.orb.get_state_vector(              \
                        self.get_departure_burn_time())[1]));
        burnAngle =                                                         \
            self.startOrbit.get_angle_in_orbital_plane(                     \
                0,                                                          \
                self.ejectionTrajectory.get_state_vector(                   \
                    self.get_departure_burn_time())[0]);
        
        self.ejectionBurnAngle = Orbit.map_angle(burnAngle-progradeAngle)
        if self.ejectionBurnAngle > math.pi:
            self.ejectionBurnAngle = self.ejectionBurnAngle - 2*math.pi
    
    
    def get_insertion_details(self, tol = 0.1, maxIt = 50):
        """Get insertion trajectory with burn details."""
        
        # Assume that parking orbit is circular
        mu = self.endOrbit.prim.mu
        rEnc = self.endOrbit.prim.soi  # distance from body at encounter
        
        # Get velocity of primary body and velocity needed at encounter
        vPrim = self.endOrbit.prim.orb.get_state_vector(self.startTime +    \
                                                        self.flightTime)[1];
        vTrans = self.transferOrbit.get_state_vector(self.startTime +       \
                                                      self.flightTime)[1];
        err = tol + 1
        it = 0
        roNext = self.endOrbit.a
        while abs(err) > tol:
            it = it + 1
            if it > maxIt:
                # print('Insertion burn position failed to converge')
                # print(err)
                break
            
            # Periapsis radius of insertion trajectory (also burn position)
            ro = roNext
        
            # Excess velocity at encounter at primary's sphere of influence
            vRel = vTrans - vPrim
            
            # speed before insertion burn
            vo = math.sqrt(norm(vRel)**2 + 2*(mu/ro - mu/rEnc))
            
            # insertion trajectory elements
            e = math.sqrt(1+2*(vo**2/2 - mu/ro) * ro**2 * vo**2 / mu**2)
            a = 1 / (2/ro - vo**2/mu)
            
            # Describe positions at the SOI encounter in the hyperbolic
            # insertion trajectory's orbital plane
            # true anomaly at escape
            try:
                thetaEncounter = -math.acos(1/e * (a*(1-e**2)/rEnc - 1))
            except ValueError:
                thetaEncounter = -math.acos(
                    math.copysign(1, 1/e * (a*(1-e**2)/rEnc - 1)))
            # flight path angle at encounter
            phiEncounter = math.atan(e*math.sin(thetaEncounter) /               \
                                  (1+e*math.cos(thetaEncounter)));
            # velocity vector at encounter in orbital reference bases
            vEncounter = math.sqrt(mu * (2/rEnc - 1/a)) *                    \
                np.array([math.cos(thetaEncounter + math.pi/2 - phiEncounter),  \
                          math.sin(thetaEncounter + math.pi/2 - phiEncounter),  \
                          0]);
            
            # end orbit basis vectors
            if not self.cheapEndOrb:
                Xo, Yo, Zo = self.endOrbit.get_basis_vectors()
            else:
                Xo = np.array([1,0,0])
                Yo = np.array([0,1,0])
                Zo = np.array([0,0,1])
            
            # pre-burn position and velocity vectors at periapsis, in the orbit's
            # reference bases
            if not self.cheapEndOrb:
                roVec = self.endOrbit.from_primary_to_orbit_bases(ro * Xo)
                voVec = self.endOrbit.from_primary_to_orbit_bases(vo * Yo)
            else:
                roVec = ro * Xo
                voVec = vo * Yo
            
            # Reperesent the encounter velocity in the orbital reference bases
            if not self.cheapEndOrb:
                vRel = self.endOrbit.from_primary_to_orbit_bases(vRel)
            
            # Rotate the insertion trajectory to match the desired escape velocity
            # An assumption is made that the periapsis lies in the primary body's
            # x-y plane.
            # First rotate around x-axis to match z-component
            phi = math.atan2(vRel[2], math.sqrt(abs(norm(vEncounter)**2 -       \
                                                vEncounter[0]**2 - vRel[2]**2)));
            R1 = np.array([[1,              0,              0],                 \
                           [0,      math.cos(phi),      -math.sin(phi)],        \
                           [0,      math.sin(phi),       math.cos(phi)]]);
            R1vEncounter = np.matmul(R1, vEncounter)
            
            # Then rotate around z-axis to match insertion direction
            theta = math.atan2(vRel[1], vRel[0]) - math.atan2(R1vEncounter[1],  \
                                                              R1vEncounter[0]);
            R2 = np.array([[math.cos(theta),    -math.sin(theta),   0],         \
                           [math.sin(theta),     math.cos(theta),   0],         \
                           [0,              0,              1]]);
            
            # Apply rotations
            roVec = np.matmul(R2, np.matmul(R1, roVec))
            voVec = np.matmul(R2, np.matmul(R1, voVec))
            
            # Represent periapsis state vector in primary bases
            if not self.cheapEndOrb:
                roVec = self.endOrbit.from_orbit_to_primary_bases(roVec)
                voVec = self.endOrbit.from_orbit_to_primary_bases(voVec)
            
            if self.cheapEndOrb:
                err = 0
            else:
                angleDiff = self.endOrbit.get_angle_in_orbital_plane(       \
                    self.endOrbit.get_time(0),roVec);
                roVecActual, vParkActual =                                  \
                    self.endOrbit.get_state_vector(                         \
                        self.endOrbit.get_time(angleDiff));
                prevErr = err
                err = norm(roVec) - norm(roVecActual)
                if abs(err)/abs(prevErr) > 0.9:
                    roNext = (norm(roVecActual) + ro)/2
                else:
                    roNext = norm(roVecActual)
        
        # Get burn vector
        if self.cheapEndOrb:
            Zo = np.matmul(R2, np.matmul(R1, Zo))
        
        vPark = np.cross(Zo,roVec)
        vPark = math.sqrt(mu/ro) * vPark/norm(vPark)
        
        if not self.ignoreInsertion:
            self.insertionDV = vPark - voVec;
        
        # adjust so that the mean anomaly at epoch is compatible with the 
        # transfer encounter time
        if e < 1:
            # elliptical case
            eccAnomEncounter = 2*math.atan(math.tan(thetaEncounter/2) /     \
                                  math.sqrt((1+e) / (1-e)))
            dMeanAnom = eccAnomEncounter - e*math.sin(eccAnomEncounter)
        else:
            # hyperbolic case
            hypAnomEncounter = math.copysign(                               \
                                math.acosh((math.cos(thetaEncounter)+e) /   \
                                       (1 + e*math.cos(thetaEncounter))),   \
                                           thetaEncounter);
            dMeanAnom = e*math.sinh(hypAnomEncounter) - hypAnomEncounter
        
        # Add the correct insertion trajectory and duration to the transfer
        self.insertionDT = abs(dMeanAnom * math.sqrt((abs(a))**3/mu))
        self.insertionTrajectory =                                          \
            Orbit.from_state_vector(roVec, voVec,                           \
                                    self.get_arrival_burn_time(),           \
                                    self.endOrbit.prim);
        
        # Reset end parking orbit to match departure burn timing
        if self.cheapEndOrb:
            self.endOrbit =                                                 \
                Orbit.from_state_vector(roVec,vPark,                        \
                                        self.get_arrival_burn_time(),       \
                                        self.endOrbit.prim);
    
    
    def adjust_start_orbit_mo(self):
        """Modify starting orbit to have mean anomaly at epoch matching burn.
        """
        
        burnTime = self.get_departure_burn_time()
        
        trueAnom = self.startOrbit.get_angle_in_orbital_plane(              \
            self.startOrbit.get_time(0),                                    \
            self.ejectionTrajectory.get_state_vector(burnTime)[0]);
        
        dMeanAnom = self.startOrbit.get_mean_anomaly(                       \
                                    self.startOrbit.get_time(trueAnom)) -   \
                    self.startOrbit.get_mean_anomaly(burnTime);
        
        self.startOrbit.mo =                                                \
            self.startOrbit.map_angle(self.startOrbit.mo + dMeanAnom);
    
    
    def adjust_end_orbit_mo(self):
        """Modify starting orbit to have mean anomaly at epoch matching burn.
        """
        
        burnTime = self.get_arrival_burn_time()
        
        trueAnom = self.endOrbit.get_angle_in_orbital_plane(                \
            self.endOrbit.get_time(0),                                      \
            self.insertionTrajectory.get_state_vector(burnTime)[0]);
        
        dMeanAnom = self.endOrbit.get_mean_anomaly(                         \
                                    self.endOrbit.get_time(trueAnom)) -     \
                    self.endOrbit.get_mean_anomaly(burnTime);
        
        self.endOrbit.mo =                                                  \
            self.endOrbit.map_angle(self.endOrbit.mo + dMeanAnom);
    
    
    def match_start_mean_anomaly(self, tol = 0.1, maxIt = 20):
        if self.ejectionTrajectory is None:
            self.genetic_refine()
            return
        self.startOrbit = copy(self.originalStartOrbit)
        originalStartTime = self.startTime
        it = 0
        err = tol+1
        while abs(err) > tol:
            it = it+1
            
            if it>maxIt:
                self.startTime = originalStartTime
                self.genetic_refine()
                return
            
            if it > 1:
                self.startPos = None
                self.endPos = None
                self.startTime = self.startTime + dT
                gen = self.genetic_refine()
                if gen is None:
                    break
            
            burnTime = self.get_departure_burn_time()
            orbPos = self.startOrbit.get_state_vector(burnTime)[0]
            burnPos = self.ejectionTrajectory.get_state_vector(burnTime)[0];
            
            burnTrueAnom = self.startOrbit.get_angle_in_orbital_plane(      \
                    self.startOrbit.get_time(0), burnPos);
            burnMeanAnom = self.startOrbit.get_mean_anomaly(                \
                self.startOrbit.get_time(burnTrueAnom))
            orbMeanAnom = self.startOrbit.get_mean_anomaly(burnTime)
            dMeanAnom = burnMeanAnom - orbMeanAnom
            
            while abs(dMeanAnom) > math.pi:
                dMeanAnom = dMeanAnom + math.copysign(2*math.pi, -dMeanAnom)
                
            dT = self.startOrbit.get_period() * dMeanAnom / (2*math.pi)
            err = norm(burnPos-orbPos)
        
        return
    
    
    def get_departure_burn_time(self):
        """Get the time since epoch of departure burn.
        
        Returns:
            the time in seconds at departure burn
        """
        
        return self.startTime - self.ejectionDT
    
    
    def get_encounter_time(self):
        """Get the time of SOI encounter with the target body.
        
        Returns:
            the time in seconds at target SOI encounter
        """
        
        return self.startTime + self.flightTime
    
    
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
    
    
    def genetic_refine(self, num = 10, tol = 1, maxGen = 40):
        """Genetic algorithm to find start and end positions for Transfer"""
        
        # TO DO: figure out better crossover/mutation methods, 
        # convergence for high-inclination transfers
        
        if (self.ejectionTrajectory is None) or                             \
            (self.insertionTrajectory is None):
            
            gen = 0
            maxGen = maxGen * 10
            err = self.get_error()
            while abs(err) > tol:
                gen = gen+1
                if gen > maxGen:
                    self.startPos = None
                    self.endPos = None
                    self.get_transfer_details()
                    self.convergenceFail = True
                    if not self.ejectionTrajectory is None:
                        self.adjust_start_orbit_mo()
                    if not self.insertionTrajectory is None:
                        self.adjust_end_orbit_mo()
                    return
                err = self.get_error()
            self.convergenceFail = False
            return gen
        
        startPositions, endPositions = self.get_first_generation(num)
        startPositions, endPositions, err =                                 \
            self.get_fitness(startPositions, endPositions);
        
        gen = 0
        while np.amin(err) > tol:
            gen = gen+1
            # if gen%100 == 0:
            #     print('.')
            if gen > maxGen:
                self.startPos = None
                self.endPos = None
                self.get_transfer_details()
                self.convergenceFail = True
                if not self.ejectionTrajectory is None:
                    self.adjust_start_orbit_mo()
                if not self.insertionTrajectory is None:
                    self.adjust_end_orbit_mo()
                return
            
            startPositions, endPositions = self.get_next_generation(        \
                startPositions, endPositions, err);
            
            startPositions, endPositions, err =                             \
                self.get_fitness(startPositions, endPositions);
        
        self.startPos = startPositions[0]
        self.endPos = endPositions[0]
        self.convergenceFail = False
        # if not self.ejectionTrajectory is None:
        #     self.adjust_start_orbit_mo()
        # if not self.insertionTrajectory is None:
        #     self.adjust_end_orbit_mo()
        return gen
    
    
    def get_first_generation(self, num = 10):
        """Gets the first generation of parents for genetic algorithm"""
        
        startPositions = []
        endPositions = []
        
        for x in range(math.ceil(num/2)):
            self.get_transfer_details()
            startPositions.append(self.startPos)
            endPositions.append(self.endPos)
            startMut, endMut = self.mutate(startPositions[-1],endPositions[-1])
            startPositions.append(startMut)
            endPositions.append(endMut)
        
        return startPositions, endPositions
    
    
    def get_error(self, startPos = None, endPos = None):
        """Gets error to serve as fitness for genetic algorithm."""
        
        if startPos is None:
            if self.startPos is None:
                self.startPos = self.startOrbit.get_state_vector(           \
                    self.startTime)[0]
            startPos = self.startPos
        if endPos is None:
            if self.endPos is None:
                self.endPos = self.endOrbit.get_state_vector(               \
                    self.startTime + self.flightTime)[0];
            endPos = self.endPos
        self.startPos = startPos
        self.endPos = endPos
        self.get_transfer_details()
        err = norm(self.startPos - startPos) + norm(self.endPos - endPos)
        return err
    
    
    def get_fitness(self, startPositions, endPositions):
        """Sorts population by fitness and returns array of errors"""
        
        err = []
        for x in range(len(startPositions)):
            err.append(self.get_error(startPositions[x], endPositions[x]))
        
        order = [i[0] for i in sorted(enumerate(err), key=lambda x:x[1])]
        startPositions = [startPositions[i] for i in order]
        endPositions = [endPositions[i] for i in order]
        err = [err[i] for i in order]
        
        return startPositions, endPositions, err
    
    
    def get_next_generation(self, startPositions, endPositions, err):
        """ Gets the next generation for the genetic algorithm"""
        
        err = np.array(err)
        fitness = 1/err;
        
        probs = []
        for fit in fitness:
            probs.append((fit)/sum(fitness))
        # probs = 1 - probs;
        probs = np.cumsum(probs)
        
        nextStartPositions = [startPositions[0], startPositions[1]]
        nextEndPositions = [endPositions[0], endPositions[1]]
        for x in range(len(startPositions)-2):
            val1 = np.random.rand()
            for y, prob in enumerate(probs):
                if val1 < prob:
                    p1 = y
                    break
                p1 = len(probs)-1
            val2 = np.random.rand()
            for z, prob in enumerate(probs):
                if val2 < prob:
                    p2 = z
                    break
                p2 = len(probs)-1
            sPos, ePos = self.crossover(
                [startPositions[p1], startPositions[p2]],
                [endPositions[p1], endPositions[p2]],
                [err[p1], err[p2]]);
            val3 = np.random.rand()
            if val3 < 0.25:
                sPos, ePos = self.mutate(sPos, ePos)
            nextStartPositions.append(sPos)
            nextEndPositions.append(ePos)
        
        return nextStartPositions, nextEndPositions
    
    
    def crossover(self, starts, ends, errs):
        """Combines positions with random weighted average."""
        
        startBodyPos = self.startOrbit.prim.orb.get_state_vector(           \
            self.startTime)[0];
        endBodyPos = self.endOrbit.prim.orb.get_state_vector(               \
            self.startTime + self.flightTime)[0];
        
        startSOI = self.startOrbit.prim.soi
        endSOI = self.endOrbit.prim.soi
        
        minErrIndex = np.argmin(errs)
        maxErrIndex = np.argmax(errs)
        errRatio = errs[minErrIndex]/errs[maxErrIndex]
        # errRatio = 0.5
        
        startRatio = gauss(1, errRatio)
        endRatio = gauss(1, errRatio)
        
        start = starts[minErrIndex]*startRatio +                            \
                starts[maxErrIndex]*(1-startRatio);
        start = start - startBodyPos
        start = start/norm(start) * startSOI
        start = start + startBodyPos
        
        end = ends[minErrIndex]*endRatio +                                  \
              ends[maxErrIndex]*(1-endRatio);
        end = end - endBodyPos
        end = end/norm(end) * endSOI
        end = end + endBodyPos
        
        randomVal = np.random.rand()
        if randomVal < 0.1:
            self.startPos = start
            self.endPos = end
            self.get_transfer_details()
            start = self.startPos
            end = self.endPos
        
        return start, end
    
    
    def mutate(self, start, end):
        """Modifies positions by randomly changing spherical angles."""
        
        startBodyPos = self.startOrbit.prim.orb.get_state_vector(           \
            self.startTime)[0];
        endBodyPos = self.endOrbit.prim.orb.get_state_vector(               \
            self.startTime + self.flightTime)[0];
        
        err = self.get_error(start, end)
        
        start = start - startBodyPos
        end = end - endBodyPos
        
        # angle = err/(2*self.startOrbit.prim.soi)*180 * math.pi/180
        # vecStart = np.array([np.random.rand()-0.5, np.random.rand()-0.5])
        # vecStart = vecStart/norm(vecStart) * angle
        # vecEnd = np.array([np.random.rand()-0.5, np.random.rand()-0.5])
        # vecEnd = vecEnd/norm(vecEnd) * angle
        
        # thetaStart =  math.atan2(start[1],start[0])
        # phiStart = math.acos(start[2]/norm(start))
        # thetaEnd = math.atan2(end[1], end[0])
        # phiEnd =  math.acos(end[2]/norm(end))
        
        # thetaStart = thetaStart + vecStart[0]
        # phiStart = phiStart + vecStart[1]
        # thetaEnd = thetaEnd + vecEnd[0]
        # phiEnd = phiEnd + vecEnd[1]
        
        # start = norm(start) * np.array(
        #     [math.sin(phiStart)*math.cos(thetaStart),
        #       math.sin(phiStart)*math.sin(thetaStart),
        #       math.cos(phiStart)])
        # end = norm(end) * np.array(
        #     [math.sin(phiEnd)*math.cos(thetaEnd),
        #       math.sin(phiEnd)*math.sin(thetaEnd),
        #       math.cos(phiEnd)])
        
        
        angle = err/(2*self.startOrbit.prim.soi) * math.pi
        start = start + self.startOrbit.get_basis_vectors()[2] *            \
                        math.sin(gauss(0,angle)) * self.startOrbit.prim.soi;
        end = end + self.endOrbit.get_basis_vectors()[2] *                  \
                    math.sin(gauss(0,angle)) * self.endOrbit.prim.soi;
        
        start = start/norm(start) * self.startOrbit.prim.soi
        end = end/norm(end) * self.endOrbit.prim.soi
        
        start = start + startBodyPos
        end = end + endBodyPos
        
        return start, end
