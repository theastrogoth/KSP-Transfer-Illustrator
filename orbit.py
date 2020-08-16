import math
import numpy as np
from numpy.linalg import norm
from copy import copy
from body import Body

class Orbit:
    """Orbital trajectory defined by Keplerian elements and a primary body.
    
    Attributes:
        a (float): semimajor axis (meters)
        ecc (float): eccentricity (-/-)
        inc (float): inclination (radians)
        argp (float): argument of the periapsis (radians)
        lan (float): longitude of the ascending node (radians)
        mo (float): mean anomaly (radians) at epoch, t=0 seconds.
        prim (Body): the body at the node (in the middle) of the orbit
        timeshift (float): relative time in seconds of previous epoch
        period (float): length of an orbital period (s)
        X (array): first basis vector
        Y (array): second basis vector
        Z (array): third basis vector (normal to orbital plane)
        
    """
    
    def __init__(self, a=None, ecc=None, inc=None, 
                 argp=None, lan=None, mo=None, prim=None, timeShift=None):
        
        # These attributes will be filled in via methods
        self.period = None
        self.X = None
        self.Y = None
        self.Z = None
        
        # These attributes are set at initialization from input
        self.a = a
        self.ecc = ecc
        self.inc = inc
        self.argp = argp
        self.lan = lan
        self.prim = prim
        
        if timeShift is None:
            self.mo = mo
        else:
            self.mo = self.map_angle(mo - timeShift*2*math.pi               \
                                     / self.get_period())
        
    @classmethod
    def from_state_vector(cls,pos,vel,t,primaryBody):
        """Generates an Orbit object from position and velocity vectors.
       
        Args:
            pos (array): three dimensional position vector (meters)
            vel (array): three dimensional velocity vector (meters)
            t (float): time at position occurrence (seconds)
            primaryBody (Body): the body at the node of the orbit
        """
        
        # Copy input vectors to prevent mutation
        r = copy(pos)
        drdt = copy(vel)
        # Shape vectors to row vectors for numpy math
        r.shape = (3,)
        drdt.shape = (3,)
        
        # If the position vector is at the origin, throw an exception
        if norm(r) == 0:
            raise Exception('invalid position');
        
        # Calculate the semimajor axis
        a = 1 / (2/norm(r) - ((norm(drdt))**2)/primaryBody.mu)
        
        # Calculate the orbital momentum vector (perpendicular to plane)
        h = np.cross(r,drdt)
        
        # Calculate the inclination using the momentum vector
        inc = math.acos(h[2]/norm(h))
        
        # Calculate eccentricity vector (points from apoapsis to periapsis)
        # The magnitude of the vector is the orbit's eccentricity
        eccVec = np.cross(drdt,h)/primaryBody.mu - r/norm(r)
        ecc = norm(eccVec)
        
        # If orbit is circular, set periapsis at reference direction
        if norm(eccVec) == 0:
            eccVec = np.array([1,0,0])
        
        # Calculate the vector pointing toward the ascending node
        n = np.cross(np.array([0,0,1]),h)
        
        # If orbit is non-inclined, set ascending node at periapsis
        if norm(n) == 0:
            n = eccVec
        
        # Calculate the longitude of the ascending node
        if n[1]>=0:
            lan = math.acos(n[0]/norm(n))
        else:
            lan = 2*math.pi -                                               \
                math.acos(n[0]/norm(n))
        
        # Calculate the argument of the periapsis
        if (n == eccVec).all():
            argp = 0
        elif eccVec[2] >= 0:
            try:
                argp = math.acos(np.dot(n,eccVec)/(norm(n)*norm(eccVec)))
            except ValueError:
                argp = math.acos(math.copysign(1, np.dot(n,eccVec)))
        else:
            try:
                argp = 2*math.pi -                                          \
                    math.acos(np.dot(n,eccVec) / (norm(n)*norm(eccVec)))
            except ValueError:
                argp = math.acos(math.copysign(1, np.dot(n,eccVec)))
        
        # Calculate the true anomaly of the position at time t
        if np.dot(r,drdt) >= 0:
            try:
                nu = math.acos(np.dot(r,eccVec) / (norm(r)*norm(eccVec)))
            except ValueError:
                nu = math.acos(math.copysign(1, np.dot(r,eccVec)))
        else:
            try:
                nu = 2*math.pi -                                            \
                    math.acos(np.dot(r,eccVec) / (norm(r)*norm(eccVec)))
            except ValueError:
                nu = math.acos(math.copysign(1, np.dot(r,eccVec)))
        
        # Calculate the mean anomaly when t=0.
        # Parabolic case
        if ecc == 1:
            raise Exception('parabolic case (e=1) not implemented')
            
        # Elliptical case
        elif ecc < 1:
            # Get eccentric anomaly at time t
            eccAnom = 2*math.atan(math.tan(nu/2) /                          \
                                  math.sqrt((1+ecc) / (1-ecc)))
            
            # Get mean anomaly at time t via Kepler's equation
            meanAnom = eccAnom - ecc*math.sin(eccAnom)
            
            # Subtract angle traveled in time t to get mean anomaly at t=0
            mo = meanAnom - t*math.sqrt(primaryBody.mu/a**3)
            
            # Adjust angle to fall with 0 and 2pi radians
            mo = Orbit.map_angle(mo)
        
        # Hyperbolic case
        elif ecc > 1:
            # Get hyperbolic anomaly at time t
            try:
                hypAnom = math.acosh((ecc+math.cos(nu)) /                   \
                                     (1+ecc*math.cos(nu)))
            except ValueError:
                hypAnom = 2*math.atanh(math.sqrt((ecc-1)/(ecc+1))           \
                                       * math.tan(nu/2));
            if nu < 0:
                hypAnom = -hypAnom
            
            # Get mean anomaly at time t via Kepler's equation
            meanAnom = ecc * math.sinh(hypAnom) - hypAnom
            
            # Subtract angle traveled in time t to get mean anomaly at t=0
            mo = meanAnom - t*math.sqrt(primaryBody.mu/-a**3)
        
        else:
            raise Exception('invalid eccentricity')
        # Return an Orbit with the calculated Keplerian elements
        return cls(a,ecc,inc,argp,lan,mo,primaryBody);
    
    
    @staticmethod
    def map_angle(theta):
        """Maps an angle to the range [0 2*pi).
        
        Args:
            theta (float): an angle in radians
        
        Returns:
            an angle between 0 and 2*pi
        """
        
        if theta < 0 or theta >= 2*math.pi:
            twoPiRatio = theta/(2*math.pi)
            return theta - 2*math.pi*math.floor(twoPiRatio);
        else:
            return theta;
    
    
    @staticmethod
    def solve_Keplers(meanAnom, ecc, tol = 1E-12, maxIt = 1000):
        """Solves the inverse Kepler's equation for eccentric anomaly.
        
        Args:
            meanAnom (float): mean anomaly (radians)
            ecc (float): eccentricity (-/-)
            tol (float): error tolerance for iteration
            maxIt (int): maximum number of iterations before terminating
            
        Returns:
            eccentric or hyperbolic anomaly (radians) as appropriate
        """
        
        # Parabolic case
        if ecc == 1:
            raise Exception('parabolic case (e=1) not implemented')
            
        # Elliptical case
        elif ecc < 1:
            # Set first guess before iterating
            if ecc < 0.08:
                eccAnom = meanAnom
            else:
                eccAnom = math.pi
            
            # Initialize previous guess at a value guaranteed not to pass
            eccAnomPrev = (meanAnom+math.pi)*4
            
            # Iterate using Newton's method to find eccentric anomaly
            it = 0
            while abs(eccAnom-eccAnomPrev) > tol:
                it = it+1
                if it > maxIt:
                    break
                eccAnomPrev = eccAnom
                eccAnom = eccAnom -                                         \
                    (eccAnom - ecc*math.sin(eccAnom) - meanAnom)/           \
                        (1-ecc*math.cos(eccAnom))
            return eccAnom;
        
        # Hyperbolic case
        else:
            # Set first guess before iterating
            if abs(meanAnom) > 4*math.pi:
                hypAnom = math.copysign(4*math.pi, meanAnom)
            else:
                hypAnom = meanAnom
            
            # Initialize previous guess at a value guaranteed not to pass
            hypAnomPrev = (meanAnom+math.pi)*4
            
            # Iterate using Newton's method to find hyperbolic anomaly
            it = 0
            while abs(hypAnom-hypAnomPrev) > tol:
                it = it+1
                if it > maxIt:
                    break
                hypAnomPrev = hypAnom
                try:
                    hypAnom = hypAnom +                                     \
                        (meanAnom - ecc*math.sinh(hypAnom) + hypAnom)/      \
                            (ecc*math.cosh(hypAnom) - 1);
                except OverflowError:
                    hypAnom = math.pi
            return hypAnom
        
        # If this point is reached, the input must be invalid
        raise Exception('invalid eccentricity')
    
    
    @staticmethod
    def rotate_to_bases(vec, basis1, basis2, invert = False):
        """Rotates a vector to a new set of bases.
        
        Args:
            vec (array): three dimensional vector
            basis1 (array): three dimensional vector
            basis2 (array): three dimensional vector
            invert (bool): if true, inverse rotation is performed
            
        Returns:
            a new three dimensional vector in the new set of bases
        """
        
        # Copy input vectors to prevent mutation
        r = copy(vec)
        X = copy(basis1)
        Y = copy(basis2)
        
        # Ensure vectors are orthonormal
        Y = Y - np.dot(X,Y)/norm(X)**2
        X = X / norm(X)
        Y = Y / norm(Y)
        Z = np.cross(X,Y)
        Z = Z / norm(Z)
        
        # # Calcualte Tiat-Bryan angles NO NEED, SINCE WE KNOW THE NEXT BASES
        # psi = math.atan2(X[1], X[0])
        # theta = math.atan2(-X[2], math.sqrt(1-X[2]**2))
        # phi = math.atan2(Y[2], math.sqrt(1-Y[2]**2))
        
        # psi = math.atan2(X[1], X[0])
        # theta = math.atan2(-X[2], math.sqrt(1-X[2]**2))
        # phi = math.atan2(Y[2], math.sqrt(1-Y[2]**2))
        
        # c1 = math.cos(psi);     s1 = math.sin(psi);
        # c2 = math.cos(theta);   s2 = math.sin(theta);
        # c3 = math.cos(phi);     s3 = math.sin(phi);
        
        # # Construct rotation matrix
        # rot = np.array([[c1*c2, c1*s2*s3-c3*s1, s1*s3+c1*c3*s2],            \
        #                [c2*s1, c1*c3+s1*s2*s3, c3*s1*s2-c1*s3],             \
        #                [-s2,   c2*s3,          c2*c3]]);
        
        # Transformation Matrix
        rot = np.array([X,Y,Z])
        
        # Apply the rotation (or inverse of the rotation)
        if invert:
            r = np.matmul(np.transpose(rot),r)
        else:
            r = np.matmul(rot,r)
        return r
    
    
    def get_period(self):
        """Returns the sidereal period of the orbit.
        
        Returns:
            orbital period (seconds)
        """
        if self.period is None:
            self.period = 2*math.pi * math.sqrt((abs(self.a)**3)/self.prim.mu)
        return self.period
    
    
    def get_mean_anomaly(self, t):
        """Returns the mean anomaly of the orbital trajectory at time t.
        
        Args:
            t (float): time (seconds)
        
        Returns:
            mean anomaly (radians)
        """
        
        if t == 0:
            meanAnom = self.mo
        else:
            meanAnom = self.mo + t / (self.get_period()/(2*math.pi))
            if self.ecc < 1:
                meanAnom = self.map_angle(meanAnom)
        return meanAnom
    
    
    def get_true_anomaly(self, t):
        """Returns the true anomaly of the orbital trajectory at time t.
        
        Args:
            t (float): time (seconds)
        
        Returns:
            true anomaly (radians)
        """
        
        # Get the mean anomaly
        meanAnom = self.get_mean_anomaly(t)
        
        # Parabolic case
        if self.ecc == 1:
            raise Exception('parabolic case (e=1) not implemented')
        
        # Elliptical case
        elif self.ecc < 1:
            # Solve Kepler's equation to get eccentric anomaly
            eccAnom = self.solve_Keplers(meanAnom, self.ecc)
            
            # Calculate true anomaly
            return 2*math.atan2(math.sqrt(1+self.ecc)*math.sin(eccAnom/2),  \
                                math.sqrt(1-self.ecc)*math.cos(eccAnom/2))
        
        # Hyperbolic case
        else:
            # Solve Kepler's equation to get hyperbolic anomaly
            hypAnom = self.solve_Keplers(meanAnom, self.ecc)
            
            # Calculate position vector in orbit's reference frame
            pos = -np.array(                                                \
                [self.a*(self.ecc-math.cosh(hypAnom)),                      \
                 self.a*math.sqrt(self.ecc**2 -1)*math.sinh(hypAnom),       \
                 0])
    
            # get true anomaly by finding angle of position in plane
            return math.atan2(pos[1],pos[0])
        
        # If this position is reached, the input must be invalid
        raise Exception('invalid eccentricity')
    
    
    def get_time(self, trueAnom, tMin=0):
        """Returns the time when the orbit has the input true anomaly.
        
        Args:
            trueAnom (float): true anomaly (radians)
            tMin (float): earliest possible time (seconds) for elliptical case
        
        Returns:
            first time (s) after tMin where the orbit has the input trueAnom
        """
        
        period = self.get_period()
        
        # Parabolic case
        if self.ecc == 1:
            raise Exception('parabolic case (e=1) not implemented')
        
        # Elliptical case
        elif self.ecc < 1:
            # Calculate the eccentric anomaly
            eccAnom = 2*math.atan(math.tan(trueAnom/2) /                    \
                                  math.sqrt((1+self.ecc) / (1-self.ecc)))
            
            # Calculate mean anomaly via Keplers equation
            meanAnom = eccAnom - self.ecc*math.sin(eccAnom)
            meanAnom = self.map_angle(meanAnom)
            
            # Calculate time passed since epoch (t=0)
            t = (meanAnom-self.mo)*period/(2*math.pi)
            
            # Adjust to move beyond minimum time
            diffRatio = (tMin-t)/period
            return t + math.ceil(diffRatio)*period
        
        # Hyperbolic case
        else:
            # Calculate hyperbolic anomaly
            try:
                hypAnom = math.acosh((self.ecc+math.cos(trueAnom)) /        \
                                     (1+self.ecc*math.cos(trueAnom)));
            except ValueError:
                hypAnom = 2*math.atanh(math.sqrt((self.ecc-1)/(self.ecc+1)) \
                                       * math.tan(trueAnom/2));
                    
            if trueAnom < 0:
                hypAnom = -hypAnom
            
            # Calculate mean anomaly via Keplers equation
            meanAnom = self.ecc*math.sinh(hypAnom) - hypAnom
            
            # Calculate time difference from epoch
            # Since hyperbolic orbits are not periodic, tMin is ignored
            return (meanAnom-self.mo)*period/(2*math.pi)
    
    
    def get_state_vector(self, t):
        """Returns the position and velocity vectors of the orbit at time t.
        
        Args:
            t (float): time (seconds)
        
        Returns:
            The position (m) and velocity (m/s) vectors in a list  
        """
       
        # If the orbit is for the system root, make it stationary at origin
        if self.a is None:
            r = np.array([0,0,0])
            drdt = np.array([0,0,0])
            return r, drdt;
        
        # Get true anomaly at time t
        nu = self.get_true_anomaly(t)
        
        # Get magnitude of position vector (meters)
        rMag = self.a*(1-self.ecc**2) / (1+self.ecc*math.cos(nu))
        
        # Get position in orbital frame (periapsis on +x axis)
        o = rMag * np.array([[math.cos(nu)],                                \
                             [math.sin(nu)],                                \
                             [0]])
        
        # Get flight path angle (radians)
        phi = math.atan(self.ecc*math.sin(nu) / (1+self.ecc*math.cos(nu)))
        
        # Get magnitude of velocity vector (meters/second)
        drdtMag = math.sqrt(self.prim.mu*(2/rMag - 1/self.a))
        
        # Get velocity vector in orbital frame
        dodt = drdtMag * np.array([[math.cos(nu+math.pi/2 - phi)],          \
                                   [math.sin(nu+math.pi/2 - phi)],          \
                                   [0]])
        
        # Set up rotation matrices to transform to primary reference frame
        # Rotation around z-axis to match longitude of ascending node
        R1 = np.array([[math.cos(-self.lan), -math.sin(-self.lan), 0],      \
                       [math.sin(-self.lan), math.cos(-self.lan), 0],       \
                       [0, 0, 1]])
        
        # Rotation around x-axis to match inclination
        R2 = np.array([[1, 0, 0],                                           \
                       [0, math.cos(-self.inc), -math.sin(-self.inc)],      \
                       [0, math.sin(-self.inc), math.cos(-self.inc)]])
        
        # Rotation around z-axis to match argument of the periapsis
        R3 = np.array([[math.cos(-self.argp), -math.sin(-self.argp), 0],    \
                       [math.sin(-self.argp), math.cos(-self.argp), 0],
                       [0, 0, 1]])
        
        # Apply rotations to position and velocity vectors
        R = np.transpose(np.matmul(R3, np.matmul(R2,R1)))
        r = np.matmul(R,o)
        drdt = np.matmul(R,dodt)
        r.shape = (3,)
        drdt.shape = (3,)
        return r, drdt
    
    
    def get_basis_vectors(self):
        """Returns the basis vectors for the reference plane of an orbit.
        
        Returns:
            X (array): a vector pointing toward the periapsis
            Y (array): a vector pointing with velocity at periapsis
            Z (array): a vector normal to the orbital plane
        
        """
        
        # Return standard bases if the orbit is for the solar system root
        if self.a is None:
            X = np.array([1, 0, 0])
            Y = np.array([0, 1, 0])
            Z = np.array([0, 0, 1])
            return X, Y, Z
        elif not self.X is None:
            return self.X, self.Y, self.Z
        else:
            # Get state vector at periapsis
            rp, vp = self.get_state_vector(self.get_time(0))
            
            # Get vector normal to the orbital plane
            Z = np.cross(rp,vp)
            Z = Z / norm(Z)
            
            # Set first basis vector based on celestial longitude
            X = np.array([1, 0, 0])     # assumed celestial longitude
            X = X - np.dot(X,Z) * Z/norm(Z)**2
            if norm(X) < 1E-15:
                X = np.array([0, math.copysign(1,Z[0]), 0])
                Z = np.array([math.copysign(1,Z[0]), 0, 0])
            else:
                X = X / norm(X)
            
            # Determine second basis through cross product of the others
            Y = np.cross(Z,X)
            Y = Y / norm(Y)
            
            self.X = X
            self.Y = Y
            self.Z = Z
            
            return X, Y, Z
    
    
    def from_primary_to_orbit_bases(self, vec):
        """Returns a copy of the input vector in the orbit's reference bases.
        
        Args:
            vec (array): a 3D vector in the primary body's reference bases.
        
        Returns:
            The input vector represented in the orbit's reference bases.
        """
        
        # Copy input vector to avoid mutation
        rPrim = copy(vec)
        
        # Get basis vectors for the orbit
        X, Y, Z = self.get_basis_vectors()
        
        # Perform rotation
        return self.rotate_to_bases(rPrim, X, Y, False)
    
    
    def from_orbit_to_primary_bases(self, vec):
        """Returns a copy of the input vector in the primary body's bases.
        
        Args:
            vec (array): a 3D vector in the orbit's reference bases.
        
        Returns:
            The input vector represented in the primary body's reference bases.
        """
        
        # Copy input vector to avoid mutation
        rPrim = copy(vec)
        
        # Get basis vectors for the orbit
        X, Y, Z = self.get_basis_vectors()
        
        # Perform rotation
        return self.rotate_to_bases(rPrim, X, Y, True)
    
    
    def get_positions(self, startTime = None, endTime = None, 
                      num = 101, times = None):
        """ Returns an array of position vectors.
            
        Args:
            startTime (float): the earliest time (s)
            endTIme (float): the latest time (s)
            num (int): the total number of positions, sampled evenly between 
            the start and end
            times (array): an array of all the times to be sampled
                            (ignores all other inputs)
                
        Returns:
            An array of position vectors at the specified times
        """
        # initialize positions array and evenly sample times
        positions = np.empty((0,3))
        velocities = np.empty((0,3))
        if times is None:
            times = np.linspace(startTime, endTime, num)
        
        # get the position vector for each time and add it to array
        for t in times:
            pos, vel = self.get_state_vector(t)
            # parent = self.prim
            # while not parent.orb.prim == parent:
            #     pos = parent.orb.from_orbit_to_primary_bases(pos)
            #     parent = parent.orb.prim
            positions = np.append(positions, [pos], axis=0)
            velocities = np.append(velocities, [vel], axis=0)
        return positions, velocities
    
    
    def get_angle_in_orbital_plane(self, t, vec):
        """Returns the angle between the position at time t and the input.
        
        Args:
            t (float): time (seconds)
            vec (array): a 3D vector in the primary body's reference bases.
        
        Returns:
            The angle between the position and the vector projected to the 
            orbital plane.
        """
        
        # Get position vector at the specified time
        r = self.get_state_vector(t)[0]
        
        # Represent position and input vectors in the orbital reference frame
        rPlane = self.from_primary_to_orbit_bases(r)
        vecPlane = self.from_primary_to_orbit_bases(vec)
        
        # Get angle between the vectors projected to the orbital plane
        thetaRPlane = math.atan2(rPlane[1], rPlane[0])
        thetaVecPlane = math.atan2(vecPlane[1], vecPlane[0])
        
        return self.map_angle(thetaVecPlane - thetaRPlane)
    
    def __str__(self):
        string = '  Semi-major axis: ' + "{:.2f}".format(self.a) + ' m\n' +  \
            '  Eccentricity: ' + "{:.6f}".format(self.ecc) + '\n'           \
            '  Inclination: '+"{:.6f}".format(self.inc*180/math.pi) + '°\n'+\
            '  Argument of Periapsis: ' +                                   \
                "{:.6f}".format(self.argp*180/math.pi) + '°\n' +            \
            '  Longitude of Ascending Node: ' +                             \
                "{:.6f}".format(self.lan*180/math.pi) + '°\n' +             \
            '  Mean Anomaly at Epoch: ' +                                   \
                "{:.6f}".format(self.mo) + ' radians\n' +                   \
            '  Primary Body: ' + self.prim.name;
        return string
