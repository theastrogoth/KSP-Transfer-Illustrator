import orbit

class Body:
    """Celestial body defined by its physical characeteristics and its orbit.
    
    Attributes:
        name (string): the body's name
        eqr (float): equatorial radius (meters)
        mu (float): standard gravitational parameter (m^3/s^2)
        soi (float): sphere of influence for patched conics (meters)
        orb (Orbit): Keplerian orbit
        satellites (list): list of other bodies orbiting this one
        color (tuple): 3 ints for and RGB color for plotting
        
    """
    
    def __init__(self, name  = None, eqr = None, mu = None, soi = None, 
                 orb = None, ref = None, satellites = None, 
                 color = (255,255,255)):
        self.name = name
        self.eqr = eqr
        self.mu = mu
        self.ref = ref
        self.color = color
        
        if orb is None:
            self.orb = orbit.Orbit(prim = self)
        else:
            self.orb = orb
            self.add_to_primary()
            
        if (soi is None) and not (orb is None):
            self.soi = orb.a * (mu/orb.prim.mu)**(2/5)
        else:
            self.soi = soi
        
        if satellites is None:
            self.satellites = []
        else:
            self.satellites = satellites
    
    def add_to_primary(self):
        """Adds the body to its primary's list of satellites."""
        if not self.orb.prim is None:
            if not(self in self.orb.prim.satellites):
                self.orb.prim.satellites.append(self)
    
    def rescale(self, factor):
        if not self.orb.a is None:
            self.orb.a = self.orb.a*factor
            self.soi = self.orb.a * (self.mu/self.orb.prim.mu)**(2/5)
    
    def resize(self, factor):
        surfGrav = self.mu/self.eqr**2
        
        self.eqr = self.eqr*factor
        self.mu = surfGrav * self.eqr**2
        
        if not self.orb.a is None:
            self.soi = self.orb.a * (self.mu/self.orb.prim.mu)**(2/5)

