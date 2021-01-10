import orbit
import numpy as np

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
                 rotPeriod = None, rotIni = None,
                 orb = None, ref = None, satellites = None, 
                 color = (255,255,255)):
        self.name = name
        self.eqr = eqr
        self.mu = mu
        self.rotPeriod = rotPeriod
        self.rotIni = rotIni
        self.ref = ref
        self.color = color
        
        if orb is None:
            self.orb = orbit.Orbit(prim = self)
        else:
            self.orb = orb
            self.add_to_primary()
            
        if (soi is None) and not (orb is None):
            self.set_soi(orb.a, mu, orb.prim.mu)
        else:
            self.soi = soi
        
        if satellites is None:
            self.satellites = []
        else:
            self.satellites = satellites
    
    def set_soi(self, sma, mu, muPrim):
        self.soi = sma * (mu/muPrim)**(2/5)
    
    def add_to_primary(self):
        """Adds the body to its primary's list of satellites."""
        if not self.orb.prim is None:
            if not(self in self.orb.prim.satellites):
                self.orb.prim.satellites.append(self)
    
    def remove_from_primary(self):
        """Removes body from list of satellites of the primary body."""
        if not self.orb.prim is None:
            satList = [sat.name for sat in self.orb.prim.satellites]
            if self.name in satList:
                delIndex = satList.index(self.name)
                del self.orb.prim.satellites[delIndex]
                
    
    def sort_satellites(self):
        """Orders satellites by semimajor axis (increasing)."""
        sats = self.satellites
        SMAs = [sat.orb.a for sat in sats]
        idxs = np.argsort(SMAs)
        newSatNames = []
        newSMAs = []
        for idx in idxs:
            newSatNames.append(sats[idx].name)
            newSMAs.append(SMAs[idx])
        
        # use name as tiebreaker
        ii = 0
        while ii < len(newSatNames):
            a = newSMAs[ii]
            duplicates = [jj for jj, sma in enumerate(newSMAs) if sma==a]
            if len(duplicates) > 1:
                names = [newSatNames[jj] for jj in duplicates]
                metaIdxs = np.argsort(names)
                for kk, metaIdx in enumerate(metaIdxs):
                    newSatNames[ii+kk] = names[metaIdx]
            ii = ii+len(duplicates)
        
        newSats = []
        for name in newSatNames:
            sat = [bd for bd in sats if bd.name==name][0]
            newSats.append(sat)
        
        
        self.satellites = newSats
        return
    
    
    def rescale(self, factor):
        """Multiplies the SMA of the body's orbit by a given factor."""
        if not self.orb.a is None:
            self.orb.a = self.orb.a*factor
            self.soi = self.orb.a * (self.mu/self.orb.prim.mu)**(2/5)
    
    def resize(self, factor):
        """Mutliplies the radius and gravity param of a body by a given factor.
        """
        surfGrav = self.mu/self.eqr**2
        
        self.eqr = self.eqr*factor
        self.mu = surfGrav * self.eqr**2
        
        if not self.orb.a is None:
            self.soi = self.orb.a * (self.mu/self.orb.prim.mu)**(2/5)
    


