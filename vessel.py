from orbit import Orbit
from body import Body

class Vessel:
    """Spacecraft defined by its orbit and maneuver nodes.
    
    Attributes:
        name (string): the vessel's name
        orb (Orbit): Keplerian orbit
        maneuverNodes (list): list  of lists containing prograde, normal, and 
            radial delta V and universal time
        color (tuple): 3 ints for and RGB color for plotting
        
    """
    
    def __init__(self, name, orb, nodes = None, color = (255,255,255)):
        self.name = name
        self.orb = orb
        self.maneuverNodes = nodes
        self.color = color
    
    def add_maneuver_node(self, prograde, normal, radial, time):
        self.maneuverNodes.append([prograde, normal, radial, time])
    
    def add_to_primary(self):
        """Adds the body to its primary's list of satellites."""
        if not self.orb.prim is None:
            if not(self in self.orb.prim.satellites):
                self.orb.prim.satellites.append(self)
