# -*- coding: utf-8 -*-

import math
import numpy as np
from orbit import Orbit
from body import Body
from transfer import Transfer

class PorkchopTable:
    """Table of delta v values for transfers between the specified orbits.
    
    Attributes:
        startOrbit (Orbit): orbit prior to departure burns
        endOrbit (Orbit): orbit following arrival burns
        transferType (string): specifies whether the transfer is ballistic,
            has a plane change maneuver, or is the "cheaper" of the two
        ignoreInsertion (bool): if true, arrival burn is ignored.
        fixedEndTime (float): if given as input, the target position 
            will be specified by the position of the end orbit at this 
            time (s). Used for setting up successive transfers.
        minStartTime (float): earliest time at the beginning of transfer 
            trajectory (s)
        minStartTime (float): latest time at the beginning of transfer 
            trajectory (s)
        minFlightTime (float): shortest duration of transfer trajectory (s)
        maxFlightTime (float): longest duration of transfer trajectory (s)
        startTimeSize (int): number of samples taken on the start time axis
        flightTimeSize (int): number of samples taken on the flight time axis
        deltaV: a table of values with the sum of the magnitue of all burn 
            maneuvers (m/s) at each choice of start and flight times
    
    """
    
    def __init__(self, startOrbit, endOrbit, transferType = 'ballistic',
                 ignoreInsertion = False, fixedEndTime = None, 
                 minStartTime = 0, maxStartTime = None, 
                 minFlightTime = None, maxFlightTime = None,
                 startTimeSize = 101, flightTimeSize = 101):
        
        self.startOrbit = startOrbit
        self.endOrbit = endOrbit
        self.transferType = transferType
        self.ignoreInsertion = ignoreInsertion
        self.fixedEndTime = fixedEndTime
        self.minStartTime = minStartTime
        self.startTimeSize = startTimeSize
        self.flightTimeSize = flightTimeSize
        
        if (endOrbit.prim in startOrbit.prim.satellites or                  \
            startOrbit.prim == endOrbit.prim):
                startPeriod = startOrbit.get_period()
        else:
            startPeriod = startOrbit.prim.orb.get_period()
            
        if (startOrbit.prim in endOrbit.prim.satellites or                  \
            endOrbit.prim == startOrbit.prim):
                endPeriod = endOrbit.get_period()
        else:
            endPeriod = endOrbit.prim.orb.get_period()
        
        if maxStartTime is None:
            self.maxStartTime = minStartTime + 2 * min(startPeriod,endPeriod)
        else:
            self.maxStartTime = maxStartTime
            
        if minFlightTime is None:
            self.minFlightTime = math.sqrt((startPeriod+endPeriod)**2)/8
        else:
            self.minFlightTime = minFlightTime
            
        if maxFlightTime is None:
            self.maxFlightTime = self.minFlightTime * 4
        else:
            self.maxFlightTime = maxFlightTime
            
        self.startTimes = np.linspace(self.minStartTime,                    \
                                      self.maxStartTime,                    \
                                      num = self.startTimeSize);
        self.flightTimes = np.linspace(self.minFlightTime,                  \
                                       self.maxFlightTime,                  \
                                       num = self.flightTimeSize)
        
        # The attributes defined here will be filled in with methods
        self.deltaV = None
        
        
        
        self.fill_table()
    
    
    def fill_table(self):
        """Calculates the delta v for each choice of start and flight time."""
        
        deltaVTable = np.zeros((self.flightTimeSize,self.startTimeSize))
        
        for xx, flightTime in enumerate(self.flightTimes, start=0):
            for yy, startTime in enumerate(self.startTimes, start=0):
                deltaVTable[xx][yy] =                                       \
                    self.get_chosen_transfer(startTime, flightTime)         \
                        .get_total_delta_V()
                
        self.deltaV = deltaVTable
    
    
    # def get_ejection_delta_v(self, vRel):
    #     ejectionInc = vRel
    
    
    def get_best_transfer(self):
        """Returns the transfer with the lowest delta V among sampled points."""
        
        minDV = np.nanmin(self.deltaV)
        index = np.where(self.deltaV == minDV)
        
        startTime = self.startTimes[index[1][0]]
        flightTime = self.flightTimes[index[0][0]]
        
        return self.get_chosen_transfer(startTime, flightTime)
    
    def get_chosen_transfer (self, startTime, flightTime):
        """Returns the transfer with the specified start and flight times.
        
        Arguments:
            startTime (float): time in seconds since epoch of transfer start
            flightTIme (float): time in seconds of transfer duration
            
        Returns:
            The transfer at with the specified start and flight times
        """
        
        if self.transferType == 'ballistic':
            return Transfer(self.startOrbit, self.endOrbit, startTime,      \
                            flightTime, False, self.ignoreInsertion,        \
                            self.fixedEndTime);
        
        elif self.transferType == 'plane change':
            return Transfer(self.startOrbit, self.endOrbit, startTime,      \
                            flightTime, True, self.ignoreInsertion,         \
                            self.fixedEndTime);
        
        elif self.transferType == 'optimal':
            btr = Transfer(self.startOrbit, self.endOrbit,                  \
                                   startTime, flightTime,                   \
                                   False, self.ignoreInsertion,             \
                                   self.fixedEndTime);
            ptr = Transfer(self.startOrbit, self.endOrbit,                  \
                                   startTime, flightTime,                   \
                                   True, self.ignoreInsertion,              \
                                   self.fixedEndTime);
            bdv = btr.get_total_delta_V()
            pdv = ptr.get_total_delta_V()
                    
            if bdv <= pdv:
                return btr
            else:
                return ptr
        
        else:
            raise Exception('uncrecognized transfer type')