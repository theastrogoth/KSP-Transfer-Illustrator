# -*- coding: utf-8 -*-

import math
import numpy as np
from orbit import Orbit
from body import Body
from transfer import Transfer

class PorkchopTable:
    """
    
    """
    
    def __init__(self, startOrbit, endOrbit, transferType = 'ballistic',
                 ignoreInsertion = False, fixedEndTime = None, 
                 minStartTime = 0, maxStartTime = None, 
                 minFlightTime = None, maxFlightTime = None,
                 startTimeSize = 101, flightTimeSize = 101):
        """
        
        """
        
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
        """
        
        """
        deltaVTable = np.zeros((self.flightTimeSize,self.startTimeSize))
        
        for xx, flightTime in enumerate(self.flightTimes, start=0):
            for yy, startTime in enumerate(self.startTimes, start=0):
                deltaVTable[xx][yy] =                                       \
                    self.get_chosen_transfer(startTime, flightTime)         \
                        .get_total_delta_V()
                
        self.deltaV = deltaVTable
    
    
    def get_best_transfer(self):
        """
        
        """
        
        minDV = np.nanmin(self.deltaV)
        index = np.where(self.deltaV == minDV)
        
        startTime = self.startTimes[index[1][0]]
        flightTime = self.flightTimes[index[0][0]]
        
        return self.get_chosen_transfer(startTime, flightTime)
    
    def get_chosen_transfer (self, startTime, flightTime):
        """
        
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