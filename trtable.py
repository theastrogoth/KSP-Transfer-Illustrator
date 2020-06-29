# -*- coding: utf-8 -*-

import math
import numpy as np
from numpy import linalg as LA
from copy import copy
from orbit import Orbit
from body import Body
from transfer import Transfer

class d:
    """
    
    """
    
    def __init__(self, startOrbit, endOrbit, transferType = 'ballistic', 
                 startInclined = False, ignoreInsertion = False,
                 endTime = None, minStartTime = 0, maxStartTime = None, 
                 minFlightTime = None, maxFlightTime = None,
                 startTimeSize = 101, flightTimeSize = 101):
        """
        
        """
        
        self.startOrbit = startOrbit
        self.endOrbit = endOrbit
        self.transferType = transferType
        self.startInclined = startInclined
        self.ignoreInsertion = ignoreInsertion
        self.endTime = endTime
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
            self.maxStartTime = minStartTime + 2 * startPeriod
        else:
            self.maxStartTime = maxStartTime
            
        if minFlightTime is None:
            self.minFlightTime = math.sqrt((startPeriod+endPeriod)**2)/12
        else:
            self.minFlightTime = minFlightTime
            
        if maxFlightTime is None:
            self.maxFlightTime = self.minFlightTime * 4
        else:
            self.maxFlightTime = maxFlightTime
        
        # The attributes defined here will be filled in with methods
        self.deltaV = None
        
        self.fill_table()
    
    
    def fill_table(self):
        """
        
        """
        
        startTimes = np.linspace(self.minStartTime,                         \
                                 self.maxStartTime,                         \
                                 num = self.startTimeSize);
        flightTimes = np.linspace(self.minFlightTime,                       \
                                  self.maxFlightTime,                       \
                                  num = self.flightTimeSize)
        
        deltaVTable = np.zeros((self.startTimeSize,self.flightTimeSize))
        
        for xx, startTime in enumerate(startTimes, start=0):
            for yy, flightTime in enumerate(flightTimes, start=0):
                deltaVTable[xx][yy] =                                       \
                    self.get_chosen_transfer(startTime, flightTime)         \
                        .get_total_delta_V()
                
        self.deltaV = deltaVTable
    
    
    def get_best_transfer(self):
        minDV = np.nanmin(self.deltaV)
        index = np.where(self.deltaV == minDV)
        
        startTimes = np.linspace(self.minStartTime,                         \
                                 self.maxStartTime,                         \
                                 num = self.startTimeSize);
        flightTimes = np.linspace(self.minFlightTime,                       \
                                  self.maxFlightTime,                       \
                                  num = self.flightTimeSize);
            
        startTime = startTimes[index[0][0]]
        flightTime = flightTimes[index[1][0]]
        
        return self.get_chosen_transfer(startTime, flightTime)
    
    def get_chosen_transfer (self, startTime, flightTime):
        """
        
        """
        if self.transferType == 'ballistic':
            return Transfer(self.startOrbit, self.endOrbit, startTime,      \
                            flightTime, False, self.startInclined,          \
                            self.ignoreInsertion, self.endTime);
        
        elif self.transferType == 'plane change':
            return Transfer(self.startOrbit, self.endOrbit, startTime,      \
                            flightTime, True, self.startInclined,           \
                            self.ignoreInsertion, self.endTime);
        
        elif self.transferType == 'optimal':
            btr = Transfer(self.startOrbit, self.endOrbit,                  \
                                   startTime, flightTime,                   \
                                   False, self.startInclined,               \
                                   self.ignoreInsertion,                    \
                                   self.endTime);
            ptr = Transfer(self.startOrbit, self.endOrbit,                  \
                                   startTime, flightTime,                   \
                                   True, self.startInclined,                \
                                   self.ignoreInsertion,                    \
                                   self.endTime);
            bdv = btr.get_total_delta_V()
            pdv = ptr.get_total_delta_V()
                    
            if bdv <= pdv:
                return btr
            else:
                return ptr
        
        else:
            raise Exception('uncrecognized transfer type')