import math
import numpy as np
from numpy.linalg import norm
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
        cheapStartOrb (bool): if true, the only parameter of the starting
            park orbit used is the semimajor axis
        cheapEndOrb (bool): if true, the only parameter of the ending park
            orbit used is the semimajor axis
        minStartTime (float): earliest time at the beginning of transfer 
            trajectory (s)
        minStartTime (float): latest time at the beginning of transfer 
            trajectory (s)
        minFlightTime (float): shortest duration of transfer trajectory (s)
        maxFlightTime (float): longest duration of transfer trajectory (s)
        startTimeSize (int): number of samples taken on the start time axis
        flightTimeSize (int): number of samples taken on the flight time axis
        startTimes (floats): list holding all start times sampled (s)
        flightTimes (floats): list hold all flight times sampled (s)
        deltaV: a table of values with the sum of the magnitue of all burn 
            maneuvers (m/s) at each choice of start and flight times
    
    """
    
    def __init__(self, startOrbit, endOrbit, transferType = 'ballistic',
                 ignoreInsertion = False,
                 cheapStartOrb = False, cheapEndOrb = True,
                 minStartTime = 0, maxStartTime = None, 
                 minFlightTime = None, maxFlightTime = None,
                 startTimeSize = 25, flightTimeSize = 25):
        
        self.startOrbit = startOrbit
        self.endOrbit = endOrbit
        self.transferType = transferType
        self.ignoreInsertion = ignoreInsertion
        self.cheapStartOrb = cheapStartOrb
        self.cheapEndOrb = cheapEndOrb
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
        if self.minFlightTime > self.maxFlightTime:
            self.minFlightTime = self.maxFlightTime/4
            
        self.startTimes = np.linspace(self.minStartTime,                    \
                                      self.maxStartTime,                    \
                                      num = self.startTimeSize);
        self.flightTimes = np.linspace(self.minFlightTime,                  \
                                       self.maxFlightTime,                  \
                                       num = self.flightTimeSize)
        
        # The attributes defined here will be filled in with methods
        self.totalDeltaV = None
        self.ejectionDeltaV = None
        self.insertionDeltaV = None
        
        # Fill in the empty attributes
        self.fill_table()
    
    
    def fill_table(self):
        """Calculates the delta v for each choice of start and flight time."""
        
        # f = open('gens.csv','w')
        
        totalDeltaVTable = np.zeros((self.flightTimeSize,self.startTimeSize))
        ejectDeltaVTable = np.zeros((self.flightTimeSize,self.startTimeSize))
        insertDeltaVTable = np.zeros((self.flightTimeSize,self.startTimeSize))
        
        
        for xx, flightTime in enumerate(self.flightTimes, start=0):
            for yy, startTime in enumerate(self.startTimes, start=0):
                # trs, gen = self.get_chosen_transfer(startTime, flightTime)
                # f.write(str(startTime))
                # f.write(',')
                # f.write(str(flightTime))
                # f.write(',')
                # f.write(str(gen))
                # f.write(',')
                # f.write(str(trs.convergenceFail))
                # f.write(',')
                # f.write('\n')
                trs = self.get_chosen_transfer(startTime, flightTime)
                totalDeltaVTable[xx][yy] = trs.get_total_delta_v()
                ejectDeltaVTable[xx][yy] = norm(trs.ejectionDV)
                insertDeltaVTable[xx][yy] = norm(trs.insertionDV)
        # f.close()
        self.totalDeltaV = totalDeltaVTable
        self.ejectionDeltaV = ejectDeltaVTable
        self.insertionDeltaV = insertDeltaVTable
    
    
    def get_best_transfer(self):
        """Returns the transfer with the lowest delta V among sampled points."""
        
        minDV = np.nanmin(self.totalDeltaV)
        index = np.where(self.totalDeltaV == minDV)
        
        startTime = self.startTimes[index[1][0]]
        flightTime = self.flightTimes[index[0][0]]
        
        return self.get_chosen_transfer(startTime, flightTime) # [0]
    
    def get_chosen_transfer (self, startTime, flightTime):
        """Returns the transfer with the specified start and flight times.
        
        Arguments:
            startTime (float): time in seconds since epoch of transfer start
            flightTIme (float): time in seconds of transfer duration
            
        Returns:
            The transfer at with the specified start and flight times
        """
        
        if self.transferType == 'ballistic':
            trs = Transfer(self.startOrbit, self.endOrbit, startTime,       \
                            flightTime, False, self.ignoreInsertion,        \
                            self.cheapStartOrb, self.cheapEndOrb);
        
        elif self.transferType == 'plane change':
            trs = Transfer(self.startOrbit, self.endOrbit, startTime,       \
                            flightTime, True, self.ignoreInsertion,         \
                            self.cheapStartOrb, self.cheapEndOrb);
        
        elif self.transferType == 'optimal':
            btr = Transfer(self.startOrbit, self.endOrbit,                  \
                                   startTime, flightTime,                   \
                                   False, self.ignoreInsertion,             \
                                   self.cheapStartOrb, self.cheapEndOrb);
            ptr = Transfer(self.startOrbit, self.endOrbit,                  \
                                   startTime, flightTime,                   \
                                   True, self.ignoreInsertion,              \
                                   self.cheapStartOrb, self.cheapEndOrb);
            bdv = btr.get_total_delta_v()
            pdv = ptr.get_total_delta_v()
                    
            if bdv <= pdv:
                trs = btr
            else:
                trs = ptr
        
        else:
            raise Exception('uncrecognized transfer type')
        
        # gen = trs.genetic_refine()
        return trs # , gen
