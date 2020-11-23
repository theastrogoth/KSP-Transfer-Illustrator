import math
from math import asin, acos, atan
from random import gauss
import numpy as np
from numpy import sin, cos, tan, dot, cross
from numpy.linalg import norm
from orbit import Orbit
from body import Body
from transfer import Transfer
from copy import copy
from scipy.optimize import dual_annealing, differential_evolution

def rodrigues_rotate(vec, axis, angle):
    
    axis = axis/norm(axis)
    return vec*cos(angle) + cross(axis, vec)*sin(angle) +                   \
           axis*dot(axis, vec)*(1-cos(angle))

def flyby_height_objective(periapsis, *params):
    
    aIn, aOut, deltaAngle, body = params
    
    eccIn = 1 - periapsis/aIn
    eccOut = 1 - periapsis/aOut
    
    inOrb = Orbit(aIn, eccIn, 0,0,0,0,0,body)
    outOrb = Orbit(aOut, eccOut, 0,0,0,0,0,body)
    
    deltaIn = inOrb.get_flight_angle_at_soi()
    deltaOut = outOrb.get_flight_angle_at_soi()
    
    return abs(deltaAngle - deltaIn - deltaOut)

def get_flyby_orbit(a, ecc, body, inOutDir, nDir, isOut=False):
    
    periapsis = a*(1-ecc)
    
    tempOrb = Orbit(a, ecc, 0,0,0,0,0,body)
    delta = tempOrb.get_flight_angle_at_soi(isOut)
    
    rpDir = rodrigues_rotate(inOutDir, nDir, delta)
    vpDir = cross(nDir, rpDir)
    
    rp = rpDir * periapsis
    vp = vpDir * math.sqrt(body.mu * (2/periapsis - 1/a))
    
    return Orbit.from_state_vector(rp, vp, 0, body)
    

def get_flyby_orbits(vIn, vOut, body):
    
    inDir = vIn/norm(vIn)
    outDir = vOut/norm(vOut)
    nDir = np.cross(inDir, outDir)
    nDir = nDir/norm(nDir)
    
    aIn = 1 / (2/body.soi - norm(vIn)**2/body.mu)
    aOut = 1 / (2/body.soi - norm(vOut)**2/body.mu)
    
    deltaAngle = math.acos(np.dot(inDir, outDir))
    params = (aIn, aOut, deltaAngle, body)
    
    res = dual_annealing(func=flyby_height_objective, args = params,        \
                         bounds=[(body.eqr, body.soi)], maxiter=100)
        
    periapsis = res.x
    eccIn = 1 - periapsis/aIn
    eccOut = 1- periapsis/aOut
    
    inOrb = get_flyby_orbit(aIn, eccIn, body, inDir, nDir, False)
    outOrb = get_flyby_orbit(aOut, eccOut, body, outDir, nDir, True)
    
    return inOrb, outOrb

def flyby_delta_v_objective(t, *params):
    
    startOrbit, flybyBody, endOrbit, ignoreInsertion, refine = params
    
    transfer1, inOrb, outOrb, transfer2 =                                   \
        calculate_flyby(t, startOrbit, flybyBody, endOrbit, ignoreInsertion, 
                        refine)
    
    vPeriIn = inOrb.get_state_vector(0)[1]
    vPeriOut = outOrb.get_state_vector(0)[1]
    flybyDV = vPeriOut - vPeriIn
    
    totalDV = norm(transfer1.ejectionDV) + norm(flybyDV) + norm(transfer2.insertionDV)
    return totalDV

def calculate_flyby(t, *params):
    
    startTime = t[0]
    flightTime1 = t[1]
    flightTime2 = t[2]
    
    startOrbit, flybyBody, endOrbit, ignoreInsertion, refine = params
    
    transfer1 = Transfer(startOrbit, flybyBody.orb, startTime, flightTime1,
                         cheapEndOrb=False)
    transfer2 = Transfer(flybyBody.orb, endOrbit, startTime + flightTime1, 
                         flightTime2, cheapEndOrb=False, 
                         ignoreInsertion=ignoreInsertion)
    
    vIn = -transfer1.insertionDV
    vOut = transfer2.ejectionDV
    
    inOrb, outOrb = get_flyby_orbits(vIn, vOut, flybyBody)
    
    if refine:
        count = 0
        while count < 10:
            inNu = -inOrb.get_nu_at_soi()
            inTime = inOrb.get_time(inNu)
            outNu = outOrb.get_nu_at_soi()
            outTime = outOrb.get_time(outNu)
            
            inPos = flybyBody.orb.get_state_vector(startTime+flightTime1)[0] +      \
                    inOrb.get_state_vector(inTime)[0]
            outPos = flybyBody.orb.get_state_vector(startTime+flightTime1+outTime-inTime)[0] +      \
                    outOrb.get_state_vector(outTime)[0]
            
            transfer1 = Transfer(startOrbit, flybyBody.orb, startTime, flightTime1,
                         cheapEndOrb=False, endPos = inPos)
            transfer2 = Transfer(flybyBody.orb, endOrbit, startTime + flightTime1, 
                         flightTime2, cheapEndOrb=False, startPos = outPos,
                         ignoreInsertion=ignoreInsertion)
            transfer1.match_start_mean_anomaly()
            transfer2.genetic_refine()
            
            vIn = -transfer1.insertionDV
            vOut = transfer2.ejectionDV
            inOrb, outOrb = get_flyby_orbits(vIn, vOut, flybyBody)
            
            count = count + 1
        
        inOrb.epoch = startTime+flightTime1-inTime
        inOrb.epoch = startTime+flightTime1-inTime
    
    return transfer1, inOrb, outOrb, transfer2

class Flyby:
    
    def __init__(self, startOrbit, flybyBody, endOrbit,
                 minStartTime=None, maxStartTime=None,
                 minFT1=None, maxFT1=None, minFT2=None, maxFT2=None,
                 ignoreInsertion = True):
        
        if startOrbit.prim.is_system_center():
            startPeriod = startOrbit.get_period()
        else:
            startPeriod = startOrbit.prim.orb.get_period()
        flybyPeriod = flybyBody.orb.get_period()
        if endOrbit.prim.is_system_center():
            endPeriod = endOrbit.get_period()
        else:
            endPeriod = endOrbit.prim.orb.get_period()
        
        self.startOrbit = startOrbit
        self.endOrbit = endOrbit
        self.flybyBody = flybyBody
        if minStartTime is None:
            self.minStartTime = 0
        else:
            self.minStartTime = minStartTime
        if maxStartTime is None:
            self.maxStartTime = self.minStartTime  +                        \
                startPeriod * 10;
        else:
            self.maxStartTime = maxStartTime
            
        midPeriod1 = ((startPeriod**(2/3) + flybyPeriod**(2/3))/2)**(3/2)
        midPeriod2 = ((flybyPeriod**(2/3) + endPeriod**(2/3))/2)**(3/2)
        
        if minFT1 is None:
            self.minFT1 = midPeriod1/25
        else:
            self.minFT1 = minFT1
        if maxFT1 is None:
            self.maxFT1 = midPeriod1*5
        else:
            self.maxFT1 = maxFT1
        if minFT2 is None:
            self.minFT2 = midPeriod2/25
        else:
            self.minFT2 = minFT2
        if maxFT2 is None:
            self.maxFT2 = midPeriod2*5
        else:
            self.maxFT2 = maxFT2
            
        self.ignoreInsertion = ignoreInsertion
        
        # filled in by methods
        self.transfer1 = None
        self.transfer2 = None
        self.inOrb = None
        self.outOrb = None
        self.flybyDV = None
        self.bestTimes = None
        
        # first pass at optimization
        self.optimize_flyby()
    
    def optimize_flyby(self, guess = None):
        
        if guess is None:
            guess = self.bestTimes
        
        params = (self.startOrbit, self.flybyBody, self.endOrbit, 
                  self.ignoreInsertion, False)
        
        bounds = [(self.minStartTime, self.maxStartTime),                   \
                  (self.minFT1, self.maxFT1),                               \
                  (self.minFT2, self.maxFT2)];
        
        if guess is None:
            res = differential_evolution(func=flyby_delta_v_objective,
                                         args = params,
                                         bounds = bounds, disp=True, popsize=4, 
                                         strategy='best2bin', maxiter = 20,
                                         mutation=(0.5, 1), recombination=0.7,
                                         polish=True)
        else:
            init = [guess]
            for ii in range(11):
                nextGuess = []
                for jj, param in enumerate(guess):
                    boundSize = bounds[jj][1]-bounds[jj][0]
                    newParam = param + boundSize*gauss(0,0.5)
                    nextGuess.append(newParam)
                init.append(nextGuess)
            
            res = differential_evolution(func=flyby_delta_v_objective,
                                 args = params,
                                 bounds = bounds, disp=True, popsize=4, 
                                 strategy='best2bin', maxiter = 20,
                                 mutation=(0.5, 1), recombination=0.7,
                                 polish=True, init = init)
        
        transfer1, inOrb, outOrb, transfer2 =                               \
            calculate_flyby(res.x, self.startOrbit,                         \
                            self.flybyBody, self.endOrbit, 
                            self.ignoreInsertion, True);
        
        vPeriIn = inOrb.get_state_vector(0)[1]
        vPeriOut = outOrb.get_state_vector(0)[1]
        flybyDV = vPeriOut - vPeriIn
        
        self.transfer1 = transfer1
        self.inOrb = inOrb
        self.outOrb = outOrb
        self.transfer2 = transfer2
        self.flybyDV = flybyDV
        self.bestTimes = res.x
        
        return

# #%%
# # testing

# import jsonpickle
# infile = open('kerbol_system.json','r')
# kerbol_system = jsonpickle.decode(infile.read())
# infile.close

# startBody = kerbol_system[4]
# flybyBody = kerbol_system[2]
# endBody = kerbol_system[7]

# startOrbit = Orbit(startBody.eqr+100000, 0,0,0,0,0,0, startBody)
# flybyOrbit = Orbit(flybyBody.eqr+100000, 0,0,0,0,0,0, flybyBody)
# endOrbit = Orbit(endBody.eqr+100000, 0,0,0,0,0,0, endBody)

# testFlyby = Flyby(startOrbit, flybyBody, endOrbit)
# testFlyby.optimize_flyby()

# #%%
# t = np.array([1.72605416e+08, 5.25832607e+06, 6.82303613e+06])
# testFlyby.optimize_flyby(t)