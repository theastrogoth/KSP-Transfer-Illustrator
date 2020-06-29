# -*- coding: utf-8 -*-

# test script
import time
import pickle
import math
import numpy as np
from numpy.linalg import norm
from copy import copy
from orbit import Orbit
from body import Body
from transfer import Transfer
from prktable import PorkchopTable

infile = open('kerbol_system','rb')
kerbol_system = pickle.load(infile)
infile.close

startName = 'Moho'
endName = 'Kerbin'

startBody = [x for x in kerbol_system if x.name == startName][0]
endBody = [x for x in kerbol_system if x.name == endName][0]

trType = 'ballistic'
print('\n\n',trType,'transfer from',startBody.name,'to',endBody.name)

sOrb = Orbit(100000 + startBody.eqr, 0, 0, 0, 0, 0, startBody)
eOrb = Orbit(100000 + endBody.eqr, 0, 0, 0, 0, 0, endBody)

t0 = time.time()
testTable = PorkchopTable(sOrb,eOrb, transferType = trType)
t1 = time.time()

bestTransfer = testTable.get_best_transfer()
minDV = bestTransfer.get_total_delta_V()
departDay = math.floor(bestTransfer.get_departure_burn_time()/(3600*6))
arriveDay = math.floor(bestTransfer.get_arrival_burn_time()/(3600*6))
if bestTransfer.planeChange is True:
    planeChangeDay = math.floor(minDV/(3600*6))

print('\nMinimum delta v: ', bestTransfer.get_total_delta_V(), 'm/s')
print('Departure day:  day', departDay)
if bestTransfer.planeChange is True:
    print('Plane change day:  day', planeChangeDay)
print('Arrival day:  day', arriveDay)
print('\nTime elapsed during table generation: ', t1-t0, 's')

#%%
startBody.orb.get_basis_vectors()