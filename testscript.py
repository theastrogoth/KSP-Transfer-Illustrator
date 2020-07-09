# -*- coding: utf-8 -*-

# test script
import time
import jsonpickle
import math
import numpy as np
from numpy.linalg import norm
from copy import copy
from orbit import Orbit
from body import Body
from transfer import Transfer
from prktable import PorkchopTable

infile = open('kerbol_system.json','r')
kerbol_system = jsonpickle.decode(infile.read())
infile.close

startName = 'Jool'
endName = 'Kerbin'

startBody = [x for x in kerbol_system if x.name == startName][0]
endBody = [x for x in kerbol_system if x.name == endName][0]

trType = 'ballistic'
print('\n\n',trType,'transfer from',startBody.name,'to',endBody.name)

sOrb = Orbit(100000 + startBody.eqr, 0, 0, 0, 0, 0, startBody)
eOrb = Orbit(100000 + endBody.eqr, 0, 0, 0, 0, 0, endBody)

t0 = time.time()
porkTable = PorkchopTable(sOrb,eOrb, trType, False, None, 0, None, None, None, 31, 31)
t1 = time.time()

# #%%
# bestTransfer = porkTable.get_best_transfer()
# minDV = bestTransfer.get_total_delta_V()
# departDay = math.floor(bestTransfer.get_departure_burn_time()/(3600*6))
# arriveDay = math.floor(bestTransfer.get_arrival_burn_time()/(3600*6))
# if bestTransfer.planeChange is True:
#     planeChangeDay = math.floor((bestTransfer.startTime +                    \
#                                 bestTransfer.planeChangeDT)/(3600*6));

# print('\nMinimum delta v: ', bestTransfer.get_total_delta_V(), 'm/s')
# print('Departure day:  day', departDay)
# if bestTransfer.planeChange is True:
#     print('Plane change day:  day', planeChangeDay)
# print('Arrival day:  day', arriveDay)
# print('\nTime elapsed during table generation: ', t1-t0, 's')

# #%%
# t0 = t0 = time.time()
# for xx in range(1000):
#     testOrbit = Orbit.from_state_vector(np.array([-4.85410069e+09,2.26952608e+09,-8.12986030e+07]),
#                                         np.array([889.74802429,-416.47936992,14.89479467]),
#                                         5432790.187884141,
#                                         startBody)
# t1 = time.time()
# print('\nAverage time elapsed during orbit calculation: ', (t1-t0)/1000, 's')

# #%%
# t0 = t0 = time.time()
# for xx in range(1000):
#     testPos, testVel = bestTransfer.startOrbit.get_state_vector(5432790.187884141)
# t1 = time.time()
# print('\nAverage time elapsed during state vector calculation: ', (t1-t0)/1000, 's')

# #%%
# t0 = t0 = time.time()
# for xx in range(1000):
#     testPos = bestTransfer.startOrbit.from_primary_to_orbit_bases(np.array([1,1,1]))
# t1 = time.time()
# print('\nAverage time elapsed during coordinate change calculation: ', (t1-t0)/1000, 's')


# #%%
# t0 = t0 = time.time()
# for xx in range(1000):
#     testTransfer = porkTable.get_best_transfer()
# t1 = time.time()
# print('\nAverage time elapsed during transfer calculation: ', (t1-t0)/1000, 's')
# t0 = t0 = time.time()
# for xx in range(1000):
#     testTransfer.get_ejection_details()
# t1 = time.time()
# print('\nAverage time elapsed during ejection calculation: ', (t1-t0)/1000, 's')