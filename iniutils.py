"""Parser for KSPTOT bodies file"""
from math import pi
from orbit import Orbit
from body import Body
import re

def dicts_from_ini_file(ini, path=True):
    
    if path:
        data = open(ini, "r").read()
    else:
        data = ini
    
    # removes semicolon comments
    data = re.sub('(?m);.*', '', data)
    # removes all tabs and cursor marks
    data = data.replace("    ", "\t")   # remove large numbers of spaces
    data = data.replace("=", " =")      # ensure there is a space before each =
    data = data.replace("\t", "")       # remove tabs
    data = data.replace("\r", "")       # remove cursor mark
    # replace multiple newlines or spaces with singles
    while "\n\n" in data:
        data = data.replace("\n\n", "\n")
    while "  " in data:
        data = data.replace("  ", " ")
    # add newline at end
    if not data[-1] == '\n':
        data = data + '\n'
    
    dict_list = []
    # key_read contains the start and end index of the key being read
    key_read = [0, None]
    # value_read contains the start index of the value being read
    value_read = None
    trigger = set(("\n", "]", "="))
    for index, char in enumerate(data):
        # check if the char is one of the chars which leads to an action
        # this is an optimisation only
        if char in trigger:
            if char == "]":
                dict_list.append(dict())
                key_read = [index+2, None]
                value_read = None
            elif char == "=":
                key_read[1] = index - 1
                value_read = index + 2
            elif char == "\n":
                if not data[index - 1] == "\n":
                    if not value_read is None and not "]" in data[index-2:index]:
                        key = data[key_read[0]:key_read[1]]
                        value = data[value_read:index]
                        dict_list[-1][key] = value
                        key_read = [index + 1, None]
                else:
                    key_read = [index + 1, None]
            else:
                pass
        
    return dict_list


def dict_to_body(body_dict, system_list):
    
    name = body_dict['name']
    mu = float(body_dict['gm']) * 1E9
    eqr = float(body_dict['radius']) * 1000
    rotPeriod = float(body_dict['rotperiod'])
    rotIni = float(body_dict['rotini']) * pi/180
    
    a = float(body_dict['sma']) * 1000 
    ecc = float(body_dict['ecc'])
    inc = float(body_dict['inc']) * pi/180
    argp = float(body_dict['arg']) * pi/180
    lan = float(body_dict['raan']) * pi/180
    mo = float(body_dict['mean']) * pi/180
    epoch = float(body_dict['epoch'])
    
    primName = body_dict['parent']
    if primName == '':
        orb = None
    else:
        prim = [bd for bd in system_list if bd.name == primName][0]
        orb = Orbit(a, ecc, inc, argp, lan, mo, epoch, prim)
    
    try:
        colorStr = body_dict['color']
        color = eval(colorStr)
    except:
        color = (255, 255, 255)
    
    ref = int(body_dict['id'])
    
    return Body(name, eqr, mu, None, rotPeriod, rotIni, orb, ref, None, color)

def dicts_to_system(system_dicts):
    system = []
    idx = 0
    while len(system_dicts)>0:
        if idx >= len(system_dicts):
            idx = 0
        try:
            system.append(dict_to_body(system_dicts[idx], system))
            del system_dicts[idx]
        except IndexError:
            idx = idx + 1
    
    return system

def add_body_and_satellites(system, body):
    system.append(body)
    if not body.satellites is None:
        for sat in body.satellites:
            add_body_and_satellites(system, sat)

def sort_system(unsortedSystem):
    # sort Bodies by their orbits' SMAs
    sortedSystem = []
    for bd in unsortedSystem:
        bd.sort_satellites()
    
    try:
        sun = [bd for bd in unsortedSystem if bd.name == 'Sun'][0]
    except:
        sun = [bd for bd in unsortedSystem if bd.orb.prim.name == bd.name][0]
        
    add_body_and_satellites(sortedSystem, sun)
    return sortedSystem

def ini_to_system(ini, path=True):
    system_dicts = dicts_from_ini_file(ini, path)
    system = dicts_to_system(system_dicts)
    return sort_system(system)

def system_to_ini(system, path=None):
    if path is None:
        path = "bodies.ini"
    f = open(path, "w")
    for bd in system:
        f.write('['+bd.name+']\n')
        if bd.orb.prim.name == bd.name:
            f.write('epoch = 0\n')
            f.write('sma = 0\n')
            f.write('ecc = 0\n')
            f.write('inc = 0\n')
            f.write('raan = 0\n')
            f.write('arg = 0\n')
            f.write('mean = 0\n')
            f.write('gm = ' + str(bd.mu / 1E9) + '\n')
            f.write('radius = ' + str(bd.eqr / 1000) + '\n')
            f.write('rotperiod = ' + str(bd.rotPeriod) + '\n')
            f.write('rotini = ' + str(bd.rotIni*180/pi) + '\n')
            f.write('parent = \n')
            f.write('parentID = -1\n')
            f.write('name = ' + bd.name + '\n')
            f.write('id = ' + str(bd.ref) + '\n')
        else:
            f.write('epoch = ' + str(bd.orb.epoch)+'\n')
            f.write('sma = ' + str(bd.orb.a / 1000)+'\n')
            f.write('ecc = ' + str(bd.orb.ecc)+'\n')
            f.write('inc = ' + str(bd.orb.inc *180/pi)+'\n')
            f.write('raan = ' + str(bd.orb.lan *180/pi)+'\n')
            f.write('arg = ' + str(bd.orb.argp *180/pi)+'\n')
            f.write('mean = ' + str(bd.orb.mo *180/pi)+'\n')
            f.write('gm = ' + str(bd.mu / 1E9) + '\n')
            f.write('radius = ' + str(bd.eqr / 1000) + '\n')
            f.write('rotperiod = ' + str(bd.rotPeriod) + '\n')
            f.write('rotini = ' + str(bd.rotIni*180/pi) + '\n')
            f.write('parent = ' + bd.orb.prim.name + '\n')
            f.write('parentID = ' + str(bd.orb.prim.ref) + '\n' )
            f.write('name = ' + bd.name + '\n')
            f.write('id = ' + str(bd.ref) + '\n')
        
        f.write('color = ' + str(bd.color) + '\n')
        f.write('\n')
        
    f.close()
    
    return


