# Porkchop and 3D orbit plotting utility functions

import plotly.graph_objects as go
import jsonpickle
import math
import numpy as np
from numpy.linalg import norm
from orbit import Orbit
from body import Body
from transfer import Transfer
from prktable import PorkchopTable
from imageutils import image_colormap, map_url, get_pixel_values

#%% misc functions

def seconds_to_dates(times, dateFormat):
    
    times = np.array(times)
    
    day = dateFormat['day']     # hours per day
    year = dateFormat['year']   # days per year
    
    ydhms = np.stack((np.floor(times/(3600*day*year))+1,                    \
                      np.floor(times%(3600*day*year)/(day*3600)+1),         \
                      np.floor((times%(3600*day))/3600),                    \
                      np.floor(times%3600/60),                              \
                      times%60),                                            \
                     axis=1);
    return ydhms 

def dates_to_seconds(ydhms, dateFormat):
    
    ydhms = np.array(ydhms)
    
    day = dateFormat['day']
    year = dateFormat['year']
    
    ydhms = np.transpose(ydhms)
    
    seconds = ((((ydhms[0]-1) * year +                                      \
                    (ydhms[1]-1)) * day +                                   \
                        ydhms[2]) * 60 +                                    \
                           ydhms[3]) * 60 +                                 \
                                ydhms[4]
    
    seconds = np.transpose(seconds)
    return seconds

def seconds_to_date_string(time, dateFormat):
    
    date = seconds_to_dates([time], dateFormat)[0]
    
    dateString = 'Year ' + "{:d}".format(int(date[0])) +                    \
                 ', Day ' + "{:d}".format(int(date[1]))+ ', '               \
                 "{:02d}".format(int(date[2])) + ':' +                      \
                 "{:02d}".format(int(date[3])) + ':' +                      \
                 "{:02d}".format(int(date[4]));
    return dateString

def seconds_to_days(time, dateFormat):
    
    day = dateFormat['day']
    return time/day/3600

def get_burn_components(burnDV, position, velocity, normalDir = None):
    
    if normalDir is None:
        normalDir = np.cross(position, velocity)
        normalDir = normalDir/norm(normalDir)
    
    progradeDir = velocity/norm(velocity)
    radialDir = np.cross(progradeDir, normalDir)
    radialDir = radialDir/norm(radialDir)
    
    prograde = np.dot(burnDV, progradeDir)
    normal =   np.dot(burnDV, normalDir)
    radial =   np.dot(burnDV, radialDir)
    
    return prograde, normal, radial

def burn_components_string(burnDV, position, velocity, normalDir = None):
    
    prograde, normal, radial =                                              \
        get_burn_components(burnDV, position, velocity, normalDir);
    
    threshold = 0.05
    
    burnString = ''
    if prograde > threshold:
        burnString = burnString +                                           \
            "{:.2f}".format(prograde) + ' m/s prograde';
    elif prograde < -threshold:
        burnString = burnString +                                           \
            "{:.2f}".format(-prograde) + ' m/s retrograde';
    if normal > threshold:
        burnString = burnString + ', ' +                                \
            "{:.2f}".format(normal) + ' m/s normal';
    if normal < -threshold:
        burnString = burnString + ', ' +                                \
            "{:.2f}".format(-normal) + ' m/s anti-normal';
    if radial > threshold:
        burnString = burnString + ', ' +                                \
            "{:.2f}".format(radial) + ' m/s radial';
    elif radial < -threshold:
        burnString = burnString + ', ' +                                \
            "{:.2f}".format(-radial) + ' m/s anti-radial';
    
    return burnString

def burn_components_to_absolute(prograde, normal, radial, position, velocity,
                                normalDir = None):
    if normalDir is None:
        normalDir = np.cross(position, velocity)
        normalDir = normalDir/norm(normalDir)
    
    progradeDir = velocity/norm(velocity)
    radialDir = np.cross(progradeDir, normalDir)
    radialDir = radialDir/norm(radialDir)
    
    burnDV = prograde*progradeDir + normal*normalDir + radial*radialDir
    return burnDV

def cartesian_to_spherical(xyz):
    spherical = xyz*0
    xy = xyz[:,0]**2 + xyz[:,1]**2
    spherical[:,0] = np.sqrt(xy + xyz[:,2]**2)
    spherical[:,1] = np.arctan2(xyz[:,1], xyz[:,0])
    spherical[:,2] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    return spherical

def spherical_to_cartesian(spherical):
    xyz = spherical*0
    xyz[:,0] = spherical[:,0] * np.cos(spherical[:,1]) * np.sin(spherical[:,2])
    xyz[:,1] = spherical[:,0] * np.sin(spherical[:,1]) * np.sin(spherical[:,2])
    xyz[:,2] = spherical[:,0] * np.cos(spherical[:,2])
    return xyz

def lat_lon_to_spherical(lon, lat, radius):
    
    # xs = radius * np.cos(lon) * np.transpose(np.cos(lat))
    # ys = radius * np.sin(lon) * np.transpose(np.cos(lat))
    # zs = radius * lon/lon * np.transpose(np.sin(lat))
    # return np.reshape(xs, xs.size), np.reshape(ys, ys.size), np.reshape(zs, zs.size)
    
    lat, lon = np.meshgrid(lat, lon)
    
    xs = radius * np.cos(lon) * np.cos(lat)
    ys = radius * np.sin(lon) * np.cos(lat)
    zs = radius * np.sin(lat)
    
    return xs, ys, zs


#%% porkchop plot functions

def add_lines(figure, x, y, minX, maxX, minY, maxY, color = 'black'):
    """ Adds lines indicating a point at (x,y) """
    
    figure.add_trace(go.Scatter(
        name = 'tr-x',
        x = [x,x],
        y = [minY, maxY],
        mode = 'lines',
        line = dict(
            color = color
            ),
        hoverinfo = "skip",
        showlegend = False,
        ))
    
    figure.add_trace(go.Scatter(
        name = 'tr-y',
        x = [minX, maxX],
        y = [y, y],
        mode = 'lines',
        line = dict(
            color = color
            ),
        hoverinfo = "skip",
        showlegend = False,
        ))

def add_marker(figure, x, y, symbol = 'x', color = 'black', size = 12):
    """ Adds a marker indicating a point at (x,y) """
    
    
    
    figure.add_trace(go.Scatter(
        name = 'tr-x',
        x = [x],
        y = [y],
        mode = 'markers',
        marker = dict(
            color = color,
            symbol = symbol,
            size = size
            ),
        hoverinfo = "skip",
        showlegend = False,
        ))

#%% trajectory plot functions

def fade_color(color, div = 2):
    """Divides each element of the tuple by the specified number."""
    
    return tuple(math.floor(c/div) for c in color)

def add_orbit(figure, orb, startTime, endTime=None, numPts=201,
              dateFormat=None, apses=False, nodes=False, fullPeriod=True,
              color=(255,255,255), name='', style='solid', fade=True,):
    
    if fade:
        fadedColor = fade_color(color,3)
    else:
        fadedColor = color
    
    period = orb.get_period()
    if fullPeriod and (endTime is None):
        endTime = startTime + period
    
    # start and end mean anomalies
    mStart = orb.get_mean_anomaly(startTime)
    if (period < (endTime-startTime) or fullPeriod) and (orb.ecc < 1):
        mEnd = mStart + 2*math.pi
    else:
        mEnd = mStart + 2*math.pi/period * (endTime-startTime)
    
    # get points clustered around apoapsis and periapsis
    if orb.ecc < 1:
        a = mStart - (mStart%math.pi)
        b = a + math.pi
    else:
        if mStart < 0:
            a = mStart
            b = 0
        else:
            a = 0
            b = mEnd
    # orbit crosses two apo/peri-apses
    if ((mEnd >= b + math.pi) and (orb.ecc < 1)):
        c = b + math.pi
        d = c + math.pi
        n = math.ceil(math.pi/(mEnd-mStart)*numPts)
        kStart = n*math.acos((a+b-2*mStart)/(b-a))/math.pi
        kEnd =   n*math.acos((c+d-2*mEnd)/(d-c))/math.pi
        ks1 = [[kStart]]
        ks1.append([*range(math.ceil(kStart), n)])
        ks1 = [k for sublist in ks1 for k in sublist]
        meanAnoms1 =                                                        \
            [0.5*(a+b) + 0.5*(b-a) *math.cos((n-k)/n*math.pi) for k in ks1];
        ks2 = [*range(0,n)]
        meanAnoms2 =                                                        \
            [0.5*(b+c) + 0.5*(c-b) *math.cos((n-k)/n*math.pi) for k in ks2];
        ks3 = [*range(0,math.ceil(kEnd))]
        ks3.append(kEnd)
        meanAnoms3 =                                                        \
            [0.5*(c+d) + 0.5*(d-c) *math.cos((n-k)/n*math.pi) for k in ks3];
        meanAnoms = np.append(meanAnoms1, meanAnoms2)
        meanAnoms = np.append(meanAnoms, meanAnoms3)
        times = startTime + period/(2*math.pi) * (meanAnoms - mStart)
    # orbit crosses one apo/peri-apsis
    elif mEnd > b:
        if orb.ecc < 1:
            c = b + math.pi
            n1 = math.ceil(math.pi/(mEnd-mStart)*numPts)
            n2 = n1
        else:
            c = mEnd
            n1 = math.ceil(abs(mStart/(mEnd-mStart))*numPts)
            n2 = math.ceil(abs(mEnd/(mEnd-mStart))*numPts)
        kStart = n1*math.acos((a+b-2*mStart)/(b-a))/math.pi
        kEnd =   n2*math.acos((b+c-2*mEnd)/(c-b))/math.pi
        ks1 = [[kStart]]
        ks1.append([*range(math.ceil(kStart), n1)])
        ks1 = [k for sublist in ks1 for k in sublist]
        meanAnoms1 =                                                        \
            [0.5*(a+b) + 0.5*(b-a) *math.cos((n1-k)/n1*math.pi) for k in ks1];
        ks2 = [*range(0,math.ceil(kEnd))]
        ks2.append(kEnd)
        meanAnoms2 =                                                        \
            [0.5*(b+c) + 0.5*(c-b) *math.cos((n2-k)/n2*math.pi) for k in ks2];
        meanAnoms = np.append(meanAnoms1, meanAnoms2)
        times = startTime + period/(2*math.pi) * (meanAnoms - mStart)
    # orbit crosses no apo/peri-apses
    else:
        if orb.ecc < 1:
            n = math.ceil(2*math.pi/(mEnd-mStart)*numPts)
        else:
            n = numPts
        kStart = n*math.acos((a+b-2*mStart)/(b-a))/math.pi
        kEnd = n*math.acos((a+b-2*mEnd)/(b-a))/math.pi
        ks = [[kStart]]
        ks.append([*range(math.ceil(kStart), math.ceil(kEnd))])
        ks.append([kEnd])
        ks = [k for sublist in ks for k in sublist]
        meanAnoms = np.array(                                               \
            [0.5*(a+b) + 0.5*(b-a) *math.cos((n-k)/n*math.pi) for k in ks]);
        meanAnoms = meanAnoms.flatten()
        times = startTime + period/(2*math.pi) * (meanAnoms - mStart)
    
    pos, vel = orb.get_positions(times = times)
    pos = np.transpose(pos)
    vel = np.transpose(vel)
    
    if orb.ecc<1:
        for ii, m in enumerate(meanAnoms):
            meanAnoms[ii] = Orbit.map_angle(m)
    
    if not dateFormat is None:
        day = dateFormat['day']
        year = dateFormat['year']
        
        cData = np.stack((norm(pos, axis = 0)/1000,
                          norm(vel, axis = 0),
                          np.floor(times/(3600*day*year))+1,
                          np.floor(times%(3600*day*year)/(day*3600)+1),
                          np.floor((times%(3600*day))/3600),
                          np.floor(((times%(3600*day))%3600)/60),
                          np.floor(((times%(3600*day))%3600)%60),
                          times,
                          meanAnoms),
                          axis=1);
        hoverLabel = "r = %{customdata[0]:.3e} km" + "<br>" +\
                     "v = %{customdata[1]:.3e} m/s" + "<br>" + "<br>" +\
                     "Year %{customdata[2]:.0f}, " +\
                     "Day %{customdata[3]:.0f} " +\
                     "%{customdata[4]:0>2d}" + ":" +\
                     "%{customdata[5]:0>2d}" + ":" +\
                     "%{customdata[6]:0>2d}" + "<br>" +\
                     "UT: %{customdata[7]:.3f} s" + "<br>" +\
                     "Mean Anomaly: %{customdata[8]:.5f} rad" + "<br>" + "<br>" +\
                     "Semi-major Axis = " + "{:.0f}".format(orb.a) + " m" + "<br>" +\
                     "Eccentricity = " + "{:.4f}".format(orb.ecc) + "<br>" +\
                     "Inclination = " + "{:.4f}".format(orb.inc*180/math.pi) + "°" + "<br>" +\
                     "Argument of the Periapsis = " + "{:.4f}".format(orb.argp*180/math.pi) + "°" + "<br>" +\
                     "Longitude of Ascending Node = " + "{:.4f}".format(orb.lan*180/math.pi) + "°" + "<br>" +\
                     "Mean Anomaly at Epoch = " + "{:.4f}".format(orb.mo) + " rad" + "<br>" +\
                     "Epoch = " + "{:.2f}".format(orb.epoch) + " s"
    else:
        cData = norm(pos, axis = 0)/1000
        hoverLabel = "r = %{customdata:.4e} km"
    
    figure.add_trace(go.Scatter3d(
        x = pos[0],
        y = pos[1],
        z = pos[2],
        customdata = cData,
        mode = "lines",
        line = dict(
            color = times,
            colorscale = [[0, 'rgb'+str(fadedColor)], 
                      [1, 'rgb'+str(color)]],
            dash = style
            ),
        hovertemplate = hoverLabel,
        name = name,
        showlegend = False,
        ))
    
    if nodes:
        add_nodes(figure, orb)
    if apses:
        add_apses(figure, orb)

def add_apses(figure, orb, size = 4, color = (0,0,255)):
    
    if orb.ecc == 0:
        return
    
    periPos = orb.get_state_vector(orb.get_time(0))[0]
    
    if orb.ecc < 1:
        apoPos  = orb.get_state_vector(orb.get_time(math.pi))[0]
        x = [apoPos[0],periPos[0]]
        y = [apoPos[1],periPos[1]]
        z = [apoPos[2],periPos[2]]
    else:
        x = [periPos[0]]
        y = [periPos[1]]
        z = [periPos[2]]
    
    figure.add_trace(go.Scatter3d(
                                  x = x,
                                  y = y,
                                  z = z,
                                  mode = "markers",
                                  marker = dict(
                                      color = 'rgb'+str(color),
                                      symbol = 'circle',
                                      size = size),
                                  showlegend = False,
                                  hoverinfo = 'skip',
                                  ))

def add_nodes(figure, orb, size = 4, color = (0,255,0)):
    
    if (orb.inc == 0) or (orb.ecc > 1):
        return
    
    ascTrueAnom  = -orb.argp
    ascPos  = orb.get_state_vector(orb.get_time(ascTrueAnom))[0]
    
    # if orb.ecc < 1:
    descTrueAnom = ascTrueAnom + math.pi
    descPos = orb.get_state_vector(orb.get_time(descTrueAnom))[0]
    x = [ascPos[0],descPos[0]]
    y = [ascPos[1],descPos[1]]
    z = [ascPos[2],descPos[2]]
    # else:
    #     x = [ascPos[0]]
    #     y = [ascPos[1]]
    #     z = [ascPos[2]]
    
    figure.add_trace(go.Scatter3d(
                                  x = x,
                                  y = y,
                                  z = z,
                                  mode = "markers",
                                  marker = dict(
                                      color = 'rgb'+str(color),
                                      symbol = 'diamond',
                                      size = size),
                                  showlegend = False,
                                  hoverinfo = 'skip',
                                  ))

def add_primary(figure, bd, surf = True, lat = None, lon = None, pix = None):
    
    fadedColor = fade_color(bd.color)
    
    size = 10
    phi = np.linspace(0, 2*math.pi, size)
    theta = np.linspace(-math.pi/2, math.pi/2, size)
    phi, theta = np.meshgrid(phi, theta)
    
    x = bd.eqr * np.cos(theta) * np.sin(phi)
    y = bd.eqr * np.cos(theta) * np.cos(phi)
    z = bd.eqr * np.sin(theta)
    
    if not lat is None:
        xs, ys, zs = lat_lon_to_spherical(lon, lat, bd.eqr)
        cmap, mapDict, pix = image_colormap(pix, rounded=True)
        val = []
        for pick in pix:
            val.append(mapDict[pick])
        # val.reverse()
        
        figure.add_trace(go.Surface(
                                    x=xs,
                                    y=ys,
                                    z=-zs,
                                    colorscale=cmap,
                                    surfacecolor=np.transpose(np.reshape(val,xs.shape)),
                                    showscale=False,
                                    name = bd.name,
                                    showlegend = False,
                                    hovertemplate = "Central body"
                                    ))
    elif surf:
        figure.add_trace(go.Mesh3d(
                                    x = np.ndarray.flatten(x),
                                    y = np.ndarray.flatten(y),
                                    z = np.ndarray.flatten(z),
                                    alphahull = 0,
                                    color = 'rgb'+str(fadedColor),
                                    opacity = 0.5,
                                    name = bd.name,
                                    showlegend = False,
                                    hovertemplate = "Central body"
                                    # hoverinfo = 'skip'
                                        ))
    else:
        figure.add_trace(go.Scatter3d(
                                      x = np.array([0]),
                                      y = np.array([0]),
                                      z = np.array([0]),
                                      mode = "markers",
                                      marker = dict(
                                          color = 'rgb'+str(bd.color),
                                          symbol = "circle"),
                                      name = bd.name,
                                      showlegend = False,
                                      hovertemplate = "Central body"
                                      ))

def add_body(figure, bd, time, surf=True, pos=None, size=8, symbol='circle'):
    
    if pos is None:
        pos = bd.orb.get_state_vector(time)[0]
    fadedColor = fade_color(bd.color)
    
    if surf:
        size = 10
        phi = np.linspace(0, 2*math.pi, size)
        theta = np.linspace(-math.pi/2, math.pi/2, size)
        phi, theta = np.meshgrid(phi, theta)
        
        x = bd.eqr * np.cos(theta) * np.sin(phi) + pos[0]
        y = bd.eqr * np.cos(theta) * np.cos(phi) + pos[1]
        z = bd.eqr * np.sin(theta) + pos[2]
        
        figure.add_trace(go.Mesh3d(
                                    x = np.ndarray.flatten(x),
                                    y = np.ndarray.flatten(y),
                                    z = np.ndarray.flatten(z),
                                    alphahull = 0,
                                    color = 'rgb'+str(fadedColor),
                                    opacity = 0.5,
                                    showlegend = False,
                                    hoverinfo = 'skip'
                                        ))
    else:    
        figure.add_trace(go.Scatter3d(
                                      x = [pos[0]],
                                      y = [pos[1]],
                                      z = [pos[2]],
                                      mode = "markers",
                                      marker = dict(
                                          color = 'rgb'+str(bd.color),
                                          symbol = symbol,
                                          size = size),
                                      showlegend = False,
                                      hoverinfo = 'skip',
                                      ))

def add_soi(figure, bd, time, pos=None):
    
    if pos is None:
        pos = bd.orb.get_state_vector(time)[0]
    fadedColor = fade_color(bd.color)
    
    size = 10
    phi = np.linspace(0, 2*math.pi, size)
    theta = np.linspace(-math.pi/2, math.pi/2, size)
    phi, theta = np.meshgrid(phi, theta)
    
    x = bd.soi * np.cos(theta) * np.sin(phi) + pos[0]
    y = bd.soi * np.cos(theta) * np.cos(phi) + pos[1]
    z = bd.soi * np.sin(theta) + pos[2]
    
    figure.add_trace(go.Mesh3d(
                                x = np.ndarray.flatten(x),
                                y = np.ndarray.flatten(y),
                                z = np.ndarray.flatten(z),
                                alphahull = 0,
                                color = 'rgb'+str(fadedColor),
                                opacity = 0.1,
                                showlegend = False,
                                hoverinfo = 'skip'
                                    ))

def add_burn_arrow(figure, burnDV, burnTime, startOrbit, dateFormat = None,
                       scale=1/2, name = 'Burn', color = (255,0,0), 
                       burnDVisAbsolute = True):
    
    burnPos, preBurnVel = startOrbit.get_state_vector(burnTime)
    
    if burnDVisAbsolute:
        prograde, normal, radial = get_burn_components(burnDV, burnPos, preBurnVel)
    else:
        prograde = burnDV[0]
        normal = burnDV[1]
        radial = burnDV[2]
        
        burnDV =  burn_components_to_absolute(prograde, normal, radial, 
                                              burnPos, preBurnVel)
    
    arrowLength = abs((1-startOrbit.ecc)*startOrbit.a)*scale
    arrowHeadPos = burnPos + burnDV/norm(burnDV)*arrowLength
    figure.add_trace(go.Scatter3d(
        x = [burnPos[0], arrowHeadPos[0]],
        y = [burnPos[1], arrowHeadPos[1]],
        z = [burnPos[2], arrowHeadPos[2]],
        mode = "lines",
        line = dict(
            color = 'rgb'+str(color)),
        hoverinfo = 'skip',
        showlegend = False,
        ))
    figure.add_trace(go.Scatter3d(
        x = [burnPos[0]], 
        y = [burnPos[1]], 
        z = [burnPos[2]],
        mode = "markers",
        marker = dict(
            color = 'rgb'+str(color),
            symbol = 'circle-open',
            size = 5
            ),
        hoverinfo = 'skip',
        showlegend = False,
        ))
    
    if dateFormat is None:
        dateFormat = dict(day = 6, year = 426)
    
    dateString = seconds_to_date_string(burnTime, dateFormat)
    
    direction = burnDV/norm(burnDV) * arrowLength/2
    figure.add_trace(go.Cone(
        x = [arrowHeadPos[0]],
        y = [arrowHeadPos[1]],
        z = [arrowHeadPos[2]],
        u = [direction[0]],
        v = [direction[1]],
        w = [direction[2]],
        colorscale = [[0, 'rgb'+str(color)],
                      [1, 'rgb'+str(color)]],
        showscale = False,
        showlegend = False,
        name = name,
        hovertemplate =                                                     \
                     "{:.2f}".format(prograde) + 'm/s prograde' + "<br>" +  \
                     "{:.2f}".format(normal) + 'm/s normal' + "<br>" +      \
                     "{:.2f}".format(radial) + 'm/s radial'+"<br>"+"<br>"+  \
                     dateString
        ))
    
def add_transfer_phase_angle(figure, transfer, r = None):
    
    if r is None:
        r = 1.5*transfer.transferOrbit.a
    
    if transfer.ejectionBurnAngle is None:
        startAngle = -transfer.startOrbit                                   \
        .get_angle_in_orbital_plane(                                        \
            transfer.get_departure_burn_time(),                             \
            np.array([1,0,0]));
        endAngle = startAngle + transfer.phaseAngle
        arcAngles = np.linspace(startAngle, endAngle, 101)
        arcPos = r * np.array([
            np.cos(arcAngles),
            np.sin(arcAngles),
            0*arcAngles])
        arcPos = transfer.startOrbit.from_orbit_to_primary_bases(arcPos)
    else:
        startAngle = -transfer.startOrbit.prim.orb                          \
        .get_angle_in_orbital_plane(                                        \
            transfer.get_departure_burn_time(),                             \
            np.array([1,0,0]));
        endAngle = startAngle + transfer.phaseAngle
        arcAngles = np.linspace(startAngle, endAngle, 101)
        arcPos = r * np.array([
            np.cos(arcAngles),
            np.sin(arcAngles),
            0*arcAngles])
        arcPos = transfer.startOrbit.prim.orb                               \
            .from_orbit_to_primary_bases(arcPos)
    
    figure.add_trace(go.Scatter3d(
        x = arcPos[0],
        y = arcPos[1],
        z = arcPos[2],
        mode = "lines",
        line = dict(
            color = 'white',
            dash = 'dash'
            ),
        hovertemplate = "Phase Angle: " + 
                        '%.2f' % (transfer.phaseAngle/math.pi*180) +        \
                        "°" + "<extra></extra>",
        name = "phaseAngle",
        showlegend = False,
        ))
    
    midArcX = arcPos[0][int(len(arcPos[0])/2)]
    midArcY = arcPos[1][int(len(arcPos[1])/2)]
    midArcZ = arcPos[2][int(len(arcPos[2])/2)]
    
    figure.update_layout(
        scene = dict(
            annotations=[dict(
                showarrow=False,
                x=midArcX*1.1,
                y=midArcY*1.1,
                z=midArcZ*1.1,
                text= '%.2f' % (transfer.phaseAngle/math.pi*180) + "°",
                font = dict(color = 'white')
            )]
            )
        )

def add_ejection_angle(figure, transfer, r = None):
    
    if not (transfer.ejectionBurnAngle is None):
        if r is None:
            r =  1.5*transfer.startOrbit.a
        rBurn = transfer.startOrbit.from_primary_to_orbit_bases(
            transfer.ejectionTrajectory.get_state_vector(
                transfer.get_departure_burn_time())[0])
        startAngle = math.atan2(rBurn[1],rBurn[0])
        endAngle = startAngle - transfer.ejectionBurnAngle
        arcAngles =  np.linspace(startAngle, endAngle, 101)
        arcPos = r * np.array([
            np.cos(arcAngles),
            np.sin(arcAngles),
            0*arcAngles])
        arcPos = transfer.startOrbit.from_orbit_to_primary_bases(arcPos)
        
        figure.add_trace(go.Scatter3d(
            x = arcPos[0],
            y = arcPos[1],
            z = arcPos[2],
            mode = "lines",
            line = dict(
                color = 'white',
                dash = 'dash'
                ),
            hovertemplate = "Angle from Prograde: " + 
                            '%.2f'%(transfer.ejectionBurnAngle/math.pi*180)+\
                            "°" + "<extra></extra>",
            name = "ejBurnAngle",
            showlegend = False,
            ))
    
    midArcX = arcPos[0][int(len(arcPos[0])/2)]
    midArcY = arcPos[1][int(len(arcPos[1])/2)]
    midArcZ = arcPos[2][int(len(arcPos[2])/2)]
    
    figure.update_layout(
        scene = dict(
            annotations=[dict(
                showarrow=False,
                x=midArcX*1.1,
                y=midArcY*1.1,
                z=midArcZ*1.1,
                text= '%.2f' % (transfer.ejectionBurnAngle/math.pi*180) + "°",
                font = dict(color = 'white')
            )]
            )
        )

def add_prograde_trace(figure, body, t, interval=None, numPts = 201):
    
    if interval is None:
        interval = body.orb.get_period()/8
    startTime = t - interval
    endTime = t + interval
    
    color = body.color
    times = np.linspace(startTime, endTime, numPts)
    
    pos = np.transpose(body.orb.                                            \
                        get_positions(times = times)[0])
    pos = [dim - dim[int(len(dim)/2)] for dim in pos]
    
    figure.add_trace(go.Scatter3d(
        x = pos[0],
        y = pos[1],
        z = pos[2],
        mode = "lines",
        line = dict(
            color = times,
            colorscale = [[0, 'rgb'+str(fade_color(color,4))], 
                      [1, 'rgb'+str(color)]],
            ),
        hoverinfo = 'skip',
        name = 'prograde',
        showlegend = False,
        ))

def add_reference_line(figure, lim, style='dash'):
    
    figure.add_trace(go.Scatter3d(
        x = [0, lim],
        y = [0,0],
        z = [0,0],
        mode = 'lines',
        line = dict(
            color = 'white',
            dash = style,
            ),
        hoverinfo = 'skip',
        name = 'prograde',
        showlegend = False,
        ))

def plot_system(fig, centralBody, t, dateFormat, displays, surfTexture='Solid'):
    
    # add all orbits
    if ('orbits' in displays):
        
        if 'apses' in displays:
            apses = True
        else:
            apses = False
        if 'nodes' in displays:
            nodes = True
        else:
            nodes = False
        
        # add orbits for all satellites around the primary body
        for bd in centralBody.satellites:
            add_orbit(fig, bd.orb, t, None, 201, dateFormat,                \
                      apses = apses, nodes = nodes, color = bd.color,       \
                      name = bd.name);
    
    # add the primary body at the origin
    add_primary(fig, centralBody, False)
    if (not surfTexture == 'Solid') and ('3dSurfs' in displays):
        try:
            mapURL = map_url(centralBody.name, surfTexture+'Small')
            pix = get_pixel_values(mapURL, True)[0]
            bodyTheta = 2*np.pi/centralBody.rotPeriod * t + centralBody.rotIni
            lat = np.array([np.linspace(-np.pi/2, np.pi/2, 512)])
            lon = np.array([np.linspace(-np.pi, np.pi, 512)]) + bodyTheta
            add_primary(fig, centralBody, True, lat, lon, pix)
        except:
            add_primary(fig, centralBody, True)
    elif '3dSurfs' in displays:
        add_primary(fig, centralBody, True)
    
    # add trace for primary body's position centered at the specified time
    if not (centralBody == centralBody.orb.prim):
        add_prograde_trace(fig, centralBody, t);
    
    # finalize axis limit value
    if centralBody.soi is None:
        furthestSatellite = centralBody.satellites[-1]
        lim = 1.25 * furthestSatellite.orb.a * (1 + furthestSatellite.orb.ecc)
    else:
        lim = centralBody.soi
    
    if lim < centralBody.eqr * 5:
        lim = centralBody.eqr * 5
    
    try:
        furthestSatellite = centralBody.satellites[-1]
        soi = centralBody.soi
        satLim = furthestSatellite.orb.a * (1 + furthestSatellite.orb.ecc)
        if soi > 3*satLim:
            lim = 3*satLim
    except:
        pass
    
    # add reference direction line
    if 'ref' in displays:
        add_reference_line(fig, lim)
    
    # add body, SoI positions at specified time
    for bd in centralBody.satellites:
        add_body(fig, bd, t, False)
        if ('3dSurfs' in displays):
            add_body(fig, bd, t, True)
        if ('SoIs' in displays):
            add_soi(fig, bd, t)
    
    return lim

def set_trajectory_plot_layout(fig, lim, cameraDist=None, uirev=None):
    
    if cameraDist is None:
        cameraDist = 0.25
    if uirev is None:
        uirev = 'test'
    
    fig.update_layout(
        margin=dict(l=0, r=0, t=15, b=30),
        paper_bgcolor="rgb(30, 30, 30)",
        scene = dict(
            xaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(30, 30, 30)",
                          ticks='',
                          showticklabels=False,),
            yaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(30, 30, 30)",
                          ticks='',
                          showticklabels=False,),
            zaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(30, 30, 30)",
                          ticks='',
                          showticklabels=False,),
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            aspectratio=dict(x=1, y=1, z=1),
            camera = dict(
                eye = dict( x=cameraDist,
                            y=cameraDist,
                            z=cameraDist)
                ),
            uirevision = uirev
            )
        )

def blank_plot():
    # fig = go.Figure()
    fig = go.Figure(layout =dict(
                        xaxis = dict(visible=False),
                        yaxis = dict(visible=False)
                        )
                    )
    fig.update_layout(
        margin=dict(l=0, r=0, t=15, b=30),
        paper_bgcolor="rgb(30, 30, 30)",
        plot_bgcolor="rgb(30, 30, 30)",
        )
    return fig

#%% surface plot
def project_to_surface(orb, times):
    """Takes list of positions and times and projects to body's surface."""
    
    bd = orb.prim
    positions = orb.get_positions(times=times)[0]
    
    bodyThetas = 2*np.pi/bd.rotPeriod * np.array(times) + bd.rotIni
    sphericalPositions = cartesian_to_spherical(positions)
    
    surfaceCoords = sphericalPositions[:,1:]
    for ii in range(len(times)):
        longitude = surfaceCoords[ii,0] - Orbit.map_angle(bodyThetas[ii])
        count = 0
        while abs(longitude) > np.pi:
            longitude =  longitude - np.sign(longitude) * 2*np.pi
            count = count+1
            if count > 10:
                print('Why?')
                break
        surfaceCoords[ii,0] = longitude
    
    return surfaceCoords

def add_orbit_surface_projection(fig, orb, startTime, endTime=None, numPts=1001,
                                 name = None,
                                 color=(255, 255, 255),
                                 symbol = 'circle', markerSize = 4,
                                 borderColor = None):
    
    if color is None:
        color = (255, 255, 255)
    
    if borderColor is None:
        borderDict = dict()
    else:
        borderDict = dict(color=borderColor,
                          width=2)
    if endTime is None:
        times = [startTime]
        colorscale = [[0.0, "rgb"+str(color)],
                      [1.0, "rgb"+str(color)]]
    else:
        times = np.linspace(startTime, endTime, numPts)
        colorscale = [[0.0, "rgb"+str(color)],
                      [1.0, "rgb"+str(fade_color(color, 3))]]
    longLats = project_to_surface(orb, times)
    
    fig.add_trace(go.Scatter(x=longLats[:,0]*180/np.pi, 
                             y=longLats[:,1]*180/np.pi,
                             name = name,
                             mode='markers',
                             marker = dict(
                                 symbol = symbol,
                                 size = markerSize,
                                 line = borderDict,
                                 color = times,
                                 colorscale = colorscale),
                             hovertemplate =                                \
                                 "Longitude = %{x:.6f} °" + "<br>" +    \
                                 "Latitude = %{y:.6f} m/s" + "<br>" +  \
                                 "UT = %{marker.color:.2f} s",
                             ))

def set_surface_projection_layout(fig, mapUrl=None, uirev=None):
    
    if uirev is None:
        uirev = 'test'
    
    fig.update_layout(
            xaxis = dict(range=[-180, 180]),
            yaxis = dict(range=[-90, 90]),
            xaxis_title='Longitude (°)',
            yaxis_title='Latitude (°)',
            font=dict(
                    color="rgb(200, 200, 200)"
                    )
                )
    fig.update_layout(
        margin=dict(l=0, r=0, t=10, b=30),
        paper_bgcolor="rgba(30, 30, 30, 0)",
        plot_bgcolor="rgba(30, 30, 30, 0)",
        xaxis = dict(
            tickmode = 'array',
            tickvals = [-180, -150, -120, -90, -60, -30, 0,                 \
                        30, 60, 90, 120, 150, 180],
            ticktext = ['180° W', '150° W', '120° W', '90° W', '60° W', '30° W',
                        '0° E', '30° E', '60° E', '90° E', '120° E', '150° E', '180° E'],
            ),
        yaxis = dict(
            tickmode = 'array',
            tickvals = [-90, -60, -30, 0, 30, 60, 90],
            ticktext = ['90° S', '60° S', '30° S', '0° N', '30° N', '60° N', '90° N'],
            ),
        uirevision = uirev
        )
    fig.update_xaxes(color = "rgb(200, 200, 200)",
                     gridcolor="rgb(200, 200, 200)")
    fig.update_yaxes(color = "rgb(200, 200, 200)",
                     gridcolor="rgb(200, 200, 200)")
    
    if not mapUrl is None:
        try:
            fig.update_layout(
                images = [dict(
                    x=-180,
                    sizex=360,
                    y=90,
                    sizey=180,
                    xref="x",
                    yref="y",
                    opacity=1.0,
                    layer="below",
                    sizing="stretch",
                    source=mapUrl)
                    ])
        except:
            pass
