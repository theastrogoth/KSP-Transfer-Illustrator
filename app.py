# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import jsonpickle
import math
import numpy as np
from numpy.linalg import norm
from orbit import Orbit
from body import Body
from transfer import Transfer
from prktable import PorkchopTable
import time

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

server = app.server
app.title='KSP Transfer Illustrator'

#%% read solar system data
infile = open('kerbol_system.json','r')
kerbol_system = jsonpickle.decode(infile.read())
infile.close
infile = open('sol_system.json','r')
sol_system = jsonpickle.decode(infile.read())
infile.close


#%% porkchop plot functions

def add_lines(figure, x, y, minX, maxX, minY, maxY, color = 'black'):
    """ Adds lines and a marker indicating a point at (x,y) """
    
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

#%% trajectory plot functions

def fade_color(color, div = 2):
    """Divides each element of the tuple by the specified number."""
    
    return tuple(math.floor(c/2) for c in color)

def add_orbit(figure, orb, startTime, endTime, numPts=201, dateFormat = None,
              color = (255,255,255), name = '', style = 'solid', fade = True):
    
    if fade:
        fadedColor = fade_color(color)
    else:
        fadedColor = color
    
    period = orb.get_period()
    
    # start and end mean anomalies
    mStart = orb.get_mean_anomaly(startTime)
    if (period < (endTime-startTime) and orb.ecc < 1):
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
    maxVal = np.amax(np.absolute(pos))
    
    if not dateFormat is None:
        day = dateFormat['day']
        year = dateFormat['year']
        
        cData = np.stack((norm(pos, axis = 0)/1000,
                          norm(vel, axis = 0),
                          np.floor(times/(3600*day*year))+1,
                          np.floor(times%(3600*day*year)/(day*3600)+1),
                          np.floor((times%(3600*day))/3600),
                          np.floor(((times%(3600*day))%3600)/60),
                          np.floor(((times%(3600*day))%3600)%60)),
                          axis=1);
        hoverLabel = "r = %{customdata[0]:.3e} km" + "<br>" +\
                     "v = %{customdata[1]:.3e} m/s" + "<br>" +\
                     "Year %{customdata[2]:.0f}, " +\
                     "Day %{customdata[3]:.0f} " +\
                     "%{customdata[4]:0>2d}" + ":" +\
                     "%{customdata[5]:0>2d}" + ":" +\
                     "%{customdata[6]:0>2d}" + "<br>" + "<br>" +\
                     "Semi-major Axis = " + "{:.0f}".format(orb.a) + " m" + "<br>" +\
                     "Eccentricity = " + "{:.4f}".format(orb.ecc) + "<br>" +\
                     "Inclination = " + "{:.4f}".format(orb.inc*180/math.pi) + "°" + "<br>" +\
                     "Argument of the Periapsis = " + "{:.4f}".format(orb.argp*180/math.pi) + "°" + "<br>" +\
                     "Longitude of Ascending Node = " + "{:.4f}".format(orb.lan*180/math.pi) + "°" + "<br>" +\
                     "Mean Anomaly at Epoch = " + "{:.4f}".format(orb.mo) + " rad"
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
    
    return maxVal

def add_primary(figure, bd):
    
    fadedColor = fade_color(bd.color)
    
    size = 10
    phi = np.linspace(0, 2*math.pi, size)
    theta = np.linspace(-math.pi/2, math.pi/2, size)
    phi, theta = np.meshgrid(phi, theta)
    
    x = bd.eqr * np.cos(theta) * np.sin(phi)
    y = bd.eqr * np.cos(theta) * np.cos(phi)
    z = bd.eqr * np.sin(theta)
    
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
    
    # figure.add_trace(go.Scatter3d(
    #                               x = np.array([0]),
    #                               y = np.array([0]),
    #                               z = np.array([0]),
    #                               mode = "markers",
    #                               marker = dict(
    #                                   color = 'rgb'+str(bd.color),
    #                                   symbol = "circle"),
    #                               name = bd.name,
    #                               showlegend = False,
    #                               hovertemplate = "Central body"
    #                               ))

def add_body(figure, bd, time, surf = True):
    
    fadedColor = fade_color(bd.color)
    
    pos = bd.orb.get_state_vector(time)[0]
    
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
                                      x = pos[0],
                                      y = pos[1],
                                      z = pos[2],
                                      mode = "markers",
                                      marker = dict(
                                          color = 'rgb'+str(fadedColor),
                                          symbol = "circle"),
                                      showlegend = False,
                                      hoverinfo = 'skip'
                                      ))
        

def add_burn_arrow(figure, burnDV, burnTime, startOrbit, dateFormat = None,
                       scale=1/2, name = 'Burn', color = (255,0,0)):
    
    burnPos, preBurnVel = startOrbit.get_state_vector(burnTime)
    
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
    
    burnPos = burnPos/norm(burnPos)
    preBurnVel = preBurnVel/norm(preBurnVel)
    prograde = np.dot(preBurnVel, burnDV)
    normal = np.dot(startOrbit.get_basis_vectors()[2], burnDV);
    radial = np.dot(burnPos, burnDV)
    
    if not dateFormat is None:
        day = dateFormat['day']
        year = dateFormat['year']
    else:
        day = 6
        year = 426
    
    burnYear   = math.floor(burnTime/(3600*day*year)+1)
    burnDay    = math.floor(burnTime%(3600*day*year)/(day*3600)+1)
    burnHour   = math.floor((burnTime%(3600*day))/3600)
    burnMinute = math.floor(((burnTime%(3600*day))%3600)/60)
    burnSecond = math.floor(((burnTime%(3600*day))%3600)%60)
    
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
                     'Year ' + str(burnYear) + ', ' +                       \
                     'Day ' + str(burnDay) + ', ' +                         \
                     "{:02d}".format(burnHour) + ':' +                      \
                     "{:02d}".format(burnMinute) + ':' +                    \
                     "{:02d}".format(burnSecond)
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

def add_prograde_trace(figure, transfer, body,
                       startTime, endTime, numPts = 201):
    
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

#%% app layout

app.layout = html.Div(className='row', children=[
    html.Div(className='four columns', children=[
        dcc.Tabs(id='tabs', value='params', children=[
            dcc.Tab(
                label='Instructions',
                value='instruct',
                children = html.Div(className='ctrl-tab', children = [
                    html.H3('KSP Transfer Illustrator'),
                    dcc.Markdown('''
                                 The KSP Transer Illustrator app calculates 
                                 patched-conic transfer trajectories between 
                                 celestial bodies in Kerbal Space Program. 
                                 The app also generates interactive 3D plots 
                                 for each conic in the transfer, which can be 
                                 rotated, zoomed in/out, and panned. Plotted 
                                 orbits also contain hover-text with useful
                                 information.  
                                   
                                 Once a Porkchop plot has been generated, you 
                                 can click anywhere in the plot to view 
                                 information about the transfer with the 
                                 corresponding start time and flight duration.
                                 The 3D orbit plots will be updated when you 
                                 do this.  
                                   
                                 #### Mission Parameters  
                                 
                                 **Date Format**: Display times according to 
                                 the Kerbin or Earth calendars.  
                                   
                                 **System**: Choice of solar system.  
                                   
                                 **Starting Body**: Reference body for the 
                                 starting orbit. Any body in the system is 
                                 valid.  
                                   
                                 **Ending Body**: Reference body for the 
                                 ending orbit. Valid options include the 
                                 starting body and its satellite bodies, and 
                                 the starting body's reference body and its 
                                 satellite bodies.  
                                   
                                 **Starting altitude**: Altitude of the 
                                 starting orbit over its reference body.
                                 The starting orbit is assumed to be circular 
                                 and equatorial unless custom parameters are 
                                 entered in the Advanced Settings tab.  
                                   
                                 **Ending altitude**: Altitude of the ending 
                                 orbit over its reference body. The ending 
                                 orbit is assumed to be circular and 
                                 equatorial unless custom parameters are 
                                 are entered in the Advanced Settings tab.  
                                   
                                 **Cheapest starting orbit**: Aligns orbital 
                                 plane of the starting orbit to the ejection 
                                 trajectory to reduce Δv cost. Useful if you 
                                 haven't yet launched your craft into orbit.  
                                   
                                 **Cheapest ending orbit**: Aligns orbital 
                                 plane of the ending orbit to the insertion 
                                 trajectory to minimize Δv cost. Useful if you 
                                 are not aiming for a specific orbit at the 
                                 target body.  
                                   
                                 **No insertion burn**: Excludes a capture 
                                 burn at the target body. Useful for flybys 
                                 or aerocaptures.  
                                   
                                 **Match starting mean anomaly and epoch**: 
                                 Adjusts the start time to give a departure 
                                 burn that matches the start orbit's 
                                 position. If not selected, the mean anomaly
                                 of the start parking orbit will not match the 
                                 input under Advanced Settings.  
                                   
                                 **Transfer Type**: Choice of inclusion of a 
                                 plane change maneuver.  
                                   
                                 **Earliest Departure Year/Day**: Specifies 
                                 earlist time to search for a transfer.  
                                   
                                 #### Advanced Settings  
                                 
                                 **Departure Time and Flight Duration**: 
                                 Use these settings to customize the times to 
                                 be sampled when computing transfsers. If the 
                                 app fails to due computation timeout, reduce 
                                 the number of points sampled per axis.  
                                   
                                 **Custom Starting/Ending Orbits**: Keplerian 
                                 elements defining starting and ending orbits. 
                                 These can be copy/pasted from HyperEdit.  
                                   
                                 #### Notes
                                 
                                 For all orbit details displayed in the app, 
                                 the epoch is t=0 seconds.  
                                   
                                 The app ignores parking orbits' mean anomaly 
                                 at epoch when calculating transfers, so you 
                                 may need to adjust the start time manually to 
                                 get a transfer where your craft is at the 
                                 burn position at the right time. Hover over 
                                 the starting orbit in the 3D plot to check 
                                 its mean anomaly. For lower altitude parking 
                                 orbits this difference can usually be 
                                 ignored.  
                                   
                                 The app attempts to accurately "patch" 
                                 trajectories across SOI changes to get highly
                                 accurate trajectories, but the patching 
                                 sometimes fails. A warning message is 
                                 is displayed when this happens.  
                                   
                                 Due to the limits of Heroku (and my 
                                 knowledge of app development), this app fails 
                                 when any calculations take longer than 30s. 
                                 If this happens, try decreasing the 
                                 number of points sampled per axis in the 
                                 Advanced Settings tab.  
                                   
                                 I plan on improving this as I learn more 
                                 about astrodynamics and app deployment. 
                                 Please leave feedback about any inaccurate 
                                 trajectories, any bugs you encounter, or 
                                 suggestions you have at the KSP Forum thread 
                                 for this project or its GitHub repository:
                                 '''),
                    html.A(html.Button('KSP Forum Thread'),
                     href='https://forum.kerbalspaceprogram.com/index.php?/topic/195405-ksp-transfer-illustrator/'
                           ),
                    html.A(html.Button('Github'),
                     href='https://github.com/theastrogoth/KSP-Transfer-Illustrator/issues'
                           ),
                    ])
                ),
            dcc.Tab(
                label='Mission Parameters',
                value = 'params',
                children = html.Div([
                    html.H3('Mission Parameters'),
                    html.Div(
                        className='ctrl-name',
                        ),
                    html.Label('Date Format'),
                    dcc.RadioItems(
                        id = 'dateFormat-radio',
                        options=[
                            {'label': 'Kerbin Time (6h days, 426d years)',
                             'value': 'Kerbin'},
                            {'label': 'Earth Time (24h days, 365d years)',
                             'value': 'Earth'},
                            ],
                        value='Kerbin'
                        ),
                    html.Label('System'),
                    dcc.RadioItems(
                        id = 'system-radio',
                        options=[
                            {'label': 'Kerbol', 'value': 'Kerbol'},
                            {'label': 'Sol', 'value': 'Sol'}],
                        value='Kerbol',
                        ),
                    html.Label('Starting Body'),
                    dcc.Dropdown(
                        id = 'startingBody-dropdown',
                        value = 'Kerbin'
                        ),
                    html.Label('Ending Body'),
                    dcc.Dropdown(
                        id = 'endingBody-dropdown',
                        value = 'Duna'
                        ),
                    html.Label('Starting altitude (km)'),
                    dcc.Input(id = 'startPark-input', value=100, 
                              type='number'),
                    html.Label('Ending altitude (km)'),
                    dcc.Input(id = 'endPark-input', value=100, 
                              type='number'),
                    dcc.Checklist(
                        id = 'cheapStartOrbit-checklist',
                        value = [],
                        options=[
                            {'label': 'Cheapest starting orbit',
                             'value': 'True'},
                            ],
                        ),
                    dcc.Checklist(
                        id = 'cheapEndOrbit-checklist',
                        value = ['True'],
                        options=[
                            {'label': 'Cheapest ending orbit',
                             'value': 'True'},
                            ],
                        ),
                    dcc.Checklist(
                        id = 'noInsertion-checklist',
                        value = [],
                        options=[
                            {'label': 'No insertion burn', 'value': 'True'},
                            ],
                        ),
                    dcc.Checklist(
                        id = 'matchMo-checklist',
                        value = [],
                        options=[
                            {'label': 'Match starting mean anomaly and epoch',
                             'value': 'True'}
                            ],
                        ),
                    html.Label('Transfer Type'),
                    dcc.Dropdown(
                        id = 'trType-dropdown',
                        options=[{'label': 'Ballistic', 'value': 'ballistic'},
                                 {'label': 'Mid-course plane change', 
                                  'value': 'plane change'},
                                 {'label': 'Optimal', 'value': 'optimal'}
                                 ],
                        value='ballistic'
                        ),
                    html.Label('Earliest Departure Year'),
                    dcc.Input(id = 'earlyStartYear-input', value=1, 
                              type='number'),
                    html.Label('Earliest Departure Day'),
                    dcc.Input(id = 'earlyStartDay-input', value=1, 
                              type='number'),
                    ])
                ),
                dcc.Tab(
                    label='Advanced Settings',
                    value = 'adv', 
                    children = html.Div(className='six columns', children = [
                        html.H3('Departure Time and Flight Duration'),
                        html.Label('Earliest Departure Year'),
                        dcc.Input(id = 'earlyStartYear2-input', value=1,
                                  type='number'),
                        html.Label('Earliest Departure Day'),
                        dcc.Input(id = 'earlyStartDay2-input', value=1,
                                  type='number'),
                        html.Label('Latest Departure Year'),
                        dcc.Input(id = 'lateStartYear-input',
                                  type='number'),
                        html.Label('Latest Departure Day'),
                        dcc.Input(id = 'lateStartDay-input',
                                  type='number'),
                        html.Label('Shortest Flight Duration (days)'),
                        dcc.Input(id = 'shortFlightDays-input',
                                  type='number'),
                        html.Label('Longest Flight Duration (days)'),
                        dcc.Input(id = 'longFlightDays-input',
                                  type='number'),
                        html.Label('Number of Points Sampled per Axis'),
                        dcc.Input(id = 'numPointsSampled-input', value=25,
                                  type='number'),
                        
                        html.H3('Custom Starting Orbit'),
                        html.Label('Semi-major axis (m)'),
                        dcc.Input(id = 'starta-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Eccentricity'),
                        dcc.Input(id = 'startecc-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Inclination (°)'),
                        dcc.Input(id = 'startinc-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Argument of the Periapsis (°)'),
                        dcc.Input(id = 'startargp-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Longitude of the Ascending Node (°)'),
                        dcc.Input(id = 'startlan-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Mean anomaly at epoch (radians)'),
                        dcc.Input(id = 'startmo-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Epoch (s)'),
                        dcc.Input(id = 'startepoch-input',  
                                  type='number',
                                  value = 0),
                        
                        html.H3('Custom Ending Orbit'),
                        html.Label('Semi-major axis (m)'),
                        dcc.Input(id = 'enda-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Eccentricity'),
                        dcc.Input(id = 'endecc-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Inclination (°)'),
                        dcc.Input(id = 'endinc-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Argument of the Periapsis (°)'),
                        dcc.Input(id = 'endargp-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Longitude of the Ascending Node (°)'),
                        dcc.Input(id = 'endlan-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Mean anomaly at epoch (radians)'),
                        dcc.Input(id = 'endmo-input',  
                                  type='number',
                                  value = 0),
                        html.Label('Epoch (s)'),
                        dcc.Input(id = 'endepoch-input',  
                                  type='number',
                                  value = 0),
                        ])
                    )
            ]),
        ]),
    html.Div(className='four columns', children = [
        html.H3('Porkchop Plot'),
        dcc.Dropdown(
                        id = 'porkchopDisplay-dropdown',
                        options=[{'label': 'Total Δv', 'value': 'total'},
                                 {'label': 'Departure Δv', 'value': 'eject'},
                                 {'label': 'Arrival Δv', 'value': 'insert'}
                                 ],
                        value='total'
                        ),
        html.Button(children = 'Plot!',
                    className = 'button-primary',
                    id = 'porkchop-button',
                    n_clicks = 0
            ),
        html.Div([
            dcc.Loading(id='porkchop-loading',children=[
                dcc.Graph(
                    id='porkchop-graph',
                    figure = go.Figure(layout = dict(
                                        xaxis = dict(visible=False),
                                        yaxis = dict(visible=False))),
                    clickData = None
                    ),
                html.Div(id='porkchop-div', style={'display': 'none'}),
                html.Div(id='transfer-div', style={'display': 'none'}),
                ]),
            ]),
        html.Div([
            html.H3('Selected Transfer Details'),
            html.Div(id='failedConvergenceWarning-div',
                     style={'display': 'none'},
                     children=[
                         dcc.Markdown(id = 'failedConvergence-markdown')]),
            html.Div(id='transferDetails-div', style={'display': 'none'},
                     children=[
                         dcc.Markdown(id = 'departure-markdown'),
                         dcc.Markdown(id = 'arrival-markdown'),
                         dcc.Markdown(id = 'flightTime-markdown'),
                         dcc.Markdown(id = 'phase-markdown'),
                         dcc.Markdown(id = 'totalDV-markdown'),
                         dcc.Markdown(id = 'departureDV-markdown'),
                         dcc.Markdown(id = 'arrivalDV-markdown'),
                         ]),
            html.Div(id='planeChangeDetails-div', style={'display': 'none'},
                     children =[
                         dcc.Markdown(id = 'planeChangeDV-markdown'),
                         dcc.Markdown(id = 'planeChangeTime-markdown'),
                         ]),
            html.Div(id='ejectionDetails-div', style={'display': 'none'},
                     children =[
                         dcc.Markdown(id = 'ejectionAngle-markdown'),
                         dcc.Markdown(id = 'escapeTime-markdown'),
                         ]),
            html.Div(id='insertionDetails-div', style={'display': 'none'},
                     children =[
                         dcc.Markdown(id = 'encounterTime-markdown'),
                         ]),
            html.Div(id='orbitDetails-div', 
                     children = [
                         dcc.Markdown(id = 'transferOrbit-markdown',
                                      style={"white-space": "pre"}),
                         dcc.Markdown(id = 'transferOrbitPC-markdown',
                                      style={"white-space": "pre"}),
                         dcc.Markdown(id = 'ejectionOrbit-markdown',
                                      style={"white-space": "pre"}),
                         dcc.Markdown(id = 'insertionOrbit-markdown',
                                      style={"white-space": "pre"})
                         ])
            ])
        ]),
    html.Div(className='four columns', children = [
        html.Div([
            html.H3('Orbit Plots'),
            dcc.Loading(id='transfer-loading', type='circle', children=[
            html.Div([
                dcc.Markdown('**Transfer Trajectory**'),
                    dcc.Graph(
                        id='transfer-graph',
                        figure = go.Figure(layout = dict(
                                            xaxis = dict(visible=False),
                                            yaxis = dict(visible=False))),
                        ),
                    ]),
            html.Div(id='ejection-div', style={'display': 'none'}, children=[
                dcc.Markdown('**Ejection Trajectory**'),
                    dcc.Graph(
                        id='ejection-graph',
                        figure = go.Figure(layout = dict(
                                            xaxis = dict(visible=False),
                                            yaxis = dict(visible=False))),
                        ),
                    ]),
            html.Div(id='insertion-div', style={'display': 'none'}, children=[
                    dcc.Markdown('**Insertion Trajectory**'),
                    dcc.Graph(
                        id='insertion-graph',
                        figure = go.Figure(layout =dict(
                                            xaxis = dict(visible=False),
                                            yaxis = dict(visible=False))),
                        ),
                    ]),
                ]),
            ]),
        ]),
    # Hidden Divs to store data
    html.Div(id='dateFormat-div', style = {'display': 'none'}),
    html.Div(id='allSystems-div', style = {'display': 'none'},
             children=[
                 jsonpickle.encode(kerbol_system),
                 jsonpickle.encode(sol_system)]),
    html.Div(id='system-div', style={'display': 'none',}, 
             children=jsonpickle.encode(kerbol_system)),
    ])

@app.callback(
    Output('dateFormat-div', 'children'),
    [Input('dateFormat-radio', 'value')]
    )
def set_date_format(selected_format):
    formats = dict(Kerbin = dict(day=6, year=426),
                   Earth = dict(day=24, year = 365))
    return formats[selected_format]

@app.callback(
    [Output('system-div', 'children'),
     Output('startingBody-dropdown', 'value'),
     Output('endingBody-dropdown', 'value')],
    [Input('system-radio','value')],
    [State('allSystems-div', 'children')]
    )
def set_system(system_name, all_systems):
    if system_name == 'Kerbol':
        return all_systems[0], 'Kerbin', 'Duna'
    elif system_name == 'Sol':
        return all_systems[1], 'Earth', 'Mars'
    else:
        return dash.no_update, dash.no_update, dash.no_update

@app.callback(
    Output('startingBody-dropdown', 'options'),
    [Input('system-div', 'children')]
    )
def set_startBody_options(system_data):
    system_data_d = jsonpickle.decode(system_data)
    start_bodies = []
    for bd in system_data_d:
        start_bodies.append(bd.name) 
    return [{'label': i, 'value': i} for i in start_bodies]

@app.callback(
    Output('endingBody-dropdown', 'options'),
    [Input('startingBody-dropdown', 'value')],
    [State('system-div', 'children')]
    )
def set_endBody_options(start_body_name, system_data):
    system_data_d = jsonpickle.decode(system_data)
    sb = [x for x in system_data_d if x.name == start_body_name][0]
    end_bodies = [sb.orb.prim.name]
    if sb != sb.orb.prim:
        for sat in sb.orb.prim.satellites:
            end_bodies.append(sat.name)
    for sat in sb.satellites:
        end_bodies.append(sat.name)
    return [{'label': i, 'value': i} for i in end_bodies]

@app.callback(
    Output('starta-input', 'value'),
    [Input('startingBody-dropdown', 'value'),
     Input('startPark-input', 'value')],
    [State('system-div', 'children')],
    )
def update_start_a(start_body_name, park_alt, system_data):
    system_data_d = jsonpickle.decode(system_data)
    start_body = [x for x in system_data_d if x.name == start_body_name][0]
    return start_body.eqr + 1000*park_alt

@app.callback(
    Output('enda-input', 'value'),
    [Input('endingBody-dropdown', 'value'),
     Input('endPark-input', 'value')],
    [State('system-div', 'children')],
    )
def update_end_a(end_body_name, park_alt, system_data):
    system_data_d = jsonpickle.decode(system_data)
    end_body = [x for x in system_data_d if x.name == end_body_name][0]
    return end_body.eqr + 1000*park_alt

@app.callback(
    Output('earlyStartYear2-input', 'value'),
    [Input('earlyStartYear-input', 'value')],
    [State('earlyStartYear2-input', 'value')]
    )
def update_early_start_year2(early_start_year, prev_state):
    if early_start_year == prev_state:
        return prev_state
    else:
        return early_start_year

@app.callback(
    Output('earlyStartYear-input', 'value'),
    [Input('earlyStartYear2-input', 'value')],
    [State('earlyStartYear-input', 'value')],
    )
def update_early_start_year(early_start_year2, prev_state):
    if early_start_year2 == prev_state:
        return prev_state
    else:
        return early_start_year2

@app.callback(
    Output('earlyStartDay2-input', 'value'),
    [Input('earlyStartDay-input', 'value')],
    [State('earlyStartDay2-input', 'value')]
    )
def update_early_start_day2(early_start_day, prev_state):
    if early_start_day == prev_state:
        return prev_state
    else:
        return early_start_day

@app.callback(
    Output('earlyStartDay-input', 'value'),
    [Input('earlyStartDay2-input', 'value')],
    [State('earlyStartDay-input', 'value')]
    )
def update_early_start_day(early_start_day2, prev_state):
    if early_start_day2 == prev_state:
        return prev_state
    else:
        return early_start_day2

@app.callback(
     Output('porkchop-div','children'),
    [Input('porkchop-button','n_clicks')],
    [State('system-div','children'),
     State('dateFormat-div','children'),
     State('trType-dropdown','value'),
     State('cheapStartOrbit-checklist','value'),
     State('cheapEndOrbit-checklist','value'),
     State('noInsertion-checklist','value'),
     State('startingBody-dropdown','value'),
     State('starta-input','value'),
     State('startecc-input','value'),
     State('startinc-input','value'),
     State('startargp-input','value'),
     State('startlan-input','value'),
     State('startmo-input','value'),
     State('startepoch-input','value'),
     State('endingBody-dropdown','value'),
     State('enda-input','value'),
     State('endecc-input','value'),
     State('endinc-input','value'),
     State('endargp-input','value'),
     State('endlan-input','value'),
     State('endmo-input','value'),
     State('endepoch-input','value'),
     State('earlyStartYear-input','value'),
     State('earlyStartDay-input','value'),
     State('lateStartYear-input','value'),
     State('lateStartDay-input','value'),
     State('shortFlightDays-input','value'),
     State('longFlightDays-input','value'),
     State('numPointsSampled-input','value')]
    )
def update_porkchop_data(nClicks, system, dateFormat,
                          transferType,
                          cheapStartOrb, cheapEndOrb, noInsertion,
                          startBodyName, startA, startEcc, startInc,
                          startArgP, startLAN, startMo, startEpoch,
                          endBodyName, endA, endEcc, endInc,
                          endArgP, endLAN, endMo, endEpoch,
                          minStartYear, minStartDay, maxStartYear, maxStartDay,
                          minFlightDays, maxFlightDays, numPointsSampled):
    
    # return empty plot on page load
    if nClicks == 0:
        return dash.no_update
    
    
    # prepare system information, start & and bodies
    system = jsonpickle.decode(system)
    sBody = [x for x in system if x.name == startBodyName][0]
    eBody = [x for x in system if x.name == endBodyName][0]
    
    # prepare start and end orbit parameters
    if startA is None:
        startA = sBody.eqr + 100000
    if startEcc is None:
        startEcc = 0
    if startInc is None:
        startInc = 0
    if startArgP is None:
        startArgP = 0
    if startLAN is None:
        startLAN = 0
    if startMo is None:
        startMo = 0
    if startEpoch is None:
        startEpoch = 0
    
    if endA is None:
        endA = eBody.eqr + 100000
    if endEcc is None:
        endEcc = 0
    if endInc is None:
        endInc = 0
    if endArgP is None:
        endArgP = 0
    if endLAN is None:
        endLAN = 0
    if endMo is None:
        endMo = 0
    if endEpoch is None:
        endEpoch = 0
    
    sOrb = Orbit(startA, startEcc, startInc*math.pi/180, startArgP*math.pi/180,
                 startLAN*math.pi/180, startMo, sBody, startEpoch)
    eOrb = Orbit(endA, endEcc, endInc*math.pi/180, endArgP*math.pi/180,
                 endLAN*math.pi/180, endMo, eBody, endEpoch)
    # grab day and year formats
    day = dateFormat['day']         # hours per day
    year= dateFormat['year']        # days per year
    
    # prepare start and flight time bounds
    if not (minStartYear is None or minStartYear < 0):
        if (minStartDay is None or minStartDay < 0):
            minStartDay = 0
        minStartTime = 3600*((minStartDay-1) * day +                        \
                    (minStartYear-1) * day * year);
    else:
        minStartYear = 0
        minStartDay = 0
        minStartTime = 0
    if not (maxStartYear is None or maxStartYear < minStartYear):
        if (maxStartDay is None or maxStartDay < 0):
            maxStartDay = year
        maxStartTime = 3600*((maxStartDay-1) * day +                        \
                    (maxStartYear-1) * day * year);
        if maxStartTime < minStartTime:
            maxStartTime = None
    else:
        maxStartTime = None
    if not (minFlightDays is None or minFlightDays <= 0):
        minFlightTime = 3600*(minFlightDays * day)
    else:
        minFlightTime = None
    if not (maxFlightDays is None or maxFlightDays <= 0):
        maxFlightTime = 3600*(maxFlightDays* day)
        if not minFlightTime is None:
            if maxFlightTime < minFlightTime:
                maxFlightTime = None
    else:
        maxFlightTime = None
    
    # change orbit options from strings to bools
    if not cheapStartOrb:
        cheapStartOrb = False
    else:
        cheapStartOrb = True
    if not cheapEndOrb:
        cheapEndOrb = False
    else:
        cheapEndOrb = True
    if not noInsertion:
        noInsertion = False
    else:
        noInsertion = True
    
    # make sure number of points sampled is a valid number
    if numPointsSampled is None:
        numPointsSampled = 25
    elif numPointsSampled <2:
        numPointsSampled = 2
        
    # prepare porkchop table
    porkTable = PorkchopTable(sOrb, eOrb, transferType, noInsertion,
                              cheapStartOrb, cheapEndOrb,
                              minStartTime, maxStartTime,
                              minFlightTime, maxFlightTime,
                              numPointsSampled, numPointsSampled)
    
    return jsonpickle.encode(porkTable)

@app.callback(
    Output('transfer-div','children'),
    [Input('porkchop-div','children'),
     Input('porkchop-graph','clickData')],
    [State('dateFormat-div','children'),
     State('matchMo-checklist','value')]
    )
def update_chosen_tranfser(porkTable, clickData, dateFormat, matchMo):
    
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update
    
    porkTable = jsonpickle.decode(porkTable)
    
    # if the update comes from the porkchop data, get the best transfer
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'porkchop-div':
        transfer = porkTable.get_best_transfer()
    
    # if the update comes from user click data, get the chosen transfer
    elif ctx.triggered[0]['prop_id'].split('.')[0] == 'porkchop-graph':
        day = dateFormat['day']         # hours per day
        startDays = clickData['points'][0]['x']
        flightDays = clickData['points'][0]['y']
        startTime = startDays * 3600 * day
        flightTime = flightDays * 3600 * day
        transfer = porkTable.get_chosen_transfer(startTime, flightTime)
    
    # adjust start orbit to match burn position or vice versa
    if matchMo:
        transfer.match_start_mean_anomaly()
        if not transfer.ejectionTrajectory is None:
            transfer.adjust_end_orbit_mo()
    else:
        transfer.genetic_refine()
        if not transfer.ejectionTrajectory is None:
            transfer.adjust_start_orbit_mo()
        if not transfer.insertionTrajectory is None:
            transfer.adjust_end_orbit_mo()
    return jsonpickle.encode(transfer)

@app.callback(
    Output('porkchop-graph','figure'),
    [Input('porkchop-div','children'),
     Input('transfer-div','children'),
     Input('dateFormat-div','children'),
     Input('porkchopDisplay-dropdown','value')],
    [State('porkchop-graph','figure')]
    )
def update_porkchop_plot(porkTable, chosenTransfer, dateFormat, displayType, 
                         prevState):
    
    ctx = dash.callback_context
    if not ctx.triggered or porkTable is None:
        return prevState
    
    # load porkchop and transfer objects
    porkTable = jsonpickle.decode(porkTable)
    if chosenTransfer is None:
        chosenTransfer = porkTable.get_best_transfer()
    else:
        chosenTransfer = jsonpickle.decode(chosenTransfer)
    # grab day format
    day = dateFormat['day']         # hours per day
    
    # prep empty plot
    fig = go.Figure()
    
    # get the appropriate values according to user input
    if displayType == 'total':
        dV = porkTable.totalDeltaV
        minDV = porkTable.get_best_transfer().get_total_delta_V()
    elif displayType == 'eject':
        dV = porkTable.ejectionDeltaV
        minDV = np.amin(porkTable.ejectionDeltaV)
    elif displayType == 'insert':
        dV = porkTable.insertionDeltaV
        minDV = np.amin(porkTable.insertionDeltaV)
    else:
        dV = porkTable.totalDeltaV
        minDV = porkTable.get_best_transfer().get_total_delta_V()
        
    # prep plot data
    bs = 1.1        # exponent base for contour levels
    numlvls = 20
    lvls = minDV * bs**(np.arange(start=0,stop=numlvls,step=4))
    logdV = np.log(dV)/np.log(bs)
    logLvls = np.log(lvls)/np.log(bs)
    colorVals = logLvls
    colorLabels = [str(math.floor(bs**i)) for i in colorVals]
    
    # create porkchop plot
    fig.add_trace(
                    go.Contour(
                        z = logdV,
                        x = porkTable.startTimes/(day*3600),
                        y = porkTable.flightTimes/(day*3600),
                        contours = dict(
                            start = (logLvls[0]-0.99),
                            end = (logLvls[-1]),
                            size = 1),
                        colorscale = 'Jet',
                        colorbar = dict(
                            tickvals = colorVals,
                            ticktext = colorLabels,
                            title='Δv (m/s)',),
                        contours_coloring='heatmap',
                        customdata = bs**logdV,
                        hovertemplate = "Δv = " +
                                        "%{customdata:.2f}" +
                                        "m/s" +
                                        "<extra></extra>"
                        ))
    fig.update_xaxes(title_text='Transfer start (day #)',
                     visible = True)
    fig.update_yaxes(title_text='Transfer duration (days)',
                     visible = True)
    fig.update_layout(
        margin=dict(l=0, r=0, t=10, b=30),
        )
    add_lines(fig,
              chosenTransfer.startTime/(day*3600),
              chosenTransfer.flightTime/(day*3600),
              porkTable.minStartTime/(day*3600),
              porkTable.maxStartTime/(day*3600),
              porkTable.minFlightTime/(day*3600),
              porkTable.maxFlightTime/(day*3600),
              )
    return fig

@app.callback(
    [Output('failedConvergenceWarning-div','style'),
     Output('failedConvergence-markdown', 'children'),
     Output('transferDetails-div','style'),
     Output('departure-markdown','children'),
     Output('arrival-markdown','children'),
     Output('flightTime-markdown','children'),
     Output('phase-markdown','children'),
     Output('totalDV-markdown','children'),
     Output('departureDV-markdown','children'),
     Output('arrivalDV-markdown','children'),
     Output('planeChangeDetails-div','style'),
     Output('planeChangeDV-markdown','children'),
     Output('planeChangeTime-markdown','children'),
     Output('ejectionDetails-div','style'),
     Output('ejectionAngle-markdown','children'),
     Output('escapeTime-markdown','children'),
     Output('insertionDetails-div','style'),
     Output('encounterTime-markdown','children'),
     Output('transferOrbit-markdown','children'),
     Output('transferOrbitPC-markdown','children'),
     Output('ejectionOrbit-markdown','children'),
     Output('insertionOrbit-markdown','children')],
    [Input('transfer-div','children'),
     Input('dateFormat-div','children')]
    )
def update_transfer_details(chosenTransfer, dateFormat):
    if chosenTransfer is None:
        return dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update;
               
    chosenTransfer = jsonpickle.decode(chosenTransfer)
    # grab day and year formats
    day = dateFormat['day']         # hours per day
    year= dateFormat['year']        # days per year
    
    # Failed convergence warning
    if chosenTransfer.convergenceFail:
        convergenceFailStyle = None
        failString = '**Warning: Trajectories are not consistent across ' + \
            'SOI changes. Solutions are approximate.**';
    else:
        convergenceFailStyle = {'display': 'none'}
        failString  = ''
    
    # transfer trajectory details
    transferStyle = None
    departureTime = chosenTransfer.get_departure_burn_time()
    departureString = '**Departure:** Year ' +                              \
        str(math.floor(departureTime/(3600*day*year)+1)) +                  \
        ', Day ' +                                                          \
        str(math.floor(departureTime%(3600*day*year)/(day*3600)+1)) +       \
        ', ' +                                                              \
        "{:02d}".format(math.floor((departureTime%(3600*day))/3600)) + ':' +\
        "{:02d}".format(math.floor(((departureTime%(3600*day))%3600)/60)) + \
        ':' +                                                               \
        "{:02d}".format(math.floor(((departureTime%(3600*day))%3600)%60));
    arrivalTime = chosenTransfer.get_arrival_burn_time()
    arrivalString = '**Arrival:** Year ' +                                  \
        str(math.floor(arrivalTime/(3600*day*year)+1)) +                    \
        ', Day ' +                                                          \
        str(math.floor(arrivalTime%(3600*day*year)/(day*3600)+1)) +         \
        ', ' +                                                              \
        "{:02d}".format(math.floor((arrivalTime%(3600*day))/3600)) + ':' +  \
        "{:02d}".format(math.floor(((arrivalTime%(3600*day))%3600)/60))+':'+\
        "{:02d}".format(math.floor(((arrivalTime%(3600*day))%3600)%60));
    flightTime = chosenTransfer.get_arrival_burn_time() -                   \
        chosenTransfer.get_departure_burn_time()
    flightTimeString = '**Flight Duration:** ' +                            \
        str(math.floor(flightTime/(3600*day))+1) +                          \
            ' days, ' +                                                     \
        "{:02d}".format(math.floor((flightTime%(3600*day))/3600)) + ':' +   \
        "{:02d}".format(math.floor(((flightTime%(3600*day))%3600)/60)) +':'+\
        "{:02d}".format(math.floor(((flightTime%(3600*day))%3600)%60));
    
    phaseString = '**Phase Angle:** ' +                                     \
        "{:.2f}".format(chosenTransfer.phaseAngle*180/math.pi) + '°';
    totalDVString = '**Total Δv:** ' +                                      \
        "{:.2f}".format(chosenTransfer.get_total_delta_V()) + ' m/s';
    transferOrbitString = '**Transfer Orbit:**\n' +                         \
        str(chosenTransfer.transferOrbit);
    
    # departure burn details
    preDepartPos, preDepartVel =chosenTransfer.startOrbit.get_state_vector( \
        chosenTransfer.get_departure_burn_time())
    preDepartPos = preDepartPos/norm(preDepartPos)
    preDepartVel = preDepartVel/norm(preDepartVel)
    departPrograde = np.dot(preDepartVel, chosenTransfer.ejectionDV)
    departNormal = np.dot(chosenTransfer.startOrbit.get_basis_vectors()[2], \
                          chosenTransfer.ejectionDV);
    departRadial = np.dot(preDepartPos, chosenTransfer.ejectionDV)
    departureDVString = '**Departure Burn:** ' +                            \
        "{:.2f}".format(departPrograde) + ' m/s prograde';
    if abs(departNormal) > 0.05:
        departureDVString = departureDVString + ', ' +                      \
            "{:.2f}".format(departNormal) + ' m/s normal';
    if abs(departRadial) > 0.05:
        departureDVString = departureDVString + ', ' +                      \
            "{:.2f}".format(departRadial) + ' m/s radial';
        
    # arrival burn details
    if not chosenTransfer.ignoreInsertion:
        if chosenTransfer.endOrbit.prim == chosenTransfer.transferOrbit.prim:
            if chosenTransfer.planeChange:
                preArrivePos, preArriveVel =                                \
                    chosenTransfer.transferOrbitPC.get_state_vector(        \
                        chosenTransfer.get_arrival_burn_time());
            else:
                preArrivePos, preArriveVel =                                \
                    chosenTransfer.transferOrbit.get_state_vector(          \
                        chosenTransfer.get_arrival_burn_time());
        else:
            preArrivePos, preArriveVel =                                    \
                chosenTransfer.insertionTrajectory.get_state_vector(        \
                    chosenTransfer.get_arrival_burn_time());
        preArrivePos = preArrivePos/norm(preArrivePos)
        preArriveVel = preArriveVel/norm(preArriveVel)
        arrivePrograde = np.dot(preArriveVel, chosenTransfer.insertionDV)
        arriveNormal =np.dot(chosenTransfer.endOrbit.get_basis_vectors()[2],\
                              chosenTransfer.insertionDV)
        arriveRadial = np.dot(preArrivePos, chosenTransfer.insertionDV)
        arrivalDVString = '**Arrival Burn:** ' +                            \
            "{:.2f}".format(-arrivePrograde) + ' m/s retrograde';
        if abs(arriveNormal) > 0.05:
            arrivalDVString = arrivalDVString + ', ' +                      \
                "{:.2f}".format(arriveNormal) + ' m/s normal';
        if abs(arriveRadial) > 0.05:
            arrivalDVString = arrivalDVString + ', ' +                      \
                "{:.2f}".format(arriveRadial) + ' m/s radial';
    else:
        arrivalDVString = '**Arrival Burn:** None'
    
    # plane change details
    if chosenTransfer.planeChange is True:
        planeChangeStyle = None
        prePlaneChangePos, prePlaneChangeVel =                              \
            chosenTransfer.transferOrbit.get_state_vector(                  \
                chosenTransfer.startTime + chosenTransfer.planeChangeDT);
        prePlaneChangePos = prePlaneChangePos/norm(prePlaneChangePos)
        prePlaneChangeVel = prePlaneChangeVel/norm(prePlaneChangeVel)
        planeChangePrograde = np.dot(prePlaneChangeVel,                     \
                                     chosenTransfer.planeChangeDV);
        planeChangeNormal =                                                 \
            np.dot(chosenTransfer.transferOrbit.get_basis_vectors()[2],     \
                   chosenTransfer.planeChangeDV);
        planeChangeRadial = np.dot(prePlaneChangePos,                       \
                                   chosenTransfer.planeChangeDV);
        planeChangeDVString = '**Plane Change Burn:** ' +                   \
            "{:.2f}".format(planeChangeNormal) + ' m/s normal';
        if abs(planeChangePrograde) > 0.05:
            planeChangeDVString = planeChangeDVString + ', ' +              \
                "{:.2f}".format(planeChangePrograde) + ' m/s prograde';
        if abs(planeChangeRadial) > 0.05:
            planeChangeDVString = planeChangeDVString + ', ' +              \
                "{:.2f}".format(planeChangeRadial) + ' m/s radial';
        planeChangeTime = chosenTransfer.startTime +                        \
            chosenTransfer.planeChangeDT
        planeChangeTimeString = '**Plane Change Time:** Year ' +            \
            str(math.floor(planeChangeTime/(3600*day*year)+1)) +            \
            ', Day ' +                                                      \
            str(math.floor(planeChangeTime%(3600*day*year)/(day*3600)+1)) + \
            ', ' +                                                          \
            "{:02d}".format(math.floor((planeChangeTime%(3600*day))/3600))+ \
            ':' +                                                           \
            "{:02d}".format(math.floor(((planeChangeTime%(3600*day))%3600)/60))+':'+\
            "{:02d}".format(math.floor(((planeChangeTime%(3600*day))%3600)%60));
        transferOrbitPCString = '**Transfer Orbit (plane change):**\n' +    \
            str(chosenTransfer.transferOrbitPC);
    else:
        planeChangeStyle = {'display': 'none'}
        planeChangeDVString = ''
        planeChangeTimeString = ''
        transferOrbitPCString = ''
    
    # ejection details
    if not chosenTransfer.ejectionTrajectory is None:
        ejectionStyle = None
        ejectionAngleString = '**Ejection Angle:**\t\t' +                   \
            "{:.2f}".format(chosenTransfer.ejectionBurnAngle*180/math.pi) + \
            '° from prograde';
        escapeTime = chosenTransfer.startTime
        escapeTimeString = '**Ejection SOI Escape:**\t\tYear ' +            \
        str(math.floor(escapeTime/(3600*day*year)+1)) +                     \
        ', Day ' +                                                          \
        str(math.floor(escapeTime%(3600*day*year)/(day*3600)+1)) +          \
        ', ' +                                                              \
        "{:02d}".format(math.floor((escapeTime%(3600*day))/3600)) + ':' +   \
        "{:02d}".format(math.floor(((escapeTime%(3600*day))%3600)/60)) +':'+\
        "{:02d}".format(math.floor(((escapeTime%(3600*day))%3600)%60));
        ejectionOrbitString = '**Ejection Orbit:**\n' +                     \
            str(chosenTransfer.ejectionTrajectory)
    else:
        ejectionStyle = {'display': 'none'}
        ejectionAngleString = ''
        escapeTimeString = ''
        ejectionOrbitString = ''
    
    # insertion details
    if not chosenTransfer.insertionTrajectory is None:
        insertionStyle = None
        encounterTime = chosenTransfer.startTime + chosenTransfer.flightTime
        encounterTimeString = '**Arrival SOI Encounter:**\t\tYear ' +       \
        str(math.floor(encounterTime/(3600*day*year)+1)) +                  \
        ', Day ' +                                                          \
        str(math.floor(encounterTime%(3600*day*year)/(day*3600)+1)) +       \
        ', ' +                                                              \
        "{:02d}".format(math.floor((encounterTime%(3600*day))/3600)) + ':' +\
        "{:02d}".format(math.floor(((encounterTime%(3600*day))%3600)/60)) + \
        ':' +                                                               \
        "{:02d}".format(math.floor(((encounterTime%(3600*day))%3600)%60));
        insertionOrbitString = '**Insertion Orbit:**\n' +                   \
            str(chosenTransfer.insertionTrajectory)
    else:
        insertionStyle = {'display': 'none'}
        encounterTimeString = ''
        insertionOrbitString = ''
    
    
    return convergenceFailStyle, failString,                                \
            transferStyle, departureString, arrivalString, flightTimeString,\
            phaseString, totalDVString, departureDVString, arrivalDVString, \
            planeChangeStyle, planeChangeDVString, planeChangeTimeString,   \
            ejectionStyle, ejectionAngleString, escapeTimeString,           \
            insertionStyle, encounterTimeString,                            \
            transferOrbitString, transferOrbitPCString, ejectionOrbitString,\
            insertionOrbitString;
            
            

@app.callback(
    Output('transfer-graph', 'figure'),
    [Input('transfer-div', 'children'),
     Input('dateFormat-div', 'children')]
    )
def update_transfer_plot(chosenTransfer, dateFormat):
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    if chosenTransfer is None:
        return fig
    chosenTransfer = jsonpickle.decode(chosenTransfer)
    
    burnTime = chosenTransfer.get_departure_burn_time()
    startTime = chosenTransfer.startTime
    endTime = chosenTransfer.startTime + chosenTransfer.flightTime
    
    if chosenTransfer.planeChange is True:
        pcTime = chosenTransfer.startTime + chosenTransfer.planeChangeDT
    else:
        pcTime = endTime
    
    # record highest altitudes to use as axis limits
    # add transfer orbit
    maxVals = [add_orbit(fig, chosenTransfer.transferOrbit, startTime,      \
                         pcTime, 201, dateFormat, name = 'Transfer')];
    
    # if it exists, add the transfer orbit after plane change
    if chosenTransfer.planeChange is True:
        maxVals = np.append(maxVals,
                            add_orbit(fig, chosenTransfer.transferOrbitPC,  \
                                      pcTime, endTime, 201, dateFormat,     \
                                      name = 'Transfer (plane change)'));
        add_burn_arrow(fig, chosenTransfer.planeChangeDV, pcTime,         \
                       chosenTransfer.transferOrbit, dateFormat, scale=1/4);
    
    # if the starting orbit is around the primary body, add it
    if (chosenTransfer.startOrbit.prim == chosenTransfer.transferOrbit.prim):
        add_orbit(fig, chosenTransfer.startOrbit, startTime, endTime, 201,  \
                  dateFormat, name = 'Start');
        add_burn_arrow(fig, chosenTransfer.ejectionDV, startTime,           \
                       chosenTransfer.startOrbit, dateFormat, scale = 1/4);
    
    # if the target orbit is around the primary body, add it
    if (chosenTransfer.endOrbit.prim == chosenTransfer.transferOrbit.prim):
        add_orbit(fig, chosenTransfer.endOrbit, startTime, endTime, 201,    \
                  dateFormat, name = 'Target');
        if not chosenTransfer.ignoreInsertion:
            if chosenTransfer.planeChange:
                add_burn_arrow(fig, chosenTransfer.insertionDV, endTime,    \
                               chosenTransfer.transferOrbitPC, dateFormat,  \
                               scale = 1/4);
            else:
                add_burn_arrow(fig, chosenTransfer.insertionDV, endTime,    \
                               chosenTransfer.transferOrbit, dateFormat,    \
                               scale = 1/4);
    
    # add orbits/bodies for all satellites around the primary body
    for bd in chosenTransfer.transferOrbit.prim.satellites:
        # if bd.orb.get_period() <                                            \
        # (chosenTransfer.flightTime + chosenTransfer.ejectionDT):
        #     bdTimes = np.linspace(chosenTransfer.get_departure_burn_time(), \
        #                           chosenTransfer.get_departure_burn_time() +\
        #                           bd.orb.get_period(),                      \
        #                           501);
        # else:
        #     bdTimes = np.linspace(chosenTransfer.get_departure_burn_time(), \
        #                           chosenTransfer.startTime +                \
        #                           chosenTransfer.flightTime,                \
        #                           501);
        maxVals = np.append(maxVals,
                            add_orbit(fig, bd.orb, burnTime, endTime, 201, \
                                      dateFormat, bd.color, bd.name));
        add_body(fig, bd, chosenTransfer.get_departure_burn_time())
    
    # add the primary body at the origin
    add_primary(fig, chosenTransfer.transferOrbit.prim)
    
    # finalize axis limit value
    lim = np.amax(maxVals)*1.25
    
    # add transfer phase angle illustration
    add_transfer_phase_angle(fig, chosenTransfer,
                             1.5*chosenTransfer.transferOrbit.a)
    
    # update the plot layout with blank axes, dark grey background, etc
    fig.update_layout(
        margin=dict(l=0, r=0, t=15, b=30),
        paper_bgcolor="rgb(50, 50, 50)",
        scene = dict(
            xaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            yaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            zaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            aspectratio=dict(x=1, y=1, z=1),
            camera = dict(
                eye = dict( x=1.5*chosenTransfer.transferOrbit.a/lim,
                            y=1.5*chosenTransfer.transferOrbit.a/lim,
                            z=1.5*chosenTransfer.transferOrbit.a/lim)
                )
            )
        )
    return fig

@app.callback(
    [Output('ejection-graph', 'figure'),
     Output('ejection-div', 'style')],
    [Input('transfer-div', 'children'),
     Input('dateFormat-div', 'children')]
    )
def update_ejection_plot(chosenTransfer, dateFormat):
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    if chosenTransfer is None:
        # return fig
        return fig, dict(display = 'none')
    
    chosenTransfer = jsonpickle.decode(chosenTransfer)
    if chosenTransfer.ejectionTrajectory is None:
        # return fig
        return fig, dict(display = 'none')
    
    escTime = chosenTransfer.startTime
    burnTime = chosenTransfer.get_departure_burn_time()
    
    # add ejection trajectory
    add_orbit(fig, chosenTransfer.ejectionTrajectory, burnTime, escTime,    \
              501, dateFormat, name = 'Ejection');
    add_burn_arrow(fig, chosenTransfer.ejectionDV, burnTime,                \
                   chosenTransfer.startOrbit, dateFormat);
    
    # record highest altitudes for 3D plot limits
    maxVals = [10*chosenTransfer.startOrbit.a]
    
    # add starting orbit
    maxVals = np.append(maxVals,
                        add_orbit(fig, chosenTransfer.startOrbit,           \
                        burnTime - chosenTransfer.startOrbit.get_period()/2,\
                        burnTime + chosenTransfer.startOrbit.get_period()/2,\
                        201, dateFormat,                                    \
                        name = 'Starting Orbit', style = 'dot', fade = False)
                        );
    
    # add bodies/orbits of satellite bodies around the primary body
    for bd in chosenTransfer.ejectionTrajectory.prim.satellites:
        # if bd.orb.get_period() < chosenTransfer.ejectionDT:
        #     bdTimes = np.linspace(chosenTransfer.get_departure_burn_time(), \
        #                           chosenTransfer.get_departure_burn_time() +\
        #                           bd.orb.get_period(),                      \
        #                           501);
        # else:
        #     bdTimes = np.linspace(chosenTransfer.get_departure_burn_time(), \
        #                           chosenTransfer.startTime,                 \
        #                           501);
        maxVals = np.append(maxVals,
                            add_orbit(fig, bd.orb, burnTime, escTime, 201,  \
                                      dateFormat, bd.color, bd.name));
        add_body(fig, bd, chosenTransfer.get_departure_burn_time())
    
    # finalize value for axis limits
    lim = np.amax(maxVals)*1.25
    
    # add ejection burn angle-from-prograde illustration
    add_ejection_angle(fig, chosenTransfer)
    
    # add trace for primary body's position centered at the burn time
    add_prograde_trace(fig, chosenTransfer, chosenTransfer.startOrbit.prim, \
                       burnTime - chosenTransfer.ejectionDT,                \
                       burnTime + chosenTransfer.ejectionDT, 201);
    add_primary(fig, chosenTransfer.ejectionTrajectory.prim)
    
    # update the plot layout with blank axes, dark grey background, etc
    fig.update_layout(
        margin=dict(l=0, r=0, t=15, b=30),
        paper_bgcolor="rgb(50, 50, 50)",
        scene = dict(
            xaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            yaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            zaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            aspectratio=dict(x=1, y=1, z=1),
            camera = dict(
                eye = dict( x=(2.5*chosenTransfer.startOrbit.a/lim),
                            y=(2.5*chosenTransfer.startOrbit.a/lim),
                            z=(2.5*chosenTransfer.startOrbit.a/lim))
                )
            )
        )
    return fig, dict(display = 'block')

@app.callback(
    [Output('insertion-graph', 'figure'),
     Output('insertion-div', 'style')],
    [Input('transfer-div', 'children'),
     Input('dateFormat-div', 'children')]
    )
def update_insertion_plot(chosenTransfer, dateFormat):
    
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    if chosenTransfer is None:
        # return fig
        return fig, dict(display = 'none')
    
    chosenTransfer = jsonpickle.decode(chosenTransfer)
    if chosenTransfer.insertionTrajectory is None:
        # return fig
        return fig, dict(display = 'none')
    
    encTime = chosenTransfer.startTime + chosenTransfer.flightTime
    burnTime = chosenTransfer.get_arrival_burn_time()
    if chosenTransfer.ignoreInsertion:
        endTime = burnTime + chosenTransfer.insertionDT
    else:
        endTime = burnTime
    
    # parkTimes = np.linspace(burnTime,                                       \
    #                         chosenTransfer.endOrbit.get_period() + burnTime,\
    #                         501);
    
    # inTimes = np.linspace(encTime, burnTime, 501);
    
    # inTimesGeom = np.geomspace(-chosenTransfer.insertionDT, -1, 501) +      \
    #               burnTime + 1;
    
    # if chosenTransfer.ignoreInsertion:
    #     inTimesGeom = np.append(                                            \
    #         inTimesGeom,                                                    \
    #         np.geomspace(1,chosenTransfer.insertionDT,501) + burnTime - 1);
    
    # add insertion trajectory
    add_orbit(fig, chosenTransfer.insertionTrajectory, encTime, endTime,    \
              501, dateFormat, name = 'Insertion');
    
    # record highest altitudes for 3D plot limits
    maxVals = [10*chosenTransfer.endOrbit.a]
    
    # add ending orbit
    if not chosenTransfer.ignoreInsertion:
        maxVals=np.append(maxVals,
                          add_orbit(fig, chosenTransfer.endOrbit,           \
                          burnTime - chosenTransfer.endOrbit.get_period()/2,\
                          burnTime + chosenTransfer.endOrbit.get_period()/2,\
                                    201, dateFormat,                        \
                                    name = 'Ending Orbit', style = 'dot',   \
                                    fade = False));
        add_burn_arrow(fig, chosenTransfer.insertionDV, burnTime,           \
                   chosenTransfer.insertionTrajectory, dateFormat);
    
    # add bodies/orbits of satellite bodies around the primary body
    for bd in chosenTransfer.insertionTrajectory.prim.satellites:
        # if bd.orb.get_period() < chosenTransfer.insertionDT:
        #     bdTimes = np.linspace(encTime,                                  \
        #                           encTime + bd.orb.get_period(),            \
        #                           501);
        # else:
        #     bdTimes = np.linspace(encTime, burnTime, 501);
        maxVals = np.append(maxVals,
                            add_orbit(fig, bd.orb, encTime, endTime, 201,   \
                                      dateFormat, bd.color, bd.name));
        add_body(fig, bd, encTime)
    
    # finalize value for axis limits
    lim = np.amax(maxVals)*1.25
    if lim < chosenTransfer.startOrbit.a * 2:
        lim = chosenTransfer.startOrbit.a * 2
    
    # add trace for primary body's position centered at the burn time
    add_prograde_trace(fig, chosenTransfer, chosenTransfer.endOrbit.prim,   \
                       burnTime - chosenTransfer.insertionDT,               \
                       burnTime + chosenTransfer.insertionDT, 201);
    add_primary(fig, chosenTransfer.insertionTrajectory.prim)
    
    # update the plot layout with blank axes, dark grey background, etc
    fig.update_layout(
        margin=dict(l=0, r=0, t=15, b=30),
        paper_bgcolor="rgb(50, 50, 50)",
        scene = dict(
            xaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            yaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            zaxis = dict(nticks = 0, range=[-lim, lim],
                          showbackground=False,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            aspectratio=dict(x=1, y=1, z=1),
            camera = dict(
                eye = dict( x=(2.5*chosenTransfer.endOrbit.a/lim),
                            y=(2.5*chosenTransfer.endOrbit.a/lim),
                            z=(2.5*chosenTransfer.endOrbit.a/lim))
                )
            )
        )
    return fig, dict(display = 'block')

if __name__ == '__main__':
    app.run_server(debug=False)
