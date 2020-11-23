import os
from flask import Flask, send_from_directory
from urllib.parse import quote as urlquote

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
from plotutils import *
from sfsutils import parse_savefile
from iniutils import ini_to_system, sort_system, system_to_ini
from base64 import b64decode

import jsonpickle
import math
import numpy as np
from orbit import Orbit
from body import Body
from vessel import Vessel
from transfer import Transfer
from prktable import PorkchopTable
from flyby import Flyby


DOWNLOAD_DIRECTORY = "/tmp/app_generated_files"

if not os.path.exists(DOWNLOAD_DIRECTORY):
    os.makedirs(DOWNLOAD_DIRECTORY)

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

server = Flask(__name__)
# app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
#                 server=server)
app = dash.Dash(__name__, server=server)

app.title='KSP Transfer Illustrator'

#%% read solar system data
infile = open('kerbol_system.json','r')
kerbol_system = jsonpickle.decode(infile.read())
infile.close
infile = open('outer_planets_system.json','r')
outer_planets_system = jsonpickle.decode(infile.read())
infile.close
infile = open('sol_system.json','r')
sol_system = jsonpickle.decode(infile.read())
infile.close

#%% generate initial .ini file
filename = 'ksptiBodies.ini'
inipath = os.path.join(DOWNLOAD_DIRECTORY, filename)
inilocation = "/download/{}".format(urlquote(filename))

system_to_ini([kerbol_system[0]], inipath)

#%%

def name_options(objectList):
    nameOptions = []
    for ob in objectList:
        nameOptions.append(ob.name) 
    return [{'label': i, 'value': i} for i in nameOptions]

def plot_transfer_orbits(fig, chosenTransfer, sliderTime,
                         displays, dateFormat,
                         ignoreStart=False, ignoreEnd=False):
    
    # add the transfer orbit(s)
    burnTime = chosenTransfer.get_departure_burn_time()
    startTime = chosenTransfer.startTime
    endTime = chosenTransfer.startTime + chosenTransfer.flightTime
    
    if chosenTransfer.planeChange is True:
        pcTime = chosenTransfer.startTime + chosenTransfer.planeChangeDT
    else:
        pcTime = endTime
    
    if ('orbits' in displays):
        
        if 'apses' in displays:
            apses = True
        else:
            apses = False
        if 'nodes' in displays:
            nodes = True
        else:
            nodes = False
        
        # add transfer orbit
        add_orbit(fig, chosenTransfer.transferOrbit, startTime,             \
                  pcTime, 201, dateFormat, name = 'Transfer',               \
                  apses = False, nodes = False, fullPeriod = False);
        
        # if it exists, add the transfer orbit after plane change
        if chosenTransfer.planeChange is True:
            add_orbit(fig, chosenTransfer.transferOrbitPC,                  \
                      pcTime, endTime, 201, dateFormat,                     \
                      name = 'Transfer (plane change)',                     \
                      fullPeriod = False, apses = apses, nodes = nodes);
            if 'arrows' in displays:
                add_burn_arrow(fig, chosenTransfer.planeChangeDV, pcTime,   \
                               chosenTransfer.transferOrbit, dateFormat,    \
                               scale=1/4);
        
        # if the starting orbit is around the primary body, add it
        if (chosenTransfer.startOrbit.prim ==                               \
            chosenTransfer.transferOrbit.prim) and not ignoreStart:
            add_orbit(fig, chosenTransfer.startOrbit, startTime, endTime,
                      201, dateFormat, name = 'Start', apses = apses,       \
                      nodes = nodes);
            if 'arrows' in displays:
                add_burn_arrow(fig, chosenTransfer.ejectionDV, startTime,   \
                               chosenTransfer.startOrbit, dateFormat,       \
                               scale = 1/4);
        
        # if the target orbit is around the primary body, add it
        if (chosenTransfer.endOrbit.prim ==                                 \
            chosenTransfer.transferOrbit.prim) and not ignoreEnd:
            add_orbit(fig, chosenTransfer.endOrbit, startTime, endTime,     \
                      201, dateFormat, name = 'Target', apses = apses,      \
                      nodes = nodes);
            if 'arrows' in displays:
                if not chosenTransfer.ignoreInsertion:
                    if chosenTransfer.planeChange:
                        add_burn_arrow(fig, chosenTransfer.insertionDV,     \
                                       endTime,                             \
                                       chosenTransfer.transferOrbitPC,      \
                                       dateFormat, scale = 1/4);
                    else:
                        add_burn_arrow(fig, chosenTransfer.insertionDV,     \
                                       endTime,chosenTransfer.transferOrbit,\
                                       dateFormat, scale = 1/4);
    
    # add transfer phase angle illustration
    if 'angles' in displays:
        add_transfer_phase_angle(fig, chosenTransfer,                       \
                                 1.5*chosenTransfer.transferOrbit.a)
    
    # add marker for vessel position at slider time
    if ((sliderTime < pcTime) and (sliderTime > startTime)):
            vessel = Body('Vessel',0,0,0,chosenTransfer.transferOrbit)
            add_body(fig, vessel, sliderTime, False, size = 4, symbol = 'square')
    elif (sliderTime > pcTime) and (sliderTime < endTime):
            vessel = Body('Vessel',0,0,0,chosenTransfer.transferOrbitPC)
            add_body(fig, vessel, sliderTime, False, size = 4, symbol = 'square')

#%% download functions

@app.server.route('/download/<path:path>')
def serve_static(path):
    return send_from_directory(DOWNLOAD_DIRECTORY, path, as_attachment=True)

#%% app layout

app.layout = html.Div(id='kspti-body', children=[
  html.Div(id='header', className='row', style={'background-image': 'linear-gradient(#003236, #1E1E1E)'}, children=[
      html.Div(className='six columns', children =[
          html.H1("   KSP Transfer Illustrator",
                  style={"white-space": "pre",
                         'float' : 'left',
                         'position' : 'relative',
                         'padding-top' : 50,
                         'padding-right' : 25
                    }),
          ]),
      html.Div(className='two columns', children = [
          html.A([
            html.Img(
                src=app.get_asset_url("KSP forum logo.png"),
                style={
                    'float' : 'left',
                    'position' : 'relative',
                    'padding-top' : 50,
                    'padding-right' : 0
                })
            ], href='https://forum.kerbalspaceprogram.com/index.php?/topic/195405-web-ksp-transfer-illustrator/')
          ]),
      html.Div(className='two columns', children = [
          html.A([
            html.Img(
                src=app.get_asset_url("Github logo.png"),
                style={
                    'float' : 'left',
                    'position' : 'relative',
                    'padding-top' : 25,
                    'padding-right' : 100
                })
            ], href='https://github.com/theastrogoth/KSP-Transfer-Illustrator')
          ]),
      html.Div(className='two columns', children =[
          html.Img(
              src=app.get_asset_url("logo.png"),
              style={
                    'float' : 'left',
                    'position' : 'relative',
                    'padding-top' : 0,
                    'padding-right' : 25
                })
          ]),
      ]),
  html.Div(className='row', children=[
    html.Div(className='four columns', children=[
        dcc.Tabs(id='kspti-control-tabs', className='control-tabs', value='basic', children=[
            dcc.Tab(
                label='Basic Mission Parameters',
                value='basic',
                children = html.Div(className='control-tab', children = [
                    html.H3('Departure and Arrival'),
                    html.Label('Starting Body'),
                    dcc.Dropdown(
                        id = 'startingBody-dropdown',
                        value = 'Kerbin',
                        options = name_options(kerbol_system)
                        ),
                    html.Label('Ending Body'),
                    dcc.Dropdown(
                        id = 'endingBody-dropdown',
                        value = 'Duna',
                        ),
                    html.Label('Flyby Body'),
                    dcc.Dropdown(
                        id = 'flybyBody-dropdown',
                        value = None,
                        ),
                    html.Label('Starting altitude (km)'),
                    dcc.Input(id = 'startPark-input', value=100, 
                              type='number'),
                    html.Label('Ending altitude (km)'),
                    dcc.Input(id = 'endPark-input', value=100, 
                              type='number'),
                    html.Label('Earliest Departure Year'),
                    dcc.Input(id = 'earlyStartYear-input', value=1, 
                              type='number'),
                    html.Label('Earliest Departure Day'),
                    dcc.Input(id = 'earlyStartDay-input', value=1, 
                              type='number'),
                    html.H3('Transfer Options'),
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
                    ])),
            dcc.Tab(
                label='System Settings',
                value='system',
                children = html.Div(className='control-tab', children = [
                    dcc.Tabs(id='system-tabs', className='control-tabs', value='basic', children=[
                      dcc.Tab(
                        label='Basic Settings',
                        value='basic',
                        children = html.Div(className='control-tab', children = [
                            html.H3('System Options'),
                            dcc.RadioItems(
                                id = 'system-radio',
                                options=[
                                    {'label': 'Kerbol (Stock)', 'value': 'stock'},
                                    {'label': 'Kerbol (Outer Planets Mod)', 'value': 'opm'},
                                    {'label': 'Sol (Real Solar System)', 'value': 'rss'},
                                    {'label': 'Uploaded/Custom System', 'value': 'custom'}],
                                value='stock',
                                ),
                            html.H3('System Resizing/Rescaling'),
                            html.Label('System Resize Factor'),
                            dcc.Input(id = 'systemResize-input',  
                                      type='number',
                                      value = 1,
                                      min = 0),
                            html.Label('System Rescale Factor'),
                            dcc.Input(id = 'systemRescale-input',
                                      type='number',
                                      value = 1,
                                      min = 0),
                            html.Label('Day Length Multiplier'),
                            dcc.Input(id='systemDayScale-input',
                                      type='number',
                                      value = 1,
                                      min = 0),
                            html.H3('Load System from .ini File'),
                            dcc.Upload(
                                id='systemFile-upload',
                                className='control-upload',
                                children=html.Div([
                                    'Drag and Drop or ',
                                    html.A('Select Files')
                                ]),
                                style={
                                    'width': '100%',
                                    'height': '60px',
                                    'lineHeight': '60px',
                                    'borderWidth': '1px',
                                    'borderStyle': 'dashed',
                                    'borderRadius': '5px',
                                    'textAlign': 'center',
                                    'margin': '10px'
                                    },
                                multiple=False
                                ),
                            ])),
                      dcc.Tab(
                        label='Advanced Settings',
                        value='adv',
                        children=html.Div(
                            className='control-tab', 
                            children= [
                              html.H3('Add/Edit Body'),
                              html.Label('Body Name'),
                              dcc.Input(id = 'bodyName-input',
                                        type='text'),
                              html.Div(className='row', children=[
                                html.Div(className='five columns',children=[
                                    dcc.Markdown('**Body parameters**',
                                                 style={'padding-top' : 10}),
                                    html.Label('Equatorial Radius (m)'),
                                    dcc.Input(id = 'bodyeqr-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Gravity Parameter (m3/s2)'),
                                    dcc.Input(id = 'bodymu-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Reference ID'),
                                    dcc.Input(id = 'bodyRef-input',  
                                              type='number',
                                              step=1,
                                              value = None),
                                    dcc.Markdown('**Color (RGB)**',
                                                 style={'padding-top' : 10}),
                                    dcc.Input(id='bodyRed-input',
                                              type='number',
                                              placeholder='Red (0-255)',
                                              min=0,
                                              max=255,
                                              value=None),
                                     dcc.Input(id='bodyGreen-input',
                                              type='number',
                                              placeholder='Green (0-255)',
                                              min=0,
                                              max=255,
                                              value=None),
                                     dcc.Input(id='bodyBlue-input',
                                              type='number',
                                              placeholder='Blue (0-255)',
                                              min=0,
                                              max=255,
                                              value=None),
                                    ]),
                                html.Div(className='five columns',children=[
                                    dcc.Markdown('**Orbit parameters**',
                                                 style={'padding-top' : 10}),
                                    html.Label('Reference Body'),
                                    dcc.Dropdown(
                                        id = 'bodyPrim-dropdown',
                                        value = None,
                                        options = name_options(kerbol_system)
                                        ),
                                    html.Label('Semi-major axis (m)'),
                                    dcc.Input(id = 'bodya-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Eccentricity'),
                                    dcc.Input(id = 'bodyecc-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Inclination (°)'),
                                    dcc.Input(id = 'bodyinc-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Argument of the Periapsis (°)'),
                                    dcc.Input(id = 'bodyargp-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Longitude of the Ascending Node (°)'),
                                    dcc.Input(id = 'bodylan-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Mean anomaly at epoch (radians)'),
                                    dcc.Input(id = 'bodymo-input',  
                                              type='number',
                                              value = None),
                                    html.Label('Epoch (s)'),
                                    dcc.Input(id = 'bodyepoch-input',  
                                              type='number',
                                              value = None),
                                    ]),
                                  ]),
                            html.Div(children=[
                                html.Button(children = 'Add/Edit Body',
                                            className = 'button-primary',
                                            id = 'body-button',
                                            n_clicks = 0
                                        ),
                                    ]),
                            html.Div(style={'padding-top' : 20}, children = [
                                html.A(children=html.Button('Download INI file',className='control-button'),
                                       id='iniFile-download',
                                       download="ksptiBodies.ini", href=inilocation,
                                       target="_blank",
                                           ),
                                    ]),
                                ]),
                              ),
                            ]),
                        ]),
                    ),
            dcc.Tab(
                label='Time Settings',
                value='time',
                children = html.Div(className='control-tab', children = [
                    html.H3('Time Settings'),
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
                    html.H3('Porkchop Plot Bounds'),
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
                    ])),
            dcc.Tab(
                label='Orbit Settings',
                value='orbits',
                children = html.Div(className='control-tab', children = [
                    html.Div(className='row', children = [
                        html.Div(className='six columns', children = [
                            html.H3('Starting Orbit'),
                            dcc.Markdown(id = 'start-ref-label',
                                         children='**Reference Body**: Kerbin'),
                            html.Label('Semi-major axis (m)'),
                            dcc.Input(id = 'starta-input',  
                                      type='number',
                                      value = 700000),
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
                            ]),
                        html.Div(className='six columns', children = [
                            html.H3('Ending Orbit'),
                            dcc.Markdown(id = 'end-ref-label',
                                         children='**Reference Body**: Duna'),
                            html.Label('Semi-major axis (m)'),
                            dcc.Input(id = 'enda-input',  
                                      type='number',
                                      value = 420000),
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
                            ]),
                        ]),
                    html.H3('Load Orbits from .sfs File'),
                    dcc.Upload(
                        id='persistenceFile-upload',
                        className='control-upload',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Files')
                        ]),
                        style={
                            'width': '100%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '10px'
                            },
                        multiple=False
                        ),
                    html.Label('Select orbit to add:'),
                    dcc.Dropdown(
                        id='persistenceVessels-dropdown',
                        ),
                    html.Button(children = 'Add Starting Orbit',
                                className='control-button',
                                id = 'addStartOrbit-button',
                                n_clicks = 0
                        ),
                    html.Button(children = 'Add Ending Orbit',
                                className='control-button',
                                id = 'addEndOrbit-button',
                                n_clicks = 0
                        ),
                    ])),
            dcc.Tab(
                label='Instructions',
                value='instruct',
                children = html.Div(className='control-tab', children = [
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
                                   
                                 #### Basic Mission Parameters  
                                 
                                 **Starting Body**: Reference body for the 
                                 starting orbit.  
                                   
                                 **Ending Body**: Reference body for the 
                                 ending orbit.  
                                   
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
                                 
                                 **Earliest Departure Day/Year**: Sets the 
                                 earliest possible time of the departure burn.  
                                   
                                 **Transfer Type**: Choice of inclusion of a 
                                 plane change maneuver.  
                                   
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
                                 position. Useful to get exact timing for 
                                 burn maneuvers.  
                                   
                                 #### System Settings  
                                 
                                 **System**: Choice of solar system.  
                                   
                                 **System Resize Factor**: Multiplies 
                                 equatorial radius for every body in the 
                                 system. Also recalculates mass to match 
                                 original surface gravity.  
                                   
                                 **System Rescale Factor**: Multiplies orbit 
                                 size for every body in the system.  
                                   
                                 **Load System from .ini File**: Upload info 
                                 for your system from an .ini file generated 
                                 by KSPTOT. Necessary for most planet packs.  
                                   
                                 #### Time Settings  
                                 
                                 **Date Format**: Sets the length days and 
                                 years.  
                                   
                                 **Porkchop Plot Bounds**: 
                                 Use these settings to customize the times to 
                                 be sampled when computing transfsers. If the 
                                 app fails to due computation timeout, reduce 
                                 the number of points sampled per axis.  
                                   
                                 #### Orbit Settings  
                                 
                                 **Starting/Ending Orbits**: Keplerian 
                                 elements defining starting and ending orbits. 
                                 These can be copy/pasted from HyperEdit or 
                                 read from a savefile.  
                                   
                                 **Load Orbits from .sfs File**: Upload info 
                                 for all orbits (crafts, asteroids, etc.) from 
                                 your savefile. Space objects can be selected 
                                 by name, and their orbits can be set as the 
                                 beginning or ending orbit.
                                 
                                 #### Notes
                                   
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
                    html.A(html.Button('KSP Forum Thread',className='control-button'),
                     href='https://forum.kerbalspaceprogram.com/index.php?/topic/195405-ksp-transfer-illustrator/'
                           ),
                    html.A(html.Button('Github',className='control-button'),
                     href='https://github.com/theastrogoth/KSP-Transfer-Illustrator/issues'
                           ),
                    ])),
            ]),
        ]),
    html.Div(className='four columns', children = [
        html.Button(children = 'Plot!',
            className = 'button-primary',
            id = 'plot-button',
            n_clicks = 0
            ),
        html.Div(id='porkchop-plot-div', children=[
            html.H3('Porkchop Plot'),
            dcc.Markdown('**Display Type**'),
            dcc.Dropdown(
                        id = 'porkchopDisplay-dropdown',
                        options=[{'label': 'Total Δv', 'value': 'total'},
                                 {'label': 'Departure Δv', 'value': 'eject'},
                                 {'label': 'Arrival Δv', 'value': 'insert'}
                                 ],
                        value='total'
                        ),
            dcc.Loading(id='porkchop-loading',children=[
                dcc.Graph(
                    id='porkchop-graph',
                    figure = blank_plot(),
                    clickData = None
                    ),
                # Hidden divs under loading for porkchop generation
                html.Div(id='porkchop-div', style={'display': 'none'}),
                html.Div(id='transfer-div', style={'display': 'none'}),
                html.Div(id='depart-times-div', style={'display': 'none'}),
                html.Div(id='flight-times-div', style={'display': 'none'}),
                ]),
            ]),
        html.Div(id='transfer-details-div', children=[
            html.Div(id='failedConvergenceWarning-div',
                     style={'display': 'none'},
                     children=[
                         dcc.Markdown(id = 'failedConvergence-markdown')]),
            html.Div(id='transferDetails-div', style={'display': 'none'},
                     children=[
                         html.H3('Selected Transfer Details'),
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
            ]),
        html.Div(id='flyby-details-div', style={'display': 'none'}, children=[
            dcc.Loading(id='flyby-calcs-loading', children = [
                html.Div(style={'padding-top' : 50}, children=[
                    html.Button(children = 'Optimize',
                                className = 'button-primary',
                                id = 'flybyOptimize-button',
                                n_clicks = 0,
                                ),
                    ]),
                html.Div(id='flyby-div', style={'display': 'none'},
                         children = []),
                ])
            ]),
        ]),
    html.Div(className='four columns', children = [
        html.H3('Orbit Plots'),
        dcc.Markdown('**Display Options**'),
        dcc.Checklist(
            id = 'display-checklist',
            value = ['orbits', '3dSurfs', 'SoIs', 'arrows', 'angles'],
            options=[
                {'label': 'Orbits', 'value': 'orbits'},
                {'label': 'Body Surfaces', 'value': '3dSurfs'},
                {'label': 'Spheres of Influence', 'value': 'SoIs'},
                {'label': 'Burn Arrows', 'value': 'arrows'},
                {'label': 'Prograde/Phase Angles', 'value': 'angles'},
                {'label': 'Apses', 'value': 'apses'},
                {'label': 'Nodes', 'value': 'nodes'},
                {'label': 'Reference Direction', 'value': 'ref'},
                ],
            labelStyle={'display': 'inline-block'},
            ),
        dcc.Tabs(id='plot-tabs', className='control-tabs', value='transfer', children=[
          dcc.Tab(
              id='transfer-plot-tab',
              label='Transfer',
              value='transfer',
              children = html.Div(className='control-tab', children = [
            html.Div(id='transfer-plot-div', style={'display': 'none'}, children=[
            dcc.Loading(id='transfer-loading', type='circle', children=[
                dcc.Markdown('**Transfer Trajectory**'),
                dcc.Graph(
                    id='transfer-graph',
                    figure = blank_plot(),
                    ),
                    ]),
                dcc.Slider(
                    id='transfer-slider',
                    min=0,
                    max=1,
                    step=1,
                    marks = dict(),
                    value=0,
                    included=False,
                    updatemode='mouseup'
                    ),
                html.A(children=html.Button('Download', className='control-button'),
                       id='transferPlot-download',
                       download="TransferPlot.html", href="",
                       target="_blank"),
                ]),
            ])),
          dcc.Tab(
              id='ejection-plot-tab',
              label='Ejection',
              value='ejection',
              children = html.Div(className='control-tab', children = [
            html.Div(id='ejection-plot-div', style={'display': 'none'}, children=[
            dcc.Loading(id='ejection-loading', type='circle', children=[
                dcc.Markdown('**Ejection Trajectory**'),
                dcc.Graph(
                    id='ejection-graph',
                    figure = blank_plot(),
                    ),
                    ]),
                dcc.Slider(
                    id='ejection-slider',
                    min=0,
                    max=1,
                    step=1,
                    marks = dict(),
                    value=0,
                    included=False,
                    updatemode='mouseup'
                    ),
                html.A(children=html.Button('Download', className='control-button'),
                       id='ejectionPlot-download',
                       download="EjectionPlot.html", href="",
                       target="_blank"),
                ]),
            ])),
          dcc.Tab(
              id='insertion-plot-tab',
              label='Insertion',
              value='insertion',
              children = html.Div(className='control-tab', children = [
            html.Div(id='insertion-plot-div', style={'display': 'none'}, children=[
            dcc.Loading(id='insertion-loading', type='circle', children=[
                dcc.Markdown('**Insertion Trajectory**'),
                dcc.Graph(
                    id='insertion-graph',
                    figure = blank_plot(),
                    ),
                    ]),
                dcc.Slider(
                    id='insertion-slider',
                    min=0,
                    max=1,
                    step=1,
                    marks = dict(),
                    value=0,
                    included=False,
                    updatemode='mouseup'
                    ),
                html.A(children=html.Button('Download',className='control-button'),
                       id='insertionPlot-download',
                       download="InsertionPlot.html", href="",
                       target="_blank"),
                ]),
            ])),
          dcc.Tab(
              id='flyby-plot-tab',
              label='Flyby',
              value='flyby',
              # style={'display': 'none'},
              children = html.Div(className='control-tab', children = [
            html.Div(id='flyby-plot-div', style={'display': 'none'}, children=[
            dcc.Loading(id='flyby-loading', type='circle', children=[
                dcc.Markdown('**Flyby Trajectory**'),
                dcc.Graph(
                    id='flyby-graph',
                    figure = blank_plot(),
                    ),
                    ]),
                dcc.Slider(
                    id='flyby-slider',
                    min=0,
                    max=1,
                    step=1,
                    marks = dict(),
                    value=0,
                    included=False,
                    updatemode='mouseup'
                    ),
                html.A(children=html.Button('Download',className='control-button'),
                       id='flybyPlot-download',
                       download="FlybyPlot.html", href="",
                       target="_blank"),
                ]),
            ])),
          ]),
        ]),
    # Hidden Divs to store data
    html.Div(id='dateFormat-div', style = {'display': 'none'},
             children = []),
    html.Div(id='allSystems-div', style = {'display': 'none'},
             children=[
                 jsonpickle.encode(kerbol_system),
                 jsonpickle.encode(outer_planets_system),
                 jsonpickle.encode(sol_system),
                 jsonpickle.encode([kerbol_system[0]])]),
    html.Div(id='system-div', style={'display': 'none'}, 
             children=jsonpickle.encode(kerbol_system)),
    html.Div(id='persistenceVessels-div', style={'display': 'none'},
             children = []),
    html.Div(id='starta-div', style={'display': 'none'},
             children = [700000]),
    html.Div(id='enda-div', style={'display': 'none'},
             children = [420000]),
    ])
  ])
#%% callbacks
@app.callback(
     Output('dateFormat-div', 'children'),
    [Input('dateFormat-radio', 'value'),
     Input('systemResize-input', 'value'),
     Input('systemRescale-input','value'),
     Input('systemDayScale-input','value')]
    )
def set_date_format(selected_format, resizeFactor, rescaleFactor, dayFactor):
    formats = dict(Kerbin = dict(day=6, year=426),
                   Earth = dict(day=24, year=365))
    
    dateFormat = formats[selected_format]
    day = dateFormat['day']
    year = dateFormat['year']
    
    day = day * dayFactor
    
    aScale = rescaleFactor
    muScale = resizeFactor**2
    year = year * math.sqrt(aScale**3/muScale) / dayFactor
    
    year = round(year)
    day = round(day)
    
    return dict(day=day, year=year)

@app.callback(
    [Output('allSystems-div','children'),
     Output('system-radio','value'),
     Output('iniFile-download','href')],
    [Input('systemFile-upload','contents'),
     Input('body-button', 'n_clicks')],
    [State('allSystems-div','children'),
     State('system-radio','options'),
     State('bodyName-input','value'),
     State('bodyeqr-input','value'),
     State('bodymu-input','value'),
     State('bodyRef-input','value'),
     State('bodyRed-input','value'),
     State('bodyGreen-input','value'),
     State('bodyBlue-input','value'),
     State('bodyPrim-dropdown','value'),
     State('bodya-input','value'),
     State('bodyecc-input','value'),
     State('bodyinc-input','value'),
     State('bodyargp-input','value'),
     State('bodylan-input','value'),
     State('bodymo-input','value'),
     State('bodyepoch-input','value'),
     State('system-div','children')],
    prevent_initial_call = True
    )
def add_edit_system(iniFile, nClicks, allSystems, radioOptions,
                    name, eqr, mu, ref, red, green, blue,
                    primName, a, ecc, inc, argp, lan, mo, epoch,
                    system):
    ctx = dash.callback_context
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'systemFile-upload':
        if iniFile is None:
            return dash.no_update, dash.no_update
        
        iniFile = iniFile.split(',')[1]
        iniFile = b64decode(iniFile).decode('utf-8')
        newSystem = ini_to_system(iniFile, False)
        
        allSystems[3] = jsonpickle.encode(newSystem)
        
        # create downloadable .ini file bodies in system
        filename = 'ksptiBodies.ini'
        path = os.path.join(DOWNLOAD_DIRECTORY, filename)
        location = "/download/{}".format(urlquote(filename))
        
        system_to_ini(newSystem, path)
        
        return allSystems, 'custom', location
    
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'body-button':
        if nClicks==0:
            return dash.no_update
        system = jsonpickle.decode(system)
        existingRefs = [bd.ref for bd in system]
        if ref in existingRefs:
            ref = int(np.amax(existingRefs) + 1)
        removeID = None
        for ii, bd in enumerate(system):
            if bd.name == name:
                removeID = ii
        if not removeID is None:
            del(system[removeID])
        prim = [bd for bd in system if bd.name==primName][0]
        if primName == name:
            newBody = Body(name, eqr, mu, None, 
                           None,
                           ref, None, (red,green,blue)
                           )
        else:
            newBody = Body(name, eqr, mu, None, 
                           Orbit(a, ecc, inc*math.pi/180, argp*math.pi/180,
                                 lan*math.pi/180, mo, epoch, prim),
                           ref, None, (red,green,blue)
                           )
        system.append(newBody)
        print([sat.name for sat in system[0].satellites])
        system = sort_system(system)
        print([sat.name for sat in system[0].satellites])
        print(prim)
        print(newBody.orb.prim)
        print(system[0])
        allSystems[3] = jsonpickle.encode(system)
        
         # create downloadable .ini file bodies in system
        filename = 'ksptiBodies.ini'
        path = os.path.join(DOWNLOAD_DIRECTORY, filename)
        location = "/download/{}".format(urlquote(filename))
        
        system_to_ini(system, path)
        
        return allSystems, 'custom', location

@app.callback(
     Output('system-div', 'children'),
    [Input('system-radio','value'),
     Input('systemResize-input','value'),
     Input('systemRescale-input','value'),
     Input('allSystems-div', 'children')],
    )
def set_system(system_name, resizeFactor, rescaleFactor, all_systems):
    
    
    if system_name == 'stock':
        system = jsonpickle.decode(all_systems[0])
    elif system_name == 'opm':
        system = jsonpickle.decode(all_systems[1])
    elif system_name == 'rss':
        system = jsonpickle.decode(all_systems[2])
    elif system_name == 'custom':
        system = jsonpickle.decode(all_systems[3])
    else:
        return dash.no_update
    
    for body in system:
        body.resize(resizeFactor)
        body.rescale(rescaleFactor)
    
    return jsonpickle.encode(system)

@app.callback(
    [Output('bodyeqr-input','value'),
     Output('bodymu-input','value'),
     Output('bodyRef-input','value'),
     Output('bodyRed-input','value'),
     Output('bodyGreen-input','value'),
     Output('bodyBlue-input','value'),
     Output('bodyPrim-dropdown','value'),
     Output('bodyPrim-dropdown','options'),
     Output('bodya-input','value'),
     Output('bodyecc-input','value'),
     Output('bodyinc-input','value'),
     Output('bodyargp-input','value'),
     Output('bodylan-input','value'),
     Output('bodymo-input','value'),
     Output('bodyepoch-input','value')],
    [Input('bodyName-input','value'),
     Input('system-div','children')]
    )
def set_edit_body_params(name, system):
    system = jsonpickle.decode(system)
    systemNames = [bd.name for bd in system]
    if not name in systemNames:
        return dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, name_options(system), dash.no_update,        \
               dash.no_update, dash.no_update, dash.no_update,              \
               dash.no_update, dash.no_update, dash.no_update
    
    matchedBody = [bd for bd in system if bd.name==name][0]
    eqr = matchedBody.eqr
    mu = matchedBody.mu
    color = matchedBody.color
    ref = matchedBody.ref
    red = color[0]
    green = color[1]
    blue = color[2]
    prim = matchedBody.orb.prim.name
    if prim == name:
        a = None
        ecc = None
        inc = None
        argp = None
        lan = None
        mo = None
        epoch = None
    else:
        a = matchedBody.orb.a
        ecc = matchedBody.orb.ecc
        inc = matchedBody.orb.inc * 180/math.pi
        argp = matchedBody.orb.argp * 180/math.pi
        lan = matchedBody.orb.lan * 180/math.pi
        mo = matchedBody.orb.mo
        epoch = matchedBody.orb.epoch
    
    return eqr, mu, ref, red, green, blue, prim, name_options(system),      \
           a, ecc, inc, argp, lan, mo, epoch

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
    [Output('endingBody-dropdown', 'options'),
     Output('flybyBody-dropdown', 'options')],
    [Input('startingBody-dropdown', 'value'),
     Input('system-div', 'children')]
    )
def set_endBody_options(start_body_name, system_data):
    system_data_d = jsonpickle.decode(system_data)
    try:
        sb = [x for x in system_data_d if x.name == start_body_name][0]
    except IndexError:
        end_bodies = [bd.name for bd in system_data_d]
        return [{'label': i, 'value': i} for i in end_bodies]
    end_bodies = [sb.orb.prim.name]
    if sb != sb.orb.prim:
        for sat in sb.orb.prim.satellites:
            end_bodies.append(sat.name)
    for sat in sb.satellites:
        end_bodies.append(sat.name)
    
    opts = [{'label': i, 'value': i} for i in end_bodies]
    
    return opts, opts

@app.callback(
     Output('start-ref-label','children'),
    [Input('startingBody-dropdown', 'value')],
    )
def update_start_body_label(startBodyName):
    return ['**Reference Body**: ' + startBodyName]

@app.callback(
     Output('end-ref-label','children'),
    [Input('endingBody-dropdown', 'value')],
    )
def update_end_body_label(endBodyName):
    return ['**Reference Body**: ' + endBodyName]

@app.callback(
     Output('starta-input', 'value'),
    [Input('startingBody-dropdown', 'value'),
     Input('startPark-input', 'value'),
     Input('starta-div','children')],
    [State('system-div', 'children'),
     State('starta-input', 'value')],
    )
def update_start_a_value(start_body_name, park_alt, start_a_div, system_data,
                         prevVal):
    
    ctx = dash.callback_context
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'starta-div':
        return start_a_div
    
    system_data_d = jsonpickle.decode(system_data)
    start_body = [x for x in system_data_d if x.name == start_body_name][0]
    
    if prevVal == start_body.eqr + 1000*park_alt:
        return dash.no_update
    
    return start_body.eqr + 1000*park_alt

@app.callback(
     Output('enda-input', 'value'),
    [Input('endingBody-dropdown', 'value'),
     Input('endPark-input', 'value'),
     Input('enda-div', 'children')],
    [State('system-div', 'children'),
     State('enda-input', 'value')],
    )
def update_end_a_value(end_body_name, park_alt, end_a_div, system_data,
                       prevVal):
    
    ctx = dash.callback_context
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'enda-div':
        return end_a_div
    
    system_data_d = jsonpickle.decode(system_data)
    end_body = [x for x in system_data_d if x.name == end_body_name][0]
    
    if prevVal == end_body.eqr + 1000*park_alt:
        return dash.no_update
    
    return end_body.eqr + 1000*park_alt

@app.callback(
     Output('startPark-input', 'value'),
    [Input('startingBody-dropdown', 'value'),
     Input('starta-input', 'value'),
     Input('starta-div','children')],
    [State('system-div', 'children'),
     State('starta-input', 'value')],
    )
def update_start_altidude(start_body_name, start_a, start_a_div, system_data,
                         prevVal):
    
    system_data_d = jsonpickle.decode(system_data)
    start_body = [x for x in system_data_d if x.name == start_body_name][0]
    
    ctx = dash.callback_context
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'starta-div':
        return (start_a_div - start_body.eqr)/1000
    
    if prevVal == (start_a - start_body.eqr)/1000:
        return dash.no_update
    
    return (start_a - start_body.eqr)/1000

@app.callback(
     Output('endPark-input', 'value'),
    [Input('endingBody-dropdown', 'value'),
     Input('enda-input', 'value'),
     Input('enda-div','children')],
    [State('system-div', 'children'),
     State('enda-input', 'value')],
    )
def update_end_altidude(end_body_name, end_a, end_a_div, system_data,
                         prevVal):
    
    system_data_d = jsonpickle.decode(system_data)
    end_body = [x for x in system_data_d if x.name == end_body_name][0]
    
    ctx = dash.callback_context
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'enda-div':
        return (end_a_div - end_body.eqr)/1000
    
    if prevVal == (end_a - end_body.eqr)/1000:
        return dash.no_update
    
    return (end_a - end_body.eqr)/1000

@app.callback(
     Output('earlyStartYear2-input', 'value'),
    [Input('earlyStartYear-input', 'value')],
    [State('earlyStartYear2-input', 'value')]
    )
def update_early_start_year2(early_start_year, prev_state):
    if early_start_year == prev_state:
        return dash.no_update
    else:
        return early_start_year

@app.callback(
     Output('earlyStartYear-input', 'value'),
    [Input('earlyStartYear2-input', 'value')],
    [State('earlyStartYear-input', 'value')],
    )
def update_early_start_year(early_start_year2, prev_state):
    if early_start_year2 == prev_state:
        return dash.no_update
    else:
        return early_start_year2

@app.callback(
     Output('earlyStartDay2-input', 'value'),
    [Input('earlyStartDay-input', 'value')],
    [State('earlyStartDay2-input', 'value')]
    )
def update_early_start_day2(early_start_day, prev_state):
    if early_start_day == prev_state:
        return dash.no_update
    else:
        return early_start_day

@app.callback(
     Output('earlyStartDay-input', 'value'),
    [Input('earlyStartDay2-input', 'value')],
    [State('earlyStartDay-input', 'value')]
    )
def update_early_start_day(early_start_day2, prev_state):
    if early_start_day2 == prev_state:
        return dash.no_update
    else:
        return early_start_day2

@app.callback(
    [Output('depart-times-div','children'),
     Output('flight-times-div','children')],
    [Input('plot-button','n_clicks'),
     Input('porkchop-graph','relayoutData')],
    [State('dateFormat-div','children'),
     State('earlyStartYear-input','value'),
     State('earlyStartDay-input','value'),
     State('lateStartYear-input','value'),
     State('lateStartDay-input','value'),
     State('shortFlightDays-input','value'),
     State('longFlightDays-input','value')],
    )
def update_flight_time_bounds(nClicks, relayout, dateFormat,
                              minStartYear, minStartDay,
                              maxStartYear, maxStartDay,
                              minFlightDays, maxFlightDays):
    
    # don't update on page load
    if nClicks == 0:
        return dash.no_update, dash.no_update
    
    # grab day and year formats
    day = dateFormat['day']         # hours per day
    year= dateFormat['year']        # days per year
    
    ctx = dash.callback_context
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'porkchop-graph':
        minStartTime = relayout['xaxis.range[0]']*3600*day
        maxStartTime = relayout['xaxis.range[1]']*3600*day
        minFlightTime = relayout['yaxis.range[0]']*3600*day
        maxFlightTime = relayout['yaxis.range[1]']*3600*day
    
    else:
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
        
    return [minStartTime, maxStartTime], [minFlightTime, maxFlightTime]

@app.callback(
    [Output('porkchop-plot-div','style'),
     Output('transfer-details-div','style'),
     Output('flyby-details-div','style')],
    [Input('plot-button','n_clicks')],
    [State('flybyBody-dropdown','value'),
     State('system-div','children')]
    )
def transfer_or_flyby(nClicks, flybyBodyName, system):
    
    system = jsonpickle.decode(system)
    
    flybyBody = [x for x in system if x.name == flybyBodyName]
    print(flybyBody)
    if not len(flybyBody) > 0:
        return None, None, {'display': 'none'}
    else:
        return {'display': 'none'}, {'display': 'none'}, None
    
    
@app.callback(
    [Output('porkchop-div','children'),
     Output('flyby-div','children')],
    [Input('depart-times-div','children'),
     Input('flight-times-div','children'),
     Input('flybyOptimize-button','n_clicks')],
    [State('system-div','children'),
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
     State('flybyBody-dropdown','value'),
     State('endingBody-dropdown','value'),
     State('enda-input','value'),
     State('endecc-input','value'),
     State('endinc-input','value'),
     State('endargp-input','value'),
     State('endlan-input','value'),
     State('endmo-input','value'),
     State('endepoch-input','value'),
     State('numPointsSampled-input','value'),
     State('flyby-div','children')]
    )
def update_porkchop_data(startTimeRange, flightTimeRange, optButtonClicks,
                         system, transferType,
                         cheapStartOrb, cheapEndOrb, noInsertion,
                         startBodyName, startA, startEcc, startInc,
                         startArgP, startLAN, startMo, startEpoch,
                         flybyBodyName,
                         endBodyName, endA, endEcc, endInc,
                         endArgP, endLAN, endMo, endEpoch,
                         numPointsSampled, prevFlyby):
    
    # return empty plot on page load
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update
    
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'flybyOptimize-button':
        if prevFlyby is None:
            return dash.no_update, dash.no_update
        else:
            flyby = jsonpickle.decode(prevFlyby)
            flyby.optimize_flyby()
            return dash.no_update, jsonpickle.encode(flyby)
    
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
                 startLAN*math.pi/180, startMo, startEpoch, sBody)
    eOrb = Orbit(endA, endEcc, endInc*math.pi/180, endArgP*math.pi/180,
                 endLAN*math.pi/180, endMo, endEpoch, eBody)
    
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
    
    # set time bounds
    minStartTime = startTimeRange[0]
    maxStartTime = startTimeRange[1]
    minFlightTime = flightTimeRange[0]
    maxFlightTime = flightTimeRange[1]
    
    flybyBody = [x for x in system if x.name == flybyBodyName]
    if len(flybyBody) > 0:
        
        flyby = Flyby(sOrb, flybyBody[0], eOrb, minStartTime, maxStartTime,
                      ignoreInsertion = noInsertion)
        
        return None, jsonpickle.encode(flyby)
    
    # prepare porkchop table
    porkTable = PorkchopTable(sOrb, eOrb, transferType, noInsertion,
                              cheapStartOrb, cheapEndOrb,
                              minStartTime, maxStartTime,
                              minFlightTime, maxFlightTime,
                              numPointsSampled, numPointsSampled)
    
    return jsonpickle.encode(porkTable), None

@app.callback(
     Output('transfer-div','children'),
    [Input('porkchop-div','children'),
     Input('porkchop-graph','clickData')],
    [State('dateFormat-div','children'),
     State('matchMo-checklist','value')]
    )
def update_chosen_tranfser(porkTable, clickData, dateFormat, matchMo):
    
    ctx = dash.callback_context
    if not ctx.triggered or porkTable is None:
        return None
    
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
        if not transfer.insertionTrajectory is None:
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
        minDV = porkTable.get_best_transfer().get_total_delta_v()
    elif displayType == 'eject':
        dV = porkTable.ejectionDeltaV
        minDV = np.amin(porkTable.ejectionDeltaV)
    elif displayType == 'insert':
        dV = porkTable.insertionDeltaV
        minDV = np.amin(porkTable.insertionDeltaV)
    else:
        dV = porkTable.totalDeltaV
        minDV = porkTable.get_best_transfer().get_total_delta_v()
        
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
                        # colorscale=[[0, "rgb(60,100,225)"],
                        #             [0.2, "rgb(115,255,243)"],
                        #             [0.4, "rgb(110,255,129)"],
                        #             [0.6, "rgb(243,255,140)"],
                        #             [0.8, "rgb(255,183,115)"],
                        #             [1, "rgb(250,85,85)"]],
                        colorscale = 'deep',
                        colorbar = dict(
                            tickvals = colorVals,
                            ticktext = colorLabels,
                            title=dict(
                                text='Δv (m/s)',
                                font=dict(
                                    color="rgb(200, 200, 200)")
                                ),
                            tickcolor="rgb(200, 200, 200)",
                            tickfont=dict(color="rgb(200, 200, 200)")),
                        # contours_coloring='heatmap',
                        customdata = bs**logdV,
                        hovertemplate = "Δv = " +
                                        "%{customdata:.2f}" +
                                        "m/s" +
                                        "<extra></extra>"
                        ))
    fig.update_xaxes(title_text='Transfer start (day #)',
                     visible = True,
                     color = "rgb(200, 200, 200)")
    fig.update_yaxes(title_text='Transfer duration (days)',
                     visible = True,
                     color = "rgb(200, 200, 200)")
    
    fig.update_layout(
        margin=dict(l=0, r=0, t=10, b=30),
        paper_bgcolor="rgb(30, 30, 30)",
        plot_bgcolor="rgb(30, 30, 30)",
        )
    add_marker(fig,
              chosenTransfer.startTime/(day*3600),
              chosenTransfer.flightTime/(day*3600),
              color = 'black'
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
    departureString = '**Departure:** ' +                                   \
                      seconds_to_date_string(departureTime, dateFormat) +   \
                    ' (UT ' + "{:.3f}".format(departureTime) + ' s)';
    
    arrivalTime = chosenTransfer.get_arrival_burn_time()
    arrivalString = '**Arrival:** ' +                                       \
                    seconds_to_date_string(arrivalTime, dateFormat) +       \
                    ' (UT ' + "{:.3f}".format(arrivalTime) + ' s)';
    
    flightTime = arrivalTime - departureTime
    flightDays = seconds_to_days(flightTime, dateFormat)
    flightTimeString = '**Flight Duration:** ' +                            \
                       "{:.2f}".format(flightDays) + ' days';
    
    phaseString = '**Phase Angle:** ' +                                     \
                  "{:.2f}".format(chosenTransfer.phaseAngle*180/math.pi)+'°';
    
    totalDVString = '**Total Δv:** ' +                                      \
                    "{:.2f}".format(chosenTransfer.get_total_delta_v()) +   \
                    ' m/s';
    
    transferOrbitString = '**Transfer Orbit:**\n' +                         \
                          str(chosenTransfer.transferOrbit);
    
    # departure burn details
    departPos, departVel =chosenTransfer.startOrbit.get_state_vector(       \
        chosenTransfer.get_departure_burn_time())
    
    departureDVString = '**Departure Burn:** ' +                            \
        burn_components_string(chosenTransfer.ejectionDV,                   \
                               departPos, departVel)
    
    # arrival burn details
    if not chosenTransfer.ignoreInsertion:
        if chosenTransfer.endOrbit.prim == chosenTransfer.transferOrbit.prim:
            if chosenTransfer.planeChange:
                arrivePos, arriveVel =                                      \
                    chosenTransfer.transferOrbitPC.get_state_vector(        \
                        chosenTransfer.get_arrival_burn_time());
            else:
                arrivePos, arriveVel =                                      \
                    chosenTransfer.transferOrbit.get_state_vector(          \
                        chosenTransfer.get_arrival_burn_time());
        else:
            arrivePos, arriveVel =                                          \
                chosenTransfer.insertionTrajectory.get_state_vector(        \
                    chosenTransfer.get_arrival_burn_time());
        
        arrivalDVString = '**Arrival Burn:** ' +                            \
        burn_components_string(chosenTransfer.insertionDV,                  \
                               arrivePos, arriveVel)
        
    else:
        arrivalDVString = '**Arrival Burn:** None'
    
    # plane change details
    if chosenTransfer.planeChange is True:
        planeChangeStyle = None
        
        planeChangePos, planeChangeVel =                                    \
            chosenTransfer.transferOrbit.get_state_vector(                  \
                chosenTransfer.startTime + chosenTransfer.planeChangeDT);
        
        planeChangeDVString = '**Plane Change Burn:** ' +                   \
            burn_components_string(chosenTransfer.planeChangeDV,            \
                                   planeChangePos, planeChangeVel)
        
        planeChangeTime = chosenTransfer.startTime +                        \
            chosenTransfer.planeChangeDT
        planeChangeTimeString = '**Plane Change Time:** ' +                 \
            seconds_to_date_string(planeChangeTime, dateFormat) +           \
                        ' (UT ' + "{:.3f}".format(planeChangeTime) + ' s)';
        
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
        
        ejectionAngleString = '**Ejection Angle:** ' +                   \
            "{:.2f}".format(chosenTransfer.ejectionBurnAngle*180/math.pi) + \
            '° from prograde';
        
        escapeTime = chosenTransfer.startTime
        escapeTimeString = '**Ejection SOI Escape:** ' +                    \
            seconds_to_date_string(escapeTime, dateFormat)+                 \
                    ' (UT ' + "{:.2f}".format(escapeTime) + ' s)';
        
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
        encounterTimeString = '**Arrival SOI Encounter:** ' +               \
            seconds_to_date_string(encounterTime, dateFormat) +             \
                    ' (UT ' + "{:.2f}".format(encounterTime) + ' s)';
        
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
    [Output('transfer-slider','min'),
     Output('transfer-slider','max'),
     Output('transfer-slider','step'),
     Output('transfer-slider','marks'),
     Output('transfer-slider','value'),
     Output('ejection-slider','min'),
     Output('ejection-slider','max'),
     Output('ejection-slider','step'),
     Output('ejection-slider','marks'),
     Output('ejection-slider','value'),
     Output('insertion-slider','min'),
     Output('insertion-slider','max'),
     Output('insertion-slider','step'),
     Output('insertion-slider','marks'),
     Output('insertion-slider','value'),
     Output('flyby-slider','min'),
     Output('flyby-slider','max'),
     Output('flyby-slider','step'),
     Output('flyby-slider','marks'),
     Output('flyby-slider','value')],
    [Input('transfer-div','children'),
     Input('flyby-div','children')]
    )
def set_slider_times(chosenTransfer, flyby):
    
    if not flyby is None:
        flyby = jsonpickle.decode(flyby)
        transfer1 = flyby.transfer1
        transfer2 = flyby.transfer2
    
    else:
        if chosenTransfer is None:
            return dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update, dash.no_update, dash.no_update,          \
                   dash.no_update;
        transfer1 = jsonpickle.decode(chosenTransfer)
        transfer2 = transfer1
    
    # set slider attributes for transfer plot
    trMinTime = transfer1.get_departure_burn_time()
    trMaxTime = transfer2.startTime + transfer2.flightTime
    trStep = 1
    trMarks = {
        trMinTime: 'Departure',
        trMaxTime: 'Target Encounter'
        }
    trValue = trMinTime
    
    # set slider attributes for ejection plot
    if not transfer1.ejectionTrajectory is None:
        ejMinTime = transfer1.get_departure_burn_time()
        ejMaxTime = transfer1.startTime
        ejStep = 1
        ejMarks = {
            ejMinTime: 'Departure Burn',
            ejMaxTime: 'SoI Escape',
            }
        ejValue = ejMinTime
    else:
        ejMinTime = 0
        ejMaxTime = 1
        ejStep = 1
        ejMarks = dict()
        ejValue = 0
    
    # set slider attributes for insertion plot
    if not transfer2.insertionTrajectory is None:
        inMinTime = transfer2.startTime + transfer2.flightTime
        inMaxTime = transfer2.get_arrival_burn_time()
        inStep = 1
        inMarks = {
            inMinTime: 'SoI Encounter',
            inMaxTime: 'Arrival Burn',
            }
        inValue = inMaxTime
    else:
        inMinTime = 0
        inMaxTime = 1
        inStep = 1
        inMarks = dict()
        inValue = 1
    
    # set slider attributes for flyby plot
    if not flyby is None:
        fbMinTime = transfer1.startTime + transfer1.flightTime
        fbMaxTime = transfer2.startTime
        fbStep = 1
        fbMarks = {
            inMinTime: 'SoI Encounter',
            inMaxTime: 'SoI Escape',
            }
        fbValue = flyby.inOrb.epoch
    else:
        fbMinTime = 0
        fbMaxTime = 1
        fbStep = 1
        fbMarks = dict()
        fbValue = 1
    
    return  trMinTime, trMaxTime, trStep, trMarks, trValue,                 \
            ejMinTime, ejMaxTime, ejStep, ejMarks, ejValue,                 \
            inMinTime, inMaxTime, inStep, inMarks, inValue,                 \
            fbMinTime, fbMaxTime, fbStep, fbMarks, fbValue;

@app.callback(
    [Output('transfer-graph', 'figure'),
     Output('transfer-plot-div','style'),
     Output('transferPlot-download', 'href')],
    [Input('transfer-div', 'children'),
     Input('flyby-div','children'),
     Input('transfer-slider', 'value'),
     Input('display-checklist', 'value'),
     Input('dateFormat-div', 'children')],
    [State('transfer-graph', 'figure')]
    )
def update_transfer_plot(chosenTransfer, flyby, sliderTime, displays,
                         dateFormat, prevFig):
    
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    
    if not flyby is None:
        flyby = jsonpickle.decode(flyby)
        transfers = [flyby.transfer1, flyby.transfer2]
        chosenTransfer = flyby.transfer1
        ignoreStart = [False, True]
        ignoreEnd = [True, False]
    
    else:
        if chosenTransfer is None:
            return prevFig, dash.no_update, ""
        chosenTransfer = jsonpickle.decode(chosenTransfer)
        transfers = [chosenTransfer]
        ignoreStart = [False]
        ignoreEnd = [False]
        
    
    # plot system at slider time
    lim = plot_system(fig, chosenTransfer.transferOrbit.prim, sliderTime,   \
                dateFormat, displays)
    
    # plot all transfer orbits
    for ii, transfer in enumerate(transfers):
        plot_transfer_orbits(fig, transfer, sliderTime, displays, dateFormat,
                             ignoreStart[ii], ignoreEnd[ii])
    
    # update the plot layout with blank axes, dark grey background, etc
    uirev = chosenTransfer.transferOrbit.prim.name
    set_trajectory_plot_layout(fig, lim,                                    \
                               1.5*chosenTransfer.transferOrbit.a/lim,      \
                               uirev);
    
    # create downloadable HTML file of plot
    filename = 'TransferPlot.html'
    path = os.path.join(DOWNLOAD_DIRECTORY, filename)
    location = "/download/{}".format(urlquote(filename))
    
    fig.write_html(path)
    
    return {'data':fig.data,'layout':fig.layout}, dict(display = 'block'),  \
           location

@app.callback(
    [Output('ejection-graph', 'figure'),
     Output('ejection-plot-div', 'style'),
     Output('ejectionPlot-download','href')],
    [Input('transfer-div', 'children'),
     Input('flyby-div', 'children'),
     Input('ejection-slider', 'value'),
     Input('display-checklist', 'value'),
     Input('dateFormat-div', 'children')]
    )
def update_ejection_plot(chosenTransfer, flyby, sliderTime,
                         displays, dateFormat):
    
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    
    if not flyby is None:
        flyby = jsonpickle.decode(flyby)
        chosenTransfer = flyby.transfer1
    
    else:
        if chosenTransfer is None:
            return fig, dict(display = 'none'), ""
        chosenTransfer = jsonpickle.decode(chosenTransfer)
        
    if chosenTransfer.ejectionTrajectory is None:
        return fig, dict(display = 'none'), ""
    
    # plot system at slider time
    lim = plot_system(fig, chosenTransfer.startOrbit.prim, sliderTime,      \
                dateFormat, displays)
    
    # plot start and ejection orbits
    escTime = chosenTransfer.startTime
    burnTime = chosenTransfer.get_departure_burn_time()
    
    if ('orbits' in displays):
        
        if 'apses' in displays:
            apses = True
        else:
            apses = False
        if 'nodes' in displays:
            nodes = True
        else:
            nodes = False
    
        # add ejection trajectory
        add_orbit(fig, chosenTransfer.ejectionTrajectory, burnTime, escTime,\
                  501, dateFormat, name = 'Ejection', fullPeriod = False);
        if 'arrows' in displays:
            add_burn_arrow(fig, chosenTransfer.ejectionDV, burnTime,        \
                           chosenTransfer.startOrbit, dateFormat);
        
        # add starting orbit
        add_orbit(fig, chosenTransfer.startOrbit,                           \
                  burnTime - chosenTransfer.startOrbit.get_period()/2,      \
                  burnTime + chosenTransfer.startOrbit.get_period()/2,      \
                  201, dateFormat, apses = apses, nodes = nodes,            \
                  name = 'Starting Orbit', style = 'dot', fade = False);
    
    # add ejection burn angle-from-prograde illustration
    if 'angles' in displays:
        add_ejection_angle(fig, chosenTransfer)
    
    # add marker for vessel position at slider time
    vessel = Body('Vessel',0,0,0,chosenTransfer.ejectionTrajectory)
    add_body(fig, vessel, sliderTime, False, size = 4, symbol = 'square')
    
    # update the plot layout with blank axes, dark grey background, etc
    uirev = chosenTransfer.startOrbit.prim.name
    set_trajectory_plot_layout(fig, lim, 3*chosenTransfer.startOrbit.a/lim, \
                               uirev)
    
    # create downloadable HTML file of plot
    filename = 'EjectionPlot.html'
    path = os.path.join(DOWNLOAD_DIRECTORY, filename)
    location = "/download/{}".format(urlquote(filename))
    
    fig.write_html(path)
    
    return {'data':fig.data,'layout':fig.layout}, dict(display = 'block'),  \
           location

@app.callback(
    [Output('insertion-graph', 'figure'),
     Output('insertion-plot-div', 'style'),
     Output('insertionPlot-download', 'href')],
    [Input('transfer-div', 'children'),
     Input('flyby-div', 'children'),
     Input('insertion-slider', 'value'),
     Input('display-checklist', 'value'),
     Input('dateFormat-div', 'children')]
    )
def update_insertion_plot(chosenTransfer, flyby, sliderTime,
                          displays, dateFormat):
    
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    
    if not flyby is None:
        flyby = jsonpickle.decode(flyby)
        chosenTransfer = flyby.transfer2
    
    else:
        if chosenTransfer is None:
            return fig, dict(display = 'none'), ""
        chosenTransfer = jsonpickle.decode(chosenTransfer)
        
    if chosenTransfer.insertionTrajectory is None:
        return fig, dict(display = 'none'), ""
    
    # plot system at slider time
    lim = plot_system(fig, chosenTransfer.endOrbit.prim, sliderTime,        \
                dateFormat, displays)
    
    # plot end and insertion orbits
    encTime = chosenTransfer.startTime + chosenTransfer.flightTime
    burnTime = chosenTransfer.get_arrival_burn_time()
    if chosenTransfer.ignoreInsertion:
        endTime = burnTime + chosenTransfer.insertionDT
    else:
        endTime = burnTime
    
    if ('orbits' in displays):
        
        if 'apses' in displays:
            apses = True
        else:
            apses = False
        if 'nodes' in displays:
            nodes = True
        else:
            nodes = False
        
        # add insertion trajectory
        add_orbit(fig, chosenTransfer.insertionTrajectory, encTime, endTime,\
                  501, dateFormat, name = 'Insertion', fullPeriod = False);
        
        # add ending orbit
        if not chosenTransfer.ignoreInsertion:
            add_orbit(fig, chosenTransfer.endOrbit,                         \
                      burnTime - chosenTransfer.endOrbit.get_period()/2,    \
                      burnTime + chosenTransfer.endOrbit.get_period()/2,    \
                      201, dateFormat, apses = apses, nodes = nodes,        \
                      name = 'Ending Orbit', style = 'dot',                 \
                      fade = False);
            if 'arrows' in displays:
                add_burn_arrow(fig, chosenTransfer.insertionDV, burnTime,   \
                           chosenTransfer.insertionTrajectory, dateFormat);
    
    # add marker for vessel position at slider time
    vessel = Body('Vessel',0,0,0,chosenTransfer.insertionTrajectory)
    add_body(fig, vessel, sliderTime, False, size = 4, symbol = 'square')
    
    # update the plot layout with blank axes, dark grey background, etc
    uirev = chosenTransfer.endOrbit.prim.name
    set_trajectory_plot_layout(fig, lim, 3*chosenTransfer.endOrbit.a/lim,   \
                               uirev)
    
    # create downloadable HTML file of plot
    filename = 'InsertionPlot.html'
    path = os.path.join(DOWNLOAD_DIRECTORY, filename)
    location = "/download/{}".format(urlquote(filename))
    
    fig.write_html(path)
    
    return {'data':fig.data,'layout':fig.layout}, dict(display = 'block'),  \
           location

@app.callback(
    [Output('persistenceVessels-div', 'children'),
     Output('persistenceVessels-dropdown', 'options'),
     Output('persistenceVessels-dropdown', 'value')],
    [Input('persistenceFile-upload', 'contents')],
    [State('system-div', 'children')],
     prevent_initial_call=True
      )
def create_orbits_from_persistence_file(persistenceFile, system):
    
    if persistenceFile is None:
        return dash.no_update, dash.no_update, dash.no_update
    
    system = jsonpickle.decode(system)
    persistenceFile = persistenceFile.split(',')[1]
    persistenceFile = b64decode(persistenceFile).decode('utf-8')
    sfsData = parse_savefile(persistenceFile, False)
    sfsVessels = sfsData['GAME']['FLIGHTSTATE']['VESSEL']
    vessels = []
    for sfsVessel in sfsVessels:
        
        name = sfsVessel['name']
        
        a = float(sfsVessel['ORBIT']['SMA'])
        ecc = float(sfsVessel['ORBIT']['ECC'])
        inc = float(sfsVessel['ORBIT']['INC'])
        argp = float(sfsVessel['ORBIT']['LPE'])
        lan = float(sfsVessel['ORBIT']['LAN'])
        mo = float(sfsVessel['ORBIT']['MNA'])
        epoch = float(sfsVessel['ORBIT']['EPH'])
        if 'IDENT' in list(sfsVessel['ORBIT'].keys()):
            primName = sfsVessel['ORBIT']['IDENT']
            primName = primName.replace('Squad/','')
            prim = [bd for bd in system if bd.name == primName][0]
        else:
            primRef = float(sfsVessel['ORBIT']['REF'])
            prim = [bd for bd in system if bd.ref == primRef][0]
        orb = Orbit(a, ecc, inc, argp, lan, mo, epoch, prim)
        
        if not a==0:
            vessels.append(Vessel(name, orb))
    
    vesselOptions = name_options(vessels)
    
    return jsonpickle.encode(vessels), vesselOptions, vessels[0].name

@app.callback(
    [Output('startingBody-dropdown','value'),
     Output('starta-div','children'),
     Output('startecc-input','value'),
     Output('startinc-input','value'),
     Output('startargp-input','value'),
     Output('startlan-input','value'),
     Output('startmo-input','value'),
     Output('startepoch-input','value')],
    [Input('addStartOrbit-button', 'n_clicks')],
    [State('persistenceVessels-div','children'),
     State('persistenceVessels-dropdown','value')]
    )
def add_start_orbit(nClicks, persistenceVessels, addVesselName):
    # don't update on page load
    if nClicks == 0:
        return dash.no_update
    
    persistenceVessels = jsonpickle.decode(persistenceVessels)
    vessel = [vs for vs in persistenceVessels if vs.name == addVesselName][0]
    orb = vessel.orb
    
    return orb.prim.name, orb.a, orb.ecc, orb.inc, orb.argp, orb.lan,       \
           orb.mo, orb.epoch

@app.callback(
    [Output('endingBody-dropdown','value'),
     Output('enda-div','children'),
     Output('endecc-input','value'),
     Output('endinc-input','value'),
     Output('endargp-input','value'),
     Output('endlan-input','value'),
     Output('endmo-input','value'),
     Output('endepoch-input','value')],
    [Input('addEndOrbit-button', 'n_clicks')],
    [State('persistenceVessels-div','children'),
     State('persistenceVessels-dropdown','value')]
    )
def add_end_orbit(nClicks, persistenceVessels, addVesselName):
    # don't update on page load
    if nClicks == 0:
        return dash.no_update
    
    persistenceVessels = jsonpickle.decode(persistenceVessels)
    vessel = [vs for vs in persistenceVessels if vs.name == addVesselName][0]
    orb = vessel.orb
    
    return orb.prim.name, orb.a, orb.ecc, orb.inc, orb.argp, orb.lan,       \
           orb.mo, orb.epoch

#%% run app

if __name__ == '__main__':
    app.run_server(debug=False)