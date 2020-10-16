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
from base64 import b64decode

import jsonpickle
import math
import numpy as np
from orbit import Orbit
from body import Body
from vessel import Vessel
from transfer import Transfer
from prktable import PorkchopTable


DOWNLOAD_DIRECTORY = "/tmp/app_generated_files"

if not os.path.exists(DOWNLOAD_DIRECTORY):
    os.makedirs(DOWNLOAD_DIRECTORY)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

server = Flask(__name__)
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
                server=server)

app.title='KSP Transfer Illustrator'

#%% read solar system data
infile = open('kerbol_system.json','r')
kerbol_system = jsonpickle.decode(infile.read())
infile.close
infile = open('sol_system.json','r')
sol_system = jsonpickle.decode(infile.read())
infile.close

#%%

def name_options(objectList):
    nameOptions = []
    for ob in objectList:
        nameOptions.append(ob.name) 
    return [{'label': i, 'value': i} for i in nameOptions]

#%% download functions

@app.server.route('/download/<path:path>')
def serve_static(path):
    return send_from_directory(DOWNLOAD_DIRECTORY, path, as_attachment=True)

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
                children = [
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
                        value = 'Kerbin',
                        options = name_options(kerbol_system)
                        ),
                    html.Label('Ending Body'),
                    dcc.Dropdown(
                        id = 'endingBody-dropdown',
                        value = 'Duna',
                        options=[
                            {'label': 'Kerbol', 'value': 'Kerbol'},
                            {'label': 'Moho', 'value': 'Moho'},
                            {'label': 'Eve', 'value': 'Eve'},
                            {'label': 'Kerbin', 'value': 'Kerbin'},
                            {'label': 'Mun', 'value': 'Mun'},
                            {'label': 'Minmus', 'value': 'Minmus'},
                            {'label': 'Duna', 'value': 'Duna'},
                            {'label': 'Dres', 'value': 'Dres'},
                            {'label': 'Jool', 'value': 'Jool'},
                            {'label': 'Eeloo', 'value': 'Eeloo'},
                            ],
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
                    html.H3('Load Orbits from sfs File'),
                    dcc.Upload(
                        id='persistenceFile-upload',
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
                    html.Label('Select orbit to add'),
                    dcc.Dropdown(
                        id='persistenceVessels-dropdown',
                        ),
                    html.Button(children = 'Add Starting Orbit',
                                id = 'addStartOrbit-button',
                                n_clicks = 0
                        ),
                    html.Button(children = 'Add Ending Orbit',
                                id = 'addEndOrbit-button',
                                n_clicks = 0
                        ),
                    ]),
                dcc.Tab(
                    label='Advanced Settings',
                    value = 'adv', 
                    children = [
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
                        
                        html.H3('Custom Ending Orbit'),
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
                        
                        html.H3('System Options'),
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
                        ])
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
        html.Div([
        dcc.Loading(id='transfer-loading', type='circle', children=[
            dcc.Markdown('**Transfer Trajectory**'),
            dcc.Graph(
                id='transfer-graph',
                figure = go.Figure(layout = dict(
                                    xaxis = dict(visible=False),
                                    yaxis = dict(visible=False))),
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
            html.A(children=html.Button('Download'),
                   id='transferPlot-download',
                   download="TransferPlot.html", href="",
                   target="_blank"),
            ]),
        html.Div(id='ejection-div', style={'display': 'none'}, children=[
        dcc.Loading(id='ejection-loading', type='circle', children=[
            dcc.Markdown('**Ejection Trajectory**'),
            dcc.Graph(
                id='ejection-graph',
                figure = go.Figure(layout = dict(
                                    xaxis = dict(visible=False),
                                    yaxis = dict(visible=False))),
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
            html.A(children=html.Button('Download'),
                   id='ejectionPlot-download',
                   download="EjectionPlot.html", href="",
                   target="_blank"),
            ]),
        html.Div(id='insertion-div', style={'display': 'none'}, children=[
        dcc.Loading(id='insertion-loading', type='circle', children=[
            dcc.Markdown('**Insertion Trajectory**'),
            dcc.Graph(
                id='insertion-graph',
                figure = go.Figure(layout =dict(
                                    xaxis = dict(visible=False),
                                    yaxis = dict(visible=False))),
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
            html.A(children=html.Button('Download'),
                   id='insertionPlot-download',
                   download="InsertionPlot.html", href="",
                   target="_blank"),
            ]),
        ]),
    # Hidden Divs to store data
    html.Div(id='dateFormat-div', style = {'display': 'none'},
             children = []),
    html.Div(id='allSystems-div', style = {'display': 'none'},
             children=[
                 jsonpickle.encode(kerbol_system),
                 jsonpickle.encode(sol_system)]),
    html.Div(id='system-div', style={'display': 'none'}, 
             children=jsonpickle.encode(kerbol_system)),
    html.Div(id='persistenceVessels-div', style={'display': 'none'},
             children = []),
    html.Div(id='starta-div', style={'display': 'none'},
             children = [700000]),
    html.Div(id='enda-div', style={'display': 'none'},
             children = [420000]),
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
     Output('system-div', 'children'),
    [Input('system-radio','value'),
     Input('systemResize-input','value'),
     Input('systemRescale-input','value')],
    [State('allSystems-div', 'children')]
    )
def set_system(system_name, resizeFactor, rescaleFactor, all_systems):
    if system_name == 'Kerbol':
        system = jsonpickle.decode(all_systems[0])
    elif system_name == 'Sol':
        system = jsonpickle.decode(all_systems[1])
    else:
        return dash.no_update
    
    for body in system:
        body.resize(resizeFactor)
        body.rescale(rescaleFactor)
    
    return jsonpickle.encode(system)

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
                 startLAN*math.pi/180, startMo, startEpoch, sBody)
    eOrb = Orbit(endA, endEcc, endInc*math.pi/180, endArgP*math.pi/180,
                 endLAN*math.pi/180, endMo, endEpoch, eBody)
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
    [Output('transfer-div','children'),
     Output('transfer-slider','min'),
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
     Output('insertion-slider','value')],
    [Input('porkchop-div','children'),
     Input('porkchop-graph','clickData')],
    [State('dateFormat-div','children'),
     State('matchMo-checklist','value')]
    )
def update_chosen_tranfser(porkTable, clickData, dateFormat, matchMo):
    
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update, dash.no_update,              \
            dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
            dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
            dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
            dash.no_update;
    
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
    
    # set slider attributes for transfer plot
    trMinTime = transfer.get_departure_burn_time()
    trMaxTime = transfer.startTime + transfer.flightTime
    trStep = 1
    trMarks = {
        trMinTime: 'Departure',
        trMaxTime: 'Target Encounter'
        }
    trValue = trMinTime
    
    # set slider attributes for ejection plot
    if not transfer.ejectionTrajectory is None:
        ejMinTime = transfer.get_departure_burn_time()
        ejMaxTime = transfer.startTime
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
    if not transfer.insertionTrajectory is None:
        inMinTime = transfer.startTime + transfer.flightTime
        inMaxTime = transfer.get_arrival_burn_time()
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
    
    return jsonpickle.encode(transfer),                                     \
            trMinTime, trMaxTime, trStep, trMarks, trValue,                 \
            ejMinTime, ejMaxTime, ejStep, ejMarks, ejValue,                 \
            inMinTime, inMaxTime, inStep, inMarks, inValue,

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
    add_marker(fig,
              chosenTransfer.startTime/(day*3600),
              chosenTransfer.flightTime/(day*3600)
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
                    "{:.2f}".format(chosenTransfer.get_total_delta_V()) +   \
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
    [Output('transfer-graph', 'figure'),
     Output('transferPlot-download', 'href')],
    [Input('transfer-div', 'children'),
     Input('transfer-slider', 'value'),
     Input('display-checklist', 'value'),
     Input('dateFormat-div', 'children')],
    [State('transfer-graph', 'figure')]
    )
def update_transfer_plot(chosenTransfer, sliderTime, displays,
                         dateFormat, prevFig):
    
    if chosenTransfer is None:
        return prevFig, ""
    
    chosenTransfer = jsonpickle.decode(chosenTransfer)
    
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    
    # plot system at slider time
    lim = plot_system(fig, chosenTransfer.transferOrbit.prim, sliderTime,   \
                dateFormat, displays)
    
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
            chosenTransfer.transferOrbit.prim):
            add_orbit(fig, chosenTransfer.startOrbit, startTime, endTime,
                      201, dateFormat, name = 'Start', apses = apses,       \
                      nodes = nodes);
            if 'arrows' in displays:
                add_burn_arrow(fig, chosenTransfer.ejectionDV, startTime,   \
                               chosenTransfer.startOrbit, dateFormat,       \
                               scale = 1/4);
        
        # if the target orbit is around the primary body, add it
        if (chosenTransfer.endOrbit.prim == chosenTransfer.transferOrbit.prim):
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
    if (not chosenTransfer.planeChange) or                                  \
       (sliderTime < chosenTransfer.startTime + chosenTransfer.planeChangeDT):
           vessel = Body('Vessel',0,0,0,chosenTransfer.transferOrbit)
    else:
           vessel = Body('Vessel',0,0,0,chosenTransfer.transferOrbitPC)
    add_body(fig, vessel, sliderTime, False, size = 4, symbol = 'square')
    
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
    
    return {'data':fig.data,'layout':fig.layout}, location

@app.callback(
    [Output('ejection-graph', 'figure'),
     Output('ejection-div', 'style'),
     Output('ejectionPlot-download','href')],
    [Input('transfer-div', 'children'),
     Input('ejection-slider', 'value'),
     Input('display-checklist', 'value'),
     Input('dateFormat-div', 'children')]
    )
def update_ejection_plot(chosenTransfer, sliderTime, displays, dateFormat):
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
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
     Output('insertion-div', 'style'),
     Output('insertionPlot-download', 'href')],
    [Input('transfer-div', 'children'),
     Input('insertion-slider', 'value'),
     Input('display-checklist', 'value'),
     Input('dateFormat-div', 'children')]
    )
def update_insertion_plot(chosenTransfer, sliderTime, displays, dateFormat):
    
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
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