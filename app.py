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

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

server = app.server

#%% read solar system data
infile = open('kerbol_system.json','r')
kerbol_system = jsonpickle.decode(infile.read())
infile.close

#%% prepare solar system dictionaries
start_bodies = []
for bd in kerbol_system:
    start_bodies.append(bd.name) 

end_bodies = [
              [[bd.orb.prim.name],
              [sat.name for sat in bd.orb.prim.satellites if bd != bd.orb.prim],
              [sat.name for sat in bd.satellites]] for bd in kerbol_system]

eb = []
for startBody in end_bodies:
    eb.append([item for sublist in startBody for item in sublist if item])

start_body_options = {'Kerbol': start_bodies}
#                'Sol': ['Mercury','Venus','Earth']}

end_body_options = dict()
for x,bd in enumerate(kerbol_system):
    end_body_options[bd.name] = eb[x]


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

def add_orbit(figure, orb, times, dateFormat, 
              color = (255,255,255), name = '', style = 'solid', fade = True):
    """
    
    """
    
    if fade:
        fadedColor = fade_color(color)
    else:
        fadedColor = color
    
    day = dateFormat['day']
    year = dateFormat['year']
    
    pos = np.transpose(orb.get_grandparent_positions(times = times))
    maxVal = np.amax(np.absolute(pos))
    
    figure.add_trace(go.Scatter3d(
        x = pos[0],
        y = pos[1],
        z = pos[2],
        customdata = np.stack((norm(pos, axis = 0)/1000,
                            np.floor(times/(3600*day*year))+1,
                            np.floor(times%(3600*day*year)/(day*3600)+1),
                            np.floor((times%(3600*day))/3600),
                            np.floor(((times%(3600*day))%3600)/60),
                            np.floor(((times%(3600*day))%3600)%60)),
                          axis=1),
        mode = "lines",
        line = dict(
            color = times,
            colorscale = [[0, 'rgb'+str(fadedColor)], 
                      [1, 'rgb'+str(color)]],
            dash = style
            ),
        hovertemplate = "r = %{customdata[0]:.4e} km" + "<br>" +\
                        "Year %{customdata[1]:.0f}, " +\
                            "Day %{customdata[2]:.0f} " +\
                                "%{customdata[3]:0>2d}" + ":" +\
                                    "%{customdata[4]:0>2d}" + ":" +\
                                        "%{customdata[5]:0>2d}",
        name = name,
        showlegend = False,
        ))
    
    return maxVal

def add_primary(figure, bd):
    """
    
    """
    
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
    
def add_transfer_phase_angle(figure, transfer, r = None):
    """
    
    """
    if r is None:
        r = 1.5*transfer.transferOrbit.a
    
    if transfer.ejectionBurnAngle is None:
        startAngle = -transfer.startOrbit                                   \
        .get_angle_in_orbital_plane(                                        \
            transfer.get_departure_burn_time(),                             \
            np.array([1,0,0]));
        endAngle = startAngle + transfer.phaseAngle
        numPoints = abs(math.ceil((endAngle-startAngle)/math.pi*50))
        arcAngles = np.linspace(startAngle, endAngle,numPoints)
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
        numPoints = math.ceil((endAngle-startAngle)/math.pi*50)
        arcAngles = np.linspace(startAngle, endAngle,numPoints)
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
    """
    
    """
    if not (transfer.ejectionBurnAngle is None):
        if r is None:
            r =  1.5*transfer.startOrbit.a
        rBurn = transfer.startOrbit.prim.orb.from_primary_to_orbit_bases(
            transfer.ejectionTrajectory.get_state_vector(
                transfer.get_departure_burn_time())[0])
        startAngle = math.atan2(rBurn[1],rBurn[0])
        endAngle = startAngle - transfer.ejectionBurnAngle
        numPoints = abs(math.ceil((endAngle-startAngle)/math.pi*50))
        arcAngles =  np.linspace(startAngle, endAngle, numPoints)
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

def add_prograde_trace(figure, transfer, times):
    """
    
    """
    
    color = transfer.startOrbit.prim.color
    
    pos = np.transpose(transfer.startOrbit.prim.orb.                        \
                        get_grandparent_positions(times = times))
    pos = [dim - dim[int(len(dim)/2)] for dim in pos]
    
    figure.add_trace(go.Scatter3d(
        x = pos[0],
        y = pos[1],
        z = pos[2],
        mode = "lines",
        line = dict(
            color = times,
            colorscale = [[0, 'rgb'+str(fade_color(color))], 
                      [1, 'rgb'+str(color)]],
            ),
        hoverinfo = 'skip',
        name = 'prograde',
        showlegend = False,
        ))

#%% app layout

app.layout = html.Div(className='row', children=[
    html.Div(className='four columns', children=[
        dcc.Tabs(id='tabs', value='instruct', children=[
            dcc.Tab(
                label='Instructions',
                value='instruct',
                children = html.Div(className='ctrl-tab', children = [
                    html.H3('KSP Transfer Illustrator'),
                    html.P(
                        'Use the KSP transfer illustrator to calculate '
                        'optimal departure and arrival dates for transfers '
                        'between celestial bodies, and visualize ejection '
                        'and transfer orbits. '
                    ),
                    html.P(
                        'Select the details of the transfer in the '
                        'Mission Parameters tab. Custom orbits can be '
                        'entered in the Advanced Settings tab. '
                            )
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
                        # options=[{'label': k, 'value': k}                   \
                        #          for k in start_body_options.keys()],
                        options=[
                            {'label': 'Kerbol', 'value': 'Kerbol'},],
                            # {'label': 'Sol', 'value': 'Sol'}],
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
                    html.Label('Starting parking altitude (km)'),
                    dcc.Input(id = 'startPark-input', value=100, 
                              type='number'),
                    html.Label('Ending parking altitude (km)'),
                    dcc.Input(id = 'endPark-input', value=100, 
                              type='number'),
                    dcc.Checklist(
                        id = 'noInsertion-checklist',
                        value = [],
                        options=[
                            {'label': 'No insertion burn', 'value': 'True'},
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
                        
                        html.H3('Custom Starting Orbit'),
                        html.Label('Semimajor axis (m)'),
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
                        html.Label('Semimajor axis (m)'),
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
                                      style={"white-space": "pre"})
                         ])
            ])
        ]),
    html.Div(className='four columns', children = [
        html.Div([
            html.H3('Orbit Plots'),
            dcc.Loading(id='transfer-loading', type='circle', children=[
                dcc.Markdown('**Transfer Trajectory**'),
                dcc.Graph(
                    id='transfer-graph',
                    figure = go.Figure(layout = dict(
                                        xaxis = dict(visible=False),
                                        yaxis = dict(visible=False)))
                    ),
                dcc.Markdown('**Ejection Trajectory**'),
                dcc.Graph(
                    id='ejection-graph',
                    figure = go.Figure(layout = dict(
                                        xaxis = dict(visible=False),
                                        yaxis = dict(visible=False)))
                    ),
                ]),
            ]),
        ]),
    # Hidden Divs to store data
    html.Div(id='dateFormat-div', style = {'display': 'none'}),
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
    Output('startingBody-dropdown', 'options'),
    [Input('system-div', 'children'),
     Input('system-radio', 'value')]
    )
def set_startBody_options(system_data, selected_system):
    return [{'label': i, 'value': i} 
            for i in start_body_options[selected_system]]

@app.callback(
    Output('endingBody-dropdown', 'options'),
    [Input('startingBody-dropdown', 'value')]
    )
def set_endBody_options(selected_body):
    return [{'label': i, 'value': i} for i in end_body_options[selected_body]]

@app.callback(
    Output('starta-input', 'value'),
    [Input('system-div', 'children'),
     Input('startingBody-dropdown', 'value'),
     Input('startPark-input', 'value')]
    )
def update_start_a(system_data, start_body_name, park_alt):
    system_data_d = jsonpickle.decode(system_data)
    start_body = [x for x in system_data_d if x.name == start_body_name][0]
    return start_body.eqr + 1000*park_alt

@app.callback(
    Output('enda-input', 'value'),
    [Input('system-div', 'children'),
     Input('endingBody-dropdown', 'value'),
     Input('endPark-input', 'value')]
    )
def update_end_a(system_data, end_body_name, park_alt):
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
     State('longFlightDays-input','value')]
    )
def update_porkchop_data(nClicks, system, dateFormat,
                          transferType, noInsertion,
                          startBodyName, startA, startEcc, startInc,
                          startArgP, startLAN, startMo, startEpoch,
                          endBodyName, endA, endEcc, endInc,
                          endArgP, endLAN, endMo, endEpoch,
                          minStartYear, minStartDay, maxStartYear, maxStartDay,
                          minFlightDays, maxFlightDays):
    
    # return empty plot on page load
    if nClicks == 0:
        return dash.no_update
    
    # prepare starting and ending orbits
    system = jsonpickle.decode(system)
    sBody = [x for x in system if x.name == startBodyName][0]
    sOrb = Orbit(startA, startEcc, startInc*math.pi/180, startArgP*math.pi/180,
                 startLAN*math.pi/180, startMo, sBody, startEpoch)
    eBody = [x for x in system if x.name == endBodyName][0]
    eOrb = Orbit(endA, endEcc, endInc*math.pi/180, endArgP*math.pi/180,
                 endLAN*math.pi/180, endMo, eBody, endEpoch)
    
    # grab day and year formats
    day = dateFormat['day']         # hours per day
    year= dateFormat['year']        # days per year
    
    # prepare start and flight time bounds
    minStartTime = 3600*((minStartDay-1) * day +                            \
                    (minStartYear-1) * day * year);
    if not (maxStartDay is None or maxStartYear is None):
        maxStartTime = 3600*((maxStartDay-1) * day +                        \
                    (maxStartYear-1) * day * year);
    else:
        maxStartTime = None
    if not (minFlightDays is None):
        minFlightTime = 3600*(minFlightDays * day)
    else:
        minFlightTime = None
    if not (maxFlightDays is None):
        maxFlightTime = 3600*(maxFlightDays* day)
    else:
        maxFlightTime = None
    
    # change noInsertion to a boolean
    if not noInsertion:
        noInsertion = False
    else:
        noInsertion = True
    
    # prepare porkchop table
    porkTable = PorkchopTable(sOrb, eOrb, transferType, noInsertion,
                              None, minStartTime, maxStartTime, 
                              minFlightTime, maxFlightTime, 51, 51)
    
    return jsonpickle.encode(porkTable)

@app.callback(
    Output('transfer-div','children'),
    [Input('porkchop-div','children'),
     Input('porkchop-graph','clickData')],
    [State('dateFormat-div','children')]
    )
def update_chosen_tranfser(porkTable, clickData, dateFormat):
    
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update
    
    porkTable = jsonpickle.decode(porkTable)
    
    # if the update comes from the porkchop data, get the best transfer
    if ctx.triggered[0]['prop_id'].split('.')[0] == 'porkchop-div':
        return jsonpickle.encode(porkTable.get_best_transfer())
    
    # if the update comes from user click data, get the chosen transfer
    elif ctx.triggered[0]['prop_id'].split('.')[0] == 'porkchop-graph':
        day = dateFormat['day']         # hours per day
        startDays = clickData['points'][0]['x']
        flightDays = clickData['points'][0]['y']
        startTime = startDays * 3600 * day
        flightTime = flightDays * 3600 * day
        transfer = porkTable.get_chosen_transfer (startTime, flightTime)
        return jsonpickle.encode(transfer)

@app.callback(
    Output('porkchop-graph','figure'),
    [Input('porkchop-div','children'),
     Input('transfer-div','children'),
     Input('dateFormat-div','children')],
    [State('porkchop-graph','figure')]
    )
def update_porkchop_plot(porkTable, chosenTransfer, dateFormat, prevState):
    
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
    # prep plot data
    bs = 1.1        # exponent base for contour levels
    numlvls = 20
    dV = porkTable.deltaV
    minDV = porkTable.get_best_transfer().get_total_delta_V()
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
                            title='Total Δv required (m/s)',),
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
        margin=dict(l=0, r=0, t=15, b=30),
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
    [Output('transferDetails-div','style'),
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
     Output('ejectionOrbit-markdown','children'),],
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
               dash.no_update;
               
    chosenTransfer = jsonpickle.decode(chosenTransfer)
    # grab day and year formats
    day = dateFormat['day']         # hours per day
    year= dateFormat['year']        # days per year
    
    # transfer trajectory details
    transferStyle = None
    departureTime = chosenTransfer.get_departure_burn_time()
    departureString = '**Departure:** Year ' +                              \
        str(math.floor(departureTime/(3600*day*year)+1)) +                  \
        ', Day ' +                                                          \
        str(math.floor(departureTime%(3600*day*year)/(day*3600)+1)) +       \
        ', ' +                                                              \
        str(math.floor((departureTime%(3600*day))/3600)) + ':' +            \
        str(math.floor(((departureTime%(3600*day))%3600)/60)) + ':' +       \
        str(math.floor(((departureTime%(3600*day))%3600)%60));
    arrivalTime = chosenTransfer.get_arrival_burn_time()
    arrivalString = '**Arrival:** Year ' +                                  \
        str(math.floor(arrivalTime/(3600*day*year)+1)) +                    \
        ', Day ' +                                                          \
        str(math.floor(arrivalTime%(3600*day*year)/(day*3600)+1)) +         \
        ', ' +                                                              \
        str(math.floor((arrivalTime%(3600*day))/3600)) + ':' +              \
        str(math.floor(((arrivalTime%(3600*day))%3600)/60)) + ':' +         \
        str(math.floor(((arrivalTime%(3600*day))%3600)%60));
    flightTime = chosenTransfer.get_arrival_burn_time() -                   \
        chosenTransfer.get_departure_burn_time()
    flightTimeString = '**Flight Duration:** ' +                            \
        str(math.floor(flightTime/(3600*day))+1) +                          \
            ' days, ' +                                                     \
        str(math.floor((flightTime%(3600*day))/3600)) + ':' +               \
        str(math.floor(((flightTime%(3600*day))%3600)/60)) + ':' +          \
        str(math.floor(((flightTime%(3600*day))%3600)%60));
    phaseString = '**Phase Angle:** ' +                                     \
        "{:.2f}".format(chosenTransfer.phaseAngle*180/math.pi) + '°';
    totalDVString = '**Total Δv:** ' +                                      \
        "{:.1f}".format(chosenTransfer.get_total_delta_V()) + ' m/s';
    planeDepartureDV=chosenTransfer.startOrbit.from_primary_to_orbit_bases( \
                        chosenTransfer.ejectionDV)
    departureDVString = '**Departure Burn:** ' +                            \
        "{:.1f}".format(norm([planeDepartureDV[0],planeDepartureDV[1]])) +  \
        ' m/s prograde, ' +                                                 \
        "{:.1f}".format(planeDepartureDV[2]) +                              \
        ' m/s normal';
    arrivalDVString = '**Arrival Burn:** ' +                                \
        "{:.1f}".format(chosenTransfer.insertionDV) + ' m/s';
    transferOrbitString = '**Transfer Orbit:**\n' +                         \
        str(chosenTransfer.transferOrbit);
    # plane change details
    if chosenTransfer.planeChange is True:
        planeChangeStyle = None
        planeChangeDVString = '**Plane Change Δv:** ' +                     \
            "{:.1f}".format(norm(chosenTransfer.planeChangeDV)) + ' m/s';
        planeChangeTime = chosenTransfer.startTime +                        \
            chosenTransfer.planeChangeDT
        planeChangeTimeString = '**Plane Change Time:** Year ' +            \
            str(math.floor(planeChangeTime/(3600*day*year)+1)) +            \
            ', Day ' +                                                      \
            str(math.floor(planeChangeTime%(3600*day*year)/(day*3600)+1)) + \
            ', ' +                                                          \
            str(math.floor((planeChangeTime%(3600*day))/3600)) + ':' +      \
            str(math.floor(((planeChangeTime%(3600*day))%3600)/60)) + ':' + \
            str(math.floor(((planeChangeTime%(3600*day))%3600)%60));
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
        str(math.floor((escapeTime%(3600*day))/3600)) + ':' +               \
        str(math.floor(((escapeTime%(3600*day))%3600)/60)) + ':' +          \
        str(math.floor(((escapeTime%(3600*day))%3600)%60));
        ejectionOrbitString = '**Ejection Orbit:**\n' +                     \
            str(chosenTransfer.ejectionTrajectory)
    else:
        ejectionStyle = {'display': 'none'}
        ejectionAngleString = ''
        escapeTimeString = ''
        ejectionOrbitString = ''
    # insertion details
    if not chosenTransfer.insertionDT == 0:
        insertionStyle = None
        encounterTime = chosenTransfer.startTime + chosenTransfer.flightTime
        encounterTimeString = '**Arrival SOI Encounter:**\t\tYear ' +       \
        str(math.floor(encounterTime/(3600*day*year)+1)) +                  \
        ', Day ' +                                                          \
        str(math.floor(encounterTime%(3600*day*year)/(day*3600)+1)) +       \
        ', ' +                                                              \
        str(math.floor((encounterTime%(3600*day))/3600)) + ':' +            \
        str(math.floor(((encounterTime%(3600*day))%3600)/60)) + ':' +       \
        str(math.floor(((encounterTime%(3600*day))%3600)%60));
    else:
        insertionStyle = {'display': 'none'}
        encounterTimeString = ''
    
    
    return transferStyle, departureString, arrivalString, flightTimeString, \
            phaseString, totalDVString, departureDVString, arrivalDVString, \
            planeChangeStyle, planeChangeDVString, planeChangeTimeString,   \
            ejectionStyle, ejectionAngleString, escapeTimeString,           \
            insertionStyle, encounterTimeString,                            \
            transferOrbitString, transferOrbitPCString, ejectionOrbitString;
            
            

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
    bdTimes = np.linspace(chosenTransfer.get_departure_burn_time(),         \
                          chosenTransfer.startTime +                        \
                          chosenTransfer.flightTime, 501);
    if chosenTransfer.planeChange is True:
        trTimes = np.linspace(chosenTransfer.startTime,                     \
                                chosenTransfer.startTime +                  \
                                chosenTransfer.planeChangeDT, 501);
        trPCTimes = np.linspace(chosenTransfer.startTime +                  \
                                chosenTransfer.planeChangeDT,               \
                                chosenTransfer.startTime +                  \
                                chosenTransfer.flightTime, 501);
    else:
        trTimes = np.linspace(chosenTransfer.startTime,                     \
                                chosenTransfer.startTime +                  \
                                chosenTransfer.flightTime, 501);
    
    maxVals = [add_orbit(fig, chosenTransfer.transferOrbit, trTimes,
                          dateFormat, name = 'Transfer')]
    if chosenTransfer.planeChange is True:
        maxVals = np.append(maxVals,
                            add_orbit(fig, chosenTransfer.transferOrbitPC,  \
                                      trPCTimes, dateFormat,                \
                                      name = 'Transfer (plane change)'));
    
    for bd in chosenTransfer.transferOrbit.prim.satellites:
        maxVals = np.append(maxVals,
                            add_orbit(fig, bd.orb, bdTimes, dateFormat,
                                      bd.color, bd.name)
                            )
    
    add_primary(fig, chosenTransfer.transferOrbit.prim)
    lim = np.amax(maxVals)*1.25
    add_transfer_phase_angle(fig, chosenTransfer,
                             1.5*chosenTransfer.transferOrbit.a)
    
    fig.update_layout(
        margin=dict(l=0, r=0, t=15, b=30),
        paper_bgcolor="rgb(50, 50, 50)",
        scene = dict(
            xaxis = dict(nticks = 0, range=[-lim, lim],
                          backgroundcolor="rgb(50, 50, 50)",
                          showbackground=True,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            yaxis = dict(nticks = 0, range=[-lim, lim],
                          backgroundcolor="rgb(50, 50, 50)",
                          showbackground=True,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            zaxis = dict(nticks = 0, range=[-lim, lim],
                          backgroundcolor="rgb(50, 50, 50)",
                          showbackground=True,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            camera = dict(
                eye = dict( x=1.5*chosenTransfer.transferOrbit.a/lim,
                            y=1.5*chosenTransfer.transferOrbit.a/lim,
                            z=1.5*chosenTransfer.transferOrbit.a/lim)
                )
            )
        )
    return fig

@app.callback(
    Output('ejection-graph', 'figure'),
    [Input('transfer-div', 'children'),
     Input('dateFormat-div', 'children')]
    )
def update_ejection_plot(chosen_transfer, dateFormat):
    fig = go.Figure(layout = dict(xaxis = dict(visible=False),
                                  yaxis = dict(visible=False)))
    if chosen_transfer is None:
        return fig
    
    chosen_transfer = jsonpickle.decode(chosen_transfer)
    if chosen_transfer.ejectionTrajectory is None:
        return fig
    
    parkTimes = np.linspace(chosen_transfer.get_departure_burn_time(),
                            chosen_transfer.get_departure_burn_time()
                            + chosen_transfer.startOrbit.get_period(),
                            101)
    
    ejTimes = np.linspace(chosen_transfer.get_departure_burn_time(),
                            chosen_transfer.startTime,
                            1001)
    
    maxVals =[add_orbit(fig, chosen_transfer.ejectionTrajectory, ejTimes,
                        dateFormat, name = 'Ejection')]
    
    maxVals = np.append(maxVals,
                          add_orbit(fig, chosen_transfer.startOrbit, 
                                    parkTimes, dateFormat,
                                    name = 'Starting Orbit', style = 'dot',
                                    fade = False)
                          )
    
    for bd in chosen_transfer.ejectionTrajectory.prim.satellites:
        maxVals = np.append(maxVals,
                            add_orbit(fig, bd.orb, ejTimes, dateFormat,
                                      bd.color, bd.name)
                            )
    
    lim = np.amax(maxVals)*1.25
    add_ejection_angle(fig, chosen_transfer)
    add_prograde_trace(fig, chosen_transfer, 
                       ejTimes - chosen_transfer.ejectionDT/2)
    add_primary(fig, chosen_transfer.ejectionTrajectory.prim)
    
    fig.update_layout(
        margin=dict(l=0, r=0, t=15, b=30),
        paper_bgcolor="rgb(50, 50, 50)",
        scene = dict(
            xaxis = dict(nticks = 0, range=[-lim, lim],
                          backgroundcolor="rgb(50, 50, 50)",
                          showbackground=True,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            yaxis = dict(nticks = 0, range=[-lim, lim],
                          backgroundcolor="rgb(50, 50, 50)",
                          showbackground=True,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            zaxis = dict(nticks = 0, range=[-lim, lim],
                          backgroundcolor="rgb(50, 50, 50)",
                          showbackground=True,
                          showgrid=False,
                          zerolinecolor="rgb(60, 60, 60)",
                          ticks='',
                          showticklabels=False,),
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            camera = dict(
                eye = dict( x=(1.5*chosen_transfer.startOrbit.a/lim),
                            y=(1.5*chosen_transfer.startOrbit.a/lim),
                            z=(1.5*chosen_transfer.startOrbit.a/lim))
                )
            )
        )
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)

