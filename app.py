# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import pickle
import math
import numpy as np
from numpy.linalg import norm
from copy import copy
from orbit import Orbit
from body import Body
from transfer import Transfer
from prktable import PorkchopTable
import plotly.graph_objects as go
from plotly.subplots import make_subplots
#%% initialize app?

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

#%% read solar system data
infile = open('kerbol_system','rb')
kerbol_system = pickle.load(infile)
infile.close
kerbol_system[0].orb.prim = kerbol_system[0]

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

start_body_options = {'Kerbol': start_bodies,
                'Sol': ['Mercury','Venus','Earth']}

end_body_options = dict()
for x,bd in enumerate(kerbol_system):
    end_body_options[bd.name] = eb[x]

#%% generate plot data
# user input
startName = 'Kerbin'
endName = 'Duna'
startBody = [x for x in kerbol_system if x.name == startName][0]
endBody = [x for x in kerbol_system if x.name == endName][0]
trType = 'ballistic'

# get porkchop plot table
sOrb = Orbit(100000 + startBody.eqr, 0, 0, 0, 0, 0, startBody)
eOrb = Orbit(100000 + endBody.eqr, 0, 0, 0, 0, 0, endBody)
porkTable = PorkchopTable(sOrb,eOrb, transferType = trType)
bestTrs = porkTable.get_best_transfer()

# prep plot data
bs = 1.1
numlvls = 20
st = porkTable.startTimes
ft = porkTable.flightTimes
dV = porkTable.deltaV
minDV = bestTrs.get_total_delta_V()
lvls = minDV * bs**(np.arange(start=0,stop=numlvls,step=4))

logdV = np.log(dV)/np.log(bs)
logLvls = np.log(lvls)/np.log(bs)
color_vals = logLvls
color_labels = [str(math.floor(bs**i)) for i in color_vals]

#%% make porkchop plot
fig = make_subplots(
    rows=1, cols=1,
    specs=[[{"type": "contour"}]])

fig.add_trace(
                go.Contour(
                    z = logdV,
                    x = st/(6*3600),
                    y = ft/(6*3600),
                    contours = dict(
                        start = (logLvls[0]-0.99),
                        end = (logLvls[-1]),
                        size = 1),
                    colorscale = 'Jet',
                    colorbar = dict(
                        tickvals = color_vals,
                        ticktext = color_labels,
                        title='Total Δv required (m/s)',),
                    contours_coloring='heatmap',
                    customdata = bs**logdV,
                    hovertemplate = "Δv = " +
                                    "%{customdata:.2f}" +
                                    "m/s" +
                                    "<extra></extra>"
                    ))

fig.update_xaxes(title_text='Transfer start (day #)', row=1, col=1)
fig.update_yaxes(title_text='Transfer duration (days)', row=1, col=1)
fig.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),
    )

#%% generate orbit plots

def fade_color(color, div = 2):
    """Divides each element of the tuple by the specified number."""
    
    return tuple(math.floor(c/2) for c in color)

def add_orbit(figure, row, col, orb, times, color = (255,255,255), name = '',
              style = 'solid', fade = True):
    """
    
    """
    
    if fade:
        fadedColor = fade_color(color)
    else:
        fadedColor = color
    
    pos = np.transpose(orb.get_grandparent_positions(times = times))
    maxVal = np.amax(pos)
    
    figure.add_trace(go.Scatter3d(
        x = pos[0],
        y = pos[1],
        z = pos[2],
        customdata = np.stack((norm(pos, axis = 0)/1000,
                           np.floor(times/(60**2*6*426))+1,
                           np.floor(times%(60**2*6*426)/(6*60**2)+1),
                           np.floor((times%(60**2*6))/60**2),
                           np.floor(((times%(60**2*6))%60**2)/60),
                           np.floor(((times%(60**2*6))%60**2)%60)),
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
        ),
            row = row,
            col = col)
    
    return maxVal

def add_primary(figure,row, col, bd):
    """
    
    """
    
    figure.add_trace(go.Scatter3d(
                                  x = np.array([0]),
                                  y = np.array([0]),
                                  z = np.array([0]),
                                  mode = "markers",
                                  marker = dict(color = 'rgb'+str(bd.color),),
                                  name = bd.name,
                                  showlegend = False,
                                  hovertemplate = "Central body"
                                  ),
        row = row,
        col = col
        )
    
def add_transfer_phase_angle(figure, row, col, transfer, r = None):
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
        arcAngles = np.linspace(startAngle, endAngle,101)
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
        arcAngles = np.linspace(startAngle, endAngle,101)
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
        ),
            row = row,
            col = col)
    
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

def add_ejection_angle(figure, row, col, transfer, r = None):
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
        arcAngles =  np.linspace(startAngle, endAngle,101)
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
            ),
                row = row,
                col = col)
    
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

def add_body_vectors(figure, row, col, transfer, r = None):
    """
    
    """
    

bdTimes = np.linspace(bestTrs.get_departure_burn_time(),
                       bestTrs.startTime + bestTrs.flightTime, 501)
trTimes = np.linspace(bestTrs.startTime,
                       bestTrs.startTime + bestTrs.flightTime, 501)

fig2 = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    print_grid=False)

maxVals1 = [add_orbit(fig2,1,1, bestTrs.transferOrbit, trTimes,
                     name = 'Transfer')]

for bd in bestTrs.transferOrbit.prim.satellites:
    maxVals1 = np.append(maxVals1,
                        add_orbit(fig2,1,1, bd.orb, bdTimes, bd.color,
                                  bd.name)
                        )

add_primary(fig2,1,1,bestTrs.transferOrbit.prim)
lim1 = np.amax(maxVals1)*1.25
add_transfer_phase_angle(fig2, 1, 1, bestTrs, 1.5*bestTrs.transferOrbit.a)

fig2.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),
    paper_bgcolor="rgb(50, 50, 50)",
    scene = dict(
        xaxis = dict(nticks = 0, range=[-lim1,lim1],
                     backgroundcolor="rgb(50, 50, 50)",
                     showbackground=True,
                     showgrid=False,
                     zerolinecolor="rgb(60, 60, 60)",
                     ticks='',
                     showticklabels=False,),
        yaxis = dict(nticks = 0, range=[-lim1,lim1],
                     backgroundcolor="rgb(50, 50, 50)",
                     showbackground=True,
                     showgrid=False,
                     zerolinecolor="rgb(60, 60, 60)",
                     ticks='',
                     showticklabels=False,),
        zaxis = dict(nticks = 0, range=[-lim1,lim1],
                     backgroundcolor="rgb(50, 50, 50)",
                     showbackground=True,
                     showgrid=False,
                     zerolinecolor="rgb(60, 60, 60)",
                     ticks='',
                     showticklabels=False,),
        xaxis_title='',
        yaxis_title='',
        zaxis_title='',
        # camera = dict(
        #     eye = dict(x=0,y=0,z=2*bestTrs.transferOrbit.a/lim1)
        #     )
        )
    )

# Ejection Plot
fig3 = make_subplots(rows=1, cols=1,
                    specs=[[{'is_3d': True}]],
                    print_grid=False)

parkTimes = np.linspace(bestTrs.get_departure_burn_time(),
                        bestTrs.get_departure_burn_time()
                        + bestTrs.startOrbit.get_period(),
                        101)

ejTimes = np.linspace(bestTrs.get_departure_burn_time(),
                       bestTrs.startTime,
                       1001)

maxVals2 =[ add_orbit(fig3,1,1, bestTrs.ejectionTrajectory, ejTimes,
                     name = 'Ejection')]

maxVals2 = np.append(maxVals2,
                     add_orbit(fig3,1,1, bestTrs.startOrbit, parkTimes,
                               name = 'Starting Orbit', style = 'dot',
                               fade = False)
                     )

for bd in bestTrs.ejectionTrajectory.prim.satellites:
    maxVals2 = np.append(maxVals2,
                        add_orbit(fig3,1,1, bd.orb, ejTimes, bd.color,
                                  bd.name)
                        )

add_primary(fig3,1,1,bestTrs.ejectionTrajectory.prim)
lim2 = np.amax(maxVals2)*1.25
add_ejection_angle(fig3, 1, 1, bestTrs)

fig3.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),
    paper_bgcolor="rgb(50, 50, 50)",
    scene = dict(
        xaxis = dict(nticks = 0, range=[-lim2,lim2],
                     backgroundcolor="rgb(50, 50, 50)",
                     showbackground=True,
                     showgrid=False,
                     zerolinecolor="rgb(60, 60, 60)",
                     ticks='',
                     showticklabels=False,),
        yaxis = dict(nticks = 0, range=[-lim2,lim2],
                     backgroundcolor="rgb(50, 50, 50)",
                     showbackground=True,
                     showgrid=False,
                     zerolinecolor="rgb(60, 60, 60)",
                     ticks='',
                     showticklabels=False,),
        zaxis = dict(nticks = 0, range=[-lim2,lim2],
                     backgroundcolor="rgb(50, 50, 50)",
                     showbackground=True,
                     showgrid=False,
                     zerolinecolor="rgb(60, 60, 60)",
                     ticks='',
                     showticklabels=False,),
        xaxis_title='',
        yaxis_title='',
        zaxis_title='',
        # camera = dict(
        #     eye = dict(x=0,y=0,z=(2*sOrb.a/lim2))
        #     )
        )
    )



#%% set up app
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
                        options=[{'label': k, 'value': k}                   \
                                 for k in start_body_options.keys()],
                        value='Kerbol',
                        ),
                    html.Label('Starting Body'),
                    dcc.Dropdown(
                        id = 'startingBody-dropdown',
                        ),
                    html.Label('Ending Body'),
                    dcc.Dropdown(
                        id = 'endingBody-dropdown',
                        ),
                    html.Label('Starting parking altitude (km)'),
                    dcc.Input(id = 'startPark-input', value='100', 
                              type='number'),
                    html.Label('Ending parking altitude (km)'),
                    dcc.Input(id = 'endPark-input', value='100', 
                              type='number'),
                    dcc.Checklist(
                        id = 'noInsertion-checkbox',
                        options=[
                            {'label': 'No insertion burn', 'value': 'n'},
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
                    dcc.Input(id = 'earlyStartYear-dropdown', value=1, 
                              type='number'),
                    html.Label('Earliest Departure Day'),
                    dcc.Input(id = 'earlyStartDay-dropdown', value=1, 
                              type='number'),
                    ])
                ),
                dcc.Tab(
                    label='Advanced Settings',
                    value = 'adv', 
                    children = html.Div(className='six columns', children = [
                        html.H3('Custom Starting Orbit'),
                        html.Label('Semimajor axis (m)'),
                        dcc.Input(id = 'starta-dropdown',  
                                  type='number'),
                        html.Label('Eccentricity'),
                        dcc.Input(id = 'startecc-dropdown',  
                                  type='number'),
                        html.Label('Inclination (rad)'),
                        dcc.Input(id = 'startinc-dropdown',  
                                  type='number'),
                        html.Label('Argument of the Periapsis (rad)'),
                        dcc.Input(id = 'startargp-dropdown',  
                                  type='number'),
                        html.Label('Longitude of the Ascending Node (rad)'),
                        dcc.Input(id = 'startlan-dropdown',  
                                  type='number'),
                        html.Label('Mean anomaly at epoch (rad)'),
                        dcc.Input(id = 'startmo-dropdown',  
                                  type='number'),
                        
                        html.H3('Custom Ending Orbit'),
                        html.Label('Semimajor axis (m)'),
                        dcc.Input(id = 'enda-dropdown',  
                                  type='number'),
                        html.Label('Eccentricity'),
                        dcc.Input(id = 'endecc-dropdown',  
                                  type='number'),
                        html.Label('Inclination (rad)'),
                        dcc.Input(id = 'endinc-dropdown',  
                                  type='number'),
                        html.Label('Argument of the Periapsis (rad)'),
                        dcc.Input(id = 'endargp-dropdown',  
                                  type='number'),
                        html.Label('Longitude of the Ascending Node (rad)'),
                        dcc.Input(id = 'endlan-dropdown',  
                                  type='number'),
                        html.Label('Mean anomaly at epoch (rad)'),
                        dcc.Input(id = 'endmo-dropdown',  
                                  type='number'),
                        ])
                    )
            ]),
        ]),
    html.Div(className='four columns', children = [
        html.H3('Porkchop Plot'),
        html.Button('Plot!',
                    className = 'button-primary',
                    id = 'porkchop-button',
                    type='button',
            ),
        html.Div([
            dcc.Graph(
                id='porkchop-plot',
                figure=fig
                ),
            ]),
        ]),
    html.Div(className='four columns', children = [
        html.H3('Transfer Plot'),
        html.Div([
            dcc.Graph(
                id='trs-plots',
                figure=fig2
                ),
            ]),
        html.H3('Ejection Plot'),
        html.Div([
            dcc.Graph(
                id='ej-plots',
                figure=fig3
                ),
            ]),
        ]),
    ])

@app.callback(
    Output('startingBody-dropdown', 'options'),
    [Input('system-radio', 'value')])
def set_startBody_options(selected_system):
    return [{'label': i, 'value': i} for i in start_body_options[selected_system]]

@app.callback(
    Output('endingBody-dropdown', 'options'),
    [Input('startingBody-dropdown', 'value')])
def set_endBody_options(selected_body):
    return [{'label': i, 'value': i} for i in end_body_options[selected_body]]


if __name__ == '__main__':
    app.run_server(debug=False)
