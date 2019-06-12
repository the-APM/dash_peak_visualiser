#!/usr/bin/env python3

from dash import Dash
from dash.dependencies import Input, Output, State
from scipy.signal import savgol_filter
from scipy import sparse
from scipy.sparse.linalg import spsolve
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import numpy as np
import json
import base64
import urllib
import pandas as pd

##<-fuse

def read_file(filename):
    """Reads the given file as loads its contents as json data"""
    with open(filename) as json_load:
        json_data = json.load(json_load)
        return json_data
    
    
def read_remote_file(file_url):
    """Reads a remote file and returns contents as json data.
    The given URL must point to a JSON file."""
    json_data = json.loads(urllib.request.urlopen(file_url).read())
    return json_data
    
    
def process_json_data(json_data):
    """Extracts out the number of groups and samples
    from the given json structure"""
    if not json_data:
        return []
    groups = json_data['groups']
    sample_names = []
    prep_data = {}
    prep_data_bkp ={}
    cmp_list=[]
    d2 = pd.read_csv("kibbeylab_cpd_db.csv")
    for group in groups:
        print(group['compound']['compoundId'],'a')
        if group['compound']['compoundId'] not in cmp_list:
            cmp_list.append(group['compound']['compoundId'])
        cmpd = group['compound']['compoundId'].split()
        frag = cmpd[len(cmpd)-1]
        cmpd = cmpd[0]
        peaks = group['peaks']
        expectedRt = list(d2["expectedRt"])[list(d2["compound"]).index(group['compound']['compoundId'])]
        f = {}
        a = {}
        for peak in peaks:
            rt_list = peak["eic"]["rt"]
            intensity_list = peak["eic"]["intensity"]
            rt_step = rt_list[2]-rt_list[0]
            rt_window_intensities = []
            index_list = []
            for i in range(len(rt_list)):
                if rt_list[i]<(expectedRt+rt_step) and rt_list[i]>(expectedRt-rt_step):
                    rt_window_intensities.append(intensity_list[i])
                    index_list.append(i)
            if max(rt_window_intensities)>0:
                pk_list = [index for index, value in enumerate(intensity_list) if value==max(rt_window_intensities)]
                for i in pk_list:
                    if i in index_list:
                        pk=i
                clipped_pk = intensity_list[max(0,pk-5):min(len(intensity_list)-1,pk+6)]
                clipped_rt=rt_list[max(0,pk-5):min(len(rt_list)-1,pk+6)]
                m=max(clipped_pk)
                f[peak["sampleName"]] = [clipped_pk,clipped_rt]
                if peak["sampleName"] not in sample_names:
                    sample_names.append(peak["sampleName"])
        a[frag] = f
        if cmpd.upper() not in prep_data:
            prep_data[cmpd.upper()] = {}
            prep_data_bkp[cmpd.upper()]={}
        prep_data[cmpd.upper()].update(a)
        prep_data_bkp[cmpd.upper()].update(a)

    parent_list={}
    for key in prep_data.keys():
        parent_list[key]=[]
    for i in range(len(cmp_list)):
        print(cmp_list[i],'b')
        z=cmp_list[i].split()
        if z[0].upper() in parent_list.keys():
            if parent_list[z[0].upper()]==[]:
                parent_list[z[0].upper()]=[z[len(z)-1]]
            else:
                print(parent_list[z[0].upper()][0],'c')
                prec,prod=parent_list[z[0].upper()][0].split('/')
                print(z[len(z)-1],'d')
                cur_prec,cur_prod=z[len(z)-1].split('/')
                if prec==cur_prec:
                    parent_list[z[0].upper()].append(z[len(z)-1])
                elif prec>cur_prec:
                    parent_list[z[0].upper()]=[z[len(z)-1]]
    return prep_data


class GlobalInstance(object):
    _json_data = None
    prep_data = None
    filename = ""

    @property
    def json_data(self):
        return self._json_data

    @json_data.setter
    def json_data(self, data):
        self._json_data = data
        self.prep_data = process_json_data(data)

instance = GlobalInstance()

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)
app.config['suppress_callback_exceptions']=True
app.layout = html.Div([

    dcc.Tabs(id="tabs", value='tab-1', children=[
        dcc.Tab(label="PEAK VISUALIZER", value='tab-1')
    ],
    style={'marginBottom': '6px'}),

    html.Div(id='tabs-content'),

], style={'width': '92%', 'marginLeft': '4%', 'display': 'inline-block'})

server = app.server

spacer_top_small = html.Div([], style={'marginTop': 6,})
spacer_top_large = html.Div([], style={'marginTop': 24,})
spacer_bottom_small = html.Div([], style={'marginBottom': 6})
spacer_bottom_large = html.Div([], style={'marginBottom': 24})


def eic_content(tab):
    return html.Div([

        spacer_top_small,

        dcc.Graph(
            id='eic',
            style={'display': 'inline-block'}
        ),


        spacer_top_small,

        html.Div(style={'marginBottom': '24px', 'marginTop': '24px'})
    ])

loaded_eics_content = html.Div([

    html.Div([

        html.Div([

            dcc.Upload(
                id='json-loader',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select JSON File'),
                    ' exported from El-MAVEN'
                ]),
                style={
                    'width': '100%',
                    'height': '76px',
                    'lineHeight': '76px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'marginBottom': '6px'
                },
            ),

        ])
    ], style={'width': '48%', 'display': 'inline-block'}),

    html.Div(
        id='post-upload-dropdowns', 
        style={'width': '48%', 'float': 'right', 'display': 'inline-block'}
    )

], style={'marginBottom': '18px'})


@app.callback(Output('tabs-content', 'children'),
              [Input('tabs', 'value')])
def render_content(tab):
    if tab == 'tab-1':
        return html.Div([
            loaded_eics_content,
            eic_content(tab)
        ])


def parse_contents(contents, filename, date):
    _, content_string = contents.split(',')

    utf8_content = base64.b64decode(content_string)
    try:
        if 'json' in filename:
            # Assume that the user uploaded a JSON file
            instance.json_data = json.loads(utf8_content.decode('utf-8'))
            instance.filename = filename
      ##      instance.prep_data = process_json_data(instance._json_data)
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([

        dcc.Dropdown(
            id='compound-dropdown',
            options = [{'label': compound, 'value': compound} for compound in list(instance.prep_data.keys())],
            value = list(instance.prep_data.keys())[0]
        ),

        spacer_bottom_small,

        dcc.Dropdown(
            id='fragment-dropdown',
            options = [{'label': frags, 'value': frags} for frags in list(instance.prep_data[compound-dropdown['value']].keys())],
            value = list(instance.prep_data[compound-dropdown['value']].keys())[0]
        ),
        
        spacer_bottom_small,

        dcc.Dropdown(
            id='sample-dropdown',
            options = [{'label': sample, 'value': sample} for sample in list(instance.prep_data[compound-dropdown['value']][fragment-dropdown['value']].keys())],
            value = list(instance.prep_data[compound-dropdown['value']][fragment-dropdown['value']].keys())[0]
        )])


@app.callback(Output('post-upload-dropdowns', 'children'),
              [Input('json-loader', 'contents')],
              [State('json-loader', 'filename'),
               State('json-loader', 'last_modified')])
def update_output(contents, filename, last_modified):
    if contents is not None:
        return [parse_contents(contents, filename, last_modified)]

if __name__ == '__main__':
    app.run_server(debug=True)
