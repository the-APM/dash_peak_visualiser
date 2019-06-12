# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import json,math
import pandas as pd
import numpy as np
from scipy import spatial
import plotly.graph_objs as go


data = json.load(open("Study_1_peaks_triangular.json",'rb'))
groups = data['groups']
d2 = pd.read_csv("kibbeylab_cpd_db.csv")
sample_names = []
prep_data = {}
prep_data_bkp ={}
for group in groups:
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
			##for i in range(len(clipped_pk)):clipped_pk[i]=clipped_pk[i]/m
			f[peak["sampleName"]] = [clipped_pk,clipped_rt]
			##print(max(clipped_pk),frag,peak["sampleName"])
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
cmp_list=list(d2['compound'])
for i in range(len(cmp_list)):
        z=cmp_list[i].split()
        if z[0].upper() in parent_list.keys():
                if parent_list[z[0].upper()]==[]:
                        parent_list[z[0].upper()]=[z[len(z)-1]]
                else:
                        prec,prod=parent_list[z[0].upper()][0].split('/')
                        cur_prec,cur_prod=z[len(z)-1].split('/')
                        if prec==cur_prec:
                                parent_list[z[0].upper()].append(z[len(z)-1])
                        elif prec>cur_prec:
                                parent_list[z[0].upper()]=[z[len(z)-1]]

print("Prep_data done")

def scale_peak(p):
    mx=max(p)
    if mx==0:
        return p
    for i in range(len(p)):
        p[i]/=mx
    return p

def get_parents(compound):
    return parent_list[compound]

def get_peak_area(peak,step_size):
    return np.trapz(peak,dx=step_size)

def corrected_graph(compound,frag,sample):
    c=prep_data[compound]
    parent_frags=get_parents(compound)
    first_frag,second_frag=frag.split('/')
    for i in range(len(parent_frags)):
        first_parent_frag,second_parent_frag=parent_frags[i].split('/')
        if int(second_frag)-int(second_parent_frag)>=0 and int(second_frag)-int(second_parent_frag)<=10:
            parent_frag=parent_frags[i]
            break
    step_size=c[parent_frag][sample][1][1]-c[parent_frag][sample][1][0]
    comp_area=get_peak_area(c[frag][sample][0],step_size)
    if sample not in list(c[parent_frag].keys()):
        return c[frag][sample][1],c[frag][sample][0],'rejected',comp_area,comp_area,parent_frag,0
    mx=max(c[parent_frag][sample][0])
    if mx==0:
        return c[frag][sample][1],c[frag][sample][0],'rejected',comp_area,comp_area,parent_frag,0
    gval=[index for index, value in enumerate(c[parent_frag][sample][0]) if value==0]
    gval.append(len(c[parent_frag][sample][0])-1)
    gval.insert(0,0)
    mx_ind=c[parent_frag][sample][0].index(mx)
    for i in range(len(gval)):
        if gval[i]>mx_ind:
            b=gval[i]
            a=gval[i-1]
            break
        else:
            b=gval[i]+1
            a=gval[i-1]
    if b==len(c[parent_frag][sample][1]):
        b=b-1
    parent_peak_width=b-a
    parent_rt_span=c[parent_frag][sample][1][b]-c[parent_frag][sample][1][a]
    parent_peak=[0]*len(c[parent_frag][sample][1])
    parent_peak[a:b]=c[parent_frag][sample][0][a:b]
    parent_area=np.trapz(c[parent_frag][sample][0],dx=step_size)
    parent_peak=scale_peak(parent_peak)
    parent_peak_rt=c[parent_frag][sample][1][parent_peak.index(1)]
    if frag==parent_frag:
        return c[frag][sample][1],c[frag][sample][0],'accepted',parent_area,parent_area,parent_frag,parent_area
    else:
        ll=len(c[frag][sample][0])
        gvall = [index for index, value in enumerate(c[frag][sample][0]) if value==0]
        gvall.append(len(c[frag][sample][0])-1)
        gvall.insert(0,0)
        mxl=max(c[frag][sample][0])
        mx_indl=c[frag][sample][0].index(mxl)
        peak_area=get_peak_area(c[frag][sample][0],step_size)
        for i in range(len(gvall)):
            if gvall[i]>mx_indl:
                b=gvall[i]
                a=gvall[i-1]
                break
            else:
                b=gvall[i]+1
                a=gvall[i-1]
        if b==len(c[frag][sample][1]):
            b=b-1
        rt_span_loc=c[frag][sample][1][b]-c[frag][sample][1][a]
        if rt_span_loc<=parent_rt_span:
            t=[0]*ll
            t[a:b]=c[frag][sample][0][a:b]
            actual_area=get_peak_area(t,step_size)
            ts=scale_peak(t)
            peak_rt_loc=c[frag][sample][1][ts.index(max(ts))]
            cor=0
            cost=0
            for i in range(min(len(ts),len(parent_peak))):
                if ts[i]<parent_peak[i]-0.05:
                    cor+=1
                    cost+=parent_peak[i]-ts[i]-0.05
                    ts[i]=parent_peak[i]-0.05
                elif ts[i]>parent_peak[i]+0.05:
                    cor+=1
                    cost+=ts[i]-parent_peak[i]-0.05
                    ts[i]=parent_peak[i]+0.05
            if cor<=parent_peak_width/2 or cost<0.75:
                state="accepted"
                mx=max(c[frag][sample][0])
                pk=[x*mx for x in ts]
                corrected_area=get_peak_area(pk,step_size)
            else:
                state="rejected"
                pk=c[frag][sample][0]
                corrected_area=actual_area

        else:
            return c[frag][sample][1],c[frag][sample][0],'rejected',peak_area,peak_area,parent_frag,parent_area
    return c[frag][sample][1],pk,state,actual_area,corrected_area,parent_frag,parent_area

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    dcc.Dropdown(
        id='compound-dropdown',
        options=[{'label': k, 'value': k} for k in prep_data.keys()],
        value=list(prep_data.keys())[0]
    ),

    dcc.Dropdown(id='frag-dropdown'),

    dcc.Dropdown(
        id='sample-dropdown',
        multi=True
        ),

    dcc.Dropdown(
        id='mode-dropdown',
        options=[{'label': 'normal', 'value': 'normal'},{'label':'corrected','value':'corrected'}],
        value='normal',
        multi=True
    ),

    dcc.Graph(id='plot-area')

##    dcc.Checklist(
##            id='alts',
##            options=[
##                    {'label':'Scaled','value':'SC'},
##                    {'label':'Aligned','value':'AL'}
##                    ],
##            )
])


@app.callback(
    Output('frag-dropdown', 'options'),
    [Input('compound-dropdown', 'value')])
def set_frag_options(selected_compound):
    frag_list=list(prep_data[selected_compound].keys())
    frag_list.sort()
    return [{'label': i, 'value': i} for i in frag_list]


@app.callback(
    Output('frag-dropdown', 'value'),
    [Input('frag-dropdown', 'options')])
def set_frag_value(available_options):
    return available_options[0]['value']

@app.callback(
    Output('sample-dropdown', 'options'),
    [Input('compound-dropdown', 'value'),
     Input('frag-dropdown', 'value')])
def set_sample_options(selected_compound,selected_frag):
    sample_list=list(prep_data[selected_compound][selected_frag].keys())
    sample_list.sort()
    return [{'label': i, 'value': i} for i in sample_list]


@app.callback(
    Output('sample-dropdown', 'value'),
    [Input('sample-dropdown', 'options')])
def set_sample_value(available_options):
    return [available_options[0]['value']]

@app.callback(
        Output('plot-area', 'figure'),
    [
        Input('compound-dropdown', 'value'),
        Input('frag-dropdown', 'value'),
        Input('sample-dropdown','value'),
        Input('mode-dropdown','value')
##        Input('alts','values')
    ])
def set_display_children(compound,frag,samples,mode):
    dt=[]
    dtc=[]
    alters=[]
    for f in samples:
        t=prep_data[compound][frag][f][0]
        o=prep_data[compound][frag][f][1]
        print (t,o)
##        r=0
##        for i in range(len(alters)):
##                if alters[i]=="SC":
##                        t=scale_peak(t)
##                elif alters[i]=="AL":
##                        if r==0:
##                                r=1
##                                max_peak_rt=o[t.index(max(t))]
##                        else:
##                                max_rt_loc=o[t.index(max(t))]
##                                for j in range(len(o)):
##                                        o[j]=o[j]-max_rt_loc+max_peak_rt
        dt.append(go.Scatter(
                        x = o,
                        y = t,
                        mode = 'lines',
                        name = f
                        ))
        rt,ints,state,area,cor_area,parent,parent_area=corrected_graph(compound,frag,f)
##        r==0
##        for i in range(len(alters)):
##                if alters[i]=="SC":
##                        ints=scale_peak(ints)
##                elif alters[i]=="AL":
##                        if r==0:
##                                r=1
##                                max_peak_rt=rt[ints.index(max(ints))]
##                        else:
##                                max_rt_loc=rt[ints.index(max(ints))]
##                                for j in range(len(rt)):
##                                        rt[j]=rt[j]-max_rt_loc+max_peak_rt
        dtc.append(go.Scatter(
                        x = rt,
                        y = ints,
                        mode = 'lines',
                        name = f+"+"+state
                        ))
        print("compound: ",compound)
        print("frag: ",frag)
        print("sample: ",f)
        print("state: ",state)
        print("area: ",area)
        print("cor_area: ",cor_area)
        print("parent: ",parent)
        print("parent_area: ",parent_area,"\n")
        print(dt,dtc)
    if len(mode)==1:
        if(mode[0]=='normal'):
            return go.Figure(
                data=dt,
                layout=go.Layout(
                    title='Plot',
                    showlegend=True,
                    legend=go.layout.Legend(
                        x=0,
                        y=1.0
                    ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30)
                )
            )
        else:
            return go.Figure(
                data=dtc,
                layout=go.Layout(
                    title='Plot',
                    showlegend=True,
                    legend=go.layout.Legend(
                        x=0,
                        y=1.0
                    ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30)
                )
            )
    elif(len(mode)==2):
        return go.Figure(
            data=dt+dtc,
                layout=go.Layout(
                    title='Plot',
                    showlegend=True,
                    legend=go.layout.Legend(
                        x=0,
                        y=1.0
                    ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30)
                )
            )

if __name__ == '__main__':
    app.run_server(debug=True)
