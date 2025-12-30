#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:49:09 2024

@author: sasankatalukdar
"""

############ Edit the following variables ################
basin='nio'
# basin= 'nio' or 'gulf'

start_date="2025-11-27 00:00:00"
end_date="2025-11-29 00:00:00"

#provide central coordinates for the required domain
coords=np.array([[-92.5, 22], [-95.36, 29.76]])
history_interval_d2 = 60
history_interval_d3 = 60

d3_size=[5, 3]
d2_buffer=[3, 3]
num_moves=[-99, 0]

##############################################################
###### Edits beyond this point should not be required #######
#############################################################

import pickle
import numpy as np
import pandas as pd
import os

def st2(x):
    if x<10:
        return '0'+str(x)
    else:
        return str(x)

coords_dir='/work2/09534/st37357/stampede3/IOLA/parm/iola'

dom3_yres=0.011 # resolution of the innermost domain

dom3_xres=dom3_yres*2
dom2_yres=dom3_yres*3
dom2_xres=dom3_xres*3
dom1_yres=dom2_yres*3
dom1_xres=dom2_xres*3

d3_size=np.array(d3_size)
d2_buffer=np.array(d2_buffer)

d3_max_size=np.sqrt(2)*d3_size

def evenize(x):
    x1=x.copy()
    for i in range(len(x)):
        if x[i]%2!=0:
            x1[i]+=1
    return x1.astype(int)
            
e_we_3 = evenize(1+np.round(d3_size/dom3_xres))
e_sn_3 = evenize(1+np.round(d3_size/dom3_yres))


############

sort_idx=np.argsort(coords[:,0])
sort_idy=np.argsort(coords[:,1])

dist_x=np.zeros(len(coords[:,0]))
dist_y=np.zeros(len(coords[:,1]))
for idx,i in enumerate(sort_idx):
    if idx!=0:
        dist_x[i]=coords[:,0][i]-coords[:,0][sort_idx[idx-1]]
for idy,i in enumerate(sort_idy):
    if idy!=0:
        dist_y[i]=coords[:,1][i]-coords[:,1][sort_idy[idy-1]]    

#i_parent_start_3=np.round((d2_buffer+dist_x)/dom2_xres)
#j_parent_start_3=np.round((d2_buffer+dist_y)/dom2_yres)

i_parent_start_3=np.round(d2_buffer/dom2_xres)
j_parent_start_3=np.round(d2_buffer/dom2_yres)

e_we_2=evenize(1+np.round(2*d2_buffer/dom2_xres+e_we_3*dom3_xres/dom2_xres))
e_sn_2=evenize(1+np.round(2*d2_buffer/dom2_yres+e_sn_3*dom3_yres/dom2_yres))

if basin == 'nio':
    cords=pickle.load(open(f'{coords_dir}/ind_d01_coords.pkl','rb'))
    ref_lon=80
    ref_lat=20
elif basin == 'gulf':
    cords=pickle.load(open(f'{coords_dir}/gulf_d01_coords.pkl','rb'))
    ref_lon=-78
    ref_lat=25

hlon=cords['lon']
hlat=cords['lat']

cor_lon=coords[:,0]-d3_size/2-d2_buffer
cor_lat=coords[:,1]-d3_size/2-d2_buffer

i_parent_start_2=[]
j_parent_start_2=[]
for i in range(len(cor_lon)):
    dists=np.sqrt((hlon-cor_lon[i])**2+(hlat-cor_lat[i])**2)
    i_parent_start_2.append(np.where(dists==np.min(dists))[1][0])
    j_parent_start_2.append(np.where(dists==np.min(dists))[0][0])

##################################################################

st=pd.to_datetime(start_date)
ed=pd.to_datetime(end_date)

####################################################################

max_dom = 1+2*len(cor_lon)
grid_id = str([i+1 for i in range(max_dom)])[1:-1]
parent_id = '1,'+' 1,'*len(cor_lon)+grid_id[2:2+len(cor_lon)*3]
parent_grid_ratio = '1,'+' 3,'*2*len(cor_lon)
i_parent_start = '1, '+f'{i_parent_start_2}'[1:-1]+', '+f'{i_parent_start_3.astype(int).tolist()}'[1:-1]
j_parent_start = '1, '+f'{j_parent_start_2}'[1:-1]+', '+f'{j_parent_start_3.astype(int).tolist()}'[1:-1]
if len(d3_size)==1:
    e_we = '390, '+(f'{e_we_2.tolist()}'[1:-1]+', ')*len(cor_lon)+(f'{e_we_3.tolist()}'[1:-1]+',')*len(cor_lon)
    e_sn = '780, '+(f'{e_sn_2.tolist()}'[1:-1]+', ')*len(cor_lon)+(f'{e_sn_3.tolist()}'[1:-1]+',')*len(cor_lon)
else:
    e_we = '390, '+f'{e_we_2.tolist()}'[1:-1]+', '+f'{e_we_3.tolist()}'[1:-1]+','
    e_sn = '780, '+f'{e_sn_2.tolist()}'[1:-1]+', '+f'{e_sn_3.tolist()}'[1:-1]+','

####################################################################

wps=f'''
&share
  wrf_core = "NMM",
  start_date = "{str(st).replace(' ','_')}",
  end_date = "{str(ed).replace(' ','_')}",
  max_dom = {max_dom},
  interval_seconds = 10800,
  io_form_geogrid = 2,
  nocolons = T,
/

&geogrid
  parent_id = {parent_id}
  parent_grid_ratio = {parent_grid_ratio}
  i_parent_start = {i_parent_start},
  j_parent_start = {j_parent_start},
  e_we = {e_we}
  e_sn = {e_sn}
  geog_data_res = {'"2m", '*max_dom}
  dx = {dom1_yres},
  dy = {dom1_yres},
  map_proj = "rotated_ll",
  ref_lat = {ref_lat},
  ref_lon = {ref_lon},
  geog_data_path = "./geog-data/",
  opt_geogrid_tbl_path = "./",
/

&ungrib
  out_format = "WPS",
  prefix = "FILE",
/

&metgrid
  fg_name = "FILE",
  io_form_metgrid = 2,
  opt_metgrid_tbl_path = "./",
/

&mod_levs
  press_pa = 201300, 200100, 100000, 95000, 90000, 85000, 80000, 75000, 70000, 65000, 60000, 55000, 50000, 45000, 40000, 35000, 30000, 25000, 20000, 15000, 10000, 5000, 1000, 500, 200,
/
'''

##########################################################################


inp=f'''
&time_control
  start_year = {f'{st.year}, '*max_dom}
  start_month =  {f'{st.month}, '*max_dom}
  start_day =  {f'{st.day}, '*max_dom}
  start_hour = {f'{st.hour}, '*max_dom}
  start_minute = {f'{st.minute}, '*max_dom}
  start_second = {f'{st.second}, '*max_dom}
  end_year = {f'{ed.year}, '*max_dom}
  end_month = {f'{ed.month}, '*max_dom}
  end_day = {f'{ed.day}, '*max_dom}
  end_hour = {f'{ed.hour}, '*max_dom}
  end_minute = {f'{ed.minute}, '*max_dom}
  end_second =  {f'{ed.second}, '*max_dom}
  interval_seconds = 10800,
  history_interval = 60, {f'{history_interval_d2}, '*len(cor_lon)}{f'{history_interval_d3}, '*len(cor_lon)}
  auxhist1_interval = {'60, '*max_dom}
  auxhist2_interval = {'60, '*max_dom}
  auxhist3_interval = {'60, '*max_dom}
  history_end = {'540, '*max_dom}
  auxhist2_end = {'540, '*max_dom}
  auxhist1_outname = "wrfdiag_d<domain>",
  auxhist2_outname = "wrfout_d<domain>_<date>",
  auxhist3_outname = "wrfout_d<domain>_<date>",
  frames_per_outfile = {'1, '*max_dom}
  frames_per_auxhist1 = {'1, '*max_dom}
  frames_per_auxhist2 = {'1, '*max_dom}
  frames_per_auxhist3 = {'1, '*max_dom}
  analysis = {'F, '*max_dom}
  restart = F,
  restart_interval = 36000,
  reset_simulation_start = F,
  io_form_input = 11,
  io_form_history = 11,
  io_form_restart = 11,
  io_form_boundary = 11,
  io_form_auxinput1 = 2,
  io_form_auxhist1 = 202,
  io_form_auxhist2 = 11,
  io_form_auxhist3 = 11,
  auxinput1_inname = "met_nmm.d<domain>.<date>",
  debug_level = 1,
  tg_reset_stream = 1,
  override_restart_timers = T,
  io_form_auxhist4 = 11,
  io_form_auxhist5 = 11,
  io_form_auxhist6 = 11,
  io_form_auxinput2 = 2,
  nocolons = T,
/

&fdda

/

&domains
  time_step = 30,
  time_step_fract_num = 0,
  time_step_fract_den = 1,
  max_dom = {max_dom},
  s_we = {'1, '*max_dom}
  e_we = {e_we}
  s_sn = {'1, '*max_dom}
  e_sn = {e_sn}
  s_vert = {'1, '*max_dom}
  e_vert = {'75, '*max_dom}
  dx = {dom1_yres}, {f'{dom2_yres}, '*len(cor_lon)}{f'{dom3_yres}, '*len(cor_lon)}
  dy = {dom1_yres}, {f'{dom2_yres}, '*len(cor_lon)}{f'{dom3_yres}, '*len(cor_lon)}
  grid_id = {grid_id},
  tile_sz_x = 0,
  tile_sz_y = 0,
  numtiles = 1,
  nproc_x = -1,
  nproc_y = -1,
  parent_id = {'0'+parent_id[1:]}
  parent_grid_ratio = {parent_grid_ratio}
  parent_time_step_ratio = {parent_grid_ratio}
  i_parent_start = {i_parent_start},
  j_parent_start = {j_parent_start},
  feedback = 1,
  num_moves = 0, {(f'{num_moves}'[1:-1]+', ')*len(cor_lon)}
  num_metgrid_levels = 34,
  p_top_requested = 1000.0,
  ptsgm = 15000.0,
  eta_levels = 1.0, 0.997622, 0.995078, 0.99224, 0.989036, 0.98544, 0.981451, 0.977061, 0.972249, 0.966994, 0.96128, 0.955106, 0.948462, 0.941306, 0.933562, 0.925134, 0.915937, 0.90589, 0.894913, 0.882926, 0.869842, 0.855646, 0.840183, 0.823383, 0.805217, 0.785767, 0.7651, 0.7432, 0.720133, 0.695967, 0.670867, 0.645033, 0.6187, 0.592067, 0.565333, 0.538733, 0.5125, 0.4868, 0.461767, 0.437533, 0.4142, 0.391767, 0.370233, 0.3496, 0.329867, 0.310967, 0.292867, 0.275533, 0.258933, 0.243, 0.2277, 0.213, 0.198867, 0.1853, 0.172267, 0.159733, 0.147633, 0.135967, 0.124767, 0.114033, 0.103733, 0.093867, 0.0844, 0.075333, 0.0666, 0.058267, 0.050333, 0.042833, 0.035733, 0.029, 0.0226, 0.0165, 0.010733, 0.005267, 0.0,
  use_prep_hybrid = F,
  num_metgrid_soil_levels = 4,
  corral_x = {'9, '*max_dom}
  corral_y = {'18, '*max_dom}
  smooth_option = 0,
/

&physics
  num_soil_layers = 4,
  mp_physics = {'5, '*max_dom}
  ra_lw_physics = {'4, '*max_dom}
  ra_sw_physics = {'4, '*max_dom}
  sf_sfclay_physics = {'88, '*max_dom}
  sf_surface_physics = {'2, '*max_dom}
  bl_pbl_physics = {'3, '*max_dom}
  cu_physics = {'4, '*max_dom}
  mommix = {'1.0, '*max_dom}
  var_ric = 1.0,
  coef_ric_l = 0.16,
  coef_ric_s = 0.25,
  h_diff = {'1.0, '*max_dom}
  gwd_opt = 2, {'0, '*(max_dom-1)}
  sfenth = {'0.0, '*max_dom}
  nrads = 30,{' 90,'*len(coords)}{' 270,'*len(coords)}
  nradl = 30,{' 90,'*len(coords)}{' 270,'*len(coords)}
  nphs = 2, {'6, '*(max_dom-1)}
  ncnvc = 2, {'6, '*(max_dom-1)}
  ntrack = 6,{' 6,'*len(coords)}{' 18,'*len(coords)}
  gfs_alpha = {'-1.0, '*max_dom}
  sas_pgcon = 0.55, {'0.2, '*(max_dom-1)}
  sas_mass_flux = {'0.5, '*max_dom}
  co2tf = 1,
  vortex_tracker = 2,{' 2,'*len(coords)}{' 7,'*len(coords)}
  nomove_freq = 0, {'3.0, '*(max_dom-1)}
  tg_option = 1,
  ntornado = 2,{' 6,'*len(coords)}{' 18,'*len(coords)}
  cldovrlp = 4,
  ens_cdamp = 0.2,
  ens_pblamp = 0.2,
  ens_random_seed = 99,
  ens_sasamp = 50.0,
  ensda_physics_pert = 1,
  icloud = 3,
  icoef_sf = {'6, '*max_dom}
  iwavecpl = {'0, '*max_dom}
  lcurr_sf = {'F, '*max_dom}
  pert_cd = F,
  pert_pbl = F,
  pert_sas = F,
/

&dynamics
  non_hydrostatic = {'T, '*max_dom}
  euler_adv = F,
  wp = {'0, '*max_dom}
  coac = 1.5,{' 2.0,'*len(coords)}{' 2.6,'*len(coords)}
  codamp = {'12.0, '*max_dom}
  terrain_smoothing = 2,
  dwdt_damping_lev = {'2000.0, '*max_dom}
/

&bdy_control
  spec_bdy_width = 1,
  specified = T,
/

&namelist_quilt
  poll_servers = F,
  nio_tasks_per_group = 0,
  nio_groups = 1,
/

&logging
  compute_tasks_silent = T,
  io_servers_silent = T,
  stderr_logging = 0,
/
'''

with open('namelist.wps', 'a') as f:
    f.write(wps)
with open('namelist.input', 'a') as f:
    f.write(inp)

#T1L-gi2-x9P-LpA@
'''
#minimum distance condition for separate domains
d3_buffer=2
if len(coords)>1:
    dist=np.sqrt(np.sum((coords[0]-coords[1])**2))
    dist_x=np.abs(coords[0,0]-coords[1,0])
    dist_y=np.abs(coords[0,1]-coords[1,1])

    if dist>=d3_max_size+d3_buffer:
        e_we_3 = 1+round(d3_size/dom3_xres)
        e_sn_3 = 1+round(d3_size/dom3_yres)
    else:
        e_we_3 = 1+round(((max(coords[:,0])-min(coords[:,0]))+d3_size)/dom3_xres)
        e_sn_3 = 1+round(((max(coords[:,1])-min(coords[:,1]))+d3_size)/dom3_yres)
        
    if e_we_3%2!=0:
        e_we_3+=1
    if e_sn_3%2!=0:
       e_sn_3+=1

#if dist<7, there will be only one d3
if dist<d3_max_size+d3_buffer:
    i_parent_start_3=round(d2_buffer/dom2_xres)
    j_parent_start_3=round(d2_buffer/dom2_yres)    
else:
    if coords[0,0]<coords[1,0]:
        i_parent_start_3_1=round(d2_buffer/dom2_xres)
        i_parent_start_3_2=round(i_parent_start_3_1+dist_x/dom2_xres)
    else:
        i_parent_start_3_2=round(d2_buffer/dom2_xres)
        i_parent_start_3_1=round(i_parent_start_3_2+dist_x/dom2_xres)
    if coords[0,1]<coords[1,1]:
        j_parent_start_3_1=round(d2_buffer/dom2_yres)
        j_parent_start_3_2=round(j_parent_start_3_1+dist_y/dom2_yres)
    else:
        j_parent_start_3_2=round(d2_buffer/dom2_yres)
        j_parent_start_3_1=round(j_parent_start_3_2+dist_y/dom2_yres)


if dist<d3_max_size+d3_buffer:
    e_we_2=1+round(2*d2_buffer/dom2_xres+e_we_3*dom3_xres/dom2_xres)
    e_sn_2=1+round(2*d2_buffer/dom2_yres+e_sn_3*dom3_yres/dom2_yres)
else:
    e_we_2=1+round(2*d2_buffer/dom2_xres+dist/dom2_xres+e_we_3*dom3_xres/dom2_xres) #10 deg+dist+size_of_2nd_dom
    e_sn_2=1+round(2*d2_buffer/dom2_yres+dist/dom2_yres+e_sn_3*dom3_yres/dom2_yres) 
    
if e_we_2%2!=0:
    e_we_2+=1
if e_sn_2%2!=0:
   e_sn_2+=1

min_lon=np.min(coords[:,0])
min_lat=np.min(coords[:,1])
cor_lon=min_lon-d3_size/2-d2_buffer
cor_lat=min_lat-d3_size/2-d2_buffer

if basin == 'nio':
    cords=pickle.load(open('/work2/09534/st37357/share/ind_d01_coords.pkl','rb'))
    ref_lon=80
    ref_lat=20
elif basin == 'gulf':
    cords=pickle.load(open('/work2/09534/st37357/share/gulf_d01_coords.pkl','rb'))
    ref_lon=-78
    ref_lat=25
    
hlon=cords['lon']
hlat=cords['lat']

dists=np.sqrt((hlon-cor_lon)**2+(hlat-cor_lat)**2)
i_parent_start_2=np.where(dists==np.min(dists))[1][0]
j_parent_start_2=np.where(dists==np.min(dists))[0][0]


if dist>=d3_max_size+d3_buffer:
    max_dom = 4
    grid_id = '1, 2, 3, 4'
    parent_id = '1, 1, 2, 2'
    parent_grid_ratio = '1, 3, 3, 3'
    i_parent_start = f'1, {i_parent_start_2}, {i_parent_start_3_1}, {i_parent_start_3_2}'
    j_parent_start = f'1, {j_parent_start_2}, {j_parent_start_3_1}, {j_parent_start_3_2}'
    e_we = f'390, {e_we_2}, {e_we_3}, {e_we_3}'
    e_sn = f'780, {e_sn_2}, {e_sn_3}, {e_sn_3}'

else:
    max_dom = 3
    grid_id = '1, 2, 3'
    parent_id = '1, 1, 2'
    parent_grid_ratio = '1, 3, 3'
    i_parent_start = f'1, {i_parent_start_2}, {i_parent_start_3}'
    j_parent_start = f'1, {j_parent_start_2}, {j_parent_start_3}'
    e_we = f'390, {e_we_2}, {e_we_3}'
    e_sn = f'780, {e_sn_2}, {e_sn_3}'


'''
