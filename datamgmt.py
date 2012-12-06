#
#  Database_v1.py
#  Fetch data from the pepsi sql database.
#
#  Created by Gorm Bruun Andresen on 16/12/2010.
#  Severely altered by Rolando A Rodriguez  20/07/2012.
#  Copyright (c) 2010 University of Aarhus. All rights reserved.
#

import numpy as np
from pylab import *

from sqlalchemy.engine import reflection
from sqlalchemy import create_engine
from sqlalchemy.schema import (MetaData,Table)

#cp: engine, inspector, metadata = connect2database()
def connect2database(url='postgresql://dheide:050877-4125@localhost/iset_data',localhost=True,echo=False):

	if localhost:
		## Works from anywhere else if you setup an ssh tunnel first: ssh -L5432:localhost:5432 username@pepsi.imf.au.dk
		engine = create_engine(url, echo=echo)
	else:
		## Works when you are on the IMF intranet.
		engine = create_engine(url, echo=echo)

	inspector = reflection.Inspector.from_engine(engine)

	metadata = MetaData(engine)
	
	return engine, inspector, metadata

def get_all_table_names(inspector=None):

	if inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);

	for schema_name in inspector.get_schema_names():
		print 'Tables in schema:',schema_name
		print '---------------------'
		for table_name in inspector.get_table_names(schema=schema_name):
			print table_name
		print
	
def get_datatype(table='europe_load_dt',schema='agg_avg_1hour_pdata_caps_eu2020', inspector=None):

	if inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	datatype = inspector.get_columns(table,schema=schema)
	
	return datatype
	
def get_data(table='europe_load_dt',schema='agg_avg_1hour_pdata_caps_eu2020', metadata=None):
	
	if metadata==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	data = Table(table,metadata,schema=schema,autoload=True)
	s = data.select()
	rs = s.execute() 
	
	data = list(empty(rs.rowcount))
	for i in arange(rs.rowcount):
		data[i] = list(rs.fetchone())
	
	return_data = transposed(data)
	#data = rs.fetchall()
	
	return return_data

def transposed(lists):
   if not lists: return []
   return map(lambda *row: list(row), *lists)
	
def get_data_countries(schema='agg_avg_1hour_pdata_caps_eu2020', metadata=None, inspector=None, localhost=True):

	if metadata==None or inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	#get labels
	datatype = get_datatype('countries_load_dt',schema, inspector)
	datalabels = list(empty(size(datatype)-1));
	for i in arange(0,size(datalabels)):
		datalabels[i] = datatype[i+1]['name']
		
	#get data
	load = get_data('countries_load_dt',schema,metadata)
	time = load[0];
	datetime_offset = date2num(time[0]);
	time = date2num(time) - datetime_offset #Use num2date(time+datetime_offset) to convert back to date-type objects.
	
	load = array(load[1:]);
	wind = array(get_data('countries_wind',schema,metadata)[1:])
	solar = array(get_data('countries_pv',schema,metadata)[1:])
	
	return time, load, wind, solar, datetime_offset, datalabels
	
def get_data_regions(schema='agg_avg_1hour_pdata_caps_eu2020', metadata=None, inspector=None, localhost=True):

	if metadata==None or inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	#get labels
	datatype = get_datatype('regions_load_dt',schema, inspector)
	datalabels = list(empty(size(datatype)-1));
	for i in arange(0,size(datalabels)):
		datalabels[i] = datatype[i+1]['name']
		
	#get data
	load = get_data('regions_load_dt',schema,metadata)
	time = load[0];
	datetime_offset = date2num(time[0]);
	time = date2num(time) - datetime_offset #Use num2date(time+datetime_offset) to convert back to date-type objects.
	
	load = array(load[1:]);
	wind = array(get_data('regions_wind',schema,metadata)[1:])
	solar = array(get_data('regions_pv',schema,metadata)[1:])
	
	return time, load, wind, solar, datetime_offset, datalabels

	
def plot_load_wind_solar_countries(Ndays=365, localhost=True):
	"""Test plot load, wind, and solar data for all countries."""
	
	
	time, load, wind, solar, datetime_offset, datalabels = get_data_countries(localhost=localhost);
	
	figure(1); clf();
	for i in range(load.shape[0]):
		plot(time[:Ndays*24],load[i][:Ndays*24]/sum(load[i]),label=datalabels[i])	

	xlabel('Time [hours]')
	ylabel('Load [norm.]')
	
	legend()
	
	savefig('CountriesLoad.png')
	
	figure(1); clf();
	for i in range(load.shape[0]):
		plot(time[:Ndays*24],wind[i][:Ndays*24]/sum(wind[i]),label=datalabels[i])	

	xlabel('Time [hours]')
	ylabel('Wind power [norm.]')
	
	legend()
	
	savefig('CountriesWind.png')
	
	figure(1); clf();
	for i in range(load.shape[0]):
		plot(time[:Ndays*24],solar[i][:Ndays*24]/sum(solar[i]),label=datalabels[i])	
	
	xlabel('Time [hours]')
	ylabel('Solar power [norm.]')
	
	legend()
	
	savefig('CountriesSolar.png')

######
# Convinient access to ISET country data.
# Main function:  get_ISET_country_data()
###

def get_ISET_country_data(ISO='DK',path='./data/'):
    """
    Returns data for a specific country in the ISET data set. 
    
    The country is specified using ISO two letter names. For a list of names use get_ISET_country_names().
    
    The function retrives/stores data from/in ./data as default. If the data is not available the function
     attempts to download it from pepsi. This requires ssh access: 
     ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk 
     
    Filenames are: ISET_country_ISO.npz
    
    Returns
    -------
    t: array
        hours numbered from 0 to N-1
    L: array
        load in MW
    Gw: array
        normalized wind power generation
    Gs: array
        normalized solar power generation
    datetime_offset: scalar
        num2date(datetime_offset) yields the date and hour of t=0
    datalabel: string
        country ISO
        
    Example: Returns the Danish load, wind and solar power generation
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data('DK')
    
    """
    
    if not valid_ISO(ISO):
        sys.exit("Error (43nlksd): No such country ISO ({0}). For a list of names use get_ISET_country_names().".format(ISO))
    
    filename = 'ISET_country_' + ISO + '.npz'
    
    try:
        #Load the data file if it exists:
        npzfile = np.load(path + filename)
        print 'Loaded file: ', path + filename
        sys.stdout.flush()
        
    except IOError:
        print 'Datafile does not exist:', path + filename
        print 'Trying to download data from pepsi...'
        sys.stdout.flush()
        try: 
            #Normalized data: t, l, Gw, Gs, datetime_offset, datalabels = get_data_countries(schema='norm_agg_avg_1hour_pdata_caps_eu2020',localhost=True);
            t, L, GW, GS, datetime_offset, datalabels = get_data_countries(localhost=True);
        except:
            sys.exit("Error (sdf3dz1): Could not connect to pepsi. Setup ssh access first: ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk")
            
        #Save all country files:
        for i in arange(len(datalabels)):
            ISO_ = ISET2ISO_country_codes(datalabels[i])
            filename_ =  'ISET_country_' + ISO_ + '.npz'
            np.savez(path + filename_,t=t, L=L[i], Gw=GW[i]/mean(GW[i]), Gs=GS[i]/mean(GS[i]), datetime_offset=datetime_offset, datalabel=ISO_)
            print 'Saved file: ', path + filename_
            sys.stdout.flush()
         
        #Load the relevant file now that it has been created:       
        npzfile = np.load(path + filename)
        print 'Loaded file: ', path + filename
        sys.stdout.flush()
        
    return npzfile['t'], npzfile['L'], npzfile['Gw'], npzfile['Gs'], npzfile['datetime_offset'], npzfile['datalabel']

def valid_ISO(ISO='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    return (ISO in table['ISO'])
    
def ISO2ISET_country_codes(ISO='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISET = table['ISET'][find(table['ISO']==ISO)][0]
    
    return ISET

def ISET2ISO_country_codes(ISET='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISO = table['ISO'][find(table['ISET']==ISET)][0]
    
    return ISO

def get_ISET_country_names(filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)

    print 'Index\tISO\tISET\tName'
    print '====================================='
    for i in arange(len(table)):
        print '{0}\t{1}\t{2}\t{3}'.format(i, table[i]['ISO'], table[i]['ISET'], table[i]['name'])
    sys.stdout.flush()

get_ISET_country_data()
