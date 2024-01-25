'''
Converts nav bits prn and TOW to a custom binary file format where
[[bit+prn][bit+prn][bit+prn]....[TOW]] .... [[bit+prn][bit+prn][bit+prn]....[TOW]]

There is a posibility 
'''
import struct
from tqdm import trange
import numpy as np
import h5py

chn_data = []
dt = 0.004

for chn in range(8):
    f = h5py.File("galileo_telem%d.mat" % (chn),'r')

    data = f.get('nav_symbol')
    telem_bits = [] # For converting to a NumPy array
    for i in np.array(data):
        telem_bits.append(i[0])
        
    data = f.get('TOW_at_current_symbol_ms')
    tows = []
    for i in np.array(data):
        tows.append(i[0])
        
    data = f.get('PRN')
    prn = []
    for i in np.array(data):
        prn.append(int(i[0]))
        
    prn_str = ''
    for i in set(prn):
        prn_str = prn_str + str(i)
    
    chn_data.append({'bits':telem_bits, 'prn':prn_str, 'TOW':tows})
    
    print ("Chn - %d, PRNs - [%s], total bits %d, range - %f - %f | %f - %f" % (chn, prn_str, len(telem_bits), tows[0], tows[-1], len(telem_bits)*dt, tows[-1]-tows[0]))
    
cntrl = '18'

tow_range = np.arange(chn_data[7]['TOW'][0], chn_data[7]['TOW'][-1], step=dt)
tow_min = chn_data[7]['TOW'][0]
tow_max = chn_data[7]['TOW'][-1]

data_to_write = []
#for tow in tow_range:
for chn in chn_data:
    strt = list(np.arange(tow_min, chn['TOW'][0], dt))
    end = list(np.arange(chn['TOW'][-1], tow_max, dt))
    chn['TOW'] = strt + chn['TOW'] + end
    chn['bits'] = [0]*len(strt) + chn['bits'] + [0]*len(end)

    print ("PRN - [%s], total bits %d, range - %f - %f | %f - %f" % (chn['prn'], len(chn['bits']), chn['TOW'][0], chn['TOW'][-1], len(chn['bits'])*dt, chn['TOW'][-1]-chn['TOW'] [0]))

f = open("./rx_bits_tows.dat","wb")

for i in trange(len(chn_data[7]['TOW'])):
    data = []
    for chn in chn_data:
        prn = int(chn['prn'])
        bit = int(chn['bits'][i])
        tow = chn['TOW'][i]
        
        if bit == 1:
            d = prn * 10 + bit
        else:
            d = prn * 10 + 0
            
        data.append(np.double(d))
        
    data.append(tow)
    bytes_data = struct.pack('%sd' % len(data), *data)
    f.write(bytes_data)

f.close()