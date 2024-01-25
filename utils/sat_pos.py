import numpy as np
import pandas as pd
sv = 26
for name in ["sim_pos_debug_%d" % (sv), "sat_pos_debug_%d" % (sv+32+27)]:
    sim_pos = np.fromfile("/tmp/%s" % (name),dtype=np.double())

    records = {'pos_x':[], 'pos_y':[], 'pos_z':[], 'tow':[]}

    for i in range(int(len(sim_pos)/4)):
        records['pos_x'].append(list(sim_pos[i*4:(i*4)+4])[0])
        records['pos_y'].append(list(sim_pos[i*4:(i*4)+4])[1])
        records['pos_z'].append(list(sim_pos[i*4:(i*4)+4])[2])
        records['tow'].append(list(sim_pos[i*4:(i*4)+4])[3])

    sim_pos_df = pd.DataFrame(records)

    print("PRN - [%d], total pos %d" % (18, len(records)))

    print(sim_pos_df)

    sim_pos_df.to_pickle("./%s.pkl" % (name))