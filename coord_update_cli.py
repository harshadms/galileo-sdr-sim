import math
import time
import struct
import socket
import threading
import numpy as np
from pynput import keyboard
from geopy.distance import distance
from pynput.keyboard import Key, Controller


global coord, vx, vy, exit_flag

coord = [51.5072,0.1276,40]
vx = 0
vy = 0

count = 0

time_set = False
tstart = 0

avg = 0

exit_flag = False

## Connection
HOST = "127.0.0.1"
PORT = 7533
BUFFSIZE = 1024
ADDR = (HOST, PORT)

def send_coordinate_thread():
    global exit_flag, vx, vy, coord

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    dt = 0.1

    csv_log = open("./dynamic_motion.log","a")

    while not exit_flag:
        update_coordinates(vx, vy, dt)
        data = struct.pack("ddd", coord[0], coord[1], coord[2])

        csv_log.write(f"{coord[0]}, {coord[1]}, {coord[2]}\n")

        sock.sendto(data, ADDR)
        time.sleep(dt)

    csv_log.close()

def update_coordinates(vx, vy, dt):
    global coord

    print (f"{vx}, {vy}\n")
    ang = math.degrees(math.atan2(vy, vx)) - 90
    mag = np.sqrt(vx **2 + vy ** 2)
                
            
    old_lat = coord[0]
    old_lon = coord[1]
        
    dest = distance(meters=mag * dt).destination(
        (old_lat, old_lon), bearing=ang
    )
    
    coord[0] = dest.latitude
    coord[1] = dest.longitude

def on_key_release(key):
    if key == Key.esc:
        raise KeyboardInterrupt

def on_press(key):
    global count
    global tstart
    global time_set
    global coord
    global vx, vy
    dv = 0.5

    if key == Key.right:
        vx = vx - dv
            
    elif key == Key.left:
        vx = vx + dv

    elif key == Key.up:
        vy = vy + dv

    elif key == Key.down:
        vy = vy - dv

    elif key == Key.esc:
        exit()

# Collect events until released
try:
    coord_update_thread = threading.Thread(target=send_coordinate_thread, args=())
    coord_update_thread.start()

    with keyboard.Listener(on_press=on_press, on_release=on_key_release) as listener:            
        listener.join()

except KeyboardInterrupt:
    exit_flag = True
    print ("Exiting")
