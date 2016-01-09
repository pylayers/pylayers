import pylayers.measures.switch.ni_usb_6501 as sw
import time
switch = sw.get_adapter()
#reattach=False
if not switch:
    print("No device found")

#for each measurements, set up the NI USB 6501 mode : bit 1 means write and bit 0 read
#switch.set_io_mode(0b11111111, 0b11111111, 0b00000000) 


#switch.write_port(1, 1) #select output 2 of the switch 1-4
#switch.write_port(1, 2) #select output 3 of the switch 1-4  

#switch.write_port(0, 4) #select channel 5 of the switch 1-8
#switch.write_port(0, 5) #select channel 6 of the switch 1-8


#
#MIMO case
#

# tic = time.time()

# for k in range(8):
#     print "	Transmiter side (ULA 8) : output selected  ",k
#     switch.write_port(0,k)
#     for  l in range(4):
#         print "Receiver side (ULA 4) : output selected ",l
#         switch.write_port(1,l)
#         time.sleep(1) #time waiting of the switch between antennas
# toc = time.time()
# t = toc - tic
# print "time switch (s):",t
