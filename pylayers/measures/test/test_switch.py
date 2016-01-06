import pylayers.measures.switch.ni_usb_6501 as sw
import time
switch = sw.get_adapter()
#reattach=False
if not switch:
    print("No device found")

#
#very important : to use at the beginning of the initialization of the switch
#for each measurements, set up the NI USB 6501 mode : bit 1 means write and bit 0 read
#

switch.set_io_mode(0b11111111, 0b11111111, 0b00000000) 

#
#SISO case
#

#example for use the switch 1 to 4
#switch 1 to 4 : port 1 is allowed

#switch.write_port(1, 1) #select output 2 of the switch 1-4
#switch.write_port(1, 2) #select output 3 of the switch 1-4

#example for use the switch 1 to 8
#switch 1 to 8 : port 0 is allowed    

#switch.write_port(0, 4) #select channel 5 of the switch 1-8
#switch.write_port(0, 5) #select channel 6 of the switch 1-8


#
#MIMO case
#

tic = time.time()

for k in range(8):
    print "	Transmiter : select output number ",k
    switch.write_port(0,k)
    for  l in range(4):
        print "Receiver : select output number ",l
        switch.write_port(1,l)
        time.sleep(1) #time waiting of the switch between antennas
toc = time.time()
t = toc - tic
print "Measurement time (s) with switching :",t
