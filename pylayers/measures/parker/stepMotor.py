#!/usr/bin/python
## coding=utf-8
import serial

ser=serial.Serial(port = "/dev/ttyUSB0",
                  baudrate=9600,
                  parity='N',
                  stopbits=1,
                  bytesize=8,
                  timeout=0.5)
# Direct mode
c0="1OFF\r\n"
c1="1ON\r\n"
c2="1D400\r\n"
c3="1LIMITS(3,0,0)\n"
c4="1G\n"
lc =[c0,c1,c2,c3,c4]
command = "1W(MR,4000)\r\n1W(PA,0)\r\n1R(MA)\r\n1MA\r\n1ON\r\n1LIMITS(0,1,1)\r\n1HOME(+,-10,+10,0)\r\n1GH\r\n"
#for c in lc:
#    print c
ser.write(command)
stream = ser.read()
#ser.close()
