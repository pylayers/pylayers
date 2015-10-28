#!/usr/bin/python
## coding=utf-8
"""
The ni_usb_6501 is a digital IO module for USB from National Instruments.
Unfortunately their Linux driver is excessively large (>60MB), difficult to install
and doesn't offer off-the-shelf support for python.

This python driver is based on Marc Schutz's pioneer work on c driver
(https://github.com/schuetzm/ni-usb-6501)

INSTALLATION
1. Install the latest PyUSB (at least version 1.0.a3) from http://sourcceforge.net/projects/pyusb/

2. Change the permissions of the USB device node by creating a udev rule.
   e.g. add the following line (and file) to a file in /etc/udev/rules.d/usb.rules
   SUBSYSTEM=="usb", ENV{DEVTYPE}=="usb_device", MODE="0664", GROUP="usbusers"

   This will set the owner of the device node to root:usbusers rather than root:root
   After that add user to the usbusers group for enabling access to the device.
   adduser _<user>_ usbusers
  (Make sure you have group usbusers)

...and you are good to go.

TODO
 - Counter operations
"""
import usb.core
import usb.util

ID_VENDOR = 0x3923
ID_PRODUCT = 0x718a

def get_adapter(**kwargs):
    """
    Returns NiUsb6501 handler if only single adapter is connected to PC.
    Forwards all parameters to pyusb (http://pyusb.sourceforge.net/docs/1.0/tutorial.html)
    """
    device = usb.core.find(idVendor=ID_VENDOR, idProduct=ID_PRODUCT, **kwargs)
    if not device:
        raise ValueError('Device not found')

    return NiUsb6501(device)


    """
    Returns NiUsb6501 handle for every adapter that is connected to PC.
    Forwards all parameters to pyusb (http://pyusb.sourceforge.net/docs/1.0/tutorial.html)
    """
def find_adapters(**kwargs):
    devices = usb.core.find(find_all=True, idVendor=ID_VENDOR, idProduct=ID_PRODUCT, **kwargs)
    if not devices:
        raise ValueError('Device not found')

    return [NiUsb6501(dev) for dev in devices]


class NiUsb6501:
    """
    Typical usage:
      adapter = get_adapter()
      adapter.set_io_mode(0b00000000, 0x11111111, 0x01010101) # one bit per port 1=write, 0=read
      # start calling adapter.read_port(port) and adapter.write_port(port, values)
    """
    def __init__(self, device):
        """ used only internally via get_adapter() and find_adapters() """
        self.device = device
        cfg = self.device.get_active_configuration()
        interface_number = cfg[(0,0)].bInterfaceNumber

        if self.device.is_kernel_driver_active(interface_number):
            self.device.detach_kernel_driver(interface_number)
        # set the active configuration. With no arguments, the first
        # configuration will be the active one
        self.device.set_configuration()
        # This is needed to release interface, otherwise attach_kernel_driver fails 
        # due to "Resource busy"
        usb.util.dispose_resources(self.device)

    def set_io_mode(self, port0, port1, port2):
        """
        Set mode for every IO pin. PIN modes are given in three groups (bitmasks represented by integers)
        bit = 0: read
        bit = 1: write
        """
        buf = list("\x02\x10\x00\x00\x00\x05\x00\x00\x00\x00\x05\x00\x00\x00\x00\x00")

        buf[6] = chr(port0)
        buf[7] = chr(port1)
        buf[8] = chr(port2)
        buf = ''.join(buf)

        return self.send_request(0x12, buf)

    def read_port(self, port):
        """
        Read the value from all read-mode pins from one of the 8 PIN ports
        port is 0, 1 or 2
        """
        buf = list("\x02\x10\x00\x00\x00\x03\x00\x00")

        buf[6] = chr(port)
        buf = ''.join(buf)

        response = self.send_request(0x0e, buf)

        self.packet_matches(response,
                            "\x00\x0c\x01\x00\x00\x00\x00\x02\x00\x03\x00\x00",
                            "\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\x00\xff")

        return ord(response[10])

    def write_port(self, port, value):
        """
        Write value to all write-mode pins in one of the 8 PIN ports
        port is 0, 1 or 2
        value is 8 bits represented by integer
        """
        buf = list("\x02\x10\x00\x00\x00\x03\x00\x00\x03\x00\x00\x00")

        buf[6] = chr(port)
        buf[9] = chr(value)
        buf = ''.join(buf)

        response = self.send_request(0x0f, buf)
        self.packet_matches(response,
                            "\x00\x08\x01\x00\x00\x00\x00\x02",
                            "\xff\xff\xff\xff\xff\xff\xff\xff")

        return response

    ##########################################################
    # TODO: COUNTERS ARE NOT YET IMPLEMENTED
    ##########################################################
    def read_counter(self):
        pass

    def write_counter(self):
        pass

    def start_counter(self):
        pass

    def stop_counter(self):
        pass

    ##########################################################
    # INTERNAL UTILITY FUNCTIONS
    ##########################################################
    EP_IN, EP_OUT = 0x81, 0x01
    HEADER_PACKET, HEADER_DATA = 4, 4
    INTERFACE = 0

    def send_request(self, cmd, request):
        if len(request) + self.HEADER_PACKET + self.HEADER_DATA > 255:
            raise ValueError('Request too long (%d bytes)' % (len(request) + self.HEADER_PACKET + self.HEADER_DATA))

        buf = list("\x00\x01\x00\x00\x00\x00\x01\x00")

        buf[3] = chr(self.HEADER_PACKET + self.HEADER_DATA + len(request))
        buf[5] = chr(self.HEADER_DATA + len(request))
        buf[7] = chr(cmd)

        buf = ''.join(buf) + request

        assert self.device.write(self.EP_OUT, buf, self.INTERFACE) == len(buf)

        ret = self.device.read(self.EP_IN, len(buf), self.INTERFACE)

        return ''.join([chr(x) for x in ret])[self.HEADER_PACKET:]

    def packet_matches(self, actual, expected, mask):
        if len(actual) != len(expected):
            print repr(actual)
            print repr(expected)
            print repr(mask)
            raise ValueError('Protocol error - invalid response length %d' % len(actual))

        for b, e, m in zip(actual, expected, mask):
            if (ord(b) & ord(m)) != (ord(e) & ord(m)):
                raise ValueError("""Protocol error - invalid response
                actual:   %s
                expected: %s
                mask:     %s
                """ % (repr(actual), repr(expected), repr(mask)))


#USAGE EXAMPLE
if __name__ == "__main__":
    dev = get_adapter()

    if not dev:
        raise Exception("No device found")

    dev.set_io_mode(0b11111111, 0b11111111, 0b00000000)

    dev.write_port(0, 0b10111111)
    #dev.write_port(0, 0b00000001)
    #dev.write_port(1, 0b10101010)

    print bin(dev.read_port(0))
    #print bin(dev.read_port(1))
    #print bin(dev.read_port(2))
