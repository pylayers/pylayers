import pylayers.measures.switch.ni_usb_6501 as sw

switch = sw.get_adapter()
reattach=False
if not switch:
    print("No device found")
