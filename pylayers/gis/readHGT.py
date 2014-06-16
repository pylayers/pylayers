import struct

def get_sample(filename, n, e):
    i = 1201 - int(round(n / 3, 0))
    j = int(round(e / 3, 0))
    with open(filename, "rb") as f:
        f.seek(((i - 1) * 1201 + (j - 1)) * 2)  # go to the right spot,
        buf = f.read(2)  # read two bytes and convert them:
        val = struct.unpack('>h', buf)  # ">h" is a signed two byte integer
        if not val == -32768:  # the not-a-valid-sample value
            return val
        else:
            return None

if __name__ == "__main__":
    n = 24 * 60 + 58.888
    e = 55 * 60 + 11.377
    tile = "N48W002.hgt"  # Or some magic to figure it out from position
    print get_sample(tile, n, e)
