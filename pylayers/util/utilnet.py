from   pylayers.mobility.transit.vec3 import vec3
import numpy as np
#import MySQLdb as msql


def conv_vecarr(vin,dim=2):
    """ convert vec3 => np.array and  np.array => vec3 

    Parameters
    ----------
    dim : int 
        default (2) 
    """
    if isinstance(vin,vec3):
        if dim == 2:
            return(np.array((vin[0],vin[1])))
        elif dim == 3:
            return(np.array((vin[0],vin[1],vin[2])))
        else :
             raise AttributeError('incorrect input vec3 or nd array dimension')
    elif isinstance(vin,np.ndarray): 
        if dim == 2:
            return(vec3((vin[0],vin[1],0.0)))
        elif dim == 3:
            return(vec3((vin[0],vin[1],vin[2])))
        else :
             raise AttributeError('incorrect input vec3 or nd array dimension')
    else : 
        raise TypeError('vin must be either a vec3 or a np.ndarray instance')


def str2bool(v):
  return v.lower() in ("True", "true", "1")

def create_msql(HOST='localhost',USER='root',PASSWD='sqlsql',DATABASE='piplusdb'):
    """
    Attributes
    ----------
        self : connect to a mysql database to save network
    
    Retuns 
    ------
        db,cursor: 
            database object, cursor object
            
    >>> db,cursor=create_msql()
    >>> 
    """

    HOST='localhost'
    USER='root'
    PASSWD='sqlsql'
    DATABASE='piplusdb'
    #Connection to DB
    db=msql.connect(HOST, USER, PASSWD,DATABASE)
    # We'll need a cursor to this database to excute commands
    cursor = db.cursor()
    #
    #Create tables
    #

    #~ create a new empty table for transmitters (excuted once)
    cursor.execute("CREATE TABLE IF NOT EXISTS TxNodes (TxID INT PRIMARY KEY AUTO_INCREMENT, X DOUBLE, Y DOUBLE, Z DOUBLE)")

    #~ create a new empty table for Receivers (excuted once)
    cursor.execute("CREATE TABLE IF NOT EXISTS RxNodes (RxID INT PRIMARY KEY AUTO_INCREMENT, X DOUBLE, Y DOUBLE, Z DOUBLE)")

    #~ create a new empty table for Simulations (excuted once)
    cursor.execute("CREATE TABLE IF NOT EXISTS UWBsim (SimID INT PRIMARY KEY AUTO_INCREMENT, TxID INT, RxID INT, RSSI DOUBLE, TOA DOUBLE, TOF DOUBLE, DIST DOUBLE)")


    return (db,cursor)
