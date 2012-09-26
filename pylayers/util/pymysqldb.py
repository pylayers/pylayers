import MySQLdb as msql
import pylayers.util.pyutil as pyu
import pdb

class Database:
    def __init__(self, host, user, passwd, dbname):
        self.host=host
        self.user=user
        self.passwd=passwd
        self.dbname=dbname
    
    def connectdb(self):
        #Connection to DB
        db=msql.connect(self.host, self.user, self.passwd, self.dbname)
        # We'll need a cursor to this database to excute commands
        dbc = db.cursor()
        return dbc
        
    def dropdb(self):
        #delete DB
        dbc=self.connectdb()
        dbc.execute("drop database "+self.dbname)
    


    def showtables(self):
        dbc=self.connectdb()
        dbc.execute("show tables")
        fields=dbc.fetchall()
        print "====================================="
        print "Tables in "+self.dbname
        print "====================================="
        for field in fields:
            print field[0]
        print "====================================="
            
    def describetable(self, tablename=""):
        dbc=self.connectdb()
        dbc.execute("describe "+tablename)
        fields=dbc.fetchall()
        print "====================================="
        print "Columns in "+tablename
        print "====================================="
        for field in fields:
            print field[0],"\t\t", field[1],"\t\t", field[2],"\t\t", field[3]
        print "====================================="
    
    def droptable(self, tablename):
        #delete DB
        dbc=self.connectdb()
        dbc.execute("drop table "+tablename)
    
    def insertitem1(self, tablename="tab", columns=("","",""), values=("","","")):
        """
        insert items in case the primary key is in Auto Increment mode
        """
        dbc=self.connectdb()
        sql="INSERT INTO " + tablename + "("
        for i in range(len(columns)):
            if i==max(range(len(columns))):
                sql=sql+columns[i]
            else:
                sql=sql +columns[i]+","
        sql=sql+") VALUES (" 
        for i,v in enumerate(values):
            if i==max(range(len(values))):
                if not isinstance(v,str):
                    sql=sql+str(v)
                else :
                    sql=sql+'\"' +str(v) + '\"'
            else:
                if not isinstance(v,str):
                    sql=sql + str(v)+','
                else :
                    sql=sql+'\"' +str(v) + '\"'+','
        sql=sql+")"
        dbc.execute(sql)



    def insertitem2(self, tablename="tab", values=("","","")):
        """
        insert items in case the primary key is in not in Auto Increment mode
        """
        dbc=self.connectdb()
        sql="INSERT INTO " + tablename + " Values ("
        for i,v in enumerate(values):
            if i==max(range(len(values))):
                if not isinstance(v,str):
                    sql=sql+str(v)
                else :
                    sql=sql+'\"' +str(v) + '\"'
            else:
                if not isinstance(v,str):
                    sql=sql + str(v)+','
                else :
                    sql=sql+'\"' +str(v) + '\"'+','
        sql=sql+")"
        dbc.execute(sql)

    def writemeca(self,ID,time,p,v,a):
        """
        write mecanic information into the mysql db
        """
        self.insertitem1("TruePosition",('NodeID',
                                            'Timestamp',
                                            'X',
                                            'Y',
                                            'Z',
                                            'ReferencePointID'),
                                            (eval(ID),
                                            pyu.timestamp(time),
                                            p[0],
                                            p[1],
                                            'NULL',
                                            'NULL'))
        self.insertitem1("CEASensorMeasurements",('NodeID',
                                                     'Timestamp',
                                                     'CEA_MagX',
                                                     'CEA_MagY',
                                                     'CEA_AccX',
                                                     'CEA_AccY'),
                                                     (eval(ID),
                                                      pyu.timestamp(time),
                                                      v[0],
                                                      v[1],
                                                      a[0],
                                                      a[1] ))



    def writenet(self,net,t):
        """
        write mecanic information into the mysql db
        """
        for e in net.edges_iter(data=True):
            self.insertitem1("ACOLinkMeasurements",('NodeID',
                                                'ACO_PeerID',
                                                'ACO_RSSI',
                                                'Timestamp'),
                                                (eval(e[0]),
                                                eval(e[1]),
                                                e[2]['Pr'][0],
                                                pyu.timestamp(t)))
            self.insertitem1("CEALinkMeasurements",('NodeID',
                                                'Timestamp',
                                                'CEA_PeerID',
                                                'CEA_Dist'),
                                                (eval(e[0]),
                                                pyu.timestamp(t),
                                                eval(e[1]),
                                                e[2]['d']))

    def writenode(self,ID,name,MoA):
        self.insertitem1("Nodes",('NodeID',
                            'NodeName',
                            'NodeOwner',
                            'NodeDescription',
                            'NodeOwnerID',
                            'MobileOrAnchor',
                            'TrolleyID'),
                             (eval(ID),
                             name,
                            'NULL',
                            'node description',
                            'NULL',
                            MoA,
                            'NULL'))

if __name__=="__main__":

    HOST='localhost'
    USER='root'
    PASSWD='sqlsql'
    DATABASE='TEST'


    ldp=read_mat(1,1)
    #Connection to DB
    db = DB(HOST, USER, PASSWD,DATABASE)
#    db.connectdb()
    pdb.set_trace()
    db.insertitem1("TruePosition", (1,10,13.0,-12.5,15.0,5))
