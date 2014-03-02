#-*- coding:Utf-8 -*-
#
# Developped by Mohamed Laaraiedh
# mlaaraie@univ-rennes1.fr
#
import doctest
import pyodbc as mssql
import MySQLdb as mysql


class pysql():
    """
    This class manages connection between PyLayers and MS SQL Server
    database
    """

    def __init__(self,server= 'localhost',dbname='', userid='',passwd=''):
        """
        initialize when using windows authentication
        (login and password are required)
        Parameters
        ----------
        server : string
        dbname : string
        """
        self.server = server
        self.db = dbname
        self.userid = userid
        self.passwd = passwd

    def connect_mssql(self):
        """
        connect to MS SQL Server database and return a cursor
        """
        cmd = 'DRIVER={SQL Server};SERVER='+self.server+';DATABASE='+self.db
        self.cnxn = mssql.connect(cmd)
        self.cursor = self.cnxn.cursor()

    def connect_mysql(self):
        """
        connect to MySQL database and return a cursor
        """
        self.cnxn = mysql.connect(self.host, self.user, self.passwd, self.dbname)
        self.cursor = self.cnxn.cursor()

    def disconnect(self):
        """
        close the connexion with database
        """
        self.cnxn.close()

    def gettables(self):
        """
        get tables names
        """
        cmd = 'select table_name from information_schema.tables'
        self.cursor.execute(cmd)
        table_names = self.cursor.fetchall()
        for i in range(len(table_names)):
            table_names[i]=table_names[i][0]
        return table_names

    def getcolumns(self,table):
        """
        get columns from a table
        """
        cmd = 'select column_name from information_schema.columns where table_name = \''+table+'\''
        self.cursor.execute(cmd)
        column_names = self.cursor.fetchall()
        for i in range(len(column_names)):
            column_names[i]=column_names[i][0]
        return column_names

    def runselectcmd(self,cmd):
        """
        run a select sql command cmd
        """

        self.cursor.execute(cmd)
        result = self.cursor.fetchall()
        return result

    def runinsertcmd(self,cmd):
        """
        run an insert sql command cmd
        """

        self.cursor.execute(cmd)
        self.cnxn.commit()

    def rundeletecmd(self,cmd):
        """
        run an insert sql command cmd
        """

        self.cursor.execute(cmd)
        self.cnxn.commit()



