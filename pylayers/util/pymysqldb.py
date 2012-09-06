import MySQLdb as msql
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
				sql=sql+str(v)
			else:
				sql=sql + str(v)+','
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

#		"INSERT INTO RxNodes(X,Y,Z) VALUES(%s,%s,%s)",(X,Y,Z)

if __name__ == '__main__':
	host='localhost'
	user='root'
	passwd='sqlsql'
	dbname='RTsim'
