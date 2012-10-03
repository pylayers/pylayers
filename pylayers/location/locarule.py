from pylayers.util.project import *
import networkx as nx
import ConfigParser	
import scipy.stats as sps
import pdb


#### global info PL_model
config	 = ConfigParser.ConfigParser()
config.read(basename+'/ini/EMSolver.ini')
plm_opt = dict(config.items('PL_MODEL'))
sigmaRSS	= float(plm_opt['sigmarss'])# dBm !!!!!! 			 			
f 			= float(plm_opt['f'])
RSSnp		= float(plm_opt['rssnp'])
d0		   	= float(plm_opt['d0'])
PL_method 	= plm_opt['method'] # mean, median , mode

class Take_all():
	"""
	Take TOA and Pr for any RAT.
	
	Take_all.take(net=Network.Network,RAT=None)		
	if RAT= None : All RAT are processed
	"""
	def take(self,net,RAT=None,LDP=None):
		cd = {}


		### select RAT
		if RAT == None :
			Rat = net.RAT

		elif isinstance(RAT,str):
			Rat = [RAT]

		### select LDP
		if LDP == None :
			Ldp = net.LDP

		elif isinstance(LDP,str):
			Ldp = [LDP]

		for ldp in Ldp:
			cd[ldp]={}
			for rat in Rat:
				try:
					cd[ldp][rat]=nx.get_edge_attributes(net.SubNet[rat],ldp).items()
				except: 
					pass # if the specified RAT doesn't exist in the PN


		return cd

#class Take_all_Pr():
#	"""
#	Take Pr for any RAT.
#	
#	Take_all_Pr.take(net=Network.Network,RAT=None)		
#	if RAT= None : All RAT are processed
#	"""

#	def take(self,net,RAT=None):
#		cd = {}


#		if RAT == None :
#			Rat = net.RAT

#		elif isinstance(Rat,str):
#			Rat = [RAT]

#		cd['Pr']={}
#		for rat in Rat:
#			cd['Pr'][rat]=nx.get_edge_attributes(net.SubNet[rat],ldp).items()



#class Take_all_TOA():
#	"""
#	Take TOA for any RAT.
#	
#	Take_all_Pr.take(net=Network.Network,RAT=None)		
#	if RAT= None : All RAT are processed
#	"""
#	def take(self,net,RAT=None):
#		cd = {}

#		if RAT == None :
#			Rat = net.RAT

#		elif isinstance(Rat,str):
#			Rat = [RAT]

#		cd['TOA']={}
#		for rat in Rat:
#			cd['TOA'][rat]=nx.get_edge_attributes(net.SubNet[rat],ldp).items()


#class Qconnect():
#	"""
#		take only nodes from a given RAT assuming connectivity probability
#	
#		p(n1,n2) = Q((10 np log_10(d_(1,2)/R))/sigma_sh)

#		with R = 10^((P_0-P_th)/10np) 

#		p_0 = 1/2 with d1,2 =R


#	"""
#	
#	def take(self,net,RAT=None,LDP=None):

#		if Rat == None :
#			Rat = net.RAT

#		elif isinstance(Rat,str):
#			Rat = [Rat]


#			
#			for rat in Rat:
#				d=nx.get_edge_attributes(self.net.SubNet[rat],'d')
#				P=nx.get_edge_attributes(self.net.SubNet[rat],'Pr')
#				

#		for ldp in Ldp:
#			cd[ldp]={}
#			for rat in Rat:
#				try:
#					cd[ldp][rat]=nx.get_edge_attributes(net.SubNet[rat],ldp).items()
#				except: 
#					pass # if the specified RAT doesn't exist in the PN



#			R = pow(10,()) 
#			var = 10*RSSnp*np.log10(d/R)
#			
#			#N=sps.uniform(scale=)
#			

def merge_rules(self,RAT=None,LDP=None):
	rules = {}

	for rule in self.rule:
		rules.update(rule.take(self.PN,RAT,LDP))
	return (rules)

