from pylayers.util.project import *
import pylayers.gis.osmparser as osm 
import matplotlib.pyplot as plt 

filename = datadir+'/osm/marne.osm'
#filename = datadir+'/osm/poland.osm'
#filename = datadir+'/osm/marne.osm'
#
coords,nodes,ways,relations= osm.osmparse(filename,typ='building',verbose=True)
#bdg = osm.buildingsparse(filename)
#bdg.show()
#
##indoor ={}
##nbat = 0
##
## Get the shell of the building
##
## Ajouter les nodes id impliques
##
##for b in bats.building:
##    if 'buildingpart' in bats.building[b]['tags']:
##        if  bats.building[b]['tags']['buildingpart']=='shell':
##            pshell = bats.building[b]['poly']
##            indoor[nbat]={'id':b,'shell':pshell}
##            nbat +=1
##            fig,ax =pshell.plot(fig=fig,ax=ax)
##
###
### Get the room included within the shell
###
##for bid in indoor:
##    indoor[bid]['level']={}
##    pshell = indoor[bid]['shell']
##    for b in bats.building:
##        tags = bats.building[b]['tags']
##        if b != indoor[bid]['id']:
##            if 'buildingpart' in tags:
##                try :
##                    level = tags['level']
##                except:
##                    level = 0
##                if (tags['buildingpart']=='room') | \
##                   (tags['buildingpart']=='corridor') | \
##                   (tags['buildingpart']=='hall') | \
##                   (tags['buildingpart']=='verticalpassage'):
##                    proom  = bats.building[b]['poly']
##                    if proom.within(pshell):
##                        try:
##                            indoor[bid]['level'][level].append(proom)
##                        except:
##                            indoor[bid]['level'][level]=[proom]
##
##for bid in indoor:
##    for level in indoor[bid]['level']:
##        for r in indoor[bid]['level'][level]:
##            fig,ax = r.plot(fig=fig,ax=ax)
