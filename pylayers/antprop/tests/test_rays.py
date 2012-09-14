#       S = Simul.Simul('simul8.simul')
#    gr = GrRay3D()
#    gr.choose()
#    gt = GrRayTud()
#    gt.choose()
#    sl=SlabDB()
#    sl.load(_fileslab)
#    indoor = IndoorStr(sl,_filestr)
#    indoor.load()
#    indoor.show3(True)
#
#    os.chdir(basename+'/trace')
#    filetra = FD.askopenfilename(filetypes = [("tra file","*.tra"),("All", "*")],
#                      title="Please choose a ray tracing file",
#                      initialdir=tradir)
#    _filetra=os.path.split(filetra)[1]
#    #_filetra='ceamimo2_ceamimo2_def_Tx_1_def_Rx_1.tra'
#
#    os.chdir(curdir)
#    grRay = GrRay3D()
#    grRay.load(_filetra,indoor)
#    print "grRay has been loaded"
#    grRay.show3()
#
#    os.chdir(basename+'/tud')
#    filetud = FD.askopenfilename(filetypes = [("tud file","*.tud"),("All", "*")],
#                      title="Please choose a ray tracing file",
#                      initialdir=tuddir)
#    #grRay.save("rox.tra")
#    _filetud=pyu.getshort(filetud)
#    os.chdir(curdir)
#    mat = MatDB()
#    mat.load('def.mat')
#    #print "material database has been loaded"
#    #sl   = SlabDB()
#    #sl.load('def.slab')
#    #print "slab database has been loaded"
#    grTud = GrRayTud()
#    grTud.load(_filetud,sl,mat)
#    print "grTud has been loaded"
