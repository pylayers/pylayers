Roadmap
=======
        
        + equivalent of tra2tud fully implemented in python 
                + ajout de Node_N dans Layout sur noeud diffractant 

        + equivalent of evalfield fully implemented in python 

        + Creation d'un setup.py 

        + Export au format sql - WP4 

        + Bloc 1
                 + (Layout + SetLink) --> SetSignature_{k}
        + Bloc 2
                 + (Layout + SetLink + SetSignature_{k-1}) --> SetSignature

         Algo ::
                for Sig in  SetSignature_{k-1}
                
                        valid,ray = ConstructRay2D(Sig)

                        if valid 
                                append ray 
                        else 
                                append UnsuccessfulSignature 
                        
                for Sig in UnsuccessfulSignature
                        
                         
Layout.py 
---------
        + create the function buildGv 
        + create function for displaying properly the different Layout Associated Graphs (LAG)   

GeomUtil.py 
-----------

        + Inclure la détection de convexité dans buildGv
        + Un point convexe peut voir les segments 
               + totalement 
               + partilement 
        + Si on se limite a [tail middle head] on risque de manquer des situations de visibilité

