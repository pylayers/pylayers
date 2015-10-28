import pdb
from xml.etree.ElementTree import Element
from xml.etree import ElementTree

class Matter(Element):
    """ Metal or Material

    typ : 'Me','Ma','Cs'

    Me=Matter()
    Ma

    """
    def __init__(self,Name,typ='Me',**kwargs):

        dtyp={'Me':'Metal',
              'Ma':'Material',
              'Cs':'ConductingSheet'
             }

        Element.__init__(self,dtyp[typ],Name=Name)
        Prim = Element('Primitives')
        self.append(Prim)


        if typ!='Me':
            if typ=='Cs':
                if 'thickness' not in kwargs:
                    raise NameError('Conducting Sheet must have a thickness')
                if 'conductivity' not in kwargs:
                    raise NameError('Conducting Sheet must have a conductivity')

                assert (eval(kwargs['conductivity'])>1e6), "conductivity below 1MA/Vm is not recommended"
                assert (eval(kwargs['thickness'])<500e-6), "a thickness greater than 500 um is not recommended"
                assert (eval(kwargs['thickness'])>1e-6), "a thickness lower than 1 um is not recommended"

                self.attrib.update(kwargs)

            if typ=='Ma':
                Prop = Element('Property')
                Prop.attrib.update(kwargs)
                self.append(Prop)


class Material(object):
    def __init__(self):
        pass
    def Debye(f,**kwargs):
        """ Debye model

        Parameters
        ----------

        f :
        eps_r
        kappar :
        eps_Delta :
        t_relax :

        """
        defaults = {'f' : 1000000,
                    'eps_r':1,
                    'kappar':1,
                    'eps_Delta':1,
                    't_relax':1
                   }
        for k in defaults:
            if k not in kwargs:
                kwargs[k] = defaults[k]

        eps0 = 8.85418781762e-12
        self.eps = kwargs['eps_r']-1j*kappa/(2*np.pi*f)/eps0
        term = sum(kwargs['eps_Delta']/(1+2j*np.pi*f*kwargs['t_relax']),axis=0)
        self.eps = self.eps + term

    def Lorentz(f,**kwargs):
        """
        """
    def Drude(f,**kwargs):
        """ Calculate the Drude type dispersive material constant

        Example
        -------

        >>> f = np.linespace(300e12,1100e12,201)
        >>> M = Material()
        >>> M.Drude(f=f,eps_r=3.942,kappa=7.97e3,plas

        """

#CalcDebyeMaterial.m
#CalcDrudeMaterial.m
#CalcLorentzMaterial.m
