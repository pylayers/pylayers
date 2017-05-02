#
from multiprocessing import Pool
from functools import partial

def _pickle_method(method):
	func_name = method.im_func.__name__
	obj = method.im_self
	cls = method.im_class
	if func_name.startswith('__') and not func_name.endswith('__'): #deal with mangled names
		cls_name = cls.__name__.lstrip('_')
		func_name = '_' + cls_name + func_name
	return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
	for cls in cls.__mro__:
		try:
			func = cls.__dict__[func_name]
		except KeyError:
			pass
		else:
			break
	return func.__get__(obj, cls)

import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)



class someClass(object):
	def __init__(self):
		pass
	
	def fgo_func(self, x=None):
		#can put something expensive here to verify CPU utilization
		if x is None: return 99
		return x*x

	def go_mp(self):
		pool = Pool()             
		print pool.map(self.fgo_func, range(10))

        def build(self):
    	        self.go_mp()




if __name__=='__main__':
	sc = someClass()
	sc.build()
	#x=[someClass(),someClass(),someClass()]
	#p=Pool()
	#filled_f=partial(someClass.f,x=9)
	#print p.map(filled_f,x)
	#print p.map(someClass.f,x)
