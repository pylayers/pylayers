import time
from SimPy.Simulation import Process, hold
from Transit.World import world
import os

PASSIVATE = True
DONT_PASSIVATE = False

new_state = 'new_state'
next_state = 'next_state'

class StateProcess(Process):
	def execute(self):
		while 1:
			results = self.state_generator.next()
			if results[0] is new_state:
				self.state_generator = results[1]
			elif results[0] is next_state:
				self.state_generator = getattr(self, self.state)()
			else:
				yield results

	def pace(self, seconds_to_event, *args):
		intervals = int(seconds_to_event / self.interval) + 1
		results = [intervals, seconds_to_event / intervals]
		for arg in args:
			results.append(arg / intervals)
		while intervals > 0:
			intervals -= 1
			results[0] = intervals
			yield results

class Updater(Process):
	def __init__(self, interval,sim=None):
		Process.__init__(self,sim=sim)
		self.world = world()
		self.interval = interval
		self.display_interval = 60.0

	def execute(self):
		ii = 0
		tk = self.world.tk
		next_interval = 0
		sleep_time = 0
		if os.environ.has_key('INTERVAL'):
			sleep_time = self.interval/float(os.environ['INTERVAL'])
		while 1:
			tk.main.update()
			if ii * self.interval >= next_interval:
				next_interval = next_interval + self.display_interval
				print "%f" % (ii * self.interval)
			tk.canvas.postscript(file='images-ps/%05d.ps' % ii)
			ii += 1
			if sleep_time:
				time.sleep(sleep_time)
			yield hold, self, self.interval
