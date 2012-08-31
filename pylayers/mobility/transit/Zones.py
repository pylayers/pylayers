from Transit.World import world

class Stairs:
	def __init__(self, **kw_args):
		self.__dict__ = kw_args
		self.world = world()

	def __call__(self, boid):
		on_stairs = (self.lower_left.x < boid.position.x < self.upper_right.x
					 and self.lower_left.y < boid.position.y < self.upper_right.y)
		if not on_stairs:
			if hasattr(boid, 'on_stairs'):
				boid.desired_speed = boid.max_speed
				delattr(boid, 'on_stairs')
			return
		boid.on_stairs = True
		boid.desired_speed = boid.max_speed / 2

	def draw(self):
		tk = self.world.tk
		canvas, x_, y_ = tk.canvas, tk.x_, tk.y_
		yy = self.lower_left.y
		while yy <= self.upper_right.y:
			canvas.create_line(x_(self.lower_left.x), y_(yy),
							   x_(self.upper_right.x), y_(yy), fill='white')
			yy += 0.28
