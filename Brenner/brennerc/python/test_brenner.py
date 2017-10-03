#!/usr/local/bin/python

import _brenner

class testAtom:
    def __init__(self, elemno, x, y, z):
	self.atomic_number = elemno	
	self.pos = [x, y, z]

    def position(self):
	return self.pos

    def setPosition(self, pos):
	self.pos = pos

class testObj:
    def __init__(self, atoms):
	self.starttime = 0.0
	self.timestep = 0.5
	self.atoms = atoms

    def atomList(self):
	return self.atoms

obj = testObj([testAtom(6,0.0,0.0,0.2),
	       testAtom(1,1.0,0.0,0.0),
	       testAtom(1,-1.0,0.0,0.0),
	       testAtom(1,0.0,1.0,0.0),
	       testAtom(1,0.0,-1.0,0.0)])

b = _brenner.BrennerIntegrator(obj)

for i in range(10):
    b(obj, 10)

    print 'step',i*10
    for a in obj.atomList():
	print a.atomic_number, a.position()
