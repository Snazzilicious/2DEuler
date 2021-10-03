

class SortingList:
	def __init__(self, maxSize=10):
		self.maxEntries = []
		self.maxSize = maxSize
		
	def push(self, elem, val):
		self.insert(elem,val)
		while len(self.maxEntries) > self.maxSize :
			self.maxEntries.pop(0)
	
	def insert(self, elem, val):
		for i in range(len(self.maxEntries)):
			if val < self.maxEntries[i][1] :
				self.maxEntries.insert( i, (elem,val) )
				return
				
		self.maxEntries.append( (elem,val) )
