class_name StarData extends Object

var x : float = 0.0
var y : float = 0.0
var z : float = 0.0

var d : float = 0.0
var m : float = 0.0


func _init(_source : PackedStringArray):
	if _source.size() < 5:
		push_error("Unable to construct the instance.")
	else:
		self.x = float(_source[0])
		self.y = float(_source[1])
		self.z = float(_source[2])
		self.d = float(_source[3])
		self.m = float(_source[4])
	pass
