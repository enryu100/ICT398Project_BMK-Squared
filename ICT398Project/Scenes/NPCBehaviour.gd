extends MeshInstance

# class member variables go here, for example:
# var a = 2
# var b = "textvar"
var ai = AIEmotion.new()


func _ready():
	# Called when the node is added to the scene for the first time.
	# Initialization here
	setBaseline()
	pass

func _process(delta):
	# Called every frame. Delta is time since last frame.
	# Update game logic here.
	var current_pos = Man.get_global_pos
	pass

func setBaseline():
	ai.setEmotionLevels(4.0, 2.0, 1.0, 0.0)
	
	pass
