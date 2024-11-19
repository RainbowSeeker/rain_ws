
uint64 timestamp					# time since system start (microseconds)

uint16 UAV_CMD_FORM_START   = 100
uint16 UAV_CMD_FORM_STOP    = 101

float32 param1			# Parameter 1, as defined by MAVLink uint16 VEHICLE_CMD enum.
float64 param2			# Parameter 2, as defined by MAVLink uint16 VEHICLE_CMD enum.
float64 param3			# Parameter 3, as defined by MAVLink uint16 VEHICLE_CMD enum.
uint32 command			# Command ID