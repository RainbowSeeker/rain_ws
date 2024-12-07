uint64 timestamp					# time since system start (microseconds)

# position and velocity differences between the uav and the payload,
# direction is from the uav to the payload.
float32[3] delta_position			# delta_position = pL - pi [m]
float32[3] delta_velocity			# delta_velocity = vL - vi [m/s]
float32[3] position_uav             # position of the uav
float32[3] velocity_uav             # velocity of the uav