# Rain's Workspace with ROS2
## How to use
### Multi-quadrotor sil simulation
```
source /opt/ros/humble/setup.bash
colcon build
source $PWD/install/setup.bash
ros2 launch formation mc_formation_launch.py
```
### Multi-quadrotor actual test
Each onboard computer:
1. connect wifi hotspots
```
sudo nmcli dev wifi list                    # look for nearby wifi 
sudo nmcli dev wifi connect 'wifi_name' password 'your_password'
sudo nmcli connection down 'wifi_name'      # stop connection
sudo nmcli connection modify 'wifi_name' ipv4.addresses 192.168.1.11/24 # addresses format: 192.168.1.1x, from 1 to 3. 10 is reserved for gcs.
sudo nmcli connection modify 'wifi_name' ipv4.method manual
sudo nmcli connection up 'wifi_name'        # start connection
```
2. start ros2 node
```
ros2 launch formation mc_single_prelaunch.py amc_id:=1 &
ros2 launch formation mc_single_launch.py amc_id:=1
```
or Using `form_tools.sh`:
```
./form_tools.sh --prelaunch
./form_tools.sh --launch
./form_tools.sh --kill
```