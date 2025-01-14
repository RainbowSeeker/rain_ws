# Rain's Workspace with ROS2
## Requirement
| Software        | Version    |
|-----------------|------------|
| Operating System| Ubuntu 22.04|
| ROS2            | Humble    |
| PX4             | v1.14.3   |
| Gazebo          | gz-garden |
### Installation
```
cd ~

# ROS2
wget http://fishros.com/install -O fishros && . fishros

# PX4 + Gazebo
git clone https://github.com/PX4/PX4-Autopilot.git -b v1.14.3 --recursive
bash PX4-Autopilot/Tools/setup/ubuntu.sh

# Workspace
git clone https://github.com/RainbowSeeker/rain_ws.git --recursive

# Utils
sudo apt install ros-humble-control-toolbox
git clone https://github.com/eProsima/Micro-XRCE-DDS-Agent.git
cd Micro-XRCE-DDS-Agent
mkdir build
cd build
cmake ..
make
sudo make install
sudo ldconfig /usr/local/lib/
```
## How to use
### Multi-quadrotor sil simulation
```
source /opt/ros/humble/setup.bash
colcon build
source $PWD/install/setup.bash
ros2 launch mqsls x500_3dof_mqsls_launch.py
```
### Multi-quadrotor actual test
Each onboard computer:
1. connect wifi hotspots
```
sudo nmcli dev wifi list                    # look for nearby wifi 
sudo nmcli dev wifi connect 'wifi_name' --ask
sudo nmcli connection down 'wifi_name'      # stop connection
sudo nmcli connection modify 'wifi_name' ipv4.addresses 192.168.1.11/24 # addresses format: 192.168.1.1x, from 1 to 3. 10 is reserved for gcs.
sudo nmcli connection modify 'wifi_name' ipv4.method manual
sudo nmcli connection up 'wifi_name'        # start connection
```
2. start ros2 node
```
ros2 launch mqsls x500_3dof_mqsls_prelaunch.py &
ros2 launch mqsls x500_3dof_mqsls_launch.py
```
or Using `form_tools.sh`:
```
./form_tools.sh --prelaunch
./form_tools.sh --launch
./form_tools.sh --kill
```
