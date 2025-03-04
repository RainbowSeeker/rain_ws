# README FOR FORMATION

## Directory Structure
- config: configuration files for setup.
- launch: launch files for running the formation.
- msg: message files for communication between nodes.
- include/src: source code.

## Nodes
- mc_formation_control: control node for multicopter formation.

## How to Run
1. Build the package.
```bash
colcon build --packages-select formation
source install/setup.bash
```
2. SIL Simulation
```bash
ros2 launch formation mc_formation_launch.py
```
