# 编队手册

## 系统架构
见[系统架构](./multi_uav_framework.pdf).

## 目录结构
- *config*: 存放系统默认配置文件;
- *launch*: 存放一件部署运行的 `launch` 启动脚本;
- msg: 存放自定义消息类型,用于节点间通信;
- include/src: 源码;

## 节点
- MicroXRCEAgent: px4与ros2通信的桥接`dds`节点(px4通过tel2口将数据发送给LubanCat,再由LubanCat通过`dds`协议广播到ros2网络中);
- *mc_formation_control_node*: 核心编队控制节点;
- px4_client: sil中的px4仿真节点;
- commander_node: 用于发送指令的节点(如set_ekf_origin等);
- serial_udp_bridge: px4的USB(mavlink)与主机的UDP通信桥接节点,用于主机与三个飞控的地面站通信;

## Launch
- *mc_formation_launch.py*: sil仿真启动脚本;
- *mc_single_prelaunch.py*: 真实飞行预配置脚本;
- *mc_single_launch.py*: 真实飞行启动脚本(主控制器);

## 如何使用
1. Build the package.
```bash
colcon build --packages-select formation
source install/setup.bash
```
2. Run.
```bash
# SIL Simulation
ros2 launch formation mc_formation_launch.py
# Real-world
ros2 launch formation mc_single_prelaunch.py amc_id:=1 # amc_id: 1-3
ros2 launch formation mc_single_launch.py
```
