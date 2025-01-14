from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, DeclareLaunchArgument

def generate_launch_description():
    
    node = Node(
                package='mqsls',
                executable='sim_mqsls_control_node',
                output='screen',
                shell=True,
                name='amc',
                parameters=[{
                    'cable_len': 2.0,
                    'load_mass': 1.0,
                    'uav_mass': 2.064307692307692,
                    'hover_thrust': 0.74,
                    'eso_enable': True,
                    'kq': 0.5,
                    'kw': 3.0,
                    'min_tension': 0.0,
                    'max_tension': 7.0,
                    'traj_type': 'line',
                }],
            )

    force_planner = Node(
        package='mqsls',
        executable='force_planner_node',
        output='screen',
        shell=True,
        name='force_planner',
    )
    
    return LaunchDescription([
        node,
        force_planner,
    ])