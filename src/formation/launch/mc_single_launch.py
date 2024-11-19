import os
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    # get parameters from yaml file
    param_file_path = os.path.join(get_package_share_directory('formation'), 'config', 'params.yaml')

    # declare launch arguments
    args = [DeclareLaunchArgument('amc_id', default_value='1', description='AMC ID'),]

    node = Node(
                package='formation',
                executable='mc_formation_control_node',
                output='screen',
                shell=True,
                arguments=[LaunchConfiguration('amc_id')],
                parameters=[param_file_path],
            )


    return LaunchDescription([
        *args,
        node,
    ])