import os, yaml
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, DeclareLaunchArgument
from launch.conditions import IfCondition
from launch.substitutions import LaunchConfiguration, PythonExpression
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    world_name = 'simple_pendulum.sdf'
    # world_name = 'mqsls.sdf'
    args = [DeclareLaunchArgument('world_name', default_value=world_name, description='World name'),]
    

    model_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'models')
    world_path = os.path.join(get_package_share_directory('mqsls'), 'gz', 'worlds', LaunchConfiguration('world_name'))
    plugin_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'plugins')

    gz_client = ExecuteProcess(
                    cmd=[
                        "ps -ef | grep 'gz sim' | grep -v grep | awk '{print $2}' | xargs -r kill -9",
                        '&& gz sim ' + world_path + ' --verbose 1',
                    ],
                    shell=True,
                    additional_env={'GZ_SIM_RESOURCE_PATH': model_dir,
                                    'GZ_SIM_SYSTEM_PLUGIN_PATH': plugin_dir},
                    output='screen',
                    name='gz_client',
                )
    
    node = []
    for i in range(1):
        node.append(
            Node(
                package='mqsls',
                executable='ideal_mqsls_control_node',
                output='screen',
                shell=True,
                name='ideal_mqsls_control_node_' + str(i + 1),
                arguments=[str(i + 1)],
                
            )
        )

    return LaunchDescription([
        *args,
        gz_client,
        *node,
    ])