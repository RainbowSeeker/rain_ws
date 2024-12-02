import os, yaml
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import LogInfo, ExecuteProcess, DeclareLaunchArgument, RegisterEventHandler, GroupAction
from launch.substitutions import LaunchConfiguration, PythonExpression, PathJoinSubstitution
from ament_index_python.packages import get_package_share_directory
from launch.event_handlers import OnProcessExit

def generate_launch_description():

    # world_name = 'simple_pendulum'
    world_name = 'mqsls'
    args = [DeclareLaunchArgument('world_name', default_value=world_name, description='World name'),]

    model_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'models')
    world_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'worlds')
    plugin_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'plugins')
    
    full_world_path = PythonExpression(["'", LaunchConfiguration('world_name'), ".sdf'"])
    gz_client = ExecuteProcess(
                    cmd=[
                        "ps -ef | grep 'gz sim' | grep -v grep | awk '{print $2}' | xargs -r kill -9",
                        '&& gz sim', full_world_path, '--verbose 1',
                    ],
                    shell=True,
                    additional_env={'GZ_SIM_RESOURCE_PATH': model_dir + ':' + world_dir,
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
                parameters=[{
                    'world_name': LaunchConfiguration('world_name'),
                }],
            )
        )

    # run pkgconfig script
    config_dir = os.path.join(get_package_share_directory('mqsls'), 'config')
    pkgconfig = ExecuteProcess(
        cmd=['/bin/bash', PathJoinSubstitution([config_dir, 'pkgconfig.sh']), model_dir],
        shell=True,
        on_exit=[LogInfo(msg='pkgconfig script exited')],
    )

    # register event handler
    all_actions = RegisterEventHandler(
        event_handler=OnProcessExit(
            target_action=pkgconfig,
            on_exit=[GroupAction(actions=[
                *args,
                gz_client, 
                *node,
            ])],
        ),
    )

    return LaunchDescription([
        pkgconfig,
        all_actions,
    ])