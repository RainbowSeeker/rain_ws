import os
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, DeclareLaunchArgument, RegisterEventHandler, GroupAction, TimerAction
from launch.substitutions import LaunchConfiguration, PythonExpression, PathJoinSubstitution
from ament_index_python.packages import get_package_share_directory
from launch.event_handlers import OnProcessExit

def generate_launch_description():

    world_name = 'x500_3dof_mqsls'
    args = [DeclareLaunchArgument('world_name', default_value=world_name, description='World name'),]

    model_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'models')
    world_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'worlds')
    
    full_world_path = PythonExpression(["'", LaunchConfiguration('world_name'), ".sdf'"])
    gz_sim = ExecuteProcess(
        cmd=[
            'gz sim', full_world_path, '--verbose 1',
        ],
        shell=True,
        additional_env={'GZ_SIM_RESOURCE_PATH': model_dir + ':' + world_dir},
        output='screen',
        name='gz_sim',
    )

    gz_plugin = Node(
        package='mqsls',
        executable='gz_plugin_node',
        output='screen',
        shell=True,
        name='gz_plugin',
        parameters=[{
            'world_name': LaunchConfiguration('world_name'),
        }],
    )
    
    dds_agent = ExecuteProcess(
        cmd=[
            'pgrep MicroXRCEAgent > /dev/null || MicroXRCEAgent udp4 --port 8888',
        ],
        shell=True,
        output='screen',
        name='dds_agent',
    )
    
    px4_client = []
    for i in range(3):
        px4_workdir = os.environ['PX4_HOME']
        px4_env = { 'PX4_SYS_AUTOSTART': '4001', 
                    'PX4_GZ_MODEL_NAME': 'x500_' + str(i + 1)}
        px4_action = ExecuteProcess(
            cmd=[[
                px4_workdir + '/build/px4_sitl_default/bin/px4 -i ' + str(i + 1),
            ]],
            shell=True,
            additional_env=px4_env,
            output='screen',
            name='px4_client_' + str(i + 1),
        )
        px4_client.append(
            px4_action
            # TimerAction(period=float(2), actions=[px4_action])
        )

    # run pkgconfig script
    config_dir = os.path.join(get_package_share_directory('mqsls'), 'config')
    pkgconfig = ExecuteProcess(
        cmd=['/usr/bin/bash', PathJoinSubstitution([config_dir, 'pkgconfig.sh']), model_dir],
        shell=True,
    )

    # register event handler
    all_actions = RegisterEventHandler(
        event_handler=OnProcessExit(
            target_action=pkgconfig,
            on_exit=[GroupAction(actions=[
                *args,
                gz_sim,
                gz_plugin,
                *px4_client,
                dds_agent,
            ])],
        ),
    )

    return LaunchDescription([
        pkgconfig,
        all_actions,
    ])