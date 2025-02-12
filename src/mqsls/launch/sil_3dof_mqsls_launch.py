import os
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, DeclareLaunchArgument, RegisterEventHandler, GroupAction, TimerAction
from launch.substitutions import LaunchConfiguration, PythonExpression, PathJoinSubstitution
from ament_index_python.packages import get_package_share_directory
from launch.event_handlers import OnProcessExit

def generate_launch_description():

    world_name = 'x500_3dof_mqsls'
    args = [DeclareLaunchArgument('world_name', default_value=world_name, description='World name'),
            DeclareLaunchArgument('gzsim_options', default_value='', description='Additional options to pass to gz sim'),
            DeclareLaunchArgument('eso_enable', default_value='True', description='Enable ESO'),
            DeclareLaunchArgument('pwas_enable', default_value='True', description='Enable PWAS'),
            DeclareLaunchArgument('accs_enable', default_value='True', description='Enable ACCS'),
            DeclareLaunchArgument('traj_type', default_value='circle', description='Trajectory type'),]

    model_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'models')
    world_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'worlds')

    full_world_path = PythonExpression(["'", LaunchConfiguration('world_name'), ".sdf'"])
    gz_sim = ExecuteProcess(
        cmd=[
            'gz sim --verbose=1 -r', LaunchConfiguration('gzsim_options'), full_world_path,
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
            'disturbance_enable': True,
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
    
    controller = Node(
        package='mqsls',
        executable='3dof_mqsls_controller_node',
        output='screen',
        shell=True,
        name='Controller',
        parameters=[{
            'cable_len': 2.0,
            'load_mass': 1.0,
            'uav_mass': 2.064307692307692,
            'hover_thrust': 0.74,
            'kp': 1.5,
            'kv': 2.0,
            'kq': 0.5,
            'kw': 3.0,
            'eso_enable':  LaunchConfiguration('eso_enable'),
            'pwas_enable': LaunchConfiguration('pwas_enable'),
            'accs_enable': LaunchConfiguration('accs_enable'),
            'min_tension': 2.0,
            'max_tension': 7.0,
            'traj_type': LaunchConfiguration('traj_type'), # 'line', 'circle', 'rectangle', 'lissajous'
            'lasting_time': 70,
        }],
    )

    force_planner = Node(
        package='mqsls',
        executable='force_planner_node',
        output='screen',
        shell=True,
        name='force_planner',
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
                controller,
                force_planner,
                *px4_client,
                dds_agent,
            ])],
        ),
    )

    return LaunchDescription([
        pkgconfig,
        all_actions,
    ])