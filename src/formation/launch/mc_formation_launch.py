import os, yaml
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, TimerAction
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    # get parameters from yaml file
    param_file_path = os.path.join(get_package_share_directory('formation'), 'config', 'params.yaml')

    with open(param_file_path, 'r') as file:
        yaml_context = yaml.safe_load(file)
        params = yaml_context['/**']['ros__parameters']
        origin_lat = params['origin_lat']
        origin_lon = params['origin_lon']
        origin_alt = params['origin_alt']
    
    ekf_node = []
    for i in range(3):
        ekf_node.append(
            Node(
                package='formation',
                executable='commander_node',
                output='screen',
                shell=True,
                name='ekf_node_' + str(i + 1),
                arguments=[str(i + 1), '100000', '0', '0', '0', '0', 
                           str(origin_lat), str(origin_lon), str(origin_alt)],
            )
        )
    delay_ekf = TimerAction(period=18.0, actions=[*ekf_node])

    dds_agent = ExecuteProcess(
                    cmd=[
                        'pgrep MicroXRCEAgent > /dev/null || MicroXRCEAgent udp4 --port 8888',
                    ],
                    shell=True,
                    output='screen',
                    name='dds_agent',
                )
    
    node = []
    for i in range(3):
        node.append(
            Node(
                package='formation',
                executable='mc_formation_control_node',
                output='screen',
                shell=True,
                name='amc_' + str(i + 1),
                arguments=[str(i + 1)],
                parameters=[param_file_path]
            )
        )
    delay_amc = TimerAction(period=20.0, actions=[*node])

    px4_client = []
    for i in range(3):
        px4_workdir = os.path.expanduser('~') + '/PX4-Autopilot'
        px4_env = { 'PX4_SYS_AUTOSTART': '4001', 
                    'PX4_GZ_MODEL': 'x500', 
                    'PX4_GZ_MODEL_POSE': '0,' + str(i * 2) + ',0.5'}
        if i > 0:
            px4_env['HEADLESS'] = '1'
        px4_client.append( 
            ExecuteProcess(
                cmd=[[
                    px4_workdir + '/build/px4_sitl_default/bin/px4 -i ' + str(i + 1),
                ]],
                shell=True,
                additional_env=px4_env,
                output='screen',
                name='px4_client_' + str(i + 1),
            )
        )

    return LaunchDescription([
        dds_agent,
        *px4_client,
        delay_ekf,
        delay_amc,
    ])