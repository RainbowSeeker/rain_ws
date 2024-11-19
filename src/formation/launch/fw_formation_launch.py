import os
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, TimerAction
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    # get parameters from yaml file
    param_file_path = os.path.join(get_package_share_directory('formation'), 'config', 'params.yaml')

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
                executable='fw_formation_control_node',
                output='screen',
                shell=True,
                arguments=[str(i + 1)],
                parameters=[param_file_path]
            )
        )
    delay_amc = TimerAction(period=20.0, actions=[*node])
    
    px4_client = []
    for i in range(3):
        px4_workdir = os.path.expanduser('~') + '/PX4-Autopilot'
        px4_env = { 'PX4_SYS_AUTOSTART': '4003', 
                    'PX4_GZ_MODEL': 'rc_cessna', 
                    'PX4_GZ_MODEL_POSE': '0,' + str(i * 20) + ',0.5'}
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
            )
        )

    return LaunchDescription([
        dds_agent,
        *px4_client,
        delay_amc,
    ])