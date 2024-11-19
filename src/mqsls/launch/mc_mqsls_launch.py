import os, yaml
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, TimerAction
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    model_dir = os.path.join(get_package_share_directory('mqsls'), 'gz', 'models')
    world_path = os.path.join(get_package_share_directory('mqsls'), 'gz', 'worlds', 'mqsls.sdf')
    gz_client = ExecuteProcess(
                    cmd=[
                        'gz sim ' + world_path,
                    ],
                    shell=True,
                    additional_env={'GZ_SIM_RESOURCE_PATH': model_dir},
                    output='screen',
                    name='gz_client',
                )

    # dds_agent = ExecuteProcess(
    #                 cmd=[
    #                     'pgrep MicroXRCEAgent > /dev/null || MicroXRCEAgent udp4 --port 8888',
    #                 ],
    #                 shell=True,
    #                 output='screen',
    #                 name='dds_agent',
    #             )
    
    # node = []
    # for i in range(3):
    #     node.append(
    #         Node(
    #             package='mqsls',
    #             executable='mc_mqsls_control_node',
    #             output='screen',
    #             shell=True,
    #             name='mc_mqsls_control_node_' + str(i + 1),
    #             arguments=[str(i + 1)],
    #         )
    #     )
    # delay_amc = TimerAction(period=20.0, actions=[*node])

    # px4_client = []
    # for i in range(3):
    #     px4_workdir = os.path.expanduser('~') + '/PX4-Autopilot'
    #     px4_env = { 'PX4_SYS_AUTOSTART': '4001', 
    #                 'PX4_GZ_MODEL_NAME': 'x500_' + str(i + 1), 
    #                 'PX4_GZ_WORLD': 'mqsls'}
    #     if i > 0:
    #         px4_env['HEADLESS'] = '1'
    #     px4_action = ExecuteProcess(
    #             cmd=[[
    #                 px4_workdir + '/build/px4_sitl_default/bin/px4 -i ' + str(i + 1),
    #             ]],
    #             shell=True,
    #             additional_env=px4_env,
    #             output='screen',
    #             name='px4_client_' + str(i + 1),
    #         )
    #     px4_client.append(
    #         TimerAction(period= float(i * 3), actions=[px4_action])
    #     )

    return LaunchDescription([
        gz_client,
        # dds_agent,
        # delay_amc,
        # *px4_client,
    ])