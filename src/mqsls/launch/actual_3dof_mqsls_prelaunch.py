import os, yaml
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess, DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    # get parameters from yaml file
    param_file_path = os.path.join(get_package_share_directory('mqsls'), 'config', 'params.yaml')

    with open(param_file_path, 'r') as file:
        yaml_context = yaml.safe_load(file)
        params = yaml_context['/**']['ros__parameters']
        dds_baudrate = params['dds_baudrate']
        gcs_ip = params['gcs_ip']

    # declare launch arguments
    args = [DeclareLaunchArgument('dds_baudrate', default_value=str(dds_baudrate), description='DDS baudrate'),
            DeclareLaunchArgument('gcs_ip', default_value=str(gcs_ip), description='target GCS IP address')]

    dds_agent = ExecuteProcess(
                    cmd=[
                        'pgrep MicroXRCEAgent > /dev/null || MicroXRCEAgent serial --dev /dev/ttyS0 --baudrate', 
                        LaunchConfiguration('dds_baudrate'),
                    ],
                    shell=True,
                    output='screen',
                    name='dds_agent',
                )
    
    gcs_map = Node(
                package='formation',
                executable='serial_udp_bridge',
                output='screen',
                shell=True,
                arguments=[LaunchConfiguration('gcs_ip'), '14550'],
            )
    
    mocap_bridge = Node(
                    package='mqsls',
                    executable='ros2_mocap_bridge',
                    output='screen',
                    shell=True,
                    name='mocap_bridge',
                )
    
    return LaunchDescription([
        *args,
        dds_agent,
        gcs_map,
        mocap_bridge,
    ])