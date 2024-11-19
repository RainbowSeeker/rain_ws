from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import DeclareLaunchArgument, ExecuteProcess, OpaqueFunction
from launch.substitutions import LaunchConfiguration, PythonExpression
from launch.conditions import IfCondition
from ament_index_python.packages import get_package_share_directory
import os, yaml

# get parameters from yaml file
param_file_path = os.path.join(get_package_share_directory('formation'), 'config', 'params.yaml')

def launch_setup(context, *args, **kwargs):

    amc_id = LaunchConfiguration('amc_id').perform(context)
    amc_name = '/amc_' + amc_id

    with open(param_file_path, 'r') as file:
        yaml_context = yaml.safe_load(file)
        params = yaml_context['/**']['ros__parameters']
        origin_lat = params['origin_lat']
        origin_lon = params['origin_lon']
        origin_alt = params['origin_alt']
        test_phase = yaml_context[amc_name]['ros__parameters']['test_phase']
    
    set_ekf_origin = Node(
                condition=IfCondition(PythonExpression(["'", test_phase, "' == 'formation'"])),
                package='formation',
                executable='commander_node',
                output='screen',
                shell=True,
                arguments=[LaunchConfiguration('amc_id'), '100000', '0', '0', '0', '0', 
                           str(origin_lat), str(origin_lon), str(origin_alt)],
            )
    
    return [set_ekf_origin]



def generate_launch_description():

    with open(param_file_path, 'r') as file:
        yaml_context = yaml.safe_load(file)
        params = yaml_context['/**']['ros__parameters']
        dds_baudrate = params['dds_baudrate']
        gcs_ip = params['gcs_ip']
    
    # declare launch arguments
    args = [DeclareLaunchArgument('amc_id', default_value='1', description='AMC ID'),
            DeclareLaunchArgument('dds_baudrate', default_value=str(dds_baudrate), description='DDS baudrate'),
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

    return LaunchDescription([
        *args,
        OpaqueFunction(function=launch_setup),
        dds_agent,
        gcs_map,
    ])