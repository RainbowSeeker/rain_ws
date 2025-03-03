import os, yaml
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.conditions import IfCondition
from launch.actions import ExecuteProcess, DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration, PythonExpression
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():

    # get parameters from yaml file
    param_file_path = os.path.join(get_package_share_directory('mqsls'), 'config', 'params.yaml')

    with open(param_file_path, 'r') as file:
        yaml_context = yaml.safe_load(file)
        params = yaml_context['/**']['ros__parameters']
        dds_baudrate = params['dds_baudrate']
        gcs_ip = params['gcs_ip']
        task_allocation = params['task_allocation']
        plugin_id = task_allocation['plugin']

    # declare launch arguments
    args = [DeclareLaunchArgument('amc_id', default_value='1', description='AMC ID'),
            DeclareLaunchArgument('dds_baudrate', default_value=str(dds_baudrate), description='DDS baudrate'),
            DeclareLaunchArgument('gcs_ip', default_value=str(gcs_ip), description='target GCS IP address'),
            DeclareLaunchArgument('hil_enable', default_value='False', description='HIL simulation enable'),]

    dds_agent = ExecuteProcess(
        cmd=[
            'pgrep MicroXRCEAgent > /dev/null || MicroXRCEAgent serial --dev /dev/ttyS0 --baudrate', 
            LaunchConfiguration('dds_baudrate'),
        ],
        shell=True,
        output='screen',
        name='dds_agent',
        condition=IfCondition(PythonExpression(["'", LaunchConfiguration('hil_enable'), "' == 'False'"])),
    )
    
    gcs_map = Node(
        package='formation',
        executable='serial_udp_bridge',
        output='screen',
        shell=True,
        arguments=[LaunchConfiguration('gcs_ip'), '14550'],
        condition=IfCondition(PythonExpression(["'", LaunchConfiguration('hil_enable'), "' == 'False'"])),
    )
    
    px4_actuator = Node(
        package='mqsls',
        executable='px4_actuator_node',
        output='screen',
        shell=True,
        name=['px4_actuator', LaunchConfiguration('amc_id')],
        parameters=[{
            'amc_id': LaunchConfiguration('amc_id'),
        }],
    )

    # mocap_plugin = Node(
    #     package='mqsls',
    #     executable='mocap_plugin_node',
    #     output='screen',
    #     shell=True,
    #     name='mocap_plugin',
    #     condition=IfCondition(PythonExpression(["'", LaunchConfiguration('amc_id'), "' == '", str(plugin_id), "'"])),
    # )
    
    return LaunchDescription([
        *args,
        dds_agent,
        gcs_map,
        px4_actuator,
        # mocap_plugin,
    ])