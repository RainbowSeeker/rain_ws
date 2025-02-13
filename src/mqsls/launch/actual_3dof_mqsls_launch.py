import os, yaml
from launch import LaunchDescription
from launch_ros.actions import Node
from launch.conditions import IfCondition
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration, PythonExpression
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():
    
    # get parameters from yaml file
    param_file_path = os.path.join(get_package_share_directory('mqsls'), 'config', 'params.yaml')

    with open(param_file_path, 'r') as file:
        yaml_context = yaml.safe_load(file)
        params = yaml_context['/**']['ros__parameters']
        task_allocation = params['task_allocation']
        control_id = task_allocation['control']
        optimize_id = task_allocation['optimize']

    args = [DeclareLaunchArgument('amc_id', default_value='1', description='AMC ID'),
            DeclareLaunchArgument('eso_enable', default_value='True', description='Enable ESO'),
            DeclareLaunchArgument('pwas_enable', default_value='True', description='Enable PWAS'),
            DeclareLaunchArgument('accs_enable', default_value='True', description='Enable ACCS'),
            DeclareLaunchArgument('traj_type', default_value='circle', description='Trajectory type')]

    controller = Node(
        package='mqsls',
        executable='3dof_mqsls_controller_node',
        output='screen',
        shell=True,
        name='Controller',
        parameters=[{
            'cable_len': 2.0,
            'load_mass': 0.5,
            'uav_mass': 1.6,
            'hover_thrust': 0.42,
            'kp': 0.5,
            'kv': 1.0,
            'kq': 0.5,
            'kw': 3.0,
            'eso_enable':  LaunchConfiguration('eso_enable'),
            'pwas_enable': LaunchConfiguration('pwas_enable'),
            'accs_enable': LaunchConfiguration('accs_enable'),
            'min_tension': 0.5,
            'max_tension': 5.0,
            'traj_type': LaunchConfiguration('traj_type'), # 'line', 'circle', 'rectangle', 'lissajous'
            'lasting_time': 70,
        }],
        condition=IfCondition(PythonExpression(["'", LaunchConfiguration('amc_id'), "' == '", str(control_id), "'"])),
    )

    force_planner = Node(
        package='mqsls',
        executable='force_planner_node',
        output='screen',
        shell=True,
        name='force_planner',
        condition=IfCondition(PythonExpression(["'", LaunchConfiguration('amc_id'), "' == '", str(optimize_id), "'"])),
    )
    
    return LaunchDescription([
        *args,
        controller,
        force_planner,
    ])