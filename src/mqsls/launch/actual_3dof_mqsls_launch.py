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

    args = [DeclareLaunchArgument('amc_id', default_value='1', description='AMC ID'),]

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
            'eso_enable': True,
            'kq': 0.5,
            'kw': 3.0,
            'min_tension': 0.0,
            'max_tension': 7.0,
            'traj_type': 'line',
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