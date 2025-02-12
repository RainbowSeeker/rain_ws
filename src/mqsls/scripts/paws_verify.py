import os
import re
import math
import time
import subprocess
import xml.etree.ElementTree as ET
from ament_index_python.packages import get_package_share_directory

class PawsVerify:
    def __init__(self, log_path):
        self.pwas_seq = [True, False]
        # self.pwas_seq = [False, True]
        self.disturb_list = [5.0, 3.0, 2.0, 1.0]
        self.cycle_count = 20
        self.target_pose = (0.0, 50.0, 10.0)
        self.log_path = log_path

        # clear log file
        with open(self.log_path, 'w') as f:
            f.write('') 
    
    def get_distance(self):
        try:
            result = subprocess.run(f'gz topic -e -t /world/x500_3dof_mqsls/pose/info | grep -A 6 "payload" | sed "/--/q"', 
                                    capture_output=True,
                                    text=True,
                                    timeout=2,
                                    shell=True).stdout
        except subprocess.TimeoutExpired:
            return 1000
        
        pattern = r'name: "payload"\s*id: 8\s*position\s*{\s*x:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*y:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*z:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*}'
        matches = re.findall(pattern, result)
        if matches:
            now_pose = matches[0]
            distance = math.sqrt((float(now_pose[0]) - self.target_pose[0]) ** 2 + 
                                 (float(now_pose[1]) - self.target_pose[1]) ** 2 + 
                                 (float(now_pose[2]) - self.target_pose[2]) ** 2)
            return distance
        else:
            return 1000
    
        
    def run_once(self, seq, disturb, cycle):
        disturb_vec = f'{disturb * math.cos(cycle * 2 * math.pi / self.cycle_count)} {disturb * math.sin(cycle * 2 * math.pi / self.cycle_count)} 0'

        self.dump_log(f'### PWAS: {"On" if seq else "Off"}, Disturb: {disturb}, Cycle: {cycle + 1}, Disturb Vector: {disturb_vec} ###\n')

        # modify disturbance of sdf file
        sdf_path = os.path.join(get_package_share_directory('mqsls'), 'gz', 'worlds', 'x500_3dof_mqsls.sdf')
        tree = ET.parse(sdf_path)
        root = tree.getroot()
        for linear_velocity in root.iter('linear_velocity'):
            linear_velocity.text = disturb_vec
        tree.write(sdf_path)

        # kill old process
        subprocess.run('pkill -9 --full "ros2 launch"', shell=True)

        # run the simulation
        # subprocess.run('colcon build --packages-select mqsls', shell=True)
        visual = True
        subprocess.Popen(f'ros2 launch mqsls sil_3dof_mqsls_launch.py traj_type:=line {"gzsim_options:=-s" * (not visual)} accs_enable:=False pwas_enable:={seq}',
                         shell=True,
                         stdout=open(os.devnull, 'w'),
                         stderr=open(os.devnull, 'w'),
                         text=True)

        # check status
        start_time = time.time()
        print_time = start_time + 5
        while True:
            time.sleep(1)

            # is exit ?
            if subprocess.run('ps -ef | grep "gz sim" | grep -v grep', capture_output=True, text=True, shell=True).returncode != 0:
                self.dump_log(f'Cycle {cycle + 1}: Simulation Exit\n')
                break

            # is success ?
            distance = self.get_distance()
            if distance < 1:
                self.dump_log(f'Cycle {cycle + 1}: Success\n')
                break
            elif time.time() > print_time:
                print(f'Distance: {distance:.2f}')
                print_time += 5

            # is timeout ?
            if time.time() - start_time > 60:
                self.dump_log(f'Cycle {cycle + 1}: Timeout\n')
                break
    
    def dump_log(self, content):
        print(content, end='')
        with open(self.log_path, 'a') as f:
            f.write(content)

    def main(self):
        for seq in self.pwas_seq:
            for disturb in self.disturb_list:
                for cycle in range(self.cycle_count):
                    self.run_once(seq, disturb, cycle)   

if __name__ == '__main__':
    log_path = os.path.join(os.path.dirname(__file__), 'paws_verify.log')
    pv = PawsVerify(log_path)
    pv.main()