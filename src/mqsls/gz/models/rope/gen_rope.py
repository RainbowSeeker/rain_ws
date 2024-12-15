import os

class RopeGenerator:
    def __init__(self, rope_length: float, seg_num: int, 
                mass: float = 0, collision_free: bool = False):
        self.rope_length = rope_length
        self.seg_num = seg_num
        self.seg_length = rope_length / seg_num
        self.radius = 0.01
        self.seg_mass = max(mass, 0.01) / seg_num
        self.collision_free = collision_free

        self.gripper_mass = 0.001
        self.gripper_inertial = f"""
      <inertial>
        <mass>{self.gripper_mass}</mass>
        <inertia>
          <ixx>{self.gripper_mass*0.4*self.radius*self.radius}</ixx>
          <iyy>{self.gripper_mass*0.4*self.radius*self.radius}</iyy>
          <izz>{self.gripper_mass*0.4*self.radius*self.radius}</izz>
        </inertia>
      </inertial>"""
        
        self.link_inertial = f"""
      <inertial>
        <mass>{self.seg_mass}</mass>
        <inertia>
          <ixx>{1e2*self.seg_mass / 12 * (3*self.radius*self.radius + self.seg_length * self.seg_length)}</ixx>
          <iyy>{1e2*self.seg_mass / 12 * (3*self.radius*self.radius + self.seg_length * self.seg_length)}</iyy>
          <izz>{1e2*0.5*self.seg_mass*self.radius*self.radius}</izz>
        </inertia>
      </inertial>"""
        
        self.sensor_text = """
      <sensor name="imu_sensor" type="imu">
        <always_on>1</always_on>
        <update_rate>250</update_rate>
        <visualize>true</visualize>
      </sensor>"""

    def gen_rope(self):
        result = self.head_text()
        for i in range(1, self.seg_num):
            result += self.repeated_body_text(i)
        result += self.tail_text()
        return result

    def head_text(self) -> str:
        return f"""<?xml version="1.0" encoding="UTF-8"?>
<sdf version='1.9'>
  <model name='rope'>
    <enable_wind>true</enable_wind>
    <static>false</static>
    <link name='begin_link'>
      <pose>0 0 0 0 0 0</pose>{self.gripper_inertial}
      {self.collision_free * "<!-- "} <collision name='begin_link_collision'>
        <geometry>
          <sphere>
            <radius>{self.radius}</radius>
          </sphere>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='begin_link_visual'>
        <geometry>
          <sphere>
            <radius>{self.radius}</radius>
          </sphere>
        </geometry>
      </visual>
    </link>{self.seg_link_text(0)}
    <joint name='begin_joint_2' type='universal'>
      <pose relative_to='begin_link'/>
      <child>link_0</child>
      <parent>begin_link</parent>
      <axis>
        <xyz>1 0 0</xyz>
      </axis>
      <axis2>
        <xyz>0 1 0</xyz>
      </axis2>
    </joint>
    """ 

    def tail_text(self) -> str:
        return f"""
    <link name='end_link'>
      <pose>{self.rope_length} 0 0 0 0 0</pose>{self.gripper_inertial}
      {self.collision_free * "<!-- "} <collision name='end_link_collision'>
        <geometry>
          <sphere>
            <radius>{self.radius}</radius>
          </sphere>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='end_link_visual'>
        <geometry>
          <sphere>
            <radius>{self.radius}</radius>
          </sphere>
        </geometry>
      </visual>
    </link>
    <joint name='joint_end_1' type='universal'>
      <pose relative_to='end_link'/>
      <child>end_link</child>
      <parent>link_{self.seg_num-1}</parent>
      <axis>
        <xyz>0 0 1</xyz>
      </axis>
      <axis2>
        <xyz>0 1 0</xyz>
      </axis2>
    </joint>
  </model>
</sdf>"""

    def seg_link_text(self, i: int) -> str:
        return f"""
    <link name='link_{i}'>
      <pose>{self.seg_length * i + self.seg_length / 2} 0 0 0 1.570796326794896558 0</pose>{self.link_inertial}
      {self.collision_free * "<!-- "} <collision name='link_{i}_collision'>
        <geometry>
          <cylinder>
            <length>{self.seg_length}</length>
            <radius>{self.radius}</radius>
          </cylinder>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='link_{i}_vis'>
        <geometry>
          <cylinder>
            <length>{self.seg_length}</length>
            <radius>{self.radius}</radius>
          </cylinder>
        </geometry>
        <material>
          <ambient>0 1 0 1</ambient> <!-- green -->
          <diffuse>0 1 0 1</diffuse>
          <specular>0 1 0 1</specular>
        </material>
      </visual>
      <velocity_decay>
        <linear>0.1</linear>
        <angular>0.1</angular>
      </velocity_decay>
      {(i == 1) * self.sensor_text}
    </link>"""

    def seg_joint_text(self, i: int) -> str:
        return f"""
    <joint name='joint_{i}_2' type='universal'>
      <pose relative_to='link_{i-1}'>0 0 {self.seg_length / 2} 0 0 0</pose>
      <child>link_{i}</child>
      <parent>link_{i-1}</parent>
      <axis>
        <xyz expressed_in='link_{i-1}'>1 0 0</xyz>
        <limit>
          <lower>-3.141592653589793</lower>
          <upper>3.141592653589793</upper>
          <effort>100</effort>
          <velocity>10</velocity>
        </limit>
        <dynamics>
          <damping>0.5</damping>
        </dynamics>
      </axis>
      <axis2>
        <xyz expressed_in='link_{i-1}'>0 1 0</xyz>
        <limit>
          <lower>-3.141592653589793</lower>
          <upper>3.141592653589793</upper>
          <effort>100</effort>
          <velocity>10</velocity>
        </limit>
        <dynamics>
          <damping>0.5</damping>
        </dynamics>
      </axis2>
    </joint>"""

    def repeated_body_text(self, i: int) -> str:
        return self.seg_link_text(i) + self.seg_joint_text(i)


if __name__ == '__main__':
    rope_length = 2 # m
    seg_num = 10
    rope_gen = RopeGenerator(rope_length, seg_num)

    # save to current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, 'model.sdf'), 'w') as f:
        f.write(rope_gen.gen_rope())
