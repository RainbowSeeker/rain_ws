import os

class RopeGenerator:
    def __init__(self, rope_length: float, seg_num: int, 
                mass: float = 0, collision_free: bool = True):
        self.rope_length = rope_length
        self.seg_num = seg_num
        self.seg_length = rope_length / seg_num

        self.mass_content = f"""
      <inertial>
        <mass>{max(mass, 0.000000000001) / seg_num}</mass>
        <inertia>
          <ixx>0.01</ixx>
          <iyy>0.01</iyy>
          <izz>0.01</izz>
        </inertia>
      </inertial>"""
        self.collision_free = collision_free

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
    <pose>0 0 2 0 0 0</pose>
    <static>false</static>
    <link name='begin_link'>
      <pose>0 0 0 0 0 0</pose>{self.mass_content}
      {self.collision_free * "<!-- "} <collision name='begin_link_collision'>
        <geometry>
          <sphere>
            <radius>0.003</radius>
          </sphere>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='begin_link_visual'>
        <geometry>
          <sphere>
            <radius>0.003</radius>
          </sphere>
        </geometry>
      </visual>
      <!-- <velocity_decay/> -->
    </link>
    <link name='link_0'>
      <pose>0 0 0 0 0 0</pose>{self.mass_content}
      {self.collision_free * "<!-- "} <collision name='link_0_collision'>
        <pose>{self.seg_length / 2} 0 0 0 1.570796326794896558 0</pose>
        <geometry>
          <cylinder>
            <length>{self.seg_length}</length>
            <radius>0.003</radius>
          </cylinder>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='link_0_vis'>
        <pose>{self.seg_length / 2} 0 0 0 1.570796326794896558 0</pose>
        <geometry>
          <cylinder>
            <length>{self.seg_length}</length>
            <radius>0.003</radius>
          </cylinder>
        </geometry>
      </visual>
      <!-- <velocity_decay/> -->
    </link>
    <joint name='joint_0_2' type='revolute'>
      <pose relative_to='begin_link'/>
      <child>link_0</child>
      <parent>begin_link</parent>
      <axis>
        <xyz>0 1 0</xyz>
        <limit>
          <lower>-3.14</lower>
          <upper>3.14</upper>
        </limit>
        <!--<dynamics>
          <damping>0.1</damping>
        </dynamics> -->
      </axis>
    </joint>"""

    def tail_text(self) -> str:
        return f"""
    <link name='end_link'>
      <pose>{self.rope_length} 0 0 0 0 0</pose>{self.mass_content}
      {self.collision_free * "<!-- "} <collision name='end_link_collision'>
        <geometry>
          <sphere>
            <radius>0.003</radius>
          </sphere>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='end_link_visual'>
        <geometry>
          <sphere>
            <radius>0.003</radius>
          </sphere>
        </geometry>
      </visual>
      <!-- <velocity_decay/> -->
    </link>
    <joint name='joint_end_1' type='revolute'>
      <child>end_link</child>
      <parent>link_{self.seg_num-1}</parent>
      <axis>
        <xyz>0 0 1</xyz>
        <limit>
          <lower>-3.14</lower>
          <upper>3.14</upper>
        </limit>
        <!--<dynamics>
          <damping>0.1</damping>
        </dynamics> -->
      </axis>
    </joint>
  </model>
</sdf>"""

    def repeated_body_text(self, i: int) -> str:
        return f"""
    <link name='sphere_{i}'>
      <pose>{self.seg_length * i} 0 0 0 0 0</pose>{self.mass_content}
      {self.collision_free * "<!-- "} <collision name='sphere_{i}_collision'>
        <geometry>
          <sphere>
            <radius>0.003</radius>
          </sphere>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='sphere_{i}_vis'>
        <geometry>
          <sphere>
            <radius>0.003</radius>
          </sphere>
        </geometry>
      </visual>
      <!-- <velocity_decay/> -->
    </link>
    <joint name='joint_{i}_1' type='revolute'>
      <child>sphere_{i}</child>
      <parent>link_{i-1}</parent>
      <axis>
        <xyz>0 0 1</xyz>
        <limit>
          <lower>-3.14</lower>
          <upper>3.14</upper>
        </limit>
        <!--<dynamics>
          <damping>0.1</damping>
        </dynamics> -->
      </axis>
    </joint>
    <link name='link_{i}'>
      <pose>{self.seg_length * i} 0 0 0 0 0</pose>{self.mass_content}
      {self.collision_free * "<!-- "} <collision name='link_{i}_collision'>
        <pose>{self.seg_length / 2} 0 0 0 1.570796326794896558 0</pose>
        <geometry>
          <cylinder>
            <length>{self.seg_length}</length>
            <radius>0.003</radius>
          </cylinder>
        </geometry>
      </collision> {self.collision_free * "-->"}
      <visual name='link_{i}_vis'>
        <pose>{self.seg_length / 2} 0 0 0 1.570796326794896558 0</pose>
        <geometry>
          <cylinder>
            <length>{self.seg_length}</length>
            <radius>0.003</radius>
          </cylinder>
        </geometry>
      </visual>
      <!-- <velocity_decay/> -->
    </link>
    <joint name='joint_{i}_2' type='revolute'>
      <child>link_{i}</child>
      <parent>sphere_{i}</parent>
      <axis>
        <xyz>0 1 0</xyz>
        <limit>
          <lower>-3.14</lower>
          <upper>3.14</upper>
        </limit>
        <!--<dynamics>
          <damping>0.1</damping>
        </dynamics> -->
      </axis>
    </joint>"""


if __name__ == '__main__':
    rope_length = 1 # m
    seg_num = 10
    rope_gen = RopeGenerator(rope_length, seg_num)

    # save to current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, 'model.sdf'), 'w') as f:
        f.write(rope_gen.gen_rope())
