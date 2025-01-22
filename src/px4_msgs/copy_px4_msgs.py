import os
import re

# Parse the topics from the DDS topics file
# Format:
#     type: px4_msgs::msg::xxx
def parse_dds_topics(topic_file: str, msg_dir: str):
    with open(topic_file, 'r') as f:
        lines = f.readlines()
    topics = []
    for line in lines:
        if re.match('^[^#]*px4_msgs::msg::', line):
            topics.append(line.split('::')[-1].strip())
    # detect if there are custom types
    BASE_TYPES = ['bool', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64', 'float32', 'float64']
    for topic in topics:
        with open(os.path.join(msg_dir, topic + '.msg'), 'r') as f:
            lines = f.readlines()
        for line in lines:
            clean_line = line.split('#')[0].strip()
            if clean_line == '':
                continue

            type = clean_line.split(' ')[0].split('[')[0]
            if type not in BASE_TYPES:
                if type not in topics:
                    topics.append(type)
                    # print(f'Custom type detected: {type}')
    return topics


def copy_px4_msgs(px4_home: str, yaml_file: str, dst_dir: str):
    msg_dir = os.path.join(px4_home, 'msg')

    # preprocess dds_topics.yaml
    topics = parse_dds_topics(yaml_file, msg_dir)
    
    # copy all the message files
    for topic in topics:
        src_file = os.path.join(msg_dir, topic + '.msg')
        dst_file = os.path.join(dst_dir, topic + '.msg')
        os.system(f'cp {src_file} {dst_file}')
    print('PX4 message files copied successfully!')

if __name__ == '__main__':
    px4_home = os.environ['PX4_HOME']
    if px4_home == '':
        raise Exception('PX4_HOME environment variable not set!')

    yaml_file = os.path.join(os.path.dirname(__file__), 'dds_topics.yaml')
    dst_dir = os.path.join(os.path.dirname(__file__), 'msg')
    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)
    copy_px4_msgs(px4_home, yaml_file, dst_dir)