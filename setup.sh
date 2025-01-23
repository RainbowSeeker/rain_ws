
sudo apt-get update -y --quiet
sudo DEBIAN_FRONTEND=noninteractive apt-get -y --quiet --no-install-recommends install \
                python3-pip \
                libboost-system-dev \
                libboost-filesystem-dev \
                libeigen3-dev \
                sshpass \
                ros-humble-control-toolbox

pip install empy==3.3.4 numpy lark future
