# PCL Detection

## How to install
- Only support PCL **`v1.15.0`**
```c
export PCL_TARGET_VERSION=1.15.0
wget -qO pcl.tar.gz https://github.com/PointCloudLibrary/pcl/releases/download/pcl-${PCL_TARGET_VERSION}/source.tar.gz
tar -xvf pcl.tar.gz
cd pcl && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_QT=QT5 ..
make -j$(nproc)
sudo make install
```