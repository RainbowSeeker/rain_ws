#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>

namespace mqsls {

struct DataFrame {
    // timestamp
    uint64_t timestamp;
    // payload state
    Eigen::Vector3d load_position;
    Eigen::Vector3d load_velocity;

    // uav state
    Eigen::Vector3d uav_position[3];
    Eigen::Vector3d uav_velocity[3];

    // operator<<
    friend std::ostream &operator<<(std::ostream &os, const DataFrame &frame) {
        os << frame.timestamp << ' '
           << frame.load_position.x() << ' ' << frame.load_position.y() << ' ' << frame.load_position.z() << ' '
           << frame.load_velocity.x() << ' ' << frame.load_velocity.y() << ' ' << frame.load_velocity.z() << ' ';
        for (int i = 0; i < 3; i++) {
            os << frame.uav_position[i].x() << ' ' << frame.uav_position[i].y() << ' ' << frame.uav_position[i].z() << ' '
               << frame.uav_velocity[i].x() << ' ' << frame.uav_velocity[i].y() << ' ' << frame.uav_velocity[i].z() << ' ';
        }
        os << '\n';
        return os;
    }

    static std::string header() {
        std::stringstream ss;
        ss  << "timestamp" << ' '
            << "x y z" << ' '
            << "vx vy vz" << ' ';
        for (int i = 1; i <= 3; i++) {
            ss  << "x" << i << " y" << i << " z" << i << ' '
                << "vx" << i << " vy" << i << " vz" << i << ' ';
        }

        return ss.str();
    }
};


template <typename T>
class DataRecorder {
public:
    explicit DataRecorder(const std::string &file_name, int throttle = 10) : 
        _file(file_name, std::ios::out), _throttle(throttle)
    {
        _buffer.reserve(_buffer_size);

        if (!_file.is_open()) {
            throw std::runtime_error("Failed to open file: " + file_name);
        }
        _file.rdbuf()->pubsetbuf(_buffer.data(), _buffer_size);

        // write header
        _file << T::header() << std::endl;
    };

    ~DataRecorder()
    {
        _file.flush();
        _file.close();
    }

    int push(const T &frame) {
        if (!_file.good()) {
            switch (_file.rdstate()) {
                case std::ios_base::eofbit:
                    return -EOF;
                case std::ios_base::failbit:
                case std::ios_base::badbit:
                default:
                    return -EIO;
            }
        }

        // write to file
        _file << frame;

        _write_cnt++;
        if (_write_cnt % _throttle == 0 || avail_size() < _buffer_size / 2) {
            _file.flush();
        }

        return 0;
    }

    size_t avail_size() const {
        return _buffer_size - _file.rdbuf()->in_avail();
    }

private:
    // file
    std::ofstream _file;
    int _write_cnt {0};
    int _throttle;

    // buf
    size_t _buffer_size {RECORDER_MIN_BUFFER_SIZE};
    std::vector<char> _buffer;

    static const size_t RECORDER_MIN_BUFFER_SIZE = 16384;
};


};
