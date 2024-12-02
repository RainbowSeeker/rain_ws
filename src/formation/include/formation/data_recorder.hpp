#ifndef DATA_RECORDER_HPP
#define DATA_RECORDER_HPP

#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <fstream>


template <typename T>
class DataRecorder {
public:
    explicit DataRecorder(const std::string &file_name, int throttle = 10) : 
        _file(file_name, std::ios::out), _throttle(throttle)
    {
        _buffer.reserve(_buffer_size);

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


#endif // !DATA_RECORDER_HPP