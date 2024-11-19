#include <iostream>
#include <boost/asio.hpp>
#include <boost/bind/bind.hpp>
#include <filesystem>
#include <thread>

using namespace boost::asio;

class SerialUdpBridge
{
public:
    SerialUdpBridge(const std::string &serial_port, unsigned int baud_rate, const std::string &udp_ip, unsigned short udp_port)
        : _serial_port(_io, serial_port),
        _udp_endpoint(ip::address::from_string(udp_ip), udp_port),
        _udp_socket(_io, _udp_endpoint.protocol())
    {
        _serial_port.set_option(serial_port::baud_rate(baud_rate));
        _serial_port.set_option(serial_port::flow_control(serial_port::flow_control::none));
        _serial_port.set_option(serial_port::parity(serial_port::parity::none));
        _serial_port.set_option(serial_port::stop_bits(serial_port::stop_bits::one));
        _serial_port.set_option(serial_port::character_size(8));
    }

    ~SerialUdpBridge()
    {
        _io.stop();
    }

    void start()
    {
        if (open_mavlinkv2())
        {
            throw std::runtime_error("Error opening MAVLink v2");
        }
        std::cout << "Opened MAVLink v2 successfully" << std::endl;

        do_read_serial();

        do_read_udp();

        _io.run();
    }

private:
    int open_mavlinkv2()
    {
        const char mavlinkv2[] = {(char)0xfd, 0x09, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
        int retry = 10;
        while (retry > 0)
        {
            int nwrite = _serial_port.write_some(buffer(mavlinkv2, sizeof(mavlinkv2)));
            if (nwrite == sizeof(mavlinkv2))
            {
                break;
            }
            retry--;
            usleep(100 * 1000);
        }
        if (retry == 0)
        {
            return -1;
        }
        retry = 100;
        int wait_bytes = 100;
        while (wait_bytes > 0 && retry > 0)
        {
            char buf[128];
            int nread = _serial_port.read_some(buffer(buf, sizeof(buf)));
            if (nread > 0)
            {
                wait_bytes -= nread;
            }
            retry--;
            usleep(20 * 1000);
        }
        return retry > 0 ? 0 : -1;
    }

    void do_read_serial()
    {
        _serial_port.async_read_some(buffer(_serial_data), [this](const boost::system::error_code &error, size_t bytes_transferred) {
            if (!error)
            {
                auto buf = std::make_shared<std::vector<char>>(_serial_data.begin(), _serial_data.begin() + bytes_transferred);
                _udp_socket.async_send_to(buffer(*buf), _udp_endpoint, [buf](const boost::system::error_code &error, size_t) {
                    if (error)
                    {
                        throw std::runtime_error("Error sending data to UDP socket");
                    }
                });

                do_read_serial();
            }
            else
            {
                throw std::runtime_error("Error reading serial port");
            }
        });
    }

    void do_read_udp()
    {
        _udp_socket.async_receive_from(buffer(_udp_data), _udp_endpoint, [this](const boost::system::error_code &error, size_t bytes_transferred) {
            if (!error)
            {
                auto buf = std::make_shared<std::vector<char>>(_udp_data.begin(), _udp_data.begin() + bytes_transferred);
                async_write(_serial_port, buffer(*buf), [buf](const boost::system::error_code &error, size_t) {
                    if (error)
                    {
                        throw std::runtime_error("Error writing data to serial port");
                    }
                });

                do_read_udp();
            }
            else
            {
                throw std::runtime_error("Error reading UDP socket");
            }
        });
    }

    io_service _io;
    serial_port _serial_port;
    ip::udp::endpoint _udp_endpoint;
    ip::udp::socket _udp_socket;
    std::array<char, 1024> _serial_data;
    std::array<char, 1024> _udp_data;
};

std::string find_serial_port()
{
    for (const auto & entry : std::filesystem::directory_iterator("/dev/serial/by-id/"))
    {
        auto port = entry.path().string();
        if (port.find("usb-CUAV_PX4") != std::string::npos)
        {
            return port;
        }
    }
    return "";
}

int main(int argc, const char** argv) 
{
    std::string dst_ip = argc > 1 ? argv[1] : "127.0.0.1";
    int dst_port = argc > 2 ? std::stoi(argv[2]) : 14550;

    retry:
    try
    {
        std::string serial_port = find_serial_port();
        if (serial_port.empty())
        {
            throw std::runtime_error("Serial port not found");
        }
        
        std::cout << "Found serial port: " << serial_port << std::endl;
        auto bridge = std::make_shared<SerialUdpBridge>(serial_port, 2000000, dst_ip, dst_port);
        
        bridge->start();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Runtime Error: " << e.what() << std::endl;
        
        std::this_thread::sleep_for(std::chrono::seconds(2));
        goto retry;
    }

    return 0;
}