#pragma once
#include <cerrno>
#include <boost/crc.hpp>

namespace unicore // driver for UM982 of unicore
{
	namespace msg
	{
		struct header_binary
		{
			uint8_t sync[3]; 	// 0xAA, 0x44, 0xB5
			uint8_t cpu_idle; 	// [0, 100]
			uint16_t msg_id;	// message id
			uint16_t msg_len;	// message length
			uint8_t time_ref;	// (GPST or BDST)
			uint8_t time_status;	// Time status
			uint16_t Wn;		// GPS week number
			uint32_t Ms;		// GPS millisecond
			uint32_t res;		// reserved
			uint8_t version;	// version
			uint8_t leap_sec;	// leap second
			uint16_t delay_ms;	// delay ms
		} __attribute__((packed));

		static_assert(sizeof(header_binary) == 24, "sizeof(header_binary) != 24");

		struct footer_binary
		{
			uint32_t crc;		// CRC
		} __attribute__((packed));

		enum class sol_status : uint32_t
		{
			SOL_COMPUTED = 0,
			INSUFFICIENT_OBS = 1,
			NO_CONVERGENCE = 2,
			SINGULARITY = 3,
			COV_TRACE = 4,
		};

		enum class ext_sol_status : uint8_t
		{
			NOT_VERIFY = 0,
			VERIFY = 1,
		};

		enum class pv_type : uint32_t
		{
			NONE = 0,
			FIXEDPOS = 1,
			FIXEDHEIGHT = 2,
			DOPPLER_VELOCITY = 8,
			SINGLE = 16,
			PSRDIFF = 17,
			SBAS = 18,
			L1_FLOAT = 32,
			IONOFREE_FLOAT = 33,
			NARROW_FLOAT = 34,
			L1_INT = 48,
			WIDE_INT = 49,
			NARROW_INT = 50,
			INS = 52,
			INS_PSRSP = 53,
			INS_PSRDIFF = 54,
			INS_RTKFLOAT = 55,
			INS_RTKFIXED = 56,
			PPP_CONVERGING = 68,
			PPP = 69,
		};

		struct bestnavxyz
		{
			sol_status 	status;		// solution status
			pv_type 	pos_type;	// position type
			double 		pos[3];		// x, y, z [m]
			float 		pos_std[3];	// x, y, z std [m]
			sol_status 	vel_status;	// velocity status
			pv_type 	vel_type;	// velocity type
			double 		vel[3];		// vx, vy, vz [m/s]
			float 		vel_std[3];	// vx, vy, vz std [m/s]
			char		stn_id[4];	// station id
			float		v_latency;	// velocity latency
			float		diff_age;	// differential age
			float		sol_age;	// solution age
			uint8_t		sv_used;	// satellites used
			uint8_t		sv_in_sol;	// satellites in solution
			uint8_t		gg_l1;		// GG L1
			uint8_t		res[2];		// reserved
			ext_sol_status	ext_sol;	// extended solution status
			uint8_t		res2[2];	// reserved
		} __attribute__((packed));

		static_assert(sizeof(bestnavxyz) == 112, "sizeof(bestnavxyz) != 112");

		struct uniheading
		{
			sol_status 	status;		// solution status
			pv_type 	pos_type;	// position type
			float		baseline;	// baseline [m]
			float		heading;	// heading [deg] [0, 360]
			float 		pitch;		// pitch [deg] [-90, 90]
			float 		res;		// reserved
			float 		heading_std;	// heading std [deg]
			float 		pitch_std;	// pitch std [deg]
			char 		stn_id[4];	// station id
			uint8_t 	sv;		// satellites tracked
			uint8_t 	sv_in_sol;	// satellites in solution
			uint8_t 	res1[3];	// reserved
			ext_sol_status 	ext_sol;	// extended solution status
			uint8_t 	res2[2];	// reserved
		} __attribute__((packed));

		static_assert(sizeof(uniheading) == 44, "sizeof(uniheading) != 44");

		enum id : uint16_t
		{
			MSG_ID_BESTNAVXYZ = 240,
			MSG_ID_BESTNAVXYZH = 242,
			MSG_ID_UNIHEADING = 972,
		};

		struct body_binary
		{
			header_binary 	header;
			uint8_t			payload[256];
			footer_binary 	footer;
		} __attribute__((packed));

	} // namespace msg

	int parse(uint8_t ch, msg::body_binary &msg)
	{
		enum parse_state {
			SYNC_1,
			SYNC_2,
			SYNC_3,
			CPU_IDLE,
			MSG_ID_1,
			MSG_ID_2,
			MSG_LEN_1,
			MSG_LEN_2,
			DONT_CARE,
			PAYLOAD,
			CRC_1,
			CRC_2,
			CRC_3,
			CRC_4,
		};
		static parse_state state = SYNC_1;
		static int dont_care_count = 0;
		static int payload_count = 0;

		auto reset = [&]() {
			state = SYNC_1;
			dont_care_count = 0;
			payload_count = 0;
		};

		switch (state) {
		case SYNC_1:
			if (ch == 0xAA)
				state = SYNC_2;
			break;
		case SYNC_2:
			if (ch == 0x44)
				state = SYNC_3;
			else
				reset();
			break;
		case SYNC_3:
			if (ch == 0xB5)
				state = CPU_IDLE;
			else
				reset();
			break;
		case CPU_IDLE:
			msg.header.cpu_idle = ch;
			state = MSG_ID_1;
			break;
		case MSG_ID_1:
			msg.header.msg_id = ch;
			state = MSG_ID_2;
			break;
		case MSG_ID_2:
			msg.header.msg_id |= ch << 8;
			state = MSG_LEN_1;
			break;
		case MSG_LEN_1:
			msg.header.msg_len = ch;
			state = MSG_LEN_2;
			break;
		case MSG_LEN_2:
			msg.header.msg_len |= ch << 8;
			state = DONT_CARE;
			break;
		case DONT_CARE:
			((uint8_t *)(&msg.header))[8 + dont_care_count] = ch;
			dont_care_count++;
			if (dont_care_count == sizeof(msg.header) - 8)
				state = PAYLOAD;
			break;
		case PAYLOAD:
			msg.payload[payload_count] = ch;
			payload_count++;
			if (payload_count == msg.header.msg_len)
				state = CRC_1;
			break;
		case CRC_1:
			msg.footer.crc = ch;
			state = CRC_2;
			break;
		case CRC_2:
			msg.footer.crc |= ch << 8;
			state = CRC_3;
			break;
		case CRC_3:
			msg.footer.crc |= ch << 16;
			state = CRC_4;
			break;
		case CRC_4:
			msg.footer.crc |= ch << 24;

			// // Check CRC
			// boost::crc_32_type crc;
			// crc.process_bytes(&msg, sizeof(msg.header) + msg.header.msg_len);
			// if (crc.checksum() == msg.footer.crc) {
			// 	return 0;
			// } else {
			// 	// PX4_ERR("CRC error");
			// 	return -EAGAIN;
			// }

			reset();

			return 0;
			
		}
		return -EAGAIN;
	}
} // namespace unicore
