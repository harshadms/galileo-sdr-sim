#include "../include/galileo-sdr.h"


typedef union
{
	double d_data;
	unsigned int i_data[2];
} DOUBLE_INT_UNION;

int roundi(double data)
{
	if (data >= 0)
		return (int)(data + 0.5);
	else
		return (int)(data - 0.5);
}

int roundu(double data)
{
	return (unsigned int)(data + 0.5);
}

double UnscaleDouble(double value, int scale)
{
	DOUBLE_INT_UNION data;

	data.d_data = value;
	data.i_data[1] -= (scale << 20);
	return data.d_data;
}

int UnscaleInt(double value, int scale)
{
	long long int long_value = UnscaleLong(value, scale);
	return (int)long_value;
}

unsigned int UnscaleUint(double value, int scale)
{
	unsigned long long int long_value = UnscaleULong(value, scale);
	return (unsigned int)long_value;
}

long long int UnscaleLong(double value, int scale)
{
	DOUBLE_INT_UNION *data = (DOUBLE_INT_UNION *)(&value);
	long long int fraction;
	int sign;

	sign = (data->i_data[1] & 0x80000000);
	fraction = (long long int)UnscaleULong(value, scale);
	return (sign ? -fraction : fraction);
}

unsigned long long int UnscaleULong(double value, int scale)
{
	DOUBLE_INT_UNION *data = (DOUBLE_INT_UNION *)(&value);
	unsigned long long int fraction;
	int exp, shift;

	exp = ((data->i_data[1] & 0x7ff00000) >> 20);
	data->i_data[1] &= 0xfffff;	// clear sign/exp
	fraction = *(long long *)(&value);
	if (exp == 0 && fraction == 0)
		return 0;
	fraction |= 0x10000000000000LL;	// add 1. part
	shift = 1074 - exp + scale;
	fraction += (1LL << shift);	// round
	fraction >>= (shift + 1);
	return fraction;
}

const unsigned char ConvEncodeTable[256] = {
0x00, 0x33, 0xde, 0xed, 0x7b, 0x48, 0xa5, 0x96, 0xff, 0xcc, 0x21, 0x12, 0x84, 0xb7, 0x5a, 0x69, 
0xcc, 0xff, 0x12, 0x21, 0xb7, 0x84, 0x69, 0x5a, 0x33, 0x00, 0xed, 0xde, 0x48, 0x7b, 0x96, 0xa5, 
0x21, 0x12, 0xff, 0xcc, 0x5a, 0x69, 0x84, 0xb7, 0xde, 0xed, 0x00, 0x33, 0xa5, 0x96, 0x7b, 0x48, 
0xed, 0xde, 0x33, 0x00, 0x96, 0xa5, 0x48, 0x7b, 0x12, 0x21, 0xcc, 0xff, 0x69, 0x5a, 0xb7, 0x84, 
0xb7, 0x84, 0x69, 0x5a, 0xcc, 0xff, 0x12, 0x21, 0x48, 0x7b, 0x96, 0xa5, 0x33, 0x00, 0xed, 0xde, 
0x7b, 0x48, 0xa5, 0x96, 0x00, 0x33, 0xde, 0xed, 0x84, 0xb7, 0x5a, 0x69, 0xff, 0xcc, 0x21, 0x12, 
0x96, 0xa5, 0x48, 0x7b, 0xed, 0xde, 0x33, 0x00, 0x69, 0x5a, 0xb7, 0x84, 0x12, 0x21, 0xcc, 0xff, 
0x5a, 0x69, 0x84, 0xb7, 0x21, 0x12, 0xff, 0xcc, 0xa5, 0x96, 0x7b, 0x48, 0xde, 0xed, 0x00, 0x33, 
0xcc, 0xff, 0x12, 0x21, 0xb7, 0x84, 0x69, 0x5a, 0x33, 0x00, 0xed, 0xde, 0x48, 0x7b, 0x96, 0xa5, 
0x00, 0x33, 0xde, 0xed, 0x7b, 0x48, 0xa5, 0x96, 0xff, 0xcc, 0x21, 0x12, 0x84, 0xb7, 0x5a, 0x69, 
0xed, 0xde, 0x33, 0x00, 0x96, 0xa5, 0x48, 0x7b, 0x12, 0x21, 0xcc, 0xff, 0x69, 0x5a, 0xb7, 0x84, 
0x21, 0x12, 0xff, 0xcc, 0x5a, 0x69, 0x84, 0xb7, 0xde, 0xed, 0x00, 0x33, 0xa5, 0x96, 0x7b, 0x48, 
0x7b, 0x48, 0xa5, 0x96, 0x00, 0x33, 0xde, 0xed, 0x84, 0xb7, 0x5a, 0x69, 0xff, 0xcc, 0x21, 0x12, 
0xb7, 0x84, 0x69, 0x5a, 0xcc, 0xff, 0x12, 0x21, 0x48, 0x7b, 0x96, 0xa5, 0x33, 0x00, 0xed, 0xde, 
0x5a, 0x69, 0x84, 0xb7, 0x21, 0x12, 0xff, 0xcc, 0xa5, 0x96, 0x7b, 0x48, 0xde, 0xed, 0x00, 0x33, 
0x96, 0xa5, 0x48, 0x7b, 0xed, 0xde, 0x33, 0x00, 0x69, 0x5a, 0xb7, 0x84, 0x12, 0x21, 0xcc, 0xff, 
};


// do 2bit convolution encode at the same time
// EncodeBits bit 7~2 is current state
// EncodeBits bit1 and bit0 are new incomming bits with bit1 comes first
// return value are generated encoded bit with following order
// G21 G11 G20 G10 G11 G21 G10 G20 (G11 and G21 are encoded bits using bit1, G10 and G20 are encoded bits using bit0)
// so 4MSB of return value are 4 encoded bits with G2/G1 order
// 4LSB of return value are 4 encoded bits with G1/G2 order
unsigned char ConvolutionEncode(unsigned char EncodeBits)
{
	return ConvEncodeTable[EncodeBits];
}

// Length is number of bits in BitStream to do CRC24Q encode
// input BitStream filled with MSB first and start from index 0
// if encoded bits does not fit all bits in BitStream (Length is not multiple of 32)
// 0s need to be filled first (at MSBs of BitStream[0]) to make sure last encoded bits
// in bit0 of last index of BitStream array (this is because encoding 0 into all zero
// state CRC24Q encode does not change encoder status)
unsigned int Crc24qEncode(unsigned int *BitStream, int Length)
{
	int i, ByteNum;
	unsigned int Data, crc_result = 0;

	ByteNum = ((Length + 31) / 32) * 4;
	Data = BitStream[0];
	for (i = 0; i < ByteNum; i ++)
	{
		crc_result = (crc_result << 8) ^ Crc24q[(Data >> 24) ^ (unsigned char)(crc_result >> 16)];
		Data <<= 8;
		if ((i & 3) == 3)	// move to next bit
			Data = BitStream[(i >> 2) + 1];
	}

	return crc_result & 0xffffff;
}

void encode_int_to_bits(int *page, int *offset, long value, int num_bits) {
    // Make sure the value fits within the specified number of bits
    value &= (1 << num_bits) - 1;

    // for (int i = 0; i < num_bits; ++i) {
    //     page[(num_bits - i - 1)] = ((((value >> i) & 1) == 1) ? 1:0);
    // }
    int index = *offset;

    for (int j = num_bits-1; j >= 0; j--)
    {
        page[index] = BIT_ISSET(value, j);
        index++;
    }
    *offset = index;
}

void encode_double_to_bits(int *page, int *offset, double value, int num_bits) {
	//cout << "Value: " << fixed << value << endl;
    int32_t quantized_value = (int32_t)(value); // Truncate the double to an integer
	//cout << "Value: " << quantized_value << endl;

    int index = *offset;

    for (int j = num_bits - 1; j >= 0; j--) {
        page[index] = BIT_ISSET(quantized_value, j);
        index++;
    }

    *offset = index;
}

/* Extract data in uint8_t format in bits*/
uint32_t data_extract_uint32_t(const vector<uint8_t> str_in, int offset_in, int len)
{
	uint8_t str_buff_out[4] = {0};
	uint8_t str_buff[4] = {0};
	uint8_t mask1 = 0, mask2 = 0, mask3 = 0;
	int i, j;
	int offset_index_in;
	int offset_index_out;
	int nbit_offset_in, r_nbit_offset_in;
	int nbit_offset_out, r_nbit_offset_out;
	int nbit_len;
	int len_byte;

	uint32_t out = 0;
	int offset_out = 32 - len;

	if (len > 32)
	{
		return 0;
	}

	offset_index_in = floor(offset_in / 8);
	nbit_offset_in = offset_in % 8;
	r_nbit_offset_in = 8 - nbit_offset_in;

	offset_index_out = floor(offset_out / 8);
	nbit_offset_out = offset_out % 8;
	r_nbit_offset_out = 8 - nbit_offset_out;
	nbit_len = len % 8;
	len_byte = floor(len / 8);

	for (j = 0; j < r_nbit_offset_in; j++)
		mask1 |= (uint8_t)(1 << j);

	for (j = 7; j > 7 - nbit_offset_in; j--)
		mask2 |= (uint8_t)(1 << j);

	for (i = 0; i < len_byte; i++)
	{
		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i] & mask1) << nbit_offset_in;

		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i + 1] & mask2) >> (r_nbit_offset_in);
	}

	if (r_nbit_offset_in < nbit_len)
	{
		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i] & mask1) << nbit_offset_in;

		for (j = 7; j > 7 - (nbit_len - r_nbit_offset_in); j--)
			mask3 |= (uint8_t)(1 << j);

		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i + 1] & mask3) >> r_nbit_offset_in;
	}
	else
	{
		for (j = r_nbit_offset_in - 1; j > r_nbit_offset_in - 1 - nbit_len; j--)
			mask3 |= (uint8_t)(1 << j);

		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i] & mask3) << nbit_offset_in;
	}

	if (nbit_len != 0)
		len_byte++;

	if (nbit_offset_out == 0)
	{
		for (i = 0; i < len_byte; i++)
			str_buff_out[offset_index_out + i] = str_buff[i];
	}
	else
	{
		mask1 = 0;
		mask2 = 0;

		for (j = 0; j < nbit_offset_out; j++)
			mask1 |= (uint8_t)(1 << j);

		for (j = 7; j > 7 - r_nbit_offset_out; j--)
			mask2 |= (uint8_t)(1 << j);

		for (i = 0; i < len_byte; i++)
		{
			str_buff_out[offset_index_out + i] |= (uint8_t)(str_buff[i] & mask2) >> nbit_offset_out;
			str_buff_out[offset_index_out + i + 1] = (uint8_t)(str_buff[i] & mask1) << r_nbit_offset_out;
		}
	}

	out = str_buff_out[3] | str_buff_out[2] << 8 | str_buff_out[1] << 16 | str_buff_out[0] << 24;

	return out;
}

/* Extract data in uint8_t format in bits*/
void data_extract(vector<uint8_t> *str_out, const vector<uint8_t> str_in, int offset_out, int offset_in, int len)
{
	uint8_t str_buff[MAX_CHAR] = {0};
	uint8_t mask1 = 0, mask2 = 0, mask3 = 0;
	int i, j;
	int offset_index_in;
	int offset_index_out;
	int nbit_offset_in, r_nbit_offset_in;
	int nbit_offset_out, r_nbit_offset_out;
	int nbit_len;
	int len_byte;

	offset_index_in = floor(offset_in / 8);
	nbit_offset_in = offset_in % 8;
	r_nbit_offset_in = 8 - nbit_offset_in;

	offset_index_out = floor(offset_out / 8);
	nbit_offset_out = offset_out % 8;
	r_nbit_offset_out = 8 - nbit_offset_out;
	nbit_len = len % 8;
	len_byte = floor(len / 8);

	for (j = 0; j < r_nbit_offset_in; j++)
		mask1 |= (uint8_t)(1 << j);

	for (j = 7; j > 7 - nbit_offset_in; j--)
		mask2 |= (uint8_t)(1 << j);

	for (i = 0; i < len_byte; i++)
	{
		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i] & mask1) << nbit_offset_in;

		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i + 1] & mask2) >> (r_nbit_offset_in);
	}

	if (r_nbit_offset_in < nbit_len)
	{
		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i] & mask1) << nbit_offset_in;

		for (j = 7; j > 7 - (nbit_len - r_nbit_offset_in); j--)
			mask3 |= (uint8_t)(1 << j);
		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i + 1] & mask3) >> r_nbit_offset_in;
	}
	else
	{
		for (j = r_nbit_offset_in - 1; j > r_nbit_offset_in - 1 - nbit_len; j--)
			mask3 |= (uint8_t)(1 << j);
		str_buff[i] |= (uint8_t)(str_in[offset_index_in + i] & mask3) << nbit_offset_in;
	}

	if (nbit_len != 0)
		len_byte++;

	if (nbit_offset_out == 0)
	{

		for (i = 0; i < len_byte; i++)
			str_out->at(offset_index_out + i) = str_buff[i];
	}
	else
	{
		mask1 = 0;
		mask2 = 0;

		for (j = 0; j < nbit_offset_out; j++)
			mask1 |= (uint8_t)(1 << j);

		for (j = 7; j > 7 - r_nbit_offset_out; j--)
			mask2 |= (uint8_t)(1 << j);

		for (i = 0; i < len_byte; i++)
		{
			str_out->at(offset_index_out + i) |= (uint8_t)(str_buff[i] & mask2) >> nbit_offset_out;
			str_out->at(offset_index_out + i + 1) = (uint8_t)(str_buff[i] & mask1) << r_nbit_offset_out;
		}
	}
}

/* Convert Char read in the file in uint8_t */
void chartouint8_t(vector<uint8_t> *out, char *in, int len)
{
	unsigned char tmp[2];
	int i, j;

	for (i = 0; i < len; i++)
	{
		for (j = 0; j < 2; j++)
		{
			switch (in[i * 2 + j])
			{
			case 'A':
				tmp[j] = 10;
				break;
			case 'B':
				tmp[j] = 11;
				break;
			case 'C':
				tmp[j] = 12;
				break;
			case 'D':
				tmp[j] = 13;
				break;
			case 'E':
				tmp[j] = 14;
				break;
			case 'F':
				tmp[j] = 15;
				break;
			default:
				tmp[j] = in[i * 2 + j] - 48;
				break;
			}
		}
		out->push_back(16 * tmp[0] + tmp[1]);
	}
}

/* Extract signed integers, accounting for sign bit (MSB) */
int data_extract_signedInt(int len, int data)
{
	int32_t mask = (1 << len) - 1;
	int32_t sign_bit = 1 << (len - 1);

	if (data & sign_bit)
	{
		int32_t inverted = (~data) & mask;
		return -(inverted + 1);
	}
	else
	{
		return data;
	}
}