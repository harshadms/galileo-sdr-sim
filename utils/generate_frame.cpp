#include <iostream>
#include <cmath>
#include <cstring> 
#include <vector>

using namespace std;

#define BIT_ISSET(var, n) !!((long)(var) & (1 << (n)))
#define COMPOSE_BITS(data, start, width) (((data) & ((1UL << width) - 1)) << start)


void encode_value_to_bits(uint32_t *page, int *offset, long value, int num_bits) {
    // Make sure the value fits within the specified number of bits
    value &= (1 << num_bits) - 1;
	

    // for (int i = 0; i < num_bits; ++i) {
    //     page[(num_bits - i - 1)] = ((((value >> i) & 1) == 1) ? 1:0);
    // }
    int index = *offset;

    for (int j = num_bits-1; j >= 0; j--)
    {
        page[index] = BIT_ISSET(value, j);
		cout << hex << uppercase << page[index] << endl;
        index++;        
    }
    *offset = index;
}

void encode_double_to_bits(uint32_t *page, int *offset, double value, int num_bits) {    
    int32_t quantized_value = static_cast<int32_t>(value); // Truncate the double to an integer
    cout << hex << uppercase << quantized_value << endl;
    int index = *offset;
    
    for (int j = num_bits - 1; j >= 0; j--) {
        page[index] = BIT_ISSET(quantized_value, j);
        index++;
        
    }
    *offset = index;
}

// Function to encode a double value into bits with scaling factor and offset
uint64_t encode_scaled_double(double value, double scaleFactor, int offset) {
    // Scale the value and convert to bits
    uint64_t bits = static_cast<uint64_t>(value / scaleFactor);
    return bits; // << offset;
}

// Function to encode a long value into bits with scaling factor and offset
uint64_t encode_scaled_long(long value, long scaleFactor, int offset) {
    // Scale the value and convert to bits
    uint64_t bits = static_cast<uint64_t>(value / scaleFactor);
    return bits << offset;
}

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

unsigned long long int encode_ulong(double value, int scale)
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


double encode_double(double value, int scale)
{
	DOUBLE_INT_UNION data;

	data.d_data = value;
	data.i_data[1] -= (scale << 20);
	return data.d_data;
}


long long int encode_long(double value, int scale)
{
	DOUBLE_INT_UNION *data = (DOUBLE_INT_UNION *)(&value);
	long long int fraction;
	int sign;

	sign = (data->i_data[1] & 0x80000000);
	fraction = (long long int)encode_ulong(value, scale);
	return (sign ? -fraction : fraction);
}

int encode_int(double value, int scale)
{
	long long int long_value = encode_long(value, scale);
	return (int)long_value;
}

unsigned int encode_uint(double value, int scale)
{
	unsigned long long int long_value = encode_ulong(value, scale);
	return (unsigned int)long_value;
}



// put bit in Data from MSB ot LSB into BitStream, bit order from bit(BitNumber-1) to bit(0) of Data
int AssignBits(unsigned int Data, int BitNumber, int BitStream[])
{
	int i;

	Data <<= (32 - BitNumber);
	for (i = 0; i < BitNumber; i++)
	{
		BitStream[i] = (Data & 0x80000000) ? 1 : 0;
		Data <<= 1;
	}

	return BitNumber;
}


#define POW2_M31 4.656612873077393e-10
#define POW2_M33 1.164153218269348e-10
#define POW2_M43 1.136868377216160e-13

int main() {
    // double A = 2.960013e+07;
    // double af0 = -3.754718e-04;
    // double af1 = -4.192202e-12;
    // double af2 = 0.000000e+00;
    // double aop = 9.038998e-01;
    // double bgde5a = 3.026798e-09;
    // double bgde5b = 3.492460e-09;
    // double cic = -2.421439e-08;
    // double cis = 3.725290e-09;
    // double crc = 2.540000e+02;
    // double crs = -1.687500e+00;
    // double cuc = -1.210719e-07;
    // double cus = 3.883615e-06;
    // double deltan = 3.881947e-09;
    // double ecc = 3.463253e-04;
    // double idot = -8.143196e-11;
    // double inc0 = 9.536874e-01;
    // int iodnav = 52;
    // double m0 = 7.766745e-01;
    // double n = 1.239773e-04;
    // double omg0 = 1.305452e+00;
    // double omgdot = -5.833814e-09;
    // double omgkdot = -7.292699e-05;
    // double sq1e2 = 9.999999e-01;
    // double sqrta = 5.440600e+03;
    // double toc = 31200;
    // double toe = 31200;

    double A       = 2.960030e+07;
    double af0     = -1.000991e-03;
    double af1     = -8.085976e-12;
    double af2     = 0.000000e+00;
    double aop     = -6.302550e-01;
    double bgde5a  = 4.656613e-10;
    double bgde5b  = 4.656613e-10;
    double cic     = 9.126961e-08;
    double cis     = -5.774200e-08;
    double crc     = 2.680938e+02;
    double crs     = 1.671562e+02;
    double cuc     = 7.903203e-06;
    double  cus    = 3.773719e-06;
    double deltan  = 2.967624e-09;
    double ecc     = 2.578078e-04;
    double idot    = 4.385897e-10;
    double inc0    = 9.805393e-01;
    int iodnav  = 126;
    double m0      = 3.505807e-01;
    double n       = 1.239753e-04;
    double omg0    = -2.759085e+00;
    double omgdot  = -5.827028e-09;
    double omgkdot = -7.292698e-05;
    double sq1e2   = 1.000000e+00;
    double sqrta   = 5.440616e+03;
    double toc     = 459600;
    double toe     = 459600;

    uint32_t value = 52; // Example value
    int num_bits = 10;   // Number of bits
    
    int32_t word_type = 2;

    uint32_t page[130] = {};

    int offset = 0;
    // Even or odd

    int IntValue;
    unsigned int UintValue;
    unsigned int EphData[19];

    const int WordAllocationE1[15] = {
	    2, 4, 6, 7, 8, 17, 19, 16, 0, 0, 1, 3, 5, 0, 16,
    };

    // Word 1
	EphData[0] = 0x04000000 | COMPOSE_BITS(iodnav, 16, 10) | COMPOSE_BITS(((int)toe / 60), 2, 14);
	IntValue = encode_int(m0 / M_PI, -31);
	EphData[0] |= COMPOSE_BITS(IntValue >> 30, 0, 2);
	EphData[1] = COMPOSE_BITS(IntValue, 2, 30);
	UintValue = encode_uint(ecc, -33);
	EphData[1] |= COMPOSE_BITS(UintValue >> 30, 0, 2);
	EphData[2] = COMPOSE_BITS(UintValue, 2, 30);
	UintValue = encode_uint(sqrta, -19);
	EphData[2] |= COMPOSE_BITS(UintValue >> 30, 0, 2);
	EphData[3] = COMPOSE_BITS(UintValue, 2, 30);

	// Word 2
	EphData[4] = 0x08000000 | COMPOSE_BITS(iodnav, 16, 10);
	IntValue = encode_int(omg0 / M_PI, -31);
	EphData[4] |= COMPOSE_BITS(IntValue >> 16, 0, 16);
	EphData[5] = COMPOSE_BITS(IntValue, 16, 16);
	IntValue = encode_int(inc0 / M_PI, -31);
	EphData[5] |= COMPOSE_BITS(IntValue >> 16, 0, 16);
	EphData[6] = COMPOSE_BITS(IntValue, 16, 16);
	IntValue = encode_int(aop / M_PI, -31);
	EphData[6] |= COMPOSE_BITS(IntValue >> 16, 0, 16);
	EphData[7] = COMPOSE_BITS(IntValue, 16, 16);
	IntValue = encode_int(idot / M_PI, -43);
	EphData[7] |= COMPOSE_BITS(IntValue >> 16, 2, 14);
	// Word 3
	EphData[8] = 0x0c000000 | COMPOSE_BITS(iodnav, 16, 10);
	IntValue = encode_int(omgdot / M_PI, -43);
	EphData[8] |= COMPOSE_BITS(IntValue >> 8, 0, 16);	
	EphData[9] = COMPOSE_BITS(IntValue, 24, 8);
	IntValue = encode_int(deltan / M_PI, -43);
	EphData[9] |= COMPOSE_BITS(IntValue >> 8, 8, 16);
	IntValue = encode_int(cuc, -29);
	EphData[9] |= COMPOSE_BITS(IntValue >> 8, 0, 8);
	EphData[10] = COMPOSE_BITS(IntValue, 24, 8);
	IntValue = encode_int(cus, -29);
	EphData[10] |= COMPOSE_BITS(IntValue, 8, 16);
	IntValue = encode_int(crc, -5);
	EphData[10] |= COMPOSE_BITS(IntValue >> 8, 0, 8);
	EphData[11] = COMPOSE_BITS(IntValue, 24, 8);
	IntValue = encode_int(crs, -5);
	EphData[11] |= COMPOSE_BITS(IntValue, 8, 16);
	EphData[11] |= COMPOSE_BITS(32767, 0, 8); // User Range Accuracy index 32767 for URA 15

	// Word 4
	EphData[12] = 0x10000000 | COMPOSE_BITS(iodnav, 16, 10) | COMPOSE_BITS(1, 10, 6);
	IntValue = encode_int(cic, -29);
	EphData[12] |= COMPOSE_BITS(IntValue >> 6, 0, 10);
	EphData[13] = COMPOSE_BITS(IntValue, 26, 6);
	IntValue = encode_int(cis, -29);
	EphData[13] |= COMPOSE_BITS(IntValue, 10, 16);
	UintValue = toc / 60;
	EphData[13] |= COMPOSE_BITS(UintValue >> 4, 0, 10);
	EphData[14] = COMPOSE_BITS(UintValue, 28, 4);
	IntValue = encode_int(af0, -34);
	EphData[14] |= COMPOSE_BITS(IntValue >> 3, 0, 28);
	EphData[15] = COMPOSE_BITS(IntValue, 29, 3);
	IntValue = encode_int(af1, -46);
	EphData[15] |= COMPOSE_BITS(IntValue, 8, 21);
	IntValue = encode_int(af2, -59);
	EphData[15] |= COMPOSE_BITS(IntValue, 2, 6);

	// Word 5
	EphData[16] &= 0x03ffffff; EphData[16] = 0x14000000;	// put Type=5 to 6MSB
	IntValue = encode_int(bgde5a, -32);
	EphData[17] = COMPOSE_BITS(IntValue, 7, 10);
	IntValue = encode_int(bgde5b, -32);
	EphData[17] |= COMPOSE_BITS(IntValue >> 3, 0, 7);
	EphData[18] = COMPOSE_BITS(IntValue, 29, 3);
	IntValue = 1;
	EphData[18] |= COMPOSE_BITS(IntValue >> 7, 27, 2);	// E5b HS
	EphData[18] |= COMPOSE_BITS(IntValue >> 1, 25, 2);	// E1B HS
	EphData[18] |= COMPOSE_BITS(IntValue >> 5, 24, 1);	// E5b DVS
	EphData[18] |= COMPOSE_BITS(IntValue, 23, 1);		// E1B DVS

    for (int i=0; i<4;i++)
    {
        cout << std::hex <<  uppercase << EphData[(4 * (WordAllocationE1[0] - 1)) + i];
    }

	uint64_t encodedOmg0 = encode_scaled_double(omg0 / M_PI, POW2_M31, 0);
    // //uint64_t encodedDeltN = encode_scaled_long(deltN, scaleFactorDeltN, 32);

    // //std::cout << "Encoded DeltN: " << std::hex << encodedDeltN << std::dec << std::endl;

    // encode_value_to_bits(page, &offset, 2, 8); // Word type 2
    // // IODNav
    // encode_value_to_bits(page, &offset, iodnav, 10); //omg0/POW2_M31/M_PI, inc0/POW2_M31/M_PI, aop/POW2_M31/M_PI, idot/POW2_M43/M_PI
    // // Omega0
    //encode_double_to_bits(page, &offset, omg0/(POW2_M31 * M_PI), 32);

	// for (int i=0; i < 128; i++)
	// 	cout << page[i];

	// // (POW2_M31 * PI)
    // encode_double_to_bits(page, &offset, inc0/POW2_M31/M_PI, 32);
    // encode_double_to_bits(page, &offset, aop/POW2_M31/M_PI, 32);
    // encode_double_to_bits(page, &offset, idot/POW2_M43/M_PI, 14);
    return 0;

}