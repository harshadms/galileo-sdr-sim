#include <iostream>
#include <cmath>
#include <cstring> 
#include <vector>

using namespace std;

#define BIT_ISSET(var, n) !!((long)(var) & (1 << (n)))


void encode_value_to_bits(uint32_t *page, int *offset, long value, int num_bits) {
    // Make sure the value fits within the specified number of bits
    printf("\n%ld",value);
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

void encode_double_to_bits(uint32_t *page, int *offset, double value, int num_bits) {    
    int32_t quantized_value = static_cast<int32_t>(value); // Truncate the double to an integer
    
    int index = *offset;
    
    for (int j = num_bits - 1; j >= 0; j--) {
        page[index] = BIT_ISSET(quantized_value, j);
        index++;
        
    }
    *offset = index;
}

int* edit2()
{
    int i[5] = {5};

    return &i;
}

// void edit(vector<int> *a)
// {
//     a->at(1) = 1;
//     int i = 0;
//     edit2(&i);
//     a->at(2) = i;
// }


#define POW2_M31 4.656612873077393e-10
#define POW2_M33 1.164153218269348e-10
#define POW2_M43 1.136868377216160e-13

int main() {
    double A = 2.960013e+07;
    double af0 = -3.754718e-04;
    double af1 = -4.192202e-12;
    double af2 = 0.000000e+00;
    double aop = 9.038998e-01;
    double bgde5a = 3.026798e-09;
    double bgde5b = 3.492460e-09;
    double cic = -2.421439e-08;
    double cis = 3.725290e-09;
    double crc = 2.540000e+02;
    double crs = -1.687500e+00;
    double cuc = -1.210719e-07;
    double cus = 3.883615e-06;
    double deltan = 3.881947e-09;
    double ecc = 3.463253e-04;
    double idot = -8.143196e-11;
    double inc0 = 9.536874e-01;
    int iodnav = 52;
    double m0 = 7.766745e-01;
    double n = 1.239773e-04;
    double omg0 = 1.305452e+00;
    double omgdot = -5.833814e-09;
    double omgkdot = -7.292699e-05;
    double sq1e2 = 9.999999e-01;
    double sqrta = 5.440600e+03;
    double toc = 31200;
    double toe = 31200;

    uint32_t value = 52; // Example value
    int num_bits = 10;   // Number of bits
    
    int32_t word_type = 2;

    uint32_t page[128] = {};

    int offset = 0;
    // Even or odd
    
    encode_value_to_bits(page, &offset, 2, 8); // Word type 2
    // IODNav
    encode_value_to_bits(page, &offset, iodnav, 10); //omg0/POW2_M31/PI, inc0/POW2_M31/PI, aop/POW2_M31/PI, idot/POW2_M43/PI
    // Omega0
    encode_double_to_bits(page, &offset, omg0/POW2_M31/M_PI, 32);

    encode_double_to_bits(page, &offset, inc0/POW2_M31/M_PI, 32);
    encode_double_to_bits(page, &offset, aop/POW2_M31/M_PI, 32);
    encode_double_to_bits(page, &offset, idot/POW2_M43/M_PI, 16);

    return 0;

}