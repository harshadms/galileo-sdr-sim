#include "../include/galileo-sdr.h"

/* Apply FEC, Interleaving and add Sync pattern (gal-osnma-sdr)*/
void generateFrame(int *half_page, int *frame)
{   

    int fec_tmp_sbf[240] = {0};
    // FEC Encoding
    cnv_encd(N_BIT_PAGE, half_page, fec_tmp_sbf);

    int int_tmp_sbf[240] = {0};
    // Interleaving
    for (int r = 0; r < 8; r++)
    {
        for (int c = 0; c < 30; c++)
        {
            int_tmp_sbf[r * 30 + c] = fec_tmp_sbf[c * 8 + r];
        }
    }

    // Add sync
    for (int i = 0; i < 10; i++)
        frame[i] = sync_pattern[i];
    for (int i = 0; i < 240; i++)
        frame[i + 10] = int_tmp_sbf[i];
}

int generateINavMsg(galtime_t g, channel_t *chan, ephem_t *eph, ionoutc_t *ionoutc)
{
    // Start from next odd second (Even part transmit first)
    int *page = (int *)malloc(PAGE_SIZE * sizeof(int));
    int *odd_page = (int *)malloc(120 * sizeof(int));
    int *even_page = (int *)malloc(120 * sizeof(int));
	
	for (int i =0; i < 120; i++)
		{odd_page[i] = 0;
		even_page[i] = 0;}

    int word_type_index = ((int)g.sec % 30)/2;
    int word_type = WordAllocationE1[word_type_index];

    generate_page(g, eph, ionoutc, even_page, odd_page, word_type);
    generateFrame(even_page, &page[0]);
    generateFrame(odd_page, &page[250]);

	// cout << endl;
	// for(int i = 0; i < PAGE_SIZE; i++)
	// 	cout << page[i];

	// exit(1);
	
    chan->page = page;
    return 1;
}

// Convolutional encoding (FEC) ICD 4.1.4. FEC Coding and Interleaving Parameters
void cnv_encd(long input_len, int *in_array, int *out_array)
{
	const int K_CNV_ENCD = 7;

	int m = K_CNV_ENCD - 1;
	int shift_reg[K_CNV_ENCD] = {0}; /* the encoder shift register */

	// Allocate space for the zero-padded input data array
	int *unencoded_data = (int *)malloc((input_len + m) * sizeof(int));

	if (unencoded_data == NULL)
	{
		printf("\ncnv_encd.c: Can't allocate enough memory for unencoded data!  Aborting...");
		exit(1);
	}

	// Read in the data and store it in the array
	for (int t = 0; t < input_len; t++)
		unencoded_data[t] = in_array[t];

	// Zero-pad the end of the data
	for (int t = 0; t < m; t++)
	{
		unencoded_data[input_len + t] = 0;
	}

	// Initialize the shift register
	for (int j = 0; j < K_CNV_ENCD; j++)
	{
		shift_reg[j] = 0;
	}

	/* To try to speed things up a little, the shift register will be operated
	   as a circular buffer, so it needs at least a head pointer. It doesn't
	   need a tail pointer, though, since we won't be taking anything out of
	   it--we'll just be overwriting the oldest entry with the new data. */

	int sr_head = 0; /* index to the first elt in the Shift register */

	// Initialize the channel symbol output index

	/*---- Now start the encoding process -----*/

	// Compute the upper and lower mod-two adder outputs, one bit at a time
	for (int t = 0; t < input_len; t++)
	{
		shift_reg[sr_head] = unencoded_data[t];
		int G1, nG2; /* the upper and lower xor gate outputs */
		G1 = 0;
		nG2 = 0; /// @note ICD fig.13
		for (int j = 0; j < K_CNV_ENCD; j++)
		{
			int iSr = (j + sr_head) % K_CNV_ENCD;
			G1 ^= shift_reg[iSr] & ccGenMatrix[0][j];
			nG2 ^= shift_reg[iSr] & ccGenMatrix[1][j];
		}

		// Write the upper and lower xor gate outputs as channel symbols
		out_array[2 * t] = G1;
		out_array[2 * t + 1] = 1 - nG2;

		sr_head -= 1;	 // Equivalent to shifting everything right one place
		if (sr_head < 0) // But make sure we adjust pointer moduloK_CNV_ENCDK
			sr_head = m;
	}

	// Free the dynamically allocated array
	free(unencoded_data);
}

// Simple function to shift and insert elements at the specified location (Hack)
void shift_and_insert(int *arr, int index, uint32_t newElement1, uint32_t newElement2, int arr_size) {
    if (index < 0 || index >= arr_size) {
        throw std::runtime_error("Invalid index");
    }

    for (int i = arr_size-1; i >= index + 2; i--) {
        arr[i] = arr[i - 2];
    }

    arr[index] = newElement1;
    arr[index + 1] = newElement2;
}

unsigned int Crc24qEncode(int *BitStream, int Length)
{
    int i, BitNum;
    unsigned int Data = 0, crc_result = 0;

    // Calculate the number of bits needed to store the CRC result
    BitNum = 24;
	int j;

	for (j = 0; j < Length; j++) {
		unsigned char dataByte = (char)BitStream[j];
		for (i = 0; i < 8; i++) {
			crc_result = (crc_result << 8) ^ Crc24q[(dataByte >> (7 - i)) ^ (unsigned char)(crc_result >> 16)];

		}
	}
	return (crc_result & 0xffffff);
}

// Generate page from eph - 4.3.5. I/NAV Word Types
void generate_page(galtime_t g, ephem_t *eph, ionoutc_t *ionoutc, int *even_page, int *odd_page, int word_type)
{
	int offset = 0;
	int page[240] = {};

	for (int i =0; i < 240; i++)
		page[i] = 0;
	//print_eph2(eph, eph->svid);

	int TOW;
	// first determine the current TOW and subframe number
	// g.week += g.sec*1000 / 604800000;
	// g.sec = (int)(g.sec*1000) % 604800000;
	// TOW = g.sec / 2000 * 2 - ((0 == 0) ? 1 : 0);	// Param == 0 for E1 I/NAV
	
	// if (TOW < 0)
	// 	TOW += 604800;

	TOW = (int)g.sec;

	long IntValue;
	unsigned int UintValue;

	//cout << g.week << " : " << TOW << endl;

	switch (word_type)
	{
	case 0:
		// Word type - Even page
		encode_int_to_bits(page, &offset, 0, 8); // Word type 0
		encode_int_to_bits(page, &offset, 2, 2); // Word type 0
		encode_int_to_bits(page, &offset, 0, 88);
		// WN
		encode_int_to_bits(page, &offset, g.week - 1024, 12);
		// TOW
		encode_int_to_bits(page, &offset, TOW, 20);
		break;
	case 1:
		// Word type - Even page
		encode_int_to_bits(page, &offset, 1, 8); // Word type 1
		// IODNav
		encode_int_to_bits(page, &offset, eph->iode, 10); 
		// TOE
		encode_int_to_bits(page, &offset, (int)eph->toe.sec /(60), 14);
		//cout <<  "TOE: " << eph->toe.sec << endl;
		// M0
		IntValue = UnscaleInt(eph->m0 / PI, -31);
		//cout <<  "M0: " << IntValue << " : " << eph->m0 << endl;
		encode_double_to_bits(page, &offset, IntValue, 32);
		// ECC
		UintValue = UnscaleUint(eph->ecc, -33);
		//cout <<  "ECC: " << UintValue << " : " << eph->ecc << endl;
		encode_double_to_bits(page, &offset, UintValue, 32);
		// Sqrt
		IntValue = UnscaleInt(eph->sqrta, -19);
		//cout <<  "Sqrta: " << IntValue << " : " << eph->sqrta << endl;
		encode_double_to_bits(page, &offset, IntValue, 32);
		// 2 reserved bits
		encode_int_to_bits(page, &offset, 0, 2);
		// cout <<  "PRN: " << eph->svid << endl;
		// for (int i=0;i<128;i++)
		// 	cout << page[i];

		// exit(1);
		break;
	
	case 2:
		// Word type - Even page
		encode_int_to_bits(page, &offset, 2, 8); // Word type 2
		// IODNav
		encode_int_to_bits(page, &offset, eph->iode, 10);
		// Omega0
		IntValue = UnscaleInt(eph->omg0 / PI, -31);
		encode_double_to_bits(page, &offset, IntValue, 32);
		// Inc
		IntValue = UnscaleInt(eph->inc0 / PI, -31);
		encode_double_to_bits(page, &offset, IntValue, 32);
		// AoP
		IntValue = UnscaleInt(eph->aop / PI, -31);
		encode_double_to_bits(page, &offset, IntValue, 32);
		// iDot
		IntValue = UnscaleInt(eph->idot / PI, -43);
		encode_double_to_bits(page, &offset, IntValue, 14);
		// 2 reserved bits
		encode_int_to_bits(page, &offset, 0, 2);
		break;
	
	case 3:
		// Word type - Even page
		encode_int_to_bits(page, &offset, 3, 8); // Word type 3
		// IODNav
		encode_int_to_bits(page, &offset, eph->iode, 10); 
		// Omega dot
		IntValue = UnscaleInt(eph->omgdot / PI, -43);
		encode_int_to_bits(page, &offset, IntValue, 24);
		// Deltan
		IntValue = UnscaleInt(eph->deltan / PI, -43);
		encode_int_to_bits(page, &offset, IntValue, 16);
		// Cuc
		IntValue = UnscaleInt(eph->cuc, -29);
		encode_int_to_bits(page, &offset, IntValue, 16);
		// Cus
		IntValue = UnscaleInt(eph->cus, -29);
		encode_int_to_bits(page, &offset, IntValue, 16);
		// Crc
		IntValue = UnscaleInt(eph->crc, -5);
		encode_int_to_bits(page, &offset, IntValue, 16);
		// Crs
		IntValue = UnscaleInt(eph->crs, -5);
		encode_int_to_bits(page, &offset, IntValue, 16);
		// User range accuracy index - 32767 for URA 15 SISA
		IntValue = UnscaleInt(3.12, 0);
		encode_int_to_bits(page, &offset, 32767, 8);
		break;
	
	case 4:
		// Word type - Even page
		encode_int_to_bits(page, &offset, 4, 8); // Word type 4
		// IODNav
		encode_int_to_bits(page, &offset, eph->iode, 10); 
		// SVID
		encode_int_to_bits(page, &offset, eph->svid, 6);
		// Cic
		IntValue = UnscaleInt(eph->cic, -29);
		encode_int_to_bits(page, &offset, IntValue, 16);
		// Cis
		IntValue = UnscaleInt(eph->cis, -29);
		encode_int_to_bits(page, &offset, IntValue, 16);

		// if (eph->svid == 1)
		// {cout << "cis: " << IntValue << " : " << eph->svid << endl;
		// 	exit(1);}
		// TOC
		UintValue = eph->toc.sec / 60;
		//cout << "cis: " << UintValue << " : " << eph->svid << endl;
		encode_int_to_bits(page, &offset, UintValue, 14);
		// AF0
		IntValue = UnscaleInt(eph->af0, -34);
		encode_int_to_bits(page, &offset, IntValue, 31);
		// AF1
		IntValue = UnscaleInt(eph->af1, -46);
		encode_int_to_bits(page, &offset, IntValue, 21);
		// AF2
		IntValue = UnscaleInt(eph->af2, -59);
		encode_int_to_bits(page, &offset, IntValue, 6);
		// Spare bits
		encode_int_to_bits(page, &offset, 0, 2);
		break;
	
	case 5:
		// Word type - Even page
		encode_int_to_bits(page, &offset, 5, 8); // Word type 4
		// iono ai0
		UintValue = UnscaleUint(ionoutc->ai0, -2);
		encode_double_to_bits(page, &offset, UintValue, 11);
		// iono ai1
		IntValue = UnscaleInt(ionoutc->ai1, -8);
		encode_double_to_bits(page, &offset, IntValue, 11);
		// iono ai2
		IntValue = UnscaleInt(ionoutc->ai2, -15);
		encode_double_to_bits(page, &offset, IntValue, 14);
		// Regional flags 5
		encode_int_to_bits(page, &offset, 31, 5);
		// BGD E1E5a
		IntValue = UnscaleInt(eph->bgde5a, -32);
		encode_int_to_bits(page, &offset, IntValue, 10);
		// BGD E1E5b
		IntValue = UnscaleInt(eph->bgde5b, -32);
		encode_int_to_bits(page, &offset, IntValue, 10);
		// sV Health
		encode_int_to_bits(page, &offset, eph->svhlth >> 7, 2);	// E5b HS
		encode_int_to_bits(page, &offset, eph->svhlth >> 1, 2);	// E1B HS
		encode_int_to_bits(page, &offset, eph->svhlth >> 5, 1);	// E5b DVS
		encode_int_to_bits(page, &offset, eph->svhlth, 1);		// E1B DVS
		// WN
		encode_int_to_bits(page, &offset, g.week-1024, 12);
		// TOW
		encode_int_to_bits(page, &offset, TOW, 20);
		// Spare
		encode_int_to_bits(page, &offset, 0, 3);
		break;
	
	case 6:
		// Word type - Even page
		encode_int_to_bits(page, &offset, 6, 8); // Word type 6
		// GST A0
		IntValue = UnscaleInt(ionoutc->A0, -30);
		encode_double_to_bits(page, &offset, IntValue, 32);
		// GST A1
		IntValue = UnscaleInt(ionoutc->A1, -50);
		encode_double_to_bits(page, &offset, IntValue, 24);
		// Leap seconds
		encode_int_to_bits(page, &offset, ionoutc->dtls, 8);	
		// GST ToT
		encode_int_to_bits(page, &offset, ionoutc->tot/3600.0, 8);
		encode_int_to_bits(page, &offset, ionoutc->wnt, 8);
		encode_int_to_bits(page, &offset, ionoutc->wnlsf, 8);
		encode_int_to_bits(page, &offset, ionoutc->dn, 3);
		encode_int_to_bits(page, &offset, ionoutc->dtlsf, 8);
		// TOW
		encode_int_to_bits(page, &offset, TOW, 20);
		// Spare
		encode_int_to_bits(page, &offset, 0, 3);


		break;
	
	default:
		// Word type 63 - Even page (dummy frame)
		encode_int_to_bits(page, &offset, 63, 8);
		encode_int_to_bits(page, &offset, 0, 122);
		// encode_int_to_bits(page, &offset, 0x686D73, 24);
		// encode_int_to_bits(page, &offset, 0, 48);
		break;
	}

	// 40 reserved bits
	encode_int_to_bits(page, &offset, 0, 40);
	// SAR
	encode_int_to_bits(page, &offset, 2796202, 22);
	// Spare
	encode_int_to_bits(page, &offset, 0, 2);

	// Add Odd page indicator at the 114th bit
	shift_and_insert(page, 114, 1, 0, 240);

	// CRC
	unsigned int crcq = Crc24qEncode(page, 196);
	encode_int_to_bits(page, &offset, crcq, 24);

	int ssp[3] = {4, 43, 47};

	encode_int_to_bits(page, &offset, ssp[word_type % 3], 8);

	memcpy(even_page, page, 114 * sizeof(int));
	memcpy(odd_page, page + 114, 114 * sizeof(int));

	// if (eph->svid == 1)
	// 	cout << "TOW: " << TOW << " : " << g.week << endl;
	//encode_int_to_bits()
}