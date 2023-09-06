#include <iostream>
#include <cmath>
#include "galileo-sdr.h"
#include <vector>
#include <cstring>
#include <fstream>

using namespace std;

void print_eph(ephem_t *eph, int prn)
{
	printf("\nEphemeris:\n");

	printf("\tPRN %d: A\t%e\n", prn, eph->A);
	printf("\tPRN %d: af0\t%e\n", prn, eph->af0);
	printf("\tPRN %d: af1\t%e\n", prn, eph->af1);
	printf("\tPRN %d: af2\t%e\n", prn, eph->af2);
	printf("\tPRN %d: aop\t%e\n", prn, eph->aop);
	printf("\tPRN %d: bgde5a\t%e\n", prn, eph->bgde5a);
	printf("\tPRN %d: bgde5b\t%e\n", prn, eph->bgde5b);
	printf("\tPRN %d: cic\t%e\n", prn, eph->cic);
	printf("\tPRN %d: cis\t%e\n", prn, eph->cis);
	printf("\tPRN %d: crc\t%e\n", prn, eph->crc);
	printf("\tPRN %d: crs\t%e\n", prn, eph->crs);
	printf("\tPRN %d: cuc\t%e\n ", prn, eph->cuc);
	printf("\tPRN %d: cus\t%e\n", prn, eph->cus);
	printf("\tPRN %d: deltan\t%e\n", prn, eph->deltan);
	printf("\tPRN %d: ecc\t%e\n", prn, eph->ecc);
	printf("\tPRN %d: idot\t%e\n", prn, eph->idot);
	printf("\tPRN %d: inc0\t%e\n", prn, eph->inc0);
	printf("\tPRN %d: iodnav\t%d\n", prn, eph->iode); // iodnav);
	printf("\tPRN %d: m0\t%e\n", prn, eph->m0);
	printf("\tPRN %d: n\t%e\n", prn, eph->n);
	printf("\tPRN %d: omg0\t%e\n", prn, eph->omg0);
	printf("\tPRN %d: omgdot\t%e\n", prn, eph->omgdot);
	printf("\tPRN %d: omgkdot\t%e\n", prn, eph->omgkdot);
	printf("\tPRN %d: sq1e2\t%e\n", prn, eph->sq1e2);
	printf("\tPRN %d: sqrta\t%e\n", prn, eph->sqrta);
	printf("\tPRN %d: toc\t%0.f\n", prn, eph->toc.sec);
	printf("\tPRN %d: toe\t%0.f\n", prn, eph->toe.sec);
}

void write_csv(ephem_t *eph, int prn)
{
	ofstream csv;
	csv.open("./eph1.csv", ios::app);
	csv << prn << "," << eph->time.sec << "," << eph->toc.sec << "," << eph->toe.sec << "," << eph->iode << "," << eph->deltan << "," << eph->cuc << "," << eph->cus << "," << eph->cic << "," << eph->cis << "," << eph->crc << "," << eph->crs << "," << eph->ecc << "," << eph->sqrta << "," << eph->m0 << "," << eph->omg0 << "," << eph->inc0 << "," << eph->aop << "," << eph->omgdot << "," << eph->idot << "," << eph->af0 << "," << eph->af1 << "," << eph->af2 << "," << eph->bgde5a << "," << eph->bgde5b << "," << eph->n << "," << eph->sq1e2 << "," << eph->A << "," << eph->omgkdot << ","
		<< "\n";
	csv.close();
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

vector<uint8_t> load_page(FILE *TV_ptr, int prn, bool advance_fptr, int tow, galtime_t *g)
{
	char tv_str[MAX_CHAR] = {0};
	char page[MAX_CHAR] = {0};
	unsigned long tow_prev = 0;
	int field_count;
	int str_len;
	galtime_t galtime;
	int satid;
	int k = 0;

	vector<uint8_t> nav_page;
	long new_startingPosition;

	if ((tow != tow_prev))
	{
		tow_prev = tow;
		long startingPosition = ftell(TV_ptr);

		while ((fgets(tv_str, MAX_CHAR, TV_ptr) != NULL))
		{
			/* Get len of line to go back on previous line */
			str_len = strlen(tv_str);
			k++;

			/* Extract data from test vectors */
			field_count = sscanf(tv_str, "%lf,%d,%d,%99s", &galtime.sec, &galtime.week, &satid, page);
			// std::cout << tv_str <<endl;
			//  if (galtime.sec < tow)
			//  {
			//  	advance_fptr = true;
			//  	continue;
			//  }

			// if (galtime.sec < tow && advance_fptr)
			// 	new_startingPosition = ftell(TV_ptr) - str_len;

			if (field_count == 4)
			{
				// if ( galtime.sec > tow)
				// {
				// 	/* Go back to previous line, with +1 for \n*/
				// 	fseek(TV_ptr, ((str_len+1)*-1), SEEK_CUR);
				// 	break;
				// }
				if (prn == satid)
				{
					long f;
					if (prn == 1)
						f = ftell(TV_ptr);

					*g = galtime;
					chartouint8_t(&nav_page, page, 30);
					// fseek(TV_ptr, ((str_len+1)*-1), SEEK_CUR);
					break;
				}

				if (nav_page.empty())
				{
					fprintf(stderr, "\n Returning empty \n");
				}
				// fprintf(stderr, "\n HERE %d %d\n ", tow, pages[satid-1].size());
				// fprintf(stderr, "\n Satid %d \n ", satid);
				//  for (int i=0; i<30; i++)
				//  	fprintf(stderr, "%d ", pages[satid-1][i]);
			}
		}
		// long fseekpos = ftell(TV_ptr);

		// if (!advance_fptr)
		// 	fseek(TV_ptr, startingPosition, SEEK_SET);
		// else
		// 	fseek(TV_ptr, new_startingPosition, SEEK_SET);
	}
	return nav_page;
}

/* Read each line and extract ephemeris Words 1-5 */
void extract_ephemeris(vector<uint8_t> *nav_msg_uint8, ephem_t *eph, uint8_t wordtype, int svid)
{
	vector<uint8_t> word(16, 0);
	int offsetIn, datalen, res;
	uint32_t tmp32;
	uint16_t tmp16;
	/* Extract the data from even part (14 bytes) */
	data_extract(&word, *nav_msg_uint8, 0, 2, 14 * 8);

	/* Extract the data from odd part (2 bytes) */
	data_extract(&word, *nav_msg_uint8, 14 * 8, 122, 16);

	/* Refer to Galileo_OS_SIS_ICD_v1.2 - 4.3.5 */
	switch (wordtype)
	{
	case 1:
		/* Wordtype, IODnav*/
		offsetIn = 16;

		/* t0e len */
		datalen = 14;
		eph->toe.sec = (double)data_extract_uint32_t(word, offsetIn, datalen);
		eph->toe.sec *= 60.0;
		// printf("W%d toe.tow[%d]: %.2f\n",wordtype, satid, eph->toe.tow);
		offsetIn += datalen;

		/* M0 len */
		datalen = 32;
		tmp32 = data_extract_uint32_t(word, offsetIn, datalen);
		eph->m0 = (double)((int32_t)tmp32) * POW2_M31 * PI;
		offsetIn += datalen;

		// exit(1);

		/* Eccentricity len */
		datalen = 32;
		eph->ecc = (double)data_extract_uint32_t(word, offsetIn, datalen);
		eph->ecc *= POW2_M33;
		offsetIn += datalen;

		/* A1/2 len */
		datalen = 32;
		eph->sqrta = (double)data_extract_uint32_t(word, offsetIn, datalen);
		eph->sqrta *= POW2_M19;
		break;

	case 2:

		/*Wordtype */
		offsetIn = 6;

		/* IODnav len */
		datalen = 10;
		eph->iode = (int)data_extract_uint32_t(word, offsetIn, datalen);
		offsetIn += datalen;

		/* Omega0 len */
		datalen = 32;
		tmp32 = data_extract_uint32_t(word, offsetIn, datalen);

		eph->omg0 = (double)((int32_t)tmp32) * POW2_M31 * PI;
		offsetIn += datalen;

		/* i0 ken*/
		datalen = 32;
		tmp32 = data_extract_uint32_t(word, offsetIn, datalen);
		eph->inc0 = (double)((int32_t)tmp32) * POW2_M31 * PI;
		offsetIn += datalen;

		/* Argument of perigee len=32 */
		datalen = 32;
		tmp32 = data_extract_uint32_t(word, offsetIn, datalen);
		eph->aop = (double)((int32_t)tmp32) * POW2_M31 * PI;
		offsetIn += datalen;

		/* i dot Rate of change of incluation angle */
		datalen = 14;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		res = data_extract_signedInt(datalen, res);
		eph->idot = (double)res * POW2_M43 * PI;
		break;

	case 3:
		/*Wordtype, IODnav*/
		offsetIn = 16;

		/* omega dot len*/
		datalen = 24;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		res = data_extract_signedInt(datalen, res);
		eph->omgdot = (double)res * POW2_M43 * PI;
		offsetIn += datalen;

		/* Delta n len*/
		datalen = 16;
		tmp16 = (uint16_t)data_extract_uint32_t(word, offsetIn, datalen);
		eph->deltan = (double)((int16_t)tmp16) * POW2_M43 * PI;
		offsetIn += datalen;

		/* Cuc len*/
		datalen = 16;
		tmp16 = (uint16_t)data_extract_uint32_t(word, offsetIn, datalen);
		eph->cuc = (double)(int16_t)tmp16 * POW2_M29;
		offsetIn += datalen;

		/* Cus len*/
		datalen = 16;
		tmp16 = (uint16_t)data_extract_uint32_t(word, offsetIn, datalen);
		eph->cus = (double)(int16_t)tmp16 * POW2_M29;
		offsetIn += datalen;

		/* Crc len*/
		datalen = 16;
		tmp16 = (uint16_t)data_extract_uint32_t(word, offsetIn, datalen);
		eph->crc = (double)(int16_t)tmp16 * POW2_M5;
		offsetIn += datalen;

		/* Crs len*/
		datalen = 16;
		tmp16 = (uint16_t)data_extract_uint32_t(word, offsetIn, datalen);
		eph->crs = (double)(int16_t)tmp16 * POW2_M5;

		break;

	case 4:
		/* Wordtype, IODnav, SVID*/
		offsetIn = 16;

		datalen = 6;
		eph->svid = data_extract_uint32_t(word, offsetIn, datalen);
		offsetIn += datalen;

		/* Cic len*/
		datalen = 16;
		tmp16 = (uint16_t)data_extract_uint32_t(word, offsetIn, datalen);
		eph->cic = (double)(int16_t)tmp16 * POW2_M29;
		offsetIn += datalen;

		/* Cis len*/
		datalen = 16;
		tmp16 = (uint16_t)data_extract_uint32_t(word, offsetIn, datalen);
		eph->cis = (double)(int16_t)tmp16 * POW2_M29;
		offsetIn += datalen;

		/* T0c */
		datalen = 14;
		eph->toc.sec = (double)data_extract_uint32_t(word, offsetIn, datalen);
		eph->toc.sec *= 60.0;
		offsetIn += datalen;

		/* Af0 */
		datalen = 31;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		res = data_extract_signedInt(datalen, res);
		eph->af0 = (double)res * POW2_M34;
		offsetIn += datalen;

		/* Af1 */
		datalen = 21;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		res = data_extract_signedInt(datalen, res);
		eph->af1 = (double)res * POW2_M46;
		// printf("W%d af1[%d]: %e\n",wordtype, satid, eph->af1);
		offsetIn += datalen;

		/* Af2 */
		datalen = 6;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		res = data_extract_signedInt(datalen, res);
		eph->af2 = (double)res * POW2_M59;

		break;

	case 5:
		/* offset BGD */
		offsetIn = 47;

		/* BGD(E1,E5a) E1-E5a Broadcast Group Delay */
		datalen = 10;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		res = data_extract_signedInt(datalen, res);
		eph->bgde5a = (double)res * POW2_M32;
		offsetIn += datalen;

		/* BGD(E1,E5b) E1-E5b Broadcast Group Delay */
		datalen = 10;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		res = data_extract_signedInt(datalen, res);
		eph->bgde5b = (double)res * POW2_M32;

		/* Move forward by 6 bits to get to WN and TOW - store them as time */
		offsetIn += datalen;
		datalen = 6;
		offsetIn += datalen;

		datalen = 12;

		res = data_extract_uint32_t(word, offsetIn, datalen);
		eph->toc.week = res + 1024;
		eph->toe.week = res + 1024;
		eph->time.week = res + 1024;

		offsetIn += datalen;

		datalen = 20;
		res = data_extract_uint32_t(word, offsetIn, datalen);
		eph->time.sec = res;

		/* Update the working variables */
		eph->A = pow(eph->sqrta, 2);
		eph->n = sqrt(GM_EARTH / pow(eph->A, 3)) + eph->deltan;
		eph->sq1e2 = sqrt(1.0 - pow(eph->ecc, 2));
		eph->omgkdot = eph->omgdot - OMEGA_EARTH;

		eph->PRN = svid;

		/* Eph sat complete */
		eph->vflg = true;
		break;
	}
}

vector<ephem_t> load_ephemeris(const char *fname)
{
	char str[MAX_CHAR];

	galtime_t galtime;
	int satid;
	uint32_t tmp;
	char nav_msg[MAX_CHAR] = {0};
	uint8_t wordtype;

	FILE *fp;
	fp = fopen(fname, "rt");

	// if (NULL == (fp = fopen(fname, "rt")))
	// {
	// 	printf("\nERROR: Can't open test vectors");
	// 	return 0;
	// }

	/*
		Initialize ephem vector
		This vector contains MAX_SAT vectors that hold ephemeris records for each satellite.
	*/

	int index = 0;

	vector<ephem_t> eph_vector;

	ephem_t blank_eph = {0};
	eph_vector.push_back(blank_eph);

	if (fp == NULL)
		return eph_vector;

	vector<uint8_t> nav_msg_uint8;

	/* Read the entire file and gather ephemeris data for all satellites*/
	while ((fgets(str, MAX_CHAR, fp) != NULL))
	{

		sscanf(str, "%lf,%d,%d,%99s", &galtime.sec, &galtime.week, &satid, nav_msg);

		// fprintf(stderr,"\nTime %f : ",s);
		// for(int i = 0; i < 60; i++)
		//     fprintf(stderr,"%c",nav_msg[i]);
		if (strlen(nav_msg) != 60)
		{
			continue;
		}

		/* Convert in uint8_t format */
		chartouint8_t(&nav_msg_uint8, nav_msg, 30);

		/* Start with Even and nominal page-type (0,0) */
		tmp = data_extract_uint32_t(nav_msg_uint8, 0, 2);

		if (tmp != 0)
		{
			continue;
		}

		/* Extract wordtype */
		wordtype = data_extract_uint32_t(nav_msg_uint8, 2, 6);

		// fprintf(stderr,"\n Word %d", wordtype);

		/* Updates the ephemeris of only the current satellite. At word 5 it will increment the current_eph_index value by 1 */
		extract_ephemeris(&nav_msg_uint8, &eph_vector[index], wordtype, satid);

		// print_eph(&eph_vector[satid-1][current_eph_index[satid-1]-1], satid);

		if (eph_vector[index].vflg == 1)
		{
			write_csv(&eph_vector[index], satid);
			index++;
			blank_eph = {0};
			eph_vector.push_back(blank_eph);
		}

		nav_msg_uint8.clear();
	}

	fprintf(stderr, "\nLoaded %d ephemeris entries for SV-%d", (int)eph_vector.size(), satid);
	return eph_vector;
}

/* Matches TOW and relevant ephemeris */
int epoch_matcher(double obsTime, vector<ephem_t> eph, int index)
{
	for (int i = index; i < eph.size(); i++)
	{
		double diff = fabs(obsTime - eph.at(i).toe.sec);
		if (diff < 0)
		{
			index = i - 1; // save index
			break;
		}
	}
	return index;
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