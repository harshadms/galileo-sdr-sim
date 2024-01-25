#include "../include/galileo-sdr.h"

/*! \brief Used to debug - remove */
int writecsv(int prn, double tow, double values[3], int type)
{
    FILE *csvFile = fopen("data.csv", "a");

    if (csvFile == NULL)
    {
        perror("Error opening file");
        return -1;
    }

    // Write the values to the CSV file
    switch (type)
    {
    case 0:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, satpos", prn, values[0], values[1], values[2], tow);
        break;

    case 1:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, los", prn, values[0], values[1], values[2], tow);
        break;

    case 2:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, rxpos", prn, values[0], values[1], values[2], tow);
        break;

    case 3:
        fprintf(csvFile, "%d,%.4f,%.4f,%.4f,%f, range", prn, values[0], values[1], values[2], tow);
        break;
    }

    // Add a newline to separate rows
    fprintf(csvFile, "\n");

    // Close the CSV file
    fclose(csvFile);

    return 1;
}


void print_eph2(ephem_t *eph, int prn)
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
	printf("\tPRN %d: toc\t%0.f - %d\n", prn, eph->toc.sec, eph->toc.week);
	printf("\tPRN %d: toe\t%0.f - %d\n", prn, eph->toe.sec, eph->toe.week);
	printf("\tPRN %d: flag\t%d\n", prn, eph->flag);
}
