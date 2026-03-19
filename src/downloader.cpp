#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <curl/curl.h>
#include <zlib.h>
#include "../include/downloader.h"

#define DOWNLOAD_TIMEOUT 30
#define CACHE_DIR "rinex_cache"
#define MAX_FILEPATH 512

// Callback for curl write data
static size_t write_callback(void *contents, size_t size, size_t nmemb, FILE *fp) {
    return fwrite(contents, size, nmemb, fp);
}

// Callback for curl progress
static int progress_callback(void *clientp, curl_off_t dltotal, curl_off_t dlnow,
                             curl_off_t ultotal, curl_off_t ulnow) {
    if (dltotal > 0) {
        int percent = (int)(dlnow * 100 / dltotal);
        fprintf(stderr, "\r    [%3d%%] %ld / %ld bytes", percent, dlnow, dltotal);
        fflush(stderr);
    }
    return 0;
}

// Get current GPS date
void get_gps_date(int *year, int *doy, int *gps_week, int *gps_day) {
    time_t now = time(NULL);
    struct tm *tm_info = gmtime(&now);
    
    *year = tm_info->tm_year + 1900;
    *doy = tm_info->tm_yday;
    
    // GPS epoch: January 6, 1980 = Day 6, 1980
    // Current GPS week calculation
    time_t gps_epoch = 315964800;  // Jan 6, 1980 00:00:00 UTC
    long seconds_since_epoch = now - gps_epoch;
    *gps_week = (int)(seconds_since_epoch / (7 * 86400));
    *gps_day = (int)((seconds_since_epoch % (7 * 86400)) / 86400);
}

// Ensure cache directory exists
int ensure_cache_dir(const char *cache_dir) {
    struct stat st = {0};
    if (stat(cache_dir, &st) == -1) {
        // Directory doesn't exist - create it
        char cmd[512];
        snprintf(cmd, sizeof(cmd), "mkdir -p %s", cache_dir);
        if (system(cmd) != 0) {
            fprintf(stderr, "ERROR: Failed to create cache directory: %s\n", cache_dir);
            return 1;
        }
    }
    return 0;
}

// Get file modification time
long get_file_mtime(const char *filepath) {
    struct stat st;
    if (stat(filepath, &st) == 0) {
        return st.st_mtime;
    }
    return -1;
}

// Check if file is recent
int is_file_recent(const char *filepath, long max_age_secs) {
    long mtime = get_file_mtime(filepath);
    if (mtime == -1) return 0;  // File doesn't exist
    
    if (max_age_secs == 0) return 1;  // Any age acceptable
    
    time_t now = time(NULL);
    long age = now - mtime;
    return age <= max_age_secs;
}

// Decompress gzip file
int decompress_gzip(const char *gzipped_file, const char *output_file) {
    gzFile infile = gzopen(gzipped_file, "rb");
    if (!infile) {
        fprintf(stderr, "ERROR: Failed to open gzipped file: %s\n", gzipped_file);
        return 1;
    }
    
    FILE *outfile = fopen(output_file, "wb");
    if (!outfile) {
        fprintf(stderr, "ERROR: Failed to create output file: %s\n", output_file);
        gzclose(infile);
        return 1;
    }
    
    unsigned char buffer[65536];
    int bytes;
    while ((bytes = gzread(infile, buffer, sizeof(buffer))) > 0) {
        fwrite(buffer, 1, bytes, outfile);
    }
    
    fclose(outfile);
    gzclose(infile);
    return 0;
}

// Download from CDDIS (NASA/USGS)
int download_from_cddis(int year, int doy, char *output_file) {
    char url[512];
    char temp_file[MAX_FILEPATH];
    int year_short = year % 100;
    
    // CDDIS format: ftp://cddis.eosdis.nasa.gov/pub/gps/data/daily/YYYY/DDD/brdc*d.YYn.Z
    snprintf(url, sizeof(url), 
             "ftp://cddis.eosdis.nasa.gov/pub/gps/data/daily/%04d/%03d/brdc%03d0.%02dn.Z",
             year, doy, doy, year_short);
    
    snprintf(temp_file, sizeof(temp_file), "%s/cddis_temp.Z", CACHE_DIR);
    
    fprintf(stderr, "    [*] Trying CDDIS: %s\n", url);
    
    CURL *curl = curl_easy_init();
    if (!curl) {
        fprintf(stderr, "    [-] CURL initialization failed\n");
        return 1;
    }
    
    FILE *fp = fopen(temp_file, "wb");
    if (!fp) {
        fprintf(stderr, "    [-] Cannot open temp file\n");
        curl_easy_cleanup(curl);
        return 1;
    }
    
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, DOWNLOAD_TIMEOUT);
    curl_easy_setopt(curl, CURLOPT_XFERINFOFUNCTION, progress_callback);
    curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0L);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
    
    CURLcode res = curl_easy_perform(curl);
    fclose(fp);
    
    if (res != CURLE_OK) {
        fprintf(stderr, "\n    [-] Download failed: %s\n", curl_easy_strerror(res));
        curl_easy_cleanup(curl);
        remove(temp_file);
        return 1;
    }
    
    curl_easy_cleanup(curl);
    
    // Decompress
    fprintf(stderr, "\n    [*] Decompressing...\n");
    if (decompress_z(temp_file, output_file) != 0) {
        fprintf(stderr, "    [-] Decompression failed\n");
        remove(temp_file);
        return 1;
    }
    
    remove(temp_file);
    fprintf(stderr, "    [+] Success!\n");
    return 0;
}

// Download from Garner (IGS)
int download_from_garner(int year, int doy, char *output_file) {
    char url[512];
    char temp_file[MAX_FILEPATH];
    int year_short = year % 100;
    
    snprintf(url, sizeof(url),
             "http://garner.ucsd.edu/pub/rinex/nav/%04d/%03d/brdc%03d0.%02dn.Z",
             year, doy, doy, year_short);
    
    snprintf(temp_file, sizeof(temp_file), "%s/garner_temp.Z", CACHE_DIR);
    
    fprintf(stderr, "    [*] Trying Garner (IGS): %s\n", url);
    
    CURL *curl = curl_easy_init();
    if (!curl) {
        fprintf(stderr, "    [-] CURL initialization failed\n");
        return 1;
    }
    
    FILE *fp = fopen(temp_file, "wb");
    if (!fp) {
        fprintf(stderr, "    [-] Cannot open temp file\n");
        curl_easy_cleanup(curl);
        return 1;
    }
    
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, DOWNLOAD_TIMEOUT);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
    
    CURLcode res = curl_easy_perform(curl);
    fclose(fp);
    
    if (res != CURLE_OK) {
        fprintf(stderr, "    [-] Download failed: %s\n", curl_easy_strerror(res));
        curl_easy_cleanup(curl);
        remove(temp_file);
        return 1;
    }
    
    curl_easy_cleanup(curl);
    
    fprintf(stderr, "    [*] Decompressing...\n");
    if (decompress_z(temp_file, output_file) != 0) {
        fprintf(stderr, "    [-] Decompression failed\n");
        remove(temp_file);
        return 1;
    }
    
    remove(temp_file);
    fprintf(stderr, "    [+] Success!\n");
    return 0;
}

// Download from BKG (Bundesamt für Kartographie und Geodäsie - Germany)
int download_from_bkg(int year, int doy, char *output_file) {
    char url[512];
    char temp_file[MAX_FILEPATH];
    
    // BKG provides BRDC multi-GNSS consolidated file
    // Format: https://igs.bkg.bund.de/root_ftp/IGS/BRDC/YYYY/DDD/BRDC00WRD_S_YYYYDDD0000_01D_MN.rnx.gz
    // File: BRDC00WRD_S (IGS standard broadcast ephemeris with all constellations)
    snprintf(url, sizeof(url),
             "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/%04d/%03d/BRDC00WRD_S_%04d%03d0000_01D_MN.rnx.gz",
             year, doy, year, doy);
    
    snprintf(temp_file, sizeof(temp_file), "%s/bkg_temp.gz", CACHE_DIR);
    
    fprintf(stderr, "    [*] Trying BKG (German Gov): %s\n", url);
    
    CURL *curl = curl_easy_init();
    if (!curl) {
        fprintf(stderr, "    [-] CURL initialization failed\n");
        return 1;
    }
    
    FILE *fp = fopen(temp_file, "wb");
    if (!fp) {
        fprintf(stderr, "    [-] Cannot open temp file\n");
        curl_easy_cleanup(curl);
        return 1;
    }
    
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, DOWNLOAD_TIMEOUT);
    curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
    
    CURLcode res = curl_easy_perform(curl);
    fclose(fp);
    
    if (res != CURLE_OK) {
        fprintf(stderr, "    [-] Download failed: %s\n", curl_easy_strerror(res));
        curl_easy_cleanup(curl);
        remove(temp_file);
        return 1;
    }
    
    curl_easy_cleanup(curl);
    
    fprintf(stderr, "    [*] Decompressing...\n");
    if (decompress_gzip(temp_file, output_file) != 0) {
        fprintf(stderr, "    [-] Decompression failed\n");
        remove(temp_file);
        return 1;
    }
    
    remove(temp_file);
    fprintf(stderr, "    [+] Success!\n");
    return 0;
}

// Download from CENO (ESA)
int download_from_ceno(int year, int doy, char *output_file) {
    char url[512];
    char temp_file[MAX_FILEPATH];
    
    snprintf(url, sizeof(url),
             "ftp://ceno.esa.int/GNSS/public/data/daily/%04d/%03d/",
             year, doy);
    
    snprintf(temp_file, sizeof(temp_file), "%s/ceno_temp.rnx", CACHE_DIR);
    
    fprintf(stderr, "    [*] Trying CENO (ESA): %s\n", url);
    fprintf(stderr, "    [!] CENO requires specific file listing (skipped for now)\n");
    return 1;  // Skip CENO for now (requires directory listing)
}

// Main downloader function
int download_latest_galileo_nav(char *output_file, int max_age_days) {
    int year, doy, gps_week, gps_day;
    long max_age_secs = (max_age_days > 0) ? max_age_days * 86400 : 0;
    
    // Get current date
    get_gps_date(&year, &doy, &gps_week, &gps_day);
    fprintf(stderr, "[+] Detected GPS week %d, DOY %d (%04d-%02d-%02d)\n",
            gps_week, doy, year, doy/30 + 1, doy % 30 + 1);
    
    // Ensure cache directory exists
    if (ensure_cache_dir(CACHE_DIR) != 0) {
        return 1;
    }
    
    // Check if recent cached file exists
    if (is_file_recent(output_file, max_age_secs)) {
        fprintf(stderr, "[+] Using cached navigation file: %s\n", output_file);
        return 0;
    }
    
    fprintf(stderr, "[*] Downloading latest Galileo navigation file...\n");
    
    // Try sources in order
    fprintf(stderr, "[*] Trying download sources:\n");
    
    if (download_from_cddis(year, doy, output_file) == 0) {
        return 0;
    }
    
    if (download_from_garner(year, doy, output_file) == 0) {
        return 0;
    }
    
    // Try BKG (German government - excellent reliability + European redundancy)
    if (download_from_bkg(year, doy, output_file) == 0) {
        return 0;
    }
    
    if (download_from_ceno(year, doy, output_file) == 0) {
        return 0;
    }
    
    // Check if we have a cached file to fall back to
    if (get_file_mtime(output_file) != -1) {
        fprintf(stderr, "[!] All sources failed. Using cached file (may be outdated).\n");
        return 0;  // Use cached file anyway
    }
    
    fprintf(stderr, "[-] All download sources failed and no cached file available.\n");
    return 1;
}

// Stub for decompress_z (would need uncompress utility or zlib support)
int decompress_z(const char *compressed_file, const char *output_file) {
    // Use system uncompress command
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "uncompress -c %s > %s 2>/dev/null", 
             compressed_file, output_file);
    
    int ret = system(cmd);
    if (ret != 0) {
        // Try with gzip as fallback
        snprintf(cmd, sizeof(cmd), "gunzip -c %s > %s 2>/dev/null",
                 compressed_file, output_file);
        ret = system(cmd);
    }
    
    return (ret == 0) ? 0 : 1;
}
