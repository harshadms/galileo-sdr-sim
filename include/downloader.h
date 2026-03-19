#ifndef DOWNLOADER_H
#define DOWNLOADER_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Download latest Galileo navigation (RINEX) file
 * 
 * Tries multiple sources (CDDIS, Garner, CENO) in order.
 * Caches downloaded files for reuse.
 * Falls back to cached version if download fails.
 * 
 * @param output_file    Path where downloaded file will be saved
 * @param max_age_days   Maximum acceptable file age in days (0 = any age, 7 = use if <7 days old)
 * @return               0 = success, 1 = all sources failed, 2 = file too old
 */
int download_latest_galileo_nav(char *output_file, int max_age_days);

/**
 * Get current GPS date and time
 * 
 * @param year           Output: Current year (e.g., 2026)
 * @param doy            Output: Day of year (1-366)
 * @param gps_week       Output: GPS week number (0-4095+)
 * @param gps_day        Output: Day of GPS week (0-6, where 0=Sunday)
 */
void get_gps_date(int *year, int *doy, int *gps_week, int *gps_day);

/**
 * Download from CDDIS (NASA/USGS primary source)
 * 
 * @param year           Year (e.g., 2026)
 * @param doy            Day of year (1-366)
 * @param output_file    Output file path
 * @return               0 = success, 1 = failure
 */
int download_from_cddis(int year, int doy, char *output_file);

/**
 * Download from Garner (IGS backup)
 * 
 * @param year           Year (e.g., 2026)
 * @param doy            Day of year (1-366)
 * @param output_file    Output file path
 * @return               0 = success, 1 = failure
 */
int download_from_garner(int year, int doy, char *output_file);

/**
 * Download from BKG (Bundesamt für Kartographie und Geodäsie - Germany)
 * 
 * German government source with multi-GNSS support and excellent reliability.
 * Provides both daily and near-real-time (15-min) updates.
 * 
 * @param year           Year (e.g., 2026)
 * @param doy            Day of year (1-366)
 * @param output_file    Output file path
 * @return               0 = success, 1 = failure
 */
int download_from_bkg(int year, int doy, char *output_file);

/**
 * Download from CENO (ESA backup)
 * 
 * @param year           Year (e.g., 2026)
 * @param doy            Day of year (1-366)
 * @param output_file    Output file path
 * @return               0 = success, 1 = failure
 */
int download_from_ceno(int year, int doy, char *output_file);

/**
 * Decompress gzipped RINEX file
 * 
 * @param gzipped_file   Path to .gz file
 * @param output_file    Output path for decompressed file
 * @return               0 = success, 1 = failure
 */
int decompress_gzip(const char *gzipped_file, const char *output_file);

/**
 * Decompress Z-compressed RINEX file (Uncompress)
 * 
 * @param compressed_file Path to .Z file
 * @param output_file     Output path for decompressed file
 * @return                0 = success, 1 = failure
 */
int decompress_z(const char *compressed_file, const char *output_file);

/**
 * Ensure cache directory exists and is writable
 * 
 * @param cache_dir      Cache directory path
 * @return               0 = success, 1 = failure
 */
int ensure_cache_dir(const char *cache_dir);

/**
 * Get file modification time
 * 
 * @param filepath       Path to file
 * @return               Time in seconds since epoch, or -1 if file not found
 */
long get_file_mtime(const char *filepath);

/**
 * Check if file is recent enough
 * 
 * @param filepath       Path to file
 * @param max_age_secs   Maximum acceptable age in seconds
 * @return               1 = file is recent, 0 = file is too old or not found
 */
int is_file_recent(const char *filepath, long max_age_secs);

#ifdef __cplusplus
}
#endif

#endif  // DOWNLOADER_H
