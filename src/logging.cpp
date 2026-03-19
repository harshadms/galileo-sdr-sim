#include "../include/logging.h"

log_level_t g_log_level = LOG_INFO;

void logging_init(log_level_t level) {
    g_log_level = level;
}

const char* _log_level_str(log_level_t level) {
    switch (level) {
        case LOG_DEBUG: return "DEBUG";
        case LOG_INFO:  return "INFO";
        case LOG_WARN:  return "WARN";
        case LOG_ERROR: return "ERROR";
        default:        return "UNKNOWN";
    }
}

const char* _log_color(log_level_t level) {
    switch (level) {
        case LOG_DEBUG: return "\033[36m";  // Cyan
        case LOG_INFO:  return "\033[32m";  // Green
        case LOG_WARN:  return "\033[33m";  // Yellow
        case LOG_ERROR: return "\033[31m";  // Red
        default:        return "\033[0m";   // Reset
    }
}

void _log_message(log_level_t level, const char *format, va_list args) {
    if (level < g_log_level) {
        return;
    }
    
    time_t now = time(NULL);
    struct tm *timeinfo = localtime(&now);
    char timestamp[32];
    strftime(timestamp, sizeof(timestamp), "%H:%M:%S", timeinfo);
    
    // Print with color and timestamp
    fprintf(stderr, "%s[%s] %s: \033[0m", _log_color(level), timestamp, _log_level_str(level));
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    fflush(stderr);
}

void log_debug(const char *format, ...) {
    va_list args;
    va_start(args, format);
    _log_message(LOG_DEBUG, format, args);
    va_end(args);
}

void log_info(const char *format, ...) {
    va_list args;
    va_start(args, format);
    _log_message(LOG_INFO, format, args);
    va_end(args);
}

void log_warn(const char *format, ...) {
    va_list args;
    va_start(args, format);
    _log_message(LOG_WARN, format, args);
    va_end(args);
}

void log_error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    _log_message(LOG_ERROR, format, args);
    va_end(args);
}

void log_progress(const char *format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fflush(stderr);
}
