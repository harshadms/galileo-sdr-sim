#ifndef LOGGING_H
#define LOGGING_H

#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <string.h>

typedef enum {
    LOG_DEBUG = 0,
    LOG_INFO = 1,
    LOG_WARN = 2,
    LOG_ERROR = 3,
    LOG_NONE = 4
} log_level_t;

extern log_level_t g_log_level;

// Initialize logging system
void logging_init(log_level_t level);

// Log functions with level control
void log_debug(const char *format, ...);
void log_info(const char *format, ...);
void log_warn(const char *format, ...);
void log_error(const char *format, ...);

// Progress logging (doesn't respect log level, always shown)
void log_progress(const char *format, ...);

// Internal logging function
void _log_message(log_level_t level, const char *format, va_list args);

#endif // LOGGING_H
