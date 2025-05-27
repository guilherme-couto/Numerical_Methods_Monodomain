#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "../../include/config_parser.h"
#include "../../include/core_definitions.h"

extern FILE *log_file;

// ANSI color codes (for terminal only)
#define RESET   "\033[0m"
#define BLUE    "\033[34m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define RED     "\033[31m"

#ifdef __cplusplus
extern "C" {
#endif

void init_logger(const char *filename);
void close_logger(void);
void log_simulation_header(const SimulationConfig *config, const char *config_filename);
void log_machine_info(void);
void log_message(const char *color, const char *prefix, const char *fmt, ...);
void log_device_info(FILE *log_file);

#ifdef __cplusplus
}
#endif

// Macros for logging
#define INFOMSG(fmt, ...)    log_message(BLUE,   "[i] ", fmt, ##__VA_ARGS__)
#define SUCCESSMSG(fmt, ...) log_message(GREEN,  "[+] ", fmt, ##__VA_ARGS__)
#define WARNINGMSG(fmt, ...) log_message(YELLOW, "[!] ", fmt, ##__VA_ARGS__)
#define ERRORMSG(fmt, ...)   log_message(RED,    "[x] ", fmt "\n", ##__VA_ARGS__)
#define DEBUGMSG(fmt, ...)   log_message(YELLOW, "[.] ", fmt, ##__VA_ARGS__)
#define LINEMSG()            log_message(YELLOW, "[.] ", "Line number %d in file %s\n", __LINE__, __FILE__)
#define SIMPLEMSG(fmt, ...) log_message(RESET, "", fmt "\n", ##__VA_ARGS__)

#endif // LOGGER_H
