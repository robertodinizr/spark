
#ifndef LOG_H
#define LOG_H

#define LOG_COLOR_RESET "\033[0m"
#define LOG_COLOR_RED "\033[31m"
#define LOG_COLOR_YELLOW "\033[33m"
#define LOG_COLOR_BLUE "\033[34m"

#ifdef SPARK_ENABLE_LOG_DEBUG
#include <cstdio>
#define SPARK_LOG_DEBUG(msg, ...)                                                 \
    {                                                                             \
        printf(LOG_COLOR_BLUE "[spark-debug] " msg LOG_COLOR_RESET, __VA_ARGS__); \
    }
#else
#define SPARK_LOG_DEBUG(msg, ...)
#endif

#ifdef SPARK_ENABLE_LOG_INFO
#include <cstdio>
#define SPARK_LOG_INFO(msg, ...)                  \
    {                                             \
        printf("[spark-info] " msg, __VA_ARGS__); \
    }
#else
#define SPARK_LOG_INFO(msg, ...)
#endif

#ifdef SPARK_ENABLE_LOG_WARN
#include <cstdio>
#define SPARK_LOG_WARN(msg, ...)                                                      \
    {                                                                                 \
        printf(LOG_COLOR_YELLOW "[spark-warning] " msg LOG_COLOR_RESET, __VA_ARGS__); \
    }
#else
#define SPARK_LOG_WARN(msg, ...)
#endif

#ifdef SPARK_ENABLE_LOG_ERROR
#include <cstdio>
#define SPARK_LOG_ERROR(msg, ...)                                                \
    {                                                                            \
        printf(LOG_COLOR_RED "[spark-error] " msg LOG_COLOR_RESET, __VA_ARGS__); \
    }
#else
SPARK_LOG_ERROR(msg, ...)
#endif

#endif  // LOG_H
