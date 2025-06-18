#pragma once

#include <iostream>
#include <fstream>
#include <memory>
#include <string_view>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <math.h>

const int32_t barWidth = 75;  // Total number of '=' characters for 100% progress
const float barUnit = 100.0 / static_cast<float>(barWidth);  // Each '=' represents 5%

enum LoggerLevel{
    LOG_DEBUG,
    LOG_INFO,
    LOG_WARNING,
    LOG_ERROR,
    LOG_FATAL
};

constexpr const char* SeverityTable(LoggerLevel level) {
    switch(level) {
        case LOG_DEBUG:    return "[debug]    ";
        case LOG_INFO:     return "[info]     ";
        case LOG_WARNING:  return "[warning]  ";
        case LOG_ERROR:    return "[error]    ";
        case LOG_FATAL:    return "[fatal]    ";
        default:           return "[unknown]  ";
    }
}

class Logger {
private:
    Logger() = default;

protected:
    LoggerLevel minSeverityLevel = LOG_DEBUG;
    LoggerLevel severity = LOG_INFO;

public:
    static Logger& getLogger()
    {
        static Logger instance;
        // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }
    template<typename FirstParam, typename... OtherParams>
    void log(LoggerLevel level, const FirstParam& firstParam, const OtherParams&... params) {
        if (level >= minSeverityLevel) {
            auto now = std::chrono::system_clock::now();
            std::time_t timeNow = std::chrono::system_clock::to_time_t(now);

            std::cout << std::put_time(localtime(&timeNow), "[%Y-%m-%d %H:%M:%S] ") << SeverityTable(level) << firstParam;
            if constexpr (sizeof...(params) > 0) {
                ((std::cout << params), ...);
            }
            std::cout << std::endl;
        }
    }

    void set_min_severity(LoggerLevel severity) {
        minSeverityLevel = severity;
    }

    void show_progress(int32_t percentage);

    void init_logger(LoggerLevel minLevel = LOG_DEBUG, const std::string& filePath= "");

    Logger(Logger const&) = delete;
    void operator=(Logger const&) = delete;
    ~Logger() = default;

    void complete_progres();

    void print_newline();
};
