#include "Logger.h"

#include "Logger.h"

#include <string>
#include <mutex>


void Logger::init_logger(LoggerLevel minLevel, const std::string &filePath) {
    this->minSeverityLevel = minLevel;
}

void Logger::show_progress(int32_t percentage) {

    int32_t progress = percentage / barUnit;   // Calculate how many "=" should be shown (5% per "=")

    if (percentage < 0 || percentage > 100) return;
    auto now = std::chrono::system_clock::now();
    std::time_t timeNow = std::chrono::system_clock::to_time_t(now);
    //std::cout << std::put_time(localtime(&timeNow), "\r[%Y-%m-%d %H:%M:%S] ") << SeverityTable(LOG_INFO) << "Progress: " << percentage << "%     " << std::flush;
    std::cout << std::put_time(localtime(&timeNow), "\r[%Y-%m-%d %H:%M:%S] ") << SeverityTable(LOG_INFO) << "Progress: [";

    // Print the progress bar
    for (int i = 0; i < barWidth; ++i) {
        if (i < progress) {
            std::cout << "=";  // Filled part
        } else {
            std::cout << " ";  // Empty part
        }
    }

    // Close the progress bar and show the percentage
    std::cout << "] " << percentage << "%     " << std::flush;
}

void Logger::complete_progres() {
    auto now = std::chrono::system_clock::now();
    std::time_t timeNow = std::chrono::system_clock::to_time_t(now);
    std::cout << std::put_time(localtime(&timeNow), "\r[%Y-%m-%d %H:%M:%S] ") << SeverityTable(LOG_INFO) << "Progress: [";

    // Print the progress bar
    for (int i = 0; i < barWidth; ++i) {
        std::cout << "=";// Filled part
    }

    // Close the progress bar and show the percentage
    std::cout << "] " << 100 << "%     " << std::endl;}

void Logger::print_newline() {
    std::cout << std::endl;
}
