#include "log_manager.h"
#include "config.h"

using namespace std;

log_manager& log_manager::instance()
{
    static log_manager lm;
    return lm;
}

void log_manager::initialize(filesystem::path const& log_file)
{
    log_manager& lm = instance();
    auto stdout_sink = make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto file_sink = make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.string());
    std::vector<spdlog::sink_ptr> sinks = {file_sink};
#   if _DEBUG
    sinks.push_back(stdout_sink);
#   endif
    // todo add verbose condition
    if (std::find(sinks.begin(), sinks.end(), stdout_sink) == sinks.end())
        sinks.push_back(stdout_sink);

    lm.main_logger = std::make_shared<spdlog::logger>(cfg::log::main_logger_id, begin(sinks), end(sinks));
#   if _TEST
    lm.main_logger->set_level(spdlog::level::off);
#   elif _DEBUG
    lm.main_logger->set_level(spdlog::level::debug);
#   else
    lm.main_logger->set_level(spdlog::level::info);
#   endif
    // TODO call set_pattern correctly
}