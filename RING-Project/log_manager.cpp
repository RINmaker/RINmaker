#include "log_manager.h"
#include "config.h"

using namespace std;

/*vector<string> log_manager::get_ids()
{
    log_manager& lm = instance();
    vector<string> keys;
    keys.reserve(lm.loggers.size());

    for (auto kv : lm.loggers)
    {
        keys.push_back(kv.first);
    }

    return keys;
}

shared_ptr<spdlog::logger> log_manager::get(string const& id)
{
    log_manager& lm = instance();
    if (lm.loggers.find(id) == lm.loggers.end())
    {
        lm.cerr_logger->error("you tried to retrieve logger: \"" + id + "\" but it does not exist.\nyou are now using default console logger.");
        return lm.cout_logger;
    }

    else
    {
        return lm.loggers[id];
    }
}

shared_ptr<spdlog::logger> log_manager::get_default()
{
    return instance().cout_logger;
}

void log_manager::insert(string const& id, string const& sink_name)
{
    log_manager& lm = instance();
    try
    {
        if (lm.loggers.find(id) != lm.loggers.end())
        {
            lm.cerr_logger->error("\"" + id + "\" specifies an already existing logger.");
        }

        else
        {
            auto basic_sink = make_shared<spdlog::sinks::basic_file_sink_mt>((lm.log_directory / sink_name).string());
            auto basic_logger = make_shared<spdlog::logger>(id, basic_sink);
            lm.loggers.insert({ id, basic_logger });
        }
    }

    catch (spdlog::spdlog_ex const& e)
    {
        // TODO valutare una funzione generale di exception throwing che formatti i messaggi sempre allo stesso modo
        lm.cerr_logger->error("a spdlog exception occurred while inserting \"" + id + "\" using sink \"" + sink_name + "\"\n" + e.what());
    }
}*/

void log_manager::initialize(filesystem::path const& log_directory)
{
    log_manager& lm = instance();
    auto stdout_sink = make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto file_sink = make_shared<spdlog::sinks::basic_file_sink_mt>((log_directory / cfg::log::file_logger_filename).string());
    lm.console_logger = make_shared<spdlog::logger>(cfg::log::console_logger_id, stdout_sink);
    lm.file_logger = make_shared<spdlog::logger>(cfg::log::file_logger_id, file_sink);
    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(file_sink);
#   if _DEBUG
    sinks.push_back(stdout_sink);
#   endif
    lm.main_logger = std::make_shared<spdlog::logger>(cfg::log::main_logger_id, begin(sinks), end(sinks));
#   if _TEST
    lm.main_logger->set_level(spdlog::level::off);
#   elif _DEBUG
    lm.main_logger->set_level(spdlog::level::debug);
#   else
    lm.main_logger->set_level(spdlog::level::info);
#   endif
    // TODO chiamare set_pattern opportunamente
}