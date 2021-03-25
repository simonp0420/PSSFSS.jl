module Log

export pssfss_logger, ScreenLevel, FileLevel, @logfile, @logscreen

using Logging: LogLevel, global_logger, @logmsg
using LoggingExtras: FormatLogger, TeeLogger, EarlyFilteredLogger

const ScreenLevel = LogLevel(-500)
const FileLevel = LogLevel(-400)

macro logscreen(x)
    :(@logmsg(ScreenLevel, $(esc(x))))
end

macro logfile(x)
    :(@logmsg(FileLevel, $(esc(x))))
end

function newlevels_filter(log_args)
    log_args.level ∈ (FileLevel, ScreenLevel)
end

function notnew_filter(log_args)
    log_args.level ∉ (ScreenLevel,FileLevel)
end

function plain_logger(logfile)
    isfile(logfile) && rm(logfile) 
    FormatLogger() do io, args
        open(logfile,"a") do f
            println(f, args.message)
        end
        args.level == ScreenLevel && println(args.message)
    end
end

pssfss_logger(logfile) = 
    TeeLogger(
        EarlyFilteredLogger(notnew_filter, global_logger()), 
        EarlyFilteredLogger(newlevels_filter, plain_logger(logfile))
             )

end # module