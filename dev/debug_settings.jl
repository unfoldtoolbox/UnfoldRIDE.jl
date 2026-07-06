using Logging
using LoggingExtras

# Define a filter function for the logger
function module_filter(args)
    # Check if the log message originates from this module
    r = args._module === UnfoldRIDE
    return r #log_args.meta[:module] == @__MODULE__
end

# Create a filtered logger
filtered_logger = EarlyFilteredLogger(module_filter, ConsoleLogger(stderr, Logging.Debug))

# Set the filtered logger as the global logger
global_logger(filtered_logger)

# Reset the global logger to the default ConsoleLogger
global_logger(ConsoleLogger(stderr))