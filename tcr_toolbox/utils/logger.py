import contextlib
import io
import logging
import logging.handlers
import os


def init_logger(log_path_and_file="log_file_UNKNOWN_script.log", level_msg="DEBUG"):
    should_roll_over = os.path.isfile(log_path_and_file)

    if level_msg == "DEBUG":
        level_msg = logging.DEBUG
    elif level_msg == "INFO":
        level_msg = logging.INFO
    elif level_msg == "WARNING":
        level_msg = logging.WARNING

    logging.basicConfig(filename=log_path_and_file, level=level_msg, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%m/%d/%Y %I:%M:%S", filemode="w")

    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "%m/%d/%Y %I:%M:%S")
    handler = logging.handlers.RotatingFileHandler(log_path_and_file, mode="w", backupCount=0)  # ensure older backups get deleted

    handler.setFormatter(formatter)

    log = logging.getLogger(log_path_and_file)
    log.addHandler(handler)

    if should_roll_over:  # log already exists, roll over!
        handler.doRollover()

    return log


def log_function_to_file(func, print_log_file_path, *args, **kwargs):
    # Create a string IO stream to capture the prints
    string_io = io.StringIO()
    # Redirect stdout to the string IO stream
    with contextlib.redirect_stdout(string_io):
        results = func(*args, **kwargs)
    # Get the contents of the string IO stream
    output = string_io.getvalue()
    # Write the output to a file
    with open(print_log_file_path, "w") as file:
        file.write(output)

    return results


def stream_log_function_to_file(func, print_log_file_path, *args, **kwargs):
    """Run a function while streaming stdout directly to a log file (live)."""
    with open(print_log_file_path, "w", buffering=1) as log_file:  # line-buffered
        with contextlib.redirect_stdout(log_file):  # only redirect stdout
            results = func(*args, **kwargs)
    return results
