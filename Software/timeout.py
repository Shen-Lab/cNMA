# taken from http://stackoverflow.com/questions/2281850/timeout-function-if-it-takes-too-long-to-finish

from functools import wraps
import errno
import os
import signal

# number of seconds to wait before quitting a function
RMSDtimeout = 120

class TimeoutError(Exception):
    pass

def timeout(seconds=RMSDtimeout, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)
        
        @wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator 
