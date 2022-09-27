import signal
from t import f

class TimeOut(Exception): pass

# Register an handler for the timeout
def handler(signum, frame):
    print("Forever is over!")
    raise TimeOut()

# Register the signal function handler
signal.signal(signal.SIGALRM, handler)

# Define a timeout for your function
signal.alarm(10)

f()
