import time

START_TIME = time.monotonic()
LOG_TOPICS = set()
VERSION = '3.1'

def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def log(log_message, topic=''):
    if not topic or topic in LOG_TOPICS:
        print(f'{format_runtime()} {log_message}')
