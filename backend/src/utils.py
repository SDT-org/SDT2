import psutil
from hashlib import sha256


def get_child_process_info():
    process = psutil.Process()
    info = []

    for child in process.children(recursive=True):
        try:
            cpu_percent = child.cpu_percent(interval=0.1)
            memory = child.memory_info().rss
            info.append((child.pid, cpu_percent, memory, child.nice()))
        except:
            pass

    return info

def make_doc_id(input: str) -> str:
    return sha256(input.encode()).hexdigest()
