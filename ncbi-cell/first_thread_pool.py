from multiprocessing.pool import ThreadPool
import subprocess
from time import sleep

def f(x):
    sleep(1)
    return subprocess.run(["date"], capture_output=True)

if __name__ == '__main__':
    with ThreadPool(3) as p:
        result = p.map(f, [1, 2, 3, 4, 5, 6])
