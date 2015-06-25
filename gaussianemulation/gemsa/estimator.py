from multiprocessing import Process, Queue, log_to_stderr
from gaussianemulation.gemsa.gememu import *
import logging

formatter = logging.Formatter('%(asctime)s - [%(levelname)s/%(processName)s] %(message)s')
logger = log_to_stderr()
logger.handlers[0].setFormatter(formatter)
logger.setLevel(logging.INFO)

def predict(qIn, qOut, path):
    init(path)
    logger.info("Initted")
    while True:
        x = qIn.get()
        logger.debug("Predicting:")
        logger.debug(x)
        if x is None:
            break
        qOut.put(nextpoint(x))
    qOut.put(None)
    qOut.close()

class GEMSAEmulator(object):
    """
    A gem-sa emulator, loaded from path.
    """
    def __init__(self, path):
        self.path = path
        self.qIn = Queue()
        self.qOut = Queue()
        self.estimator = Process(target=predict, args=(self.qIn, self.qOut, self.path))
        self.estimator.start()

    def predict(self, point):
        self.qIn.put(point)
        return self.qOut.get()

    def __del__(self):
        self.qIn.put(None)
        self.qOut.get()
        self.qIn.close()
        self.estimator.join()
