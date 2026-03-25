import logging
import sys


def section(self, title, char='=', width=10):
    markers = char * width
    msg = f'{markers}{title}{markers}'
    self.info(msg)


def setup_logging(level=logging.INFO, logfile=None):
    fmt = logging.Formatter("%(asctime)s|%(levelname)s| %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    root.addHandler(ch)

    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        root.addHandler(fh)

    # silence noisy libs
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    
    # add section header
    logging.Logger.section = section
