"""
Shared logging setup helpers.

Adds a custom ``Logger.section`` helper for readable stage separators and
configures console/file handlers.
"""

import logging
import sys


def _section(self: logging.Logger, title: str, char: str = "=", width: int = 10) -> None:
    """
    Log a visual section marker using ``Logger.info``.

    Args:
        self (logging.Logger): Active logger instance.
        title (str):           Section title text.
        char (str):            Repeated character used for framing.
        width (int):           Number of framing characters on each side.
    """
    markers = char * width
    msg = f'{markers}{title}{markers}'
    self.info(msg)


def setup_logging(level: int = logging.INFO, logfile: str | None = None) -> None:
    """
    Configure root logging handlers used across pipeline modules.

    Args:
        level (int):          Console log level.
        logfile (str | None): Optional path for detailed file logging.
    """
    fmt = logging.Formatter("%(asctime)s| %(message)s", datefmt="%H:%M:%S")
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
    logging.Logger.section = _section
