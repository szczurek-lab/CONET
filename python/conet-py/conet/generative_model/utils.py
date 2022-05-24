import logging
import sys
from typing import TypeVar, Callable

T = TypeVar('T')


def sample_conditionally(sampler: Callable[[], T], condition: Callable[[T], bool]) -> T:
    """
        Returns first result of call to @sampler which satisfies @condition
    """
    sample = sampler()
    while not condition(sample):
        sample = sampler()
    return sample


FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s — %(message)s")


def get_logger(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(FORMATTER)
    logger.addHandler(handler)
    logger.propagate = False
    return logger
