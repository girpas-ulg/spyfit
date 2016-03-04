"""
Processing base.
"""


class ProcessingBase(object):
    """Abstract class for retrieval processing.

    """
    def preprocess():
        raise NotImplementedError

    def postprocess():
        raise NotImplementedError

    def run():
        raise NotImplementedError
