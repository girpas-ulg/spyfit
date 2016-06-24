"""
Processing base.
"""


class ProcessingBase(object):
    """Abstract class for retrieval processing.

    """
    def preprocess(self):
        raise NotImplementedError

    def postprocess(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
