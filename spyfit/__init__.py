# -*- coding: utf-8 -*-

from .retrieval import load_dataset

try:
    from .version import version as __version__
except ImportError:  # pragma: no cover
    raise ImportError('spyfit not properly installed. If you are running from '
                      'the source directory, please instead create a new '
                      'virtual environment (using conda or virtualenv) and '
                      'then install it in-place by running: pip install -e .')
