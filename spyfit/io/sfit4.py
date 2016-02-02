# -*- coding: utf-8 -*-

"""Read (write) SFIT4 input and output ascii files."""

from __future__ import absolute_import, print_function

import datetime
import re
import os
import math

import numpy as np

from . import _sfit4out as s4out
from . import _sfit4in as s4in
