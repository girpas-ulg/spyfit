"""
Default processing classes for SFIT4.
"""

import subprocess

from .. import config


__all__ = ['SFIT4ProcessingBase',]


class SFIT4ProcessingBase(ProcessingBase):
    """
    A simple SFIT4 run with default pre- and postprcessing
    (Layer0 and Layer1).
    """
    def __init__(run_dir, bnr_spec, ctl="sfit4.ctl", hbin_input="hbin.input"):
        os.chdir(run_dir)
        self.path2exe = config.sfit4_exe_path

        def generate_spectrum(self):
            pspec_exe = os.path.join(self.path2exe, "pspec")
            subprocess.Popen(pspec_exe)

        def generate_reference_profiles(self):
            pass

        def generate_hbin(self):
            hbin_exe = os.path.join(self.path2exe, "hbin")
            subprocess.Popen(hbin_exe)

        def preprocess(self):
            pass

        def run(self):
            preprocess()

            sfit4_exe = os.path.join(self.path2exe, "sfit4")
            subprocess.Popen(sfit4_exe)

            postprocess()

        def postprocess(self):
            pass
