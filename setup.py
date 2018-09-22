import sys, os, platform
from distutils.core import setup

"""
setup script for FLORA -- Fast Long-noncoding RNA Assembly Workflow
"""

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
    print >> sys.stderr, "ERROR: FLORA requires Python 2.7"

def main():
    setup(
        name='FLORA',
        version='1.0',
        description='FLORA (Fast Long-noncoding RNA Assembly Workflow)',
        packages = ['lncmodule'],
        package_dir = {'': 'lib'},
        scripts = ['bin/filterTranscripts.py', 'bin/generateFilteredBams.py', 'bin/annotateTranscripts.py'],
        author = "Alex Shi",
        platforms = ['Linux','MacOS'],
        long_description = "FLORA is a fast workflow for assembling lncRNA transcriptome",
    )

if __name__ == "__main__":
    main()