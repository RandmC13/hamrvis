"""The independent reader of VAMPIRE's files now lives in hamrvis/verify/reference.py."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from hamrvis.verify.reference import *                   # noqa: E402,F401,F403
