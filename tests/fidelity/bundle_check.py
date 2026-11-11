"""Kept so existing job scripts work: the bundle check now lives in hamrvis/verify (hamrvis verify <bundle>)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from hamrvis.verify.check import *                       # noqa: E402,F401,F403
from hamrvis.verify.check import main                    # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
