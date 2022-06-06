import pandas as pd
import sys
sys.path.append("./src/")
from ct_analizer import get_row


get_row("   421 dG =   -150.99 NC_029256.1|-|1012891-1013311|201-221", "./example.ct", 0, energy_calc_method="")