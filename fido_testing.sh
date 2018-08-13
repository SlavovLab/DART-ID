#!/usr/bin/env python3

import sys

from dart_id.fido.BayesianNetwork import real_main, pd_main
from dart_id.fido_cpp_interface import main

#real_main(sys.argv[1:])
pd_main()

#main()
