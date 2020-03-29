import sys
import inspect

FUNC_START = lambda: print(f"**** '{inspect.currentframe().f_back.f_code.co_name}': Entering. ****")
FUNC_END = lambda: print(f"**** '{inspect.currentframe().f_back.f_code.co_name}': Leaving successfully. ****")

