#######################################################################
##### timer function ##################################################
#######################################################################

import time
# function timer, you need to introduce start and end time, as you can calculate as:
# start = time.time()
# end = time.time()
def timer(start, end):
	hours, rem = divmod(end-start, 3600)
	minutes, seconds = divmod(rem, 60)
	print("Time elapsed {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
