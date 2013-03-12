rp1_arz_traffic.so: rp1_arz_traffic.f90
	f2py -m $(basename $(notdir $@)) -c $^

rp1_hll.so: rp1_hll.f90
	f2py -m $(basename $(notdir $@)) -c $^
