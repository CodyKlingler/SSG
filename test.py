
def lll(i: int) -> int:
	return pow(3*i*i,i-2)

for i in range(3,16):
    combos = lll(i)
    increase = combos / lll(i-1)
    print("{}\t{}\t{}".format(i, increase, combos))

