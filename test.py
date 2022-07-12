import random
from statistics import covariance
import numpy as np

c = 8
r = 16


arr = []

for i in range(0,r):
    l = random.sample(range(1, c+1), c)
    arr.append(l)

arr = np.array(arr).transpose()

print(arr)

np.set_printoptions(floatmode='fixed', precision=1, suppress='false', linewidth=150)

print(np.cov(arr))




x_arr = [[0.00,	0.50,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.49,	0.00,	1.00],	
        [0.00,	0.74,	0.00,	0.49,	0.99,	0.49,	0.48,	0.49,	0.49,	0.73,	0.00,	1.00],	
        [0.92,	0.74,	0.92,	0.49,	0.99,	0.93,	0.72,	0.94,	0.49,	0.73,	0.00,	1.00],	
        [0.92,	0.96,	0.92,	0.49,	0.99,	0.93,	0.92,	0.94,	0.93,	0.95,	0.00,	1.00]]	


print(np.cov(x_arr))

print(np.cov(np.transpose(x_arr)))
