import numpy as np
import matplotlib.pyplot as plt
from Region1 import *
from Region2 import *


print(Bp1.shape)
print(Bp2.shape)

plt.figure()
plt.plot(q, Bp1)
plt.plot(q, Bp2, '--')

plt.figure()
plt.plot(q, Bq1)
plt.plot(q, Bq2, '--')




plt.show()