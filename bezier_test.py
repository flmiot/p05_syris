import quantities as q
from syris.geometry import BezierAnimationKey
import numpy as np
import matplotlib.pyplot as plt

key1 =  BezierAnimationKey(0 * q.deg, 1/25. * q.sec, [-0.3487448*q.sec, 0*q.deg], [0.4320781*q.sec, 0*q.deg])
key2 = BezierAnimationKey(180 * q.deg, 1 * q.sec, [0.5809233*q.sec, 211.0208*q.deg], [1.432078*q.sec, 180*q.deg])

count = 1500
times = np.linspace(1/25., 1.0, num = count) * q.sec
values = np.linspace(1/25., 1.0, num = count) * q.deg

values = [key1.interpolate(key2, t) for t in times]

plt.plot(times * 25, values)
plt.show()
