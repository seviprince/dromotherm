import matplotlib.pyplot as plt
"""
plt.figure(1)
plt.subplot(1, 2, 1)
plt.scatter(range(5), [x ** 2 for x in range(5)], color = 'blue')
plt.subplot(2, 2, 2)
plt.plot(range(5), color = 'red')
plt.subplot(2, 2, 4)
plt.bar(range(5), range(5), color = 'green')
plt.show()
"""

plt.subplot(121)
plt.scatter(range(5), [x ** 2 for x in range(5)], color = 'blue')
plt.subplot(222)
plt.plot(range(5), color = 'red')
plt.subplot(224)
plt.bar(range(5), range(5), color = 'green')
plt.show()