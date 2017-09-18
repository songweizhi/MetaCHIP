import matplotlib.pyplot
import pylab

x = [1,2,3,4]
y = [3,4,8,6]

m = [10,20,30,40]
n = [30,40,80,60]

matplotlib.pyplot.scatter(x,y, c='r')
matplotlib.pyplot.scatter(m,n, c='b')

matplotlib.pyplot.show()