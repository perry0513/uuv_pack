from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from prism import Prism

class Plotter:
    def __init__(self):
        fig = plt.figure()
        self.ax = fig.add_subplot(projection='3d')
        self.alpha = 0.2


    def rect_prism(self, x, y, z, w, l, h, color='g'):
        vs = Prism.get_prism_corners(x, y, z, w, l, h)
        sides = [
                [vs[0], vs[1], vs[4], vs[2]],
                [vs[0], vs[1], vs[5], vs[3]],
                [vs[0], vs[2], vs[6], vs[3]],
                [vs[3], vs[5], vs[7], vs[6]],
                [vs[1], vs[4], vs[7], vs[5]],
                [vs[2], vs[4], vs[7], vs[6]]
                ]
        self.ax.add_collection3d(Poly3DCollection(sides, facecolors=color, linewidths=1, edgecolors=color, alpha=self.alpha))

    def pyramid(self, rect, peak, color='b'):
        sides = [
                [peak, rect[0], rect[1]],
                [peak, rect[1], rect[2]],
                [peak, rect[2], rect[3]],
                [peak, rect[3], rect[0]],
                [rect[0], rect[1], rect[2], rect[3]]
                ]
        self.ax.add_collection3d(Poly3DCollection(sides, facecolors=color, linewidths=1, edgecolors=color, alpha=self.alpha))

    def point(self, x, y, z):
        self.ax.scatter(x, y, z, c='r')


    def show(self):
        plt.show()
