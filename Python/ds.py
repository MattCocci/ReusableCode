import numpy as np
import matplotlib.pyplot as plt


class DS:
  """
  Dynamical System Class:

  Inputs
  ======
  - Law of Motion (lom): a Python function; it maps
    objects in the ds
  - Starting Value (x)

  Methods
  =======
  - trajectory: Generates a tranjectory
  - plot_ts: Plot a trajectory time series
  - plot_hist: Plot a histogram of terminal values
  """

  def __init__(self, lom=None, x=None):
    self.lom = lom
    self.x   = x

  def update(self, T=1, vec=0):
    for i in range(T):
      self.x = self.lom(self.x)

  def trajectory(self, T)
    traj = []
    for i in range(T)
      traj.append(self.x)
      self.update()
    return traj

  def plot_ts(self, T, title="")
    traj = self.trajectory(T)
    plt.plot(range(T),traj,'b-',linewidth=2,alpha=0.5)
    plt.title(title)
    plt.show()

  def plot_hist(self, T, N, title="")




