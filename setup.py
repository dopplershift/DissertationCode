from distutils.core import setup
from disser import __version__

setup(name = 'disser',
      version = str(__version__),
      packages = ['disser', 'disser.plots', 'disser.atten'],
      author = 'Ryan May',
      author_email = 'rmay31@gmail.com',
      platforms = ['Linux', 'UNIX', 'Windows', 'MacOSX'],
      description = 'Code to generate plots for my dissertation',
      url = 'http://weather.ou.edu/~rmay/research.html',
      )

