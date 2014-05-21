import sys
from mock import Mock
mock = Mock()

modules = {}

try:
    import h5py
except ImportError:
    modules.update({'h5py':mock})

try:
    from tvtk.api import tvtk
    import mayavi
    from mayavi import mlab
except ImportError:
    modules.update({'tvtk':mock, 'tvtk.api': mock.module, 'traits':mock, 'traits.api':mock.module,
           'mayavi':mock, 'mayavi.tools':mock.module, 'mayavi.tools.sources':mock.module,
           'mayavi.modules':mock.module, 'mayavi.modules.streamline':mock.module})

try:
    import yt.mods
except ImportError:
    modules.update({'yt':mock, 'yt.mods':mock.module})

sys.modules.update(modules)

# Load all of the global Astropy configuration
from astropy.sphinx.conf import *

extensions += ['sphinx.ext.autosummary', 'sphinx.ext.mathjax']
numpydoc_show_class_members = False
# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'pySAC'
copyright = u'2013, Stuart Mumford'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '0.0.1a'
# The full version, including alpha/beta/rc tags.
release = '0.0.1a'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'
html_favicon = ''
