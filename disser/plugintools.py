class PluginRegistry(object):
    def __init__(self):
        self._plugins = []

    def __iter__(self):
        return iter(self._plugins)

    def register(self, plugin):
        self._plugins.append(plugin)

    def unregister(self, plugin):
        try:
            self._plugins.remove(plugin)
        except ValueError:
            pass
