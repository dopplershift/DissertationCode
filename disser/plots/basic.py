from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText


class Plot(object):
    pass

# A way to automatically generate labels for subplots through iteration
class LabelGenerator(object):
    def __init__(self, first, loc=2):
        self.first = first
        self.loc = loc

    def __iter__(self):
        item = self.first
        while True:
            # loc = 2 is upper left. TODO: Should patch matplotlib to use
            # same strings as legend
            at = AnchoredText("%s)" % item, loc=self.loc,
                prop=dict(size='large'), frameon=True)

            newitem = yield at
            if newitem:
                item = newitem
            else:
                item = self.next_item(item)

    def next_item(self, current):
        return chr(ord(current) + 1)
