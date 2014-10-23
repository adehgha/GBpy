class Col(object):
    """
    This class is defined to ouput a word or sentence in a different color
    to the standard shell.
    The colors available are:
    pink
    blue
    green
    dgrn: dark green
    yellow
    amber
    """
    def __init__(self):
        self.pink = '\033[95m'
        self.blue = '\033[94m'
        self.green = '\033[92m'
        self.dgrn = '\033[1;32m'
        self.yellow = '\033[93m'
        self.amber = '\033[91m'
        self.ENDC = '\033[0m'

    def c_prnt(self, text, color):
        if color is 'pink':
            a = self.pink
        elif color is 'blue':
            a = self.blue
        elif color is 'green':
            a = self.green
        elif color is 'dgrn':
            a = self.dgrn
        elif color is 'yel':
            a = self.yellow
        elif color is 'amber':
            a = self.amber
        else:
            raise Exception('The color you selected is not acceptable')
        print a + text + self.ENDC
