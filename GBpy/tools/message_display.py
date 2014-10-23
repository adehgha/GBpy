"""
This function displays a Message (passed as input) and gives and error
in case the matrix passed to it is not integral.`
"""

from integer_manipulations import int_check
from Col import Col


def message_display(CheckMatrix, Checknumber, Message, Precis):
    cond = int_check(CheckMatrix, Precis)
    print Checknumber, '.', Message, '-> ',
    txt = Col()
    if cond.all():
        txt.c_prnt('YES', 'yel')
    else:
        txt.c_prnt('<<<Error>>>', 'amber')
        raise Exception('Something wrong!!')