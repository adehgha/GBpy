"""
This function displays a Message (passed as input) and gives and error
in case the matrix passed to it is not integral.`
"""
# Authors: Arash Dehghan Banadaki <adehgha@ncsu.edu>, Srikanth Patala <spatala@ncsu.edu>
# Copyright (c) 2014,  Arash Dehghan Banadaki and Srikanth Patala.
# License: GNU-GPL Style.

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
