#!/usr/local/bin/python3

import sys,random

choices = '123456789'
nchoices = len(choices)

oddonly = '13579'
noddonly = len(oddonly)

def make_up_password(pwdlen):
    p = ''
    for i in range(pwdlen-1):
        p += choices[random.randint(0,nchoices-1)]
        
    # For superstitious reasons, only generate odd-numbered seeds
    p += oddonly[random.randint(0,noddonly-1)]
    
    return p

if __name__ == '__main__':
    passwd_length = 7
    if len(sys.argv) > 1:
        passwd_length = int(sys.argv[1])
    assert passwd_length > 0
    passwd = make_up_password(passwd_length)
    print(passwd)

