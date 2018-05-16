import math,sys
unit_list = zip(['bytes', 'kB', 'MB', 'GB', 'TB', 'PB'], [0, 0, 1, 2, 2, 2])
def sizeof_fmt(num):
    """Human friendly file size by Fred Cirera"""
    if num > 1:
        exponent = min(int(math.log(float(num), 1024)), len(unit_list) - 1)
        print exponent
        quotient = float(num) / 1024**exponent
        print quotient
        unit, num_decimals = unit_list[exponent]
        print unit,num_decimals
        format_string = ('%%0.%sf%%s' % (num_decimals)) % (quotient,unit)
        print format_string
        return format_string
    if num == 0:
        return '0 bytes'
    if num == 1:
        return '1 byte'

if __name__ == '__main__':
   bytes = sys.argv[1]
   
   print sizeof_fmt(bytes)


