# Copyright (C) 2023 University of Texas at El Paso
#
# Contributed by: Christoph Lauter
#
#                 and the 2023 class of CS4390/5390
#
#                 Applied Numerical Computing for Multimedia
#                 Applications.
#
# All rights reserved.
#
# NO LICENSE SPECIFIED.
#

class Utep_integer_implementation:
    def __init__(self, x):
        if isinstance(x, Utep_integer_implementation):
            self.__value = x.__value
        elif isinstance(x, str):
            self.__value = int(x)
        elif isinstance(x, int):
            self.__value = x
        else:
            raise Exception("Cannot convert {} to a UTEP integer".format(x))

    def convert_to_int(self):
        return self.__value

    def convert_to_string(self):
        return str(self.__value)
        
    def addition(self, other : Utep_integer_implementation):
        return Utep_integer_implementation(self.__value + other.__value)

    def subtraction(self, other : Utep_integer_implementation):
        return Utep_integer_implementation(self.__value - other.__value)

    def multiplication(self, other : Utep_integer_implementation):
        return Utep_integer_implementation(self.__value * other.__value)

    def division(self, other : Utep_integer_implementation):
        return Utep_integer_implementation(self.__value // other.__value)

    def modulo(self, other : Utep_integer_implementation):
        return Utep_integer_implementation(self.__value % other.__value)

    def shift_left(self, other : Utep_integer_implementation):
        return Utep_integer_implementation(self.__value << other.__value)

    def shift_right(self, other : Utep_integer_implementation):
        return Utep_integer_implementation(self.__value >> other.__value)

    def negate(self):
        return Utep_integer_implementation(-self.__value)

    def compare(self, other : Utep_integer_implementation):
        if self.__value < other.__value:
            return -1
        elif self.__value == other.__value:
            return 0
        else:
            return 1

        
class Utep_integer:
    """Python binding for the UTEP C integer functions"""

    def __wrap_value_from_string(self, x : str):
        return Utep_integer_implementation(x)

    def __wrap_value_from_int(self, x : int):
        return Utep_integer_implementation(x)

    def __wrap_value_to_string(self, x : Utep_integer_implementation):
        return x.convert_to_string()

    def __wrap_value_to_int(self, x : Utep_integer_implementation):
        return x.convert_to_int()
    
    def __wrap_value_add(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.addition(y)

    def __wrap_value_sub(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.subtraction(y)

    def __wrap_value_mul(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.multiplication(y)

    def __wrap_value_div(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.division(y)

    def __wrap_value_mod(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.modulo(y)

    def __wrap_value_lshift(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.shift_left(y)

    def __wrap_value_rshift(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.shift_right(y)
    
    def __wrap_value_cmp(self, x : Utep_integer_implementation, y : Utep_integer_implementation):
        return x.compare(y)
    
    def __wrap_value_neg(self, x : Utep_integer_implementation):
        return x.negate()
    
    def __init__(self):
        self.__value = self.__wrap_value_from_int(0)

    def __init__(self, x):
        if isinstance(x, Utep_integer):
            self.__value = x.__value
        elif isinstance(x, str):
            self.__value = self.__wrap_value_from_string(x)
        elif isinstance(x, int):
            self.__value = self.__wrap_value_from_int(x)
        else:
            raise Exception("Cannot convert {} to a UTEP integer".format(x))
        
    def __str__(self):
        return self.__wrap_value_to_string(self.__value)

    def __repr__(self):
        return "utep_integers.Utep_integer(\"{}\")".format(self.__wrap_value_to_string(self.__value))

    def __int__(self):
        return self.__wrap_value_to_int(self.__value)

    def __index__(self):
        return int(self)
    
    def __add__(self, x):
        return Utep_integer(self.__wrap_value_add(self.__value, Utep_integer(x).__value))

    def __radd__(self, x):
        return Utep_integer(self.__wrap_value_add(Utep_integer(x).__value, self.__value))
    
    def __sub__(self, x):
        return Utep_integer(self.__wrap_value_sub(self.__value, Utep_integer(x).__value))

    def __rsub__(self, x):
        return Utep_integer(self.__wrap_value_sub(Utep_integer(x).__value, self.__value))
    
    def __mul__(self, x):
        return Utep_integer(self.__wrap_value_mul(self.__value, Utep_integer(x).__value))

    def __rmul__(self, x):
        return Utep_integer(self.__wrap_value_mul(Utep_integer(x).__value, self.__value))
    
    def __floordiv__(self, x):
        return Utep_integer(self.__wrap_value_div(self.__value, Utep_integer(x).__value))

    def __rfloordiv__(self, x):
        return Utep_integer(self.__wrap_value_div(Utep_integer(x).__value, self.__value))

    def __mod__(self, x):
        return Utep_integer(self.__wrap_value_mod(self.__value, Utep_integer(x).__value))

    def __rmod__(self, x):
        return Utep_integer(self.__wrap_value_mod(Utep_integer(x).__value, self.__value))

    def __pow__(self, _y):
        y = Utep_integer(_y)
        if y < 0:
            return Utep_integer(1) / (self ** (-y))
        if y == 0:
            return Utep_integer(1)
        if y == 1:
            return self
        if y % 2 == 0:
            return (self * self) ** (y // 2)
        return self * ((self * self) ** ((y - 1) // 2))

    def __rmod__(self, x):
        return Utep_integer(x) ** self

    def __neg__(self):
        return Utep_integer(self.__wrap_value_neg(self.__value))

    def __pos__(self):
        return Utep_integer(self.__value)

    def __abs__(self):
        if self >= 0:
            return self
        else:
            return -self
    
    def __lt__(self, x):
        return (self.__wrap_value_cmp(self.__value, Utep_integer(x).__value) < 0)

    def __le__(self, x):
        return (self.__wrap_value_cmp(self.__value, Utep_integer(x).__value) <= 0)

    def __eq__(self, x):
        return (self.__wrap_value_cmp(self.__value, Utep_integer(x).__value) == 0)
    
    def __ne__(self, x):
        return (self.__wrap_value_cmp(self.__value, Utep_integer(x).__value) != 0)

    def __ge__(self, x):
        return (self.__wrap_value_cmp(self.__value, Utep_integer(x).__value) >= 0)

    def __gt__(self, x):
        return (self.__wrap_value_cmp(self.__value, Utep_integer(x).__value) > 0)

    def __lshift__(self, x):
        return Utep_integer(self.__wrap_value_lshift(self.__value, Utep_integer(x).__value))

    def __rlshift__(self, x):
        return Utep_integer(self.__wrap_value_lshift(Utep_integer(x).__value, self.__value))

    def __rshift__(self, x):
        return Utep_integer(self.__wrap_value_rshift(self.__value, Utep_integer(x).__value))

    def __rrshift__(self, x):
        return Utep_integer(self.__wrap_value_rshift(Utep_integer(x).__value, self.__value))

