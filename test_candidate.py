import functools
import numpy

def create_vector(x,y,z):
    return (x,y,z)

def scalar_product(v1,v2):
    return sum([a*b for a,b in zip(v1,v2)])


new_vector = functools.partial(create_vector, z=0)

v1 = new_vector(1,3)
v2 = new_vector(4,-2)

print(scalar_product(v1,v2))


def sum_up_to(n):
    return numpy.arange(1+n).sum()

print("sum up to 3: ",sum_up_to(3))