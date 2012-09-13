from libc.math cimport sin, cos, atan2, acos, sqrt, M_1_PI
cdef class fib_param2(object):
    cdef float a,b,c,al,be
    def __init__(self, base_point, pm=1):
        self.a=base_point[0]
        self.b=base_point[1]
        self.c=base_point[2]
        self.al=sqrt(.5*(1 + self.c))
        self.be=sqrt(.5*(1 - self.c))
        
    def __call__(self, float ph):
        cdef float th = atan2(-self.a,self.b)-ph
        cdef float w,x,y,z
        (w,x,y,z) = (self.al*cos(th), -1*self.be*cos(ph), -1*self.be*sin(ph), self.al*sin(th))
        cdef float rr = acos(w)*M_1_PI*(1.0/sqrt(1-w**2))
        return (x*rr, y*rr, z*rr)
