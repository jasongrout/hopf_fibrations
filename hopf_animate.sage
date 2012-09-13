"""
Source code for drawing and animating fibers of the Hopf fibration.  See

http://www.nilesjohnson.net/hopf.html for more information.

#*****************************************************************************
#        Copyright (C) 2011 Niles Johnson <http://www.nilesjohnson.net>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

"""

"""
Benchmark Timings:

Low quality::

    basespiral = BaseSpiral(1,(0,pi/2),(pi,0),'sph', include_final=True)
    draw_frame(basespiral.points()[5:8], n=0, high_quality=False, only_show=True)

High quality::

    basespiral = BaseSpiral(1,(0,pi/2),(pi,0),'sph', include_final=True)
    draw_frame(basespiral.points()[5:8], n=0, high_quality=True, only_show=True)


36a648e09
---------

Low: .63s
High: 20.44s

a8bcb7d
-------
Low: 0.12s
High: 2.65s



"""

def draw_frame(f,n,dir="./frames/", resolution=304, high_quality=True, only_show=False):
    """
    Draws one frame of the animation.  Inputs are a list of points on
    S^2, a frame number, and keywords.

    EXAMPLES::

        sage: a = BaseSpiral(1,(0,pi/2),(pi,0),'sph', include_final=True)
        sage: draw_frame(a.points(), n=0, high_quality=False, only_show=True)

    """
    verbose('drawing frame %08d'%n, 2)
    g = HopfFibers(resolution=resolution, high_quality=high_quality)
    g.add_plane()
    g.add_base_sphere()
    g.add_base_axes()
    g.add_fibers(f)
    if only_show:
        g.show()
    else:
        out_file = dir+'hopf-frame%08d'%n+'.png'
        g.save(filename=out_file)
    return None

def fib_param2(base_point,pm=1):
    """
    Fiber parametrization function.  Returns a parametric curve for
    the fiber over base point (a,b,c).
    """
    aa,bb,cc = base_point
    if cc^2 >= 1:
        raise ValueError("c = "+str(cc)+" must have absolute value less than 1")
    verbose('c^2 = '+str(cc^2), level=3)
    pp(ph,a,b,c,al,be,w,x,y,z,r,i,th) = i*r/sqrt(1-w^2)
    p_condensed_x(ph)=(pp.subs(i=x).subs(r=arccos(w)/pi)
                        .subs(w=al*cos(th),x=-1*be*cos(ph),y=-1*be*sin(ph),z=al*sin(th)) 
                        .subs(th=atan2(-a,b)-ph).subs(al=sqrt((1+c)/2),be=sqrt((1 - c)/2))
                        .subs(a=aa,b=bb,c=cc))
    p_condensed_y(ph)=(pp.subs(i=y).subs(r=arccos(w)/pi)
                        .subs(w=al*cos(th),x=-1*be*cos(ph),y=-1*be*sin(ph),z=al*sin(th)) 
                        .subs(th=atan2(-a,b)-ph).subs(al=sqrt((1+c)/2),be=sqrt((1 - c)/2))
                        .subs(a=aa,b=bb,c=cc))
    p_condensed_z(ph)=(pp.subs(i=z).subs(r=arccos(w)/pi)
                        .subs(w=al*cos(th),x=-1*be*cos(ph),y=-1*be*sin(ph),z=al*sin(th)) 
                        .subs(th=atan2(-a,b)-ph).subs(al=sqrt((1+c)/2),be=sqrt((1 - c)/2))
                        .subs(a=aa,b=bb,c=cc))
    p_x,p_y,p_z = (fast_callable(p_condensed_x, vars=[ph], domain=RDF),
            fast_callable(p_condensed_y, vars=[ph], domain=RDF),
            fast_callable(p_condensed_z, vars=[ph], domain=RDF))
    def f(ph):
         return p_x(ph), p_y(ph), p_z(ph)
    return f

from sage.plot.plot3d.tachyon import Tachyon
from sage.misc.misc import verbose

class HopfFibers(Tachyon):
    """
    A class for modeling a view of the Hopf fibration.
    """
    def __init__(self, fiber_coordinates=None, cam_pos=(3,0,1), light_pos=(0,0,10),
                 look_at=(1,-.2,.1), high_quality=False, 
                 fiber_radius = .01,
                 antialiasing=True, resolution=608, **kwds):
        """
        Input is a list of tuples (a,b,c), which are 2-sphere coordinates
        over which to take fibers.
        """
        self.fiber_coordinates = fiber_coordinates
        self.fiber_radius = fiber_radius
        self.bounding_rings = False # use add_* methods to add bounding rings
        self.special_fibers = False # or special fibers
        self.base_sphere = False    # or base sphere
        self.bg_plane = False          # or bg plane
        self.cam_pos = cam_pos
        self.high_quality = high_quality
        Tachyon.__init__(self, xres=resolution, yres=resolution,
                         camera_center=cam_pos, look_at=look_at,
                         antialiasing=antialiasing, **kwds)

        self.light(light_pos, .1, (.8,.8,.8))

        ## default textures;
        ## use set_tex_* methods to modify
        self.tex_plane = {'name':'plane_tex',
                          'ambient':1.0, 'diffuse':0.8, 'specular':0,
                          'opacity':1.0, 'color':(1,1,1)}
        self.tex_bounding_rings = {'name':'ringtex',
                                   'ambient':1.0, 'diffuse':0, 'specular':0,
                                   'opacity':1.0, 'color':(0,0,0)}
        self.tex_base_sphere = {'name':'basetex',
                                'ambient':.7, 'diffuse':0.9, 'specular':0,
                                'opacity':0.45, 'color':(0.2,0.2,.2)}
        self.tex_base_axes = {'name':'axestex',
                              'ambient':0.2, 'diffuse':0.7, 'specular':0.5,
                              'opacity':0.7, 'color':(0,0,1.0)}
        self.tex_special_m1 = {'name':'special_m1_tex',
                               'ambient':1, 'diffuse':0.8, 'specular':0.13,
                               'opacity':1.0, 'color':Color(.7,.6,0)}
        self.tex_special_p1 = {'name':'special_p1_tex',
                               'ambient':1, 'diffuse':0.8, 'specular':0.13,
                               'opacity':1.0, 'color':Color(.3,.3,1)}
        self.tex_fiber = {'ambient':.8, 'diffuse':0.8, 'specular':0.13,
                          'opacity':1.0}



        self.two_sphere_ctr = (1.0,-1.0,-.7)
        self.two_sphere_dot_rad = .013
        self.two_sphere_radius = 1/3

    def set_tex_plane(self,tex_params=None):
        """
        Takes a dict of texture parameters and sets the texture for
        background plane.  Returns the current settings.
        """
        if tex_params is not None:
            self.tex_plane.update(tex_params)
        return self.tex_plane
    
    def add_plane(self,pos=[(-4,-4,-4),(1,1,1)]):
        self.bg_plane = True
        self.texture(**self.tex_plane)
        self.plane(pos[0],pos[1],'plane_tex')

    def set_tex_bounding_rings(self,tex_params=None):
        if tex_params is not None:
            self.tex_bounding_rings.update(tex_params)
        return self.tex_bounding_rings

    def add_bounding_rings(self,
                           big_ring_thickness = .006,
                           little_ring_thickness = .004):
        self.bounding_rings = True
        self.texture(**self.tex_bounding_rings)
        self.ring((0,0,0),(self.cam_pos[0],self.cam_pos[1],0),1,1+big_ring_thickness,'ringtex')

    def set_tex_base_sphere(self,tex_params=None):
        if tex_params is not None:
            self.tex_base_sphere.update(tex_params)
        return self.tex_base_sphere

    def add_base_sphere(self):
        """
        Note: somehow, the coordinate axes on the base 2-sphere wound
        up not having the standard orientation.  In the pictures, the
        y-axis points to the left, the x-axis points slightly to the
        right, and toward the camera, and the z-axis points up.
        Sorry.
        """
        self.base_sphere = True
        self.texture(**self.tex_base_sphere)
        self.sphere(self.two_sphere_ctr,radius=self.two_sphere_radius,texture='basetex')

    def set_tex_base_axes(self,tex_params=None):
        if tex_params is not None:
            self.tex_base_axes.update(tex_params)
        return self.tex_base_axes

    def add_base_axes(self):
        self.texture(**self.tex_base_axes)
        axis_length=.8*self.two_sphere_radius
        (c0,c1,c2) = self.two_sphere_ctr
        axes = [lambda t: (t+c0,c1,c2),
                lambda t: (c0,t+c1,c2),
                lambda t: (c0,c1,t+c2)]
        for f in axes:
            self.parametric_plot(f, 0, axis_length, 'axestex', r=self.two_sphere_dot_rad/2)

    def set_tex_special_fiber(self,tex_params=None,c=0):
        if tex_params is not None:
            if c == -1 or c == 0:
                self.tex_special_m1.update(tex_params)
            if c == 1 or c == 0:
                self.tex_special_p1.update(tex_params)
            return self.tex_special_m1, self.tex_special_p1
        
    def add_special_fiber(self,c=-1):
        if c == -1:
            self.texture(**self.tex_special_m1)
            self.parametric_plot(lambda t: (.5*cos(t), .5*sin(t), 0), 0, 2*pi,
                                 'special_m1_tex',
                                 r=self.fiber_radius,min_depth=5,max_depth=5)
            self.sphere((i+j for (i,j) in zip((0,0,-1/3),self.two_sphere_ctr)),
                        radius=self.two_sphere_dot_rad,texture='special_m1_tex')

        if c == 1:
            self.texture(**self.tex_special_p1)
            self.parametric_plot(lambda t: (0, 0, t), -1, 1,
                                 'special_p1_tex',
                                 r=self.fiber_radius,min_depth=5,max_depth=5)
            self.sphere((i+j for (i,j) in zip((0,0,1/3),self.two_sphere_ctr)),
                        radius=self.two_sphere_dot_rad,texture='special_p1_tex')

    def set_tex_fiber(self,tex_params=None):
        if tex_params is not None:
            self.tex_fiber.update(tex_params)
        return self.tex_fiber

    def set_depth(self,f=(1,0,0)):
        """
        Return values for min_depth and max_depth, depending on input f=(a,b,c).
        """
        if not self.high_quality:
            (min_depth,max_depth) = (4,5)
        else:
            if f[2] > .99:
                min_depth = 12
            if f[2] > .98:
                min_depth = 11
            elif f[2] > .8:
                min_depth = 10
            elif f[2] > .5:
                min_depth = 9
            elif f[2] > 0:
                min_depth = 8
            else:
                min_depth = 7
            max_depth = min_depth + 2
        return (min_depth,max_depth)

            
    def add_fibers(self, fiber_coordinates=None):
        if fiber_coordinates is None:
            if self.fiber_coordinates is None:
                raise ValueError('specify fiber coordinates')
            else:
                verbose('getting fiber_coordinates from stored value')
                fiber_coordinates = self.fiber_coordinates

        
        #PI = float(pi)
        sph_rad = 1
        
        two_sphere_ctr = (1.0,-1.0,-.7)
        two_sphere_dot_rad = .013

        # start here:

        RF = RealField(40)

        #num_fibs = len(fiber_coordinates)
        #verbose('looping over {0} fiber_coordinates..'.format(num_fibs),level=1)
        fib_counter = 1
        for (a,b,c) in fiber_coordinates:
            b = RF(b)
            c = RF(c)
            if b^2 + c^2 > sph_rad^2:
                if abs(sqrt(b^2 + c^2) - sph_rad) >.000001:
                    raise ValueError("input coordinates "+str((a,b,c))+" are not on 2-sphere")
                    
            if a >= 0:
                pm = 1
            else:
                pm = -1
            a = real(pm*sqrt(sph_rad^2-b^2-c^2))

            simpleRGB = Color((a+sph_rad)/2,(b+sph_rad)/2,(c+sph_rad)/2)
            fibtex_name = 'tb:{0}:{1}:{2}'.format(pm,c,b)
            fibtex = {'name':fibtex_name, 'color':simpleRGB}
            fibtex.update(self.tex_fiber)
            verbose('texture: '+str(fibtex), level=2)
            self.texture(**fibtex)
            
            verbose('  now working on '+str(tuple('%+.3f..'%x for x in (a,b,c))), level=3)
            if c > .997:
                # special fiber over north pole
                self.add_special_fiber(1)
            elif c == RR(-1):
                # special fiber over south pole
                self.add_special_fiber(-1)
            else:
                (mindp,maxdp) = self.set_depth((a,b,c))
                verbose('  computing parametric plot of fiber', level=2)
                verbose('  (b,c,pm) = '+str(tuple('%+.3f..'%x for x in (b,c,pm))), level=3)
                # parametric function for fiber:
                # (if a < 0, flip b)
                pf = fib_param2((a,b,c),pm)
                self.parametric_plot(pf,
                                     0,2*pi,
                                     tex=fibtex_name,
                                     r=self.fiber_radius,min_depth=mindp,max_depth=maxdp)

                verbose('  plotting point on base 2-sphere', level=2)
                self.sphere((i+j for (i,j) in zip((a/3,b/3,c/3),two_sphere_ctr)),radius=two_sphere_dot_rad,texture=fibtex_name)

            fib_counter += 1
            verbose('    over'+str(tuple('%+.3f..'%x for x in (a,b,c))), level=2)
        








#################
##
## Classes for drawing lists of points
## along various curves.
##
#################


class BaseList(SageObject):
    """
    Base class for a list of points on base 2-sphere.
    """
    def __init__(self,duration, init_pos=None, final_pos=None, input_type='sph', frame_rate=30, include_final=False, base_points=None, **kwds):
        self._duration = duration
        self._init_pos = init_pos
        self._final_pos = final_pos
        self._input_type = input_type
        self._frame_rate = frame_rate
        self._include_final = include_final
        self._rotation = None
        self._transformation = None
        self._repr_extra = ""
        self._other_kwds = kwds
        self.base_points = base_points

    def _repr_(self):
        return "List of points on base 2-sphere"+self._repr_extra+"."
    def duration(self):
        return self._duration
    def init_pos(self):
        return self._init_pos
    def final_pos(self):
        return self._final_pos
    def input_type(self):
        return self._input_type
    def frame_rate(self):
        return self._frame_rate
    def num_points(self):
        """
        Number of points in list.  Points will be evenly spaced along
        curve.
        """
        if self.base_points is not None:
            return len(self.base_points)
        return floor(self.duration()*self.frame_rate())
    def irange(self):
        """
        Return list of numbers over which to iterate to compute
        self.points().
        """
        return xrange(self.num_points())
    def include_final(self):
        return self._include_final
    def final_point(self):
        """
        Included only if self.include_final() is True.  This must be
        defined in a subclass.
        """
        raise NotImplementedError('final_point must be defined in a subclass')    
    def ith_point(self,i):
        """
        Function which returns the i^th point in the list, where i is
        in self.irange().  The output should always be in cartesian
        coordinates.  This must be defined in a subclass.
        """
        if self.base_points is None:
            raise NotImplementedError('ith_point must be defined in a subclass')
        else:
            return self.base_points[i]
    def points(self):
        """
        Return list of points on base 2-sphere, evenly spaced along a
        curve.
        """
        if self.include_final():
            if self._transformation is not None:
                final = [self._transformation(self.final_point())]
            else:
                final = [self.final_point()]
        else:
            final = []
        if self._transformation is not None:
            return [self._transformation(self.ith_point(i)) for i in self.irange()] + final
        else:
            return [self.ith_point(i) for i in self.irange()] + final
        
    def set_rotation(self,v,theta,M=None,return_matrix=False,multiply=False):
        """
        Set rotation matrix for output.  Note that the axes of the
        base 2-sphere are negatively oriented, so our rotations follow
        the *left hand rule*.  Sorry.

        EXAMPLES::

            sage: c = BaseLatitude(1,pi/2,5*pi/2,po=pi/2,input_type='sph')
            sage: c.set_rotation((1,pi/2,pi/2),pi/4)
            sage: draw_frame(c.points(), n=0, high_quality=False, only_show=True)

        """
        if M is not None:
            R = M
        else:
            R = rot_mat(v,theta)
        if multiply:
            self._rotation = self._rotation*R
        else:
            self._rotation = R
        def transf(x):
            return tuple(vector(x)*self._rotation)
        self._transformation = transf
        if return_matrix:
            return self._rotation
        else:
            return None
        
    def set_transformation(self,transf,compose=False,return_transf=False):
        """
        set an arbitrary transformation to be applied to points of self.

        ``transf`` should be a function whose inputs and outputs are
        tuples on base 2-sphere, in cartesian coordinates.
        """
        if compose:
            def g(x):
                return transf(self._transformation(x))
        else:
            g = transf
        self._transformation = g
        if return_transf:
            return transf
        else:
            return None
        
        

class BaseGreatCircle(BaseList):
    """
    List of base point positions for point moving from init_pos
    to end_pos along a great circle during duration (seconds).

    input_type is 'cc' for cartesian coordinates, or 'sph' for spherical
    coordinates.

    The great circle from P to Q is parametrized by

        t \mapsto P*\cos(t) + Q'*\sin(t)

    where P and Q are in cartesian coordinates, and

        Q' = w/|w|
        w = Q - <P,Q>*P

    is perpendicular to P, along the great circle passing through P
    and Q.

    EXAMPLES::

        sage: b = BaseGreatCircle(.5,(1,0,pi/2),(1,0,0),'sph')
        sage: draw_frame(b.points(), n=0, high_quality=False, only_show=True)

    To draw entire circle, use keyword argument whole_circle=True::

        sage: b = BaseGreatCircle(1.5,(1,0,pi/2),(1,0,0),'sph',whole_circle=True)
        sage: draw_frame(b.points(), n=0, high_quality=False, only_show=True)
    """
    def __init__(self,*args,**kwds):
        BaseList.__init__(self,*args,**kwds)
        self._repr_extra = " (great circle)"

        if self.input_type() == 'sph':
            self.P = vector(sph_to_cc(*self.init_pos()))
            self.Q = vector(sph_to_cc(*self.final_pos()))
        elif input_type == 'cc':
            self.P = vector(self.init_pos())
            self.Q = vector(self.final_pos())

        if self.P.norm() != 1 or self.Q.norm() != 1:
            raise ValueError("inital and final points must be on base 2-sphere, and therefore have length 1")
        if self.P == -1*self.Q:
            raise TypeError("antipodal points do not specify a unique great circle")

        a = self.P.inner_product(self.Q)
        w = self.Q - a*self.P
        self.Qprime = w/w.norm()
        if 'whole_circle' in self._other_kwds:
            self.alpha = 2*pi
        else:
            self.alpha = arccos(a)

    def ith_point(self,i):
        alpha_i = self.alpha*i/(self.num_points())
        return self.P*RR(cos(alpha_i)) + self.Qprime*RR(sin(alpha_i))

    def final_point(self):
        return self.Q

class BaseLatitude(BaseList):
    """
    List of base point positions for point moving from init_pos
    to end_pos along a line of latitude during duration (seconds).

    input_type is 'cc' for cartesian coordinates, or 'sph' for spherical
    coordinates.

    Initial and final points should be azmuthal angles, and an
    additional argument `po` specifies the polar angle of the
    latitude line.

    EXAMPLES::

        sage: c = BaseLatitude(.5,pi/2,3*pi/2,po=pi/3,input_type='sph')
        sage: draw_frame(c.points(), n=0, high_quality=False, only_show=True)
    """
    def __init__(self,*args,**kwds):
        BaseList.__init__(self,*args,**kwds)
        self._repr_extra = " (latitude circle)"

        if 'po' in self._other_kwds:
            self.polar = self._other_kwds['po']
        else:
            raise ValueError("please specify polar angle (keyword 'po')")
        # azmuthal angles of initial and final points
        self.P_az = self.init_pos()
        self.Q_az = self.final_pos()
        

    def ith_point(self,i):
        return sph_to_cc(1,
                         self.P_az + (i/self.num_points())*(self.Q_az - self.P_az),
                         self.polar)

    def final_point(self):
        return sph_to_cc(1,self.Q_az,self.polar)

    
class BaseSpiral(BaseList):
    """
    List of base point positions for point moving from init_pos
    to final_pos along a spiral during duration (seconds).

    Inputs init_pos and final_pos should be of the form::
    
        init_pos = (az_start, po_start)
        final_pos = (az_end, po_end)

    The spiral is parametrized in spherical coordinates by:

        t |--> (1, az_start*(1-t) + az_end*t, po_start*(1-t) + t*po_end)

    Input_type must be 'sph', for spherical coordinates.

    `tmax` has default value 1

    EXAMPLES::

        sage: a = BaseSpiral(1,(0,pi/2),(pi,0),'sph', include_final=True)
        sage: draw_frame(a.points(), n=0, high_quality=False, only_show=True)

    """
    def __init__(self,*args,**kwds):
        BaseList.__init__(self,*args,**kwds)        
        self._repr_extra = " (spiral)"

        if self.input_type() == 'sph':
            self.az_start = self.init_pos()[0]
            self.po_start = self.init_pos()[1]
            self.az_end = self.final_pos()[0]
            self.po_end = self.final_pos()[1]
        elif input_type == 'cc':
            raise ValueError("input_type must be 'sph'")
        if 'tmax' not in kwds.keys():
            self.tmax = 1
        else:
            self.tmax = kwds['tmax']
        if 'opt_func' in kwds.keys():
            self.opt_func = kwds['opt_func']
        else:
            self.opt_func = None

        self.__point_list = None

    def _point_list(self, recompute=False, return_gamma_delta=False):
        """
        Compute point list, evenly spacing points along the spiral by
        parametrizing with arc length.
        """
        if self.__point_list is not None and not recompute:
            return self.__point_list
        
        if self.opt_func is not None:
            gamma, arc_length_gamma = self.opt_func
        else:
            var('t')
            az_t = (1-t)*self.az_start + t*self.az_end
            po_t = (1-t)*self.po_start + t*self.po_end
            gamma(t) = (cos(az_t)*sin(po_t), sin(az_t)*sin(po_t), cos(po_t))

            # total arc length we'll traverse
            gadot(t) = (gamma(t)[0].diff(),gamma(t)[1].diff(),gamma(t)[2].diff())
            s = lambda t: sqrt(sum(x^2 for x in gadot(t)))
            L = numerical_integral(s,0,self.tmax)
            verbose('arc length: %s'%L[0])
            arc_length_gamma = L[0]


        def delta(t0,epz):
            def out(t1):
                g = vector(gamma(t0)) - vector(gamma(t1))
                return g.norm() - epz
            return out

        # approximate distance between each pair of points on base 2-sphere
        epsilon = arc_length_gamma/self.num_points()
        verbose('epsilon: %.3f..'%epsilon)

        time_points = [0]
        output = [gamma(0)]

        for i in range(self.num_points()):
            verbose('  getting time point %s'%i)
            t0 = time_points[i]
            try:
                t2 = min(t0+.15,self.tmax)
                t1 = find_root(delta(t0,epsilon), t0, t2)
            except RuntimeError:
                verbose(' find_root failed')
                verbose('  numpoints: '+str(self.num_points()))
                verbose('  out: '+str(len(output)))
                if len(output) >= self.num_points() - 1:
                    t1 = self.tmax
                else:
                    raise RuntimeError('find_root failed')
            if t1 > self.tmax:
                t1 = (1 + t0)/2
            time_points += [t1]
            output += [tuple([x.n() for x in gamma(t1)])]

        self.__point_list = output
        if return_gamma_delta:
            return (output, gamma, delta)
        else:
            return self.__point_list


    def ith_point(self,i):
        return self._point_list()[i]
    def final_point(self):
        return sph_to_cc(1,self.az_end,self.po_end)

class BaseFlower(BaseList):
    """
    List of base point positions for point moving from init_pos
    to final_pos along a flower during duration (seconds).

    The flower is parametrized in spherical coordinates by:

        t |--> (1, (2*pi)*t + A*cos(N*(2*pi)*t) + Q, B*sin(N*(2*pi)*t) + P)

        N = 4 controls the number of petals
        A = .5 controls the fattness of the petals
        B = -pi/7 controls the height of the petals
        P = pi/2 controls the latitude of the flower
        Q = 0 controls the rotation of the flower

    """
    def __init__(self,*args,**kwds):
        BaseList.__init__(self,*args,**kwds)        
        self._repr_extra = " (flower)"

        if self.input_type() == 'sph':
            pass
        elif input_type == 'cc':
            pass

        if 'flower_params' in self._other_kwds:
            (self.N,self.A,self.B,self.P,self.Q) = self._other_kwds['flower_params']
        else:
            self.N = 4 #controls the number of petals
            self.A = .5 #controls the fattness of the petals
            self.B = -pi/7 #controls the amplitude of the polar angle range
            self.P = pi/2 #controls the latitude of the flower
            self.Q = 0 #shifts the azmuthal angle


        self.__point_list = None

    def gamma(self):
        var('t')
        az_t = (2*pi)*t + self.A*cos(self.N*(2*pi)*t) + self.Q
        po_t = self.B*sin(self.N*(2*pi)*t) + self.P
        return (cos(az_t)*sin(po_t), sin(az_t)*sin(po_t), cos(po_t))


    def _point_list(self, recompute=False, return_gamma_delta=False):
        """
        Compute point list, evenly spacing points along the spiral by
        parametrizing with arc length.
        """
        if self.__point_list is not None and not recompute:
            return self.__point_list

        gamma(t) = self.gamma()

        def delta(t0,t1):
            g = gamma(t0) - gamma(t1)
            return sqrt(sum(a^2 for a in g))

        # total arc length we'll traverse
        gadot(t) = (gamma(t)[0].diff(),gamma(t)[1].diff(),gamma(t)[2].diff())
        s = lambda t: sqrt(sum(x^2 for x in gadot(t)))
        L = numerical_integral(s,0,1)
        verbose('arc length: %s'%L[0])
        arc_length_gamma = L[0]

        # approximate distance between each pair of points on base 2-sphere
        epsilon = arc_length_gamma/self.num_points()
        verbose('epsilon: %.3f..'%epsilon)

        time_points = [0]
        output = [gamma(0)]

        for i in range(self.num_points()):
            verbose('  getting time point %s'%i, level=2)
            t0 = time_points[i]
            #try:
            t1 = find_root(delta(t0,t) == epsilon, t0, t0+.2)
            #except RuntimeError:
            #    verbose('  quitting find_root')
                #raise RuntimeError('find_root failed')
            verbose('  .. time point is %s'%t1, level=2)
            time_points += [t1]
            output += [gamma(t1).n()]

        self.__point_list = output
        if return_gamma_delta:
            return (output, gamma, delta)
        else:
            return self.__point_list


    def ith_point(self,i):
        return self._point_list()[i]
    def final_point(self):
        pass






###########
##
## Additional helper functions  
##
###########


def transform_coordinates(movement_list,u,v,w,input_type='cc', return_matrix=False):
    """
    Apply a transformation matrix to the items in ``movement_list``,
    and return the resulting list.
    """
    if input_type == 'sph':
        M = [vector(sph_to_cc(*a)) for a in [u,v,w]]
    elif input_type == 'cc':
        M = [vector(a) for a in [u,v,w]]
    else:
        raise ValueError("input_type must be 'cc' or 'sph'")
    M = Matrix(M)
    out = [vector(m)*M for m in movement_list]
    if return_matrix:
        return (out,M)
    else:
        return out

def cc_to_sph(x,y,z):
    """
    Converts cartesian coordinates (x,y,z) in the unit ball to
    spherical coordintes (r,az,po).

    The polar angle is the angle with the z-axis, and the azmuthal
    angle is the angle with the x-z plane.
    """
    r = sqrt(x^2 + y^2 + z^2)
    if r == 0:
        return (0,0,0)
    if r > 1.00001:
        raise ValueError("r should not be greater than 1")
    az = atan2(y,x) # arctan, with sign accounting for quadrant of (x,y)
    # atan2 also works correctly for x == 0
    
    po = arccos(z/r)
    try:
        return (r.n(),az.n(),po.n())
    except AttributeError:
        return (r,az,po)

def sph_to_cc(r,az,po):
    """
    convert spherical coordinates in the unit ball
        0 \leq r  \leq 1
        0 \leq az \leq 2*pi
        0 \leq po \leq pi
    to cartesian coordinates
    """
    if po > pi or po < 0:
        raise ValueError("polar angle must be between 0 and pi; you entered "+str(po))

    try:
        return tuple(t.n() for t in (r*cos(az)*sin(po), r*sin(az)*sin(po), r*cos(po)))
    except AttributeError:
        return tuple(t for t in (r*cos(az)*sin(po), r*sin(az)*sin(po), r*cos(po)))

def rot_mat(v,theta,input_type='sph'):
    """
    Return rotation matrix with axis v and angle theta.  Input type is
    'sph' or 'cc'.

    The resulting matrix should act on the *right* of vectors.

    REFERENCES::

        [1]: http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Conversion_to_rotation_matrix
    """
    if input_type == 'sph':
        k = vector(sph_to_cc(*v))
    elif input_type == 'cc':
        k = vector(v)
    else:
        raise ValueError("input_type must be 'cc' or 'sph'")
    verbose('rotation axis: (%.03f, %.03f, %.03f)'%(k[0],k[1],k[2]), level=2)
    k_cross = matrix([[0, -1*k[2], k[1]],[k[2], 0, -1*k[0]],[-1*k[1], k[0], 0]])
    R = identity_matrix(RR,3) + sin(theta)*k_cross + (1 - cos(theta))*k_cross^2

    return matrix(3,3,[x.n() for x in R.transpose()])
