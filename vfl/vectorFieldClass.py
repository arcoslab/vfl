#!/usr/bin/python
# Copyright (c) 2009 Technische Universitaet Muenchen, Informatik Lehrstuhl IX.
# Author: Federico Ruiz-Ugalde <ruizf at in.tum.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
from math import sqrt
from math import sin
from numpy import array
from numpy import inner

#Tratar de accesar un miembro (variable interna) de un funcion usando la notacion de punto, 
#suponiendo que una funcion tambien es un objecto.
#Puede servir para poder cambiar los parametros de las funciones de los 
#campos vectoriales por aparte de la clase VectorField.



class VectorField:
  def __init__(self,func=lambda pos: array([0]*16)):
    self.setFunc(func)
  def setFunc(self,func):
    self.func=func
  def getFunc(self):
    return self.func
  def getVector(self,pos):
    return self.func(pos)
  def __add__(self,other):
    return VectorField(lambda pos: self.func(pos)+other.func(pos))
  def __mul__(self,other):
    if ((other.__class__.__name__)==(self.__class__.__name__)):
      return VectorField(lambda pos: self.func(pos)*other.func(pos))
    else:
      if (self.__class__.__name__=='VectorField'):
        return VectorField(lambda pos: self.func(pos)*other)
      else:
        return VectorField(lambda pos: other.func(pos)*self)
        #I think is like this:
        #return VectorField(lambda pos: other*self.func(pos))
  def __div__(self,other):
    if ((other.__class__.__name__)==(self.__class__.__name__)):
        return VectorField(lambda pos: self.func(pos)/other.func(pos))
    else:
        if (self.__class__.__name__=='VectorField'):
          return VectorField(lambda pos: self.func(pos)/other)
        else:
          return VectorField(lambda pos: other.func(pos)/self)
          
                    
  def __sub__(self,other):
    return VectorField(lambda pos: self.func(pos)-other.func(pos))
  def normCart(self):
    return VectorField(lambda pos: self.func(pos)/array([self.funcLenght(pos)]*3+[1]*3))
  def funcLenght(self,pos):
    length=sqrt(inner(self.func(pos)[0:3],self.func(pos)[0:3]))
    if length==0:
      return 1.0
    else:
      return length

class ScalarField:
  def __init__(self,func=lambda pos: array([0]*16)):
    self.setFunc(func)
  def setFunc(self,func):
    self.func=func
  def getFunc(self):
    return self.func
  def getScalar(self,pos):
    return self.func(pos)
  def __mul__(self,other):
    if ((other.__class__.__name__)==(self.__class__.__name__)):
      return ScalarField(lambda pos: self.func(pos)*other.func(pos))
    else:
      if (self.__class__.__name__=='ScalarField'):
        return ScalarField(lambda pos: self.func(pos)*other)
      else:
        return ScalarField(lambda pos: other.func(pos)*self)
  def __div__(self,other):
    if ((other.__class__.__name__)==(self.__class__.__name__)):
        return ScalarField(lambda pos: self.func(pos)/other.func(pos))
    else:
        if (self.__class__.__name__=='ScalarField'):
          return ScalarField(lambda pos: self.func(pos)/other)
        else:
          return ScalarField(lambda pos: other.func(pos)/self)

class VectorField_test:
  def __init__(self,funcX=lambda pos,param: 0, funcY=lambda pos,param: 0, funcZ=lambda pos,param: 0, param=None):
    param = param or []
    self.setFunc(funcX,funcY,funcZ)
    self.setParam(param)
  def setFunc(self,funcX, funcY, funcZ):
    self.funcX=funcX
    self.funcY=funcY
    self.funcZ=funcZ
  def setParam(self,param):
    self.param=param
  def getParam(self):
    return self.param
  def getFunc(self):
    return self.funcX, self.funcY, self.funcZ
  def getVector(self,pos):
    return self.funcX(pos,self.param),self.funcY(pos,self.param),self.funcZ(pos,self.param)
  def __add__(self,other):
    paramNum1=len(self.param)
    paramNum2=len(other.param)
    paramNew=self.param+other.param
    return VectorField(lambda pos,param: self.funcX(pos,param[:paramNum1])+other.funcX(pos,param[paramNum1:]), lambda pos,param: self.funcY(pos,param[:paramNum1])+other.funcY(pos,param[paramNum1:]), lambda pos,param: self.funcZ(pos,param[:paramNum1])+other.funcZ(pos,param[paramNum1:]),paramNew)
  def __mul__(self,other):
    paramNum1=self.param.size()
    paramNum2=other.param.size()
    paramNew=self.param+other.param
    return VectorField(lambda pos,param: self.funcX(pos,param[:paramNum1])*other.funcX(pos,param[paramNum1:]), lambda pos,param: self.funcY(pos,param[:paramNum1])*other.funcY(pos,param[paramNum1:]), lambda pos,param: self.funcZ(pos,param[:paramNum1])*other.funcZ(pos,param[paramNum1:]),paramNew)

            
def attractorSinHib(pos,param):
#x,y,z current position
#x0,y0,z0 attraction point
#xm,ym,zm starting point
    x=pos[0]
    y=pos[1]
    z=pos[2]
    x0=param[0]
    y0=param[1]
    z0=param[2]
    xm=param[3]
    ym=param[4]
    zm=param[5]
    done=False
    k2=0.00001
    k=1.0/sqrt(1+(k2*sqrt(x**2+y**2+z**2))/(sqrt(xm**2+ym**2+zm**2)))
    #k=1.0/sqrt(1+xm)
    
    distance=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))
    #distance=sqrt((x-x0)*(x-x0))
    print "Current position: ",x, y, z
    print "Distance to the attractor: ", distance
    distancem=sqrt((x0-xm)*(x0-xm)+(y0-ym)*(y0-ym)+(z0-zm)*(z0-zm))
    #print distancem
    k=3.1416/(distancem+(distancem/50.0))
    if distance!=0:
        X=-sin(distance*k)*(x-x0)/distance
        Y=-sin(distance*k)*(y-y0)/(distance)
        Z=-sin(distance*k)*(z-z0)/distance
    else:
        X=0
        Y=0
        Z=0
    if distance>(distancem-(distancem/100.0)):
        print "Special Point: Outside of the sinusoidal attractor"
        X=-(x-x0)
        Y=-(y-y0)
        Z=-(z-z0)
        vectorSize=sqrt(X*X+Y*Y+Z*Z)
        outsideVectorSpeed=0.5
        X=outsideVectorSpeed*X/vectorSize
        Y=outsideVectorSpeed*Y/vectorSize
        Z=outsideVectorSpeed*Z/vectorSize
        #oVel=0.1
        #if X>0:
        #  X=oVel
        #else:
        #  if X!=0:
        #    X=-oVel
        #if Y>0:
        #  Y=oVel
        #else:
        #  if Y!=0:
        #    Y=-oVel
        #if Z>0:
        #  Z=oVel
        #else:
        #  if Z!=0:
        #    Z=-oVel
        #normalize vector
    #Y=sin(sqrt(x*x+y*y))
    #Y=-sin(x+y)/2
    #Z=-(z-z0)*k
    #scale=2 #For the robot
    #scale=20 #For the simulator
    #scale=0.3 # for b21
    scale=0.5
    X=X*scale
    Y=Y*scale
    Z=Z*scale
    if distance<0.02:
      X=0
      Y=0
      Z=0
      #done=True
    #vector=[X,Y,Z,done]
    vector=[X,Y,Z]
    return array(vector)
  


def attractor3(x,y,z,x0,y0,z0,xm,ym,zm): 
#x,y,z current position
#x0,y0,z0 attraction point
#xm,ym,zm starting point
    done=False
    k2=0.00001
    k=1.0/sqrt(1+(k2*sqrt(x**2+y**2+z**2))/(sqrt(xm**2+ym**2+zm**2)))
    #k=1.0/sqrt(1+xm)
    
    distance=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))
    #distance=sqrt((x-x0)*(x-x0))
    print "Current position: ",x, y, z
    print "Distance to the attractor: ", distance
    distancem=sqrt((x0-xm)*(x0-xm)+(y0-ym)*(y0-ym)+(z0-zm)*(z0-zm))
    #print distancem
    k=3.1416/(distancem+(distancem/50.0))
    if distance!=0:
        X=-sin(distance*k)*(x-x0)/distance
        Y=-sin(distance*k)*(y-y0)/(distance)
        Z=-sin(distance*k)*(z-z0)/distance
    else:
        X=0
        Y=0
        Z=0
    if distance>(distancem-(distancem/100.0)):
        print "Special Point: Outside of the sinusoidal attractor"
        X=-(x-x0)
        Y=-(y-y0)
        Z=-(z-z0)
        vectorSize=sqrt(X*X+Y*Y+Z*Z)
        outsideVectorSpeed=0.5
        X=outsideVectorSpeed*X/vectorSize
        Y=outsideVectorSpeed*Y/vectorSize
        Z=outsideVectorSpeed*Z/vectorSize
        #oVel=0.1
        #if X>0:
        #  X=oVel
        #else:
        #  if X!=0:
        #    X=-oVel
        #if Y>0:
        #  Y=oVel
        #else:
        #  if Y!=0:
        #    Y=-oVel
        #if Z>0:
        #  Z=oVel
        #else:
        #  if Z!=0:
        #    Z=-oVel
        #normalize vector
    #Y=sin(sqrt(x*x+y*y))
    #Y=-sin(x+y)/2
    #Z=-(z-z0)*k
    #scale=2 #For the robot
    #scale=20 #For the simulator
    scale=0.3 # for b21
    X=X*scale
    Y=Y*scale
    Z=Z*scale
    if distance<0.02:
      X=0
      Y=0
      Z=0
      #done=True
    vector=[X,Y,Z,done]
    return vector

def attractor2(x,y,z,x0,y0,z0,xm,ym,zm):
    k2=0.00001
    k=1.0/sqrt(1+(k2*sqrt(x**2+y**2+z**2))/(sqrt(xm**2+ym**2+zm**2)))
    #k=1.0/sqrt(1+xm)
    
    #distance=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))
    distance=sqrt((x-x0)*(x-x0))
    X=-sin(k*(x-x0))/2
    #Y=sin(sqrt(x*x+y*y))
    #Y=-sin(x+y)/2
    Y=-sin(k*(y-y0))/2
    Z=-(z-z0)*k
    vector=[X,Y]
    return vector


def attractor1(x,y,z,x0,y0,z0,xm,ym,zm):
    k2=0.000001
    k=0.10/sqrt(1+(k2*sqrt(x**2+y**2+z**2))/(sqrt(xm**2+ym**2+zm**2)))
    #k=1.0/sqrt(1+xm)
    X=-(x-x0)*k
    Y=-(y-y0)*k
    Z=-(z-z0)*k
    vector=[X,Y]
    return vector

