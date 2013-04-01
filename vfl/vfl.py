#!/usr/bin/env python
# Copyright (c) 2009 Technische Universitaet Muenchen, Informatik Lehrstuhl IX.
# Author: Federico Ruiz-Ugalde <memeruiz at gmail.com>
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
from math import pow
from math import exp
from numpy import dot
from numpy import array
from numpy import matrix
from numpy import inner
from numpy import concatenate
import PyKDL as kdl
import time
from vectorFieldClass import *

printcounter = 0
epsilon=0.0000001

class vectorFieldLibrary:
  def __init__(self):
    '''Registers all vectorfield types'''
    print "Registering vector fields"
    self.vfs=[]
    self.vfs.append(nullAttractor)
    self.vfs.append(pointAttractor)
    self.vfs.append(decayAttractor)
    self.vfs.append(planeAttractor)
    self.vfs.append(hemisphereAttractor)
    self.vfs.append(funnelAttractor)
    #add here more vector field functions
  def __getitem__(self,number):
    """returns the vector field function"""
    return self.vfs[number]
  def __contains__(self,item):
    return (item < self.vfs.__len__())

#all vectorfield classes must have a constructor a setParams a getScalar and
#a getVector function, with parameters, nothing and position as arguments
#respectively

class nullAttractor:
  def __init__(self):
    pass
  def setParams(self,param):
    pass
  def getScalar(self,pos):
    return array([1,1])
  def getVector(self,pos):
    P=array([0.0]*6)
    return P
    
class pointAttractor:
  '''Unit pointAttractor with rotation'''
  def __init__(self):
    self.frame=kdl.Frame()
    self.gframe=kdl.Frame()
    self.gp=array([0,0,0])
    self.p=array([0,0,0])
    self.dToAttractor=0
    self.slowDownDistance=0.05 #normal value 0.03
    self.tw=kdl.diff(self.frame,self.gframe)
    self.rot_factor=1.0
    self.min_rot=0.09
    self.stopDistance=0.005
    self.r_stopDistance=0.03
    self.slowDownTime=0.1
    self.slowDowninitTime=time.time()
    self.time=time.time()
    self.filter=[1.0]*5
    self.filterC=array([0.05,0.3,0.3,0.3,0.05])
    #self.filterC=array([1.0/20]*20)
  def setParams(self,param):
    '''
    0-15: homogeneous goal matrix
    16: slow down distance
    '''
    gpos=matrix([param[0:4],param[4:8],param[8:12],param[12:16]])
    self.gp=gpos[0:3,3].T.A[0]
    self.dToAttractor=length(self.p-self.gp)
    #self.rot_factor=param[16]
    self.slowDownDistance=param[16]
    for i in range(3):
      for j in range(3):
        self.gframe.M[i,j]=gpos[i,j]
    self.tw=kdl.diff(self.frame,self.gframe)
  def getScalar(self,pos):
    #self.time=time.time()
    #print "Time: ", self.time
    if self.dToAttractor<self.slowDownDistance:
      #ss=(-1/self.slowDownTime)*(self.time-self.slowDowninitTime)+1
      ss=self.dToAttractor/self.slowDownDistance
    else:
      #self.slowDowninitTime=self.time
      ss=1
    #self.filter.pop(0)
    #self.filter.append(ss)
    #filtertemp=array(self.filter)
    #ss=dot(filtertemp,self.filterC)
    
    d_rot = sqrt(self.tw.rot[0]**2+self.tw.rot[1]**2+self.tw.rot[2]**2) # feeling uneasy - don't we run into Euler singularities 
    #print "twist distance: ", d_rot, " pos distance: ", self.dToAttractor
    if self.dToAttractor>0:
      #r_weight = ((self.rot_factor*ss)/self.dToAttractor) * d_rot # (ap/dToAttractor) always smaller than 1
      r_weight = ((self.rot_factor*ss)/self.dToAttractor) # (ap/dToAttractor) always smaller than 1
    else:
      if d_rot>0:
        r_weight=self.min_rot/(d_rot)
      else:
        r_weight=0
    #minimum rotation speed
    #if r_weight*d_rot < self.min_rot:
    #  print "Minimum rotation speed!"
    #  if d_rot>0:
    #    r_weight = self.min_rot/(d_rot)
    #  else:
    #    r_weigth=0
    st=r_weight
    
    #if self.dToAttractor<self.stopDistance:
    #  ss=0
    #if d_rot<self.r_stopDistance:
    #  st=0
    #print "Scaling: ", ss, " ",st
    #st=0
    return array([ss,st]) #scalar speed, scalar twist
    
  def getVector(self,pos):
    ''' 0-15: position'''
    pos=matrix([pos[0:4],pos[4:8],pos[8:12],pos[12:16]])
    self.p=pos[0:3,3].T.A[0]
    r=pos[0:3,0:3]
    #Cartesian vector
    self.dToAttractor=length(self.p-self.gp)
    if self.dToAttractor==0:
      P=array([0,0,0])
    else:
      P=-(self.p-self.gp)/self.dToAttractor
      #todo: to improve stability we could decay near the attraction point
    #rotational twist
    for i in range(3):
      for j in range(3):
        self.frame.M[i,j]=r[i,j]
    self.tw=kdl.diff(self.frame,self.gframe)

    return concatenate((P,array([self.tw.rot[0],self.tw.rot[1],self.tw.rot[2]])),1)

class decayAttractor:
  def __init__(self):
    self.apos=array([0,0,0])
    self.aradius=0.1
    self.safeD=0.05
    self.cutD=self.aradius+self.safeD
    self.order=5
  def setParams(self,param):
    '''Repeller with bellshape decay
    Parameters:
    0-2: position of obstacle
    3: obstacle radius
    4: extra safe distance
    5: order of the decay
    '''
    self.apos=array(param[0:3])
    self.aradius=param[3]
    self.safeD=param[4]
    self.cutD=self.aradius+self.safeD
    self.order=param[5]
  def getScalar(self,pos):
    return array([1,1])
  def getVector(self,pos):
    pos=matrix([pos[0:4],pos[4:8],pos[8:12],pos[12:16]])
    p=pos[0:3,3].T.A[0]
    P=p-self.apos
    l=length(P)
    if l==0:
      P=array([0,0,1])
    else:
      P=P/l
    k=bellshape(self.cutD,self.order,l)
    P=-P*k
    return concatenate((P,array([0,0,0])),1)


class planeAttractor:
  def __init__(self):
    self.epsilon=0.0000001
    self.planePoint=array([0,0,0])
    self.normalAxis=array([0,0,1])
    self.order=1.2
    self.safeDistance=0.05
    self.d=dot(-self.normalAxis,self.planePoint)
  def setParams(self,param):
    '''parameters:
    0-2: plane point
    3-5: Normal Axis
    6: decay order
    7: decay cut distance
    '''
    self.planePoint=array(param[0:3])
    normalAxis=array(param[3:6])
    self.normalAxis=normalAxis/length(normalAxis)
    self.order=param[6]
    self.safeDistance=param[7]
    self.d=dot(-self.normalAxis,self.planePoint)
  def getScalar(self,pos):
    return array([1,1])
  def getVector(self,pos):
    pos=matrix([pos[0:4],pos[4:8],pos[8:12],pos[12:16]])
    p=pos[0:3,3].T.A[0]
    distancetoPlane=dot(self.normalAxis,p)+self.d
    if (distancetoPlane>=0):
      k=bellshape(self.epsilon+self.safeDistance,self.order,distancetoPlane)
      P=-k*self.normalAxis
    else:
      P=array([0,0,0])
    return concatenate((P,array([0,0,0])),1)

class hemisphereAttractor:
  def __init__(self):
    self.epsilon=0.0000001
    self.planePoint=array([0,0,0])
    self.normalAxis=array([0,0,1])
    self.order=1.2
    self.safeDistance=0.05
    self.d=dot(-self.normalAxis,self.planePoint)
  def setParams(self,param):
    '''parameters:
    0-2: plane point
    3-5: Normal Axis
    6: decay cut distance
    7: decay order
    '''
    self.planePoint=array(param[0:3])
    normalAxis=array(param[3:6])
    self.normalAxis=normalAxis/length(normalAxis)
    self.order=param[7]
    self.safeDistance=param[6]
    self.d=dot(-self.normalAxis,self.planePoint)
  def getScalar(self,pos):
    return array([1,1])
  def getVector(self,pos):
    pos=matrix([pos[0:4],pos[4:8],pos[8:12],pos[12:16]])
    p=pos[0:3,3].T.A[0]
    distancetoPlane=dot(self.normalAxis,p)+self.d
    if (distancetoPlane>=0):
      k=bellshape(self.epsilon+self.safeDistance,self.order,distancetoPlane)
      P=-k*self.normalAxis
    else:
      P=-self.normalAxis
    return concatenate((P,array([0,0,0])),1)


class funnelAttractor:
  def __init__(self):
    self.point=array([0,0,0])
    self.axis=array([0,0,1])
    self.aT=0.2
    self.ar=10
    self.dT=0.1
    self.dr=2
  def setParams(self,param):
    '''Params:
    0-2: position
    3-5: axis
    6: cut angle
    7: angle decay order
    8: cut distance 
    9: distance decay order
    '''
    self.point=array(param[0:3])
    self.axis=array(param[3:6])/length(array(param[3:6]))
    print "Axis: ", self.axis
    self.aT=param[6]
    self.ar=param[7]
    self.dT=param[8]
    self.dr=param[9]
  def getScalar(self,pos):
    return array([1,1])
  def getVector(self,pos):
    #goal position and approach vector
    pos=matrix([pos[0:4],pos[4:8],pos[8:12],pos[12:16]])
    p=pos[0:3,3].T.A[0]
    if length(p-self.point)==0:
      P=array([0,0,0])
    else:
      P=-(p-self.point)/length(p-self.point)
    angle=dot(-P,-self.axis)
    if angle<0:
      angle=0.0
    #angle exponential decay
    k=expDecay(1-angle,self.ar,self.aT)
    #distance exponential decay
    k2=expDecay(length(p-self.point),self.dr,self.dT)
    P=P*k*k2
    return concatenate((P,array([0,0,0])),1)


def goalObstacleVField2(goal,goalObstacle,approach,obstacle1,table):
  #Deals with the approach direction by using a vector of approach
  #approach has to be unitary
  
  #parameters goalObstacle   
  g0=goal
  o0=goalObstacle
  obstacleRadius=[0.07]
  obstacleSafeDistance=[0.005]
  obstacleOrder=[1.4]
  initialForce=[2.0] 
  obstacleOrder=[2.5]
  initialForce=[5.0]
  goalObstacleParams=obstacleRadius+obstacleSafeDistance+obstacleOrder+initialForce
  
  #parameters goal Obstacle Hyperbola shadow Repeller
  #HypB=[2] #Parameter associated with the shadow width
  bPercentage=[0.9] # Can not be bigger than 1, associated with the shadow width, bigger value more w
  yOrder=[1.1]
  yCutDistance=[0.15]
  xOrder=[2.5]
  xT=[0.5]
  goalObstacleHypParams=bPercentage+yOrder+yCutDistance+xOrder+xT

  #parameters obstacle1
  obstacleRadius=[0.12]
  obstacleSafeDistance=[0.005]
  obstacleOrder=[4]
  initialForce=[5.0] 
  obstacle1Params=obstacleRadius+obstacleSafeDistance+obstacleOrder+initialForce

  #parameters Obstacle1 Hyperbola shadow Repeller
  #HypB=[2] #Parameter associated with the shadow width
  bPercentage=[0.9] # Can not be bigger than 1, associated with the shadow width, bigger value more w
  yOrder=[1.1]
  yCutDistance=[0.12]
  xOrder=[2.5]
  xT=[0.5]
  obstacle1HypParams=bPercentage+yOrder+yCutDistance+xOrder+xT
  
  #table parameters
  tableInitForce=[10.0]
  tableHeight=[0.82]
  tableHeight=table
  tableOrder=[100.0]
  tableSafeDistance=[0.00]
  tableParams=tableInitForce+tableHeight+tableOrder+tableSafeDistance
  
  vectorFieldObstacle1=VectorField(decayRepeller1,obstacle1+obstacle1Params)
  vectorFieldObstacle1Hyp=VectorField(vectorF_hyperbola,g0+obstacle1+obstacle1HypParams)
  vectorFieldtable=VectorField(XYplanarRepeller,tableParams)
  vectorFieldgoal=VectorField(simpleAttractor1,g0)
  vectorFieldgoalObstacle=VectorField(decayRepeller1,o0+goalObstacleParams)
  vectorFieldgoalObstacleHyp=VectorField(vectorF_hyperbola,g0+o0+goalObstacleHypParams)
  vectorFieldAngle=VectorField(vectorF_Angle,g0+approach)
  vectorField3=vectorFieldgoal+vectorFieldgoalObstacle+vectorFieldgoalObstacleHyp+vectorFieldAngle*10+vectorFieldtable+vectorFieldObstacle1+vectorFieldObstacle1Hyp
# With no table repeller!!!
  #vectorField3=vectorFieldgoal+vectorFieldgoalObstacle+vectorFieldgoalObstacleHyp+vectorFieldAngle*10+vectorFieldObstacle1+vectorFieldObstacle1Hyp
  vectorField=vectorField3.norm()
  #vectorField=vectorField3
  
  return vectorField



def trapezoidWeight(startpos,goal,pos):
  maxlen=len(startpos-goal)
  curlen=len(pos-goal)
  risingEdge=0.75
  fallingEdge=0.70
  
  printnow=False
  global printcounter
  printcounter = printcounter+1
  if (printcounter >= 10):
    printcounter=0
    printnow=True


  if (printnow):
    print "Distance to the attractor: ",curlen
  perlen=curlen/maxlen
  if ((perlen>risingEdge) and (perlen<=1)):
    #rising edge
    k=(-1/fallingEdge)*(perlen-risingEdge)+1
    if (printnow):
        print "Accelerating! ", k
  elif perlen<fallingEdge:
    #falling edge
    k=(1/fallingEdge)*perlen
    #print "Deaccelerating! ", k
  elif ((perlen>=fallingEdge) and (perlen<=risingEdge)):
    #flat center
    k=1
  else:
    k=0.17 # outside of the attractor
    if (printnow):
        print "Outside of the attractor ", perlen
  #if curlen<0.02: #if distance to attractor is two centimeters, stop
    #k=0
  return k

def simpleAttractor2(pos,param):
  #attrator position:
  x0=param[0]
  y0=param[1]
  z0=param[2]
  p0=array(param)
  p=array(pos)
  #distance=sqrt(pow((x-x0),2),pow((y-y0),2),pow((z-z0),2))
  dist=sqrt(inner((p-p0),(p-p0)))
  P=-(p-p0)/sqrt(inner((p-p0),(p-p0)))
  k=1-0.05*dist
  P=P*k
  return P


def simpleRepeller1(pos,param):
  #repeller position:
  x0=param[0]
  y0=param[1]
  z0=param[2]
  p0=array(param)
  p=array(pos)
  #distance=sqrt(pow((x-x0),2),pow((y-y0),2),pow((z-z0),2))
  P=(p-p0)/sqrt(inner((p-p0),(p-p0)))
  return P

def simpleRepeller1a(pos,param):
  #repeller position:
  x0=param[0]
  y0=param[1]
  z0=param[2]
  p0=array(param)
  p=array(pos)
  #distance=sqrt(pow((x-x0),2),pow((y-y0),2),pow((z-z0),2))
  P=(p-p0)/sqrt(inner((p-p0),(p-p0)))
  dist=sqrt(inner((p-p0),(p-p0)))
  distC=1
  n=10
  k=1.0/sqrt(1.0+pow(dist/distC,2*(n)))
  P=k*P*2
  return P


def simpleRepeller2(pos,param):
  #repeller position:
  x0=param[0]
  y0=param[1]
  z0=param[2]
  p0=array(param[:3])
  dist0=param[3]
  p=array(pos)
  d=p-p0
  dist=sqrt(inner(d,d))
  if dist<=dist0:
    P=(p-p0)/sqrt(inner((p-p0),(p-p0)))
  else:
    P=array([0,0,0])
  return P

def simpleRepeller3(pos,param):
  p=array(pos)
  f1=array([5,0,0]) #focus 1
  f2=-f1
  dist1=sqrt(inner(p-f1,p-f1))
  dist2=sqrt(inner(p-f2,p-f2))
  diff=dist2-dist1
  if dist1!=0:
    k=10/dist1
    k=1-0.1*dist1
  else:
    k=1
  if diff>9:
    P=(p-f1)/sqrt(inner((p-f1),(p-f1)))
    #P=P*k
  else:
    P=array([0,0,0])
  return P

def lateralRepeller1(pos,param):
  p=array(pos)
  obs=array(param)
  goal=array([0,0,0])
  shadowVector=obs-goal
  normShadowVector=sqrt(dot(shadowVector,shadowVector))
  if normShadowVector==0:
    print "obstacle and goal can't be in the same place"
  shadowVector=shadowVector/normShadowVector
  shadowVectorScaled=dot(p,shadowVector)*shadowVector
  shadowRepel=p-shadowVectorScaled
  normShadowVector=sqrt(dot(shadowRepel,shadowRepel))
  if normShadowVector==0:
    shadowRepel=shadowRepel*0
  else:
    shadowRepel=shadowRepel/normShadowVector
  distToObstacle=sqrt(dot((p-obs),(p-obs)))
  k=bellshape(3,10,distToObstacle)
  shadowRepel=shadowRepel*k
  return shadowRepel


def XYplanarRepeller(pos,param):
  #Repellor origin is z=0
  p=array(pos)
  initialForce=param[0]
  tablePos=param[1] #table height
  order=param[2]
  safeDistance=param[3]
  #k=bellshape(tablePos,order,p[2])
  k=1
  #if p[2]<0:
  #  P=array([0,0,0])
  #else:
  #  P=k*array([0,0,initialForce])
  if p[2]<tablePos:
    P=array([0,0,initialForce])
  else:
    P=array([0,0,0])
  
  return P


def bellshape(cutDistance,order,x):
  distC=cutDistance
  n=order
  if x<epsilon:
    return 1
  else:
    return (1.0/sqrt(1.0+pow(x/distC,2*(n))))

def hyperbola(x,focus,a): #a: from the goal to the first nearer point in the curve
  b=sqrt((focus*focus)-(a*a)) # less distance between focus and a means less open hyperbola
  y=b*sqrt((pow(x,2)/pow(a,2))-1.0)
  return [x,y,y]
  
def hyperbolaDisplacedRotated(pos,goal,obstacle,b):
  obstacle2temp=obstacle-goal
  normobstacle2=sqrt(dot(obstacle2temp,obstacle2temp))
  #obstacle2=array([normobstacle2,0,0])
  pos2=pos-goal
  if normobstacle2==0:
    print "Goal can't be in the same position as the obstacle"
  hyp=(pos2)-(dot(pos2,obstacle2temp/normobstacle2)*(obstacle2temp/normobstacle2))
  #hyp=pos2
  #normx2=sqrt(dot(pos2,pos2))
  x2=dot(pos2,(obstacle2temp/normobstacle2))
  lenFocus1=len(obstacle2temp)
  a=sqrt(pow(lenFocus1,2)-pow(b,2))
  hyp2temp=hyperbola(x2,normobstacle2,a) #x2 and hyp2temp[1] are x and y in a non rotated hyperbola
  return [hyp,hyp2temp] #returns vector perpendicular to the x of the hyperbola in the global frame 
  # and the x and y pair of the hyperbola in that x
  
def hyperbolaInside(goal,obstacle,pos,b): #b can not be bigger than len(obstacle-goal)
  focus1=obstacle-goal
  focus2=-focus1
  posrel=pos-goal
  dist1=sqrt(dot(posrel-focus1,posrel-focus1))
  dist2=sqrt(dot(posrel-focus2,posrel-focus2))
  lenFocus1=len(focus1)
  #print "Length Focus 1", lenFocus1, "B: ",b
  a=sqrt(pow(lenFocus1,2)-pow(b,2))
  diff0=a*2.0
  diff=dist2-dist1
  if diff>diff0:
    k=1
  else:
    k=0
  return k
  
def len(x):
  y=sqrt(dot(x,x))
  return y
def length(x):
  y=sqrt(dot(x,x))
  return y
  
def vectorF_hyperbola(pos,param):
  #goal=array([-6,-6,0])
  goal=array(param[:3])
  #obstacle=array([2,2,0])
  obstacle=array(param[3:6])
  #bPercentage=b/len(obstacle-goal)
  #b=2
  #b=param[6]
  b=param[6]*len(obstacle-goal)
  yOrder=param[7]
  yCutDistance=param[8]
  xOrder=param[9]
  xT=param[10]
  k=hyperbolaInside(goal,obstacle,pos,b)
  if k==1:
    temp=hyperbolaDisplacedRotated(pos,goal,obstacle,b)
    perpHyptemp=temp[0]
    if len(perpHyptemp)==0:
      perpHyp=array([0,0,0])
    else:
      perpHyp=perpHyptemp/len(perpHyptemp)
    P=array([perpHyp[0],perpHyp[1],perpHyp[2]])
    Hyp=temp[1]
    y=Hyp[1]
    if y==0:
      y2=0
    else:
      y2=len(perpHyptemp)/y #Percentage of y position in hyperbola
    #print "percentage ",y2, " Y: ", len(perpHyptemp), " Ymax: ", y
    #x=len(pos-obstacle)
    #yOrder=1.1
    #yCutDistance=0.1
    k=bellshape(yCutDistance,yOrder,y2) #bellshape in Y of the hyperbola
    #cutDistance2=
    #k=k*bellshape(cutDistance2,order2,x) #bellshape for the x distance
    #T=6
    #r=8
    decayX=len(Hyp[0])-len(obstacle-goal)
    kexp=expDecay(len(decayX),xOrder,xT)
    #print "X decay: ",decayX,"Decay: ",kexp
    P=P*kexp*k
  else:
    P=array([0,0,0])
  return P
  
def expDecay(x,r,T):
  y=pow(r,-x/T)
  return y

def simpleRepeller4(pos,param):
  p=array(pos)
  f1=array([5,0,0]) #focus 1
  f2=-f1
  dist1=sqrt(inner(p-f1,p-f1))
  dist2=sqrt(inner(p-f2,p-f2))
  diff=dist2-dist1
  diff0=9
  a=diff0/2
  c=5
  x=p[0]
  temp1=((pow(x,2)/pow(a,2))-1)
  temp2=(pow(c,2)-pow(a,2))
  if temp1*temp2>=0:
    y=sqrt(temp1*temp2)
  else:
    y=1
  if dist1!=0:
    k=10/dist1
    k=1-0.1*dist1
  else:
    k=1
  #P=(p-f1)/sqrt(inner((p-f1),(p-f1)))
  #P=P*2.5*pow(10,-dist1/0.8)
  if diff>9:
    P=(p-f1)/sqrt(inner((p-f1),(p-f1)))
    P=P*2.5*pow(10,-dist1/0.8)
    if (p[0]>a):
      P2=array([0,P[1],P[2]])
      k=sqrt(inner(P2,P2))
      
      P2=P2/k
      #P=P*(1-0.07*dist1-(0.8-0.001*dist1)*k)
      print "hiperbola",y
      yper=pos[1]/y
      print yper
      ##P2=P2*(1-abs(yper))
      ##P2=P2*(1-0.055*x)
      P2=P2*2.5*pow(10,-(x-a)/7)
      P2=P2*1.2*pow(10,-(abs(yper))/0.3)
    P=P*1+P2*0
    #P=P*k
  else:
    P=array([0,0,0])
  return P


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

