import math
import random as rnd
import matplotlib.pyplot as plt
from matplotlib import animation


sx=0.6
x1=-0.5
x2=0.5
ma=1
mb=5
va=1
vb=5
dt=0.001
valpha=vb-va
vbeta=va/valpha
if ma!=mb:
	malpha=mb-ma
	mbeta=ma/malpha
	m1=malpha*(mbeta+rnd.random())
	m2=malpha*(mbeta+rnd.random())
else:
	m1=ma
	m2=ma
v1x=valpha*(vbeta+rnd.random())
v1y=0
v2x=-valpha*(vbeta+rnd.random())
v2y=0
r1=m1/(10*mb)
r2=m2/(10*mb)
sm=m1+m2
dm=m2-m1
y1=0
y2=r1
sy=1.25*(r1+r2)
ket=0.5*(m1*(v1x**2+v1y**2)+m2*(v2x**2+v2y**2))
ptx=m1*v1x+m2*v2x
pty=m1*v1y+m2*v2y
pt=math.sqrt(ptx**2+pty**2)
kea=[]
ke1a=[]
ke2a=[]
pa=[]
p1a=[]
p2a=[]
p1xa=[]
p1ya=[]
p2xa=[]
p2ya=[]
pxa=[]
pya=[]
maxlist=[]
minlist=[]
time=[]
frps=60
sec=60

fig, a=plt.subplots()
fig.tight_layout()

def run(frame):
	global x1,x2,y1,y2,v1x,v1y,v2x,v2y
	plt.clf()
	dx=x2-x1
	dy=y2-y1
	dr=math.sqrt(dx**2+dy**2)
	dxp1=x2-x1+2*sx
	dxp2=x2-x1-2*sx
	dyp1=y2-y1+2*sy
	dyp2=y2-y1-2*sy
	drxp1=math.sqrt(dxp1**2+dy**2)
	drxp2=math.sqrt(dxp2**2+dy**2)
	dryp1=math.sqrt(dx**2+dyp1**2)
	dryp2=math.sqrt(dx**2+dyp2**2)
	drxp1yp1=math.sqrt(dxp1**2+dyp1**2)
	drxp1yp2=math.sqrt(dxp1**2+dyp2**2)
	drxp2yp1=math.sqrt(dxp2**2+dyp1**2)
	drxp2yp2=math.sqrt(dxp2**2+dyp2**2)
	sr=r1+r2
	if dr>sr and drxp1>sr and drxp2>sr and dryp1>sr and dryp2>sr and drxp1yp1>sr and drxp1yp2>sr and drxp2yp1>sr and drxp2yp2>sr  :
		x1+=v1x*dt
		x2+=v2x*dt
		y1+=v1y*dt
		y2+=v2y*dt
	if dr<=sr or drxp1<=sr or drxp2<=sr or dryp1<=sr or dryp2<=sr or drxp1yp1<=sr or drxp1yp2<=sr or drxp2yp1<=sr or drxp2yp2<=sr:
		x1-=v1x*dt
		x2-=v2x*dt
		y1-=v1y*dt
		y2-=v2y*dt
		dvx=v2x-v1x
		dvy=v2y-v1y
		qa=dvx**2+dvy**2
		if dr<=sr:
			drx=x2-x1
			dry=y2-y1
		elif drxp1<=sr:
			drx=x2-x1+2*sx
			dry=y2-y1
		elif drxp2<=sr:
			drx=x2-x1-2*sx
			dry=y2-y1
		elif dryp1<=sr:
			drx=x2-x1
			dry=y2-y1+2*sy
		elif dryp2<=sr:
			drx=x2-x1
			dry=y2-y1-2*sy
		elif drxp1yp1<=sr:
			drx=x2-x1+2*sx
			dry=y2-y1+2*sy
		elif drxp1yp2<=sr:
			drx=x2-x1+2*sx
			dry=y2-y1-2*sy
		elif drxp2yp1<=sr:
			drx=x2-x1-2*sx
			dry=y2-y1+2*sy
		elif drxp2yp2<=sr:
			drx=x2-x1-2*sx
			dry=y2-y1-2*sy
		qb=2*(drx*dvx+dry*dvy)
		qc=drx**2+dry**2-sr**2
		dt1=(-qb-math.sqrt(qb**2-4*qa*qc))/(2*qa)
		dt2=dt-dt1
		x1+=v1x*dt1
		x2+=v2x*dt1
		y1+=v1y*dt1
		y2+=v2y*dt1
		if dr<=sr:
			dx=x2-x1
			dy=y2-y1
		elif drxp1<=sr:
			dx=x2-x1+2*sx
			dy=y2-y1
		elif drxp2<=sr:
			dx=x2-x1-2*sx
			dy=y2-y1
		elif dryp1<=sr:
			dx=x2-x1
			dy=y2-y1+2*sy
		elif dryp2<=sr:
			dx=x2-x1
			dy=y2-y1-2*sy
		elif drxp1yp1<=sr:
			dx=x2-x1+2*sx
			dy=y2-y1+2*sy
		elif drxp1yp2<=sr:
			dx=x2-x1+2*sx
			dy=y2-y1-2*sy
		elif drxp2yp1<=sr:
			dx=x2-x1-2*sx
			dy=y2-y1+2*sy
		elif drxp2yp2<=sr:
			dx=x2-x1-2*sx
			dy=y2-y1-2*sy
		mag=math.sqrt(dx**2+dy**2)
		ih=dx/mag
		jh=dy/mag
		tc=math.acos(ih)
		if jh<0:
			tc=-tc
		v1rx=v1x*math.cos(tc)+v1y*math.sin(tc)
		v1ry=-v1x*math.sin(tc)+v1y*math.cos(tc)
		v2rx=v2x*math.cos(tc)+v2y*math.sin(tc)
		v2ry=-v2x*math.sin(tc)+v2y*math.cos(tc)
		v1fx=(-dm/sm)*v1rx+(2*m2/sm)*v2rx
		v2fx=(2*m1/sm)*v1rx+(dm/sm)*v2rx
		v1rx=v1fx
		v2rx=v2fx
		v1x=v1rx*math.cos(-tc)+v1ry*math.sin(-tc)
		v1y=-v1rx*math.sin(-tc)+v1ry*math.cos(-tc)
		v2x=v2rx*math.cos(-tc)+v2ry*math.sin(-tc)
		v2y=-v2rx*math.sin(-tc)+v2ry*math.cos(-tc)
		x1+=v1x*dt2
		x2+=v2x*dt2
		y1+=v1y*dt2
		y2+=v2y*dt2
	if x1<-sx:
		x1+=2*sx
	if y1<-sy:
		y1+=2*sy
	if x1>sx:
		x1-=2*sx
	if y1>sy:
		y1-=2*sy
	if x2<-sx:
		x2+=2*sx
	if y2<-sy:
		y2+=2*sy
	if x2>sx:
		x2-=2*sx
	if y2>sy:
		y2-=2*sy
	plt.subplot(2,2,(1,2))
	circle = plt.Circle((x1,y1), radius=r1, fc='xkcd:red')
	plt.gca().add_patch(circle)
	if x1-r1<-sx:
		circle = plt.Circle((x1+2*sx,y1), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	if x1+r1>sx:
		circle = plt.Circle((x1-2*sx,y1), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	if y1-r1<-sy:
		circle = plt.Circle((x1,y1+2*sy), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	if y1+r1>sy:
		circle = plt.Circle((x1,y1-2*sy), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	if x1-r1<-sx and y1-r1<-sy:
		circle = plt.Circle((x1+2*sx,y1+2*sy), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	if x1-r1<-sx and y1+r1>sy:
		circle = plt.Circle((x1+2*sx,y1-2*sy), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	if x1+r1>sx and y1-r1<-sy:
		circle = plt.Circle((x1-2*sx,y1+2*sy), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	if x1+r1>sx and y1+r1>sy:
		circle = plt.Circle((x1-2*sx,y1-2*sy), radius=r1, fc='xkcd:red')
		plt.gca().add_patch(circle)
	circle = plt.Circle((x2,y2), radius=r2, fc='xkcd:cyan')
	plt.gca().add_patch(circle)
	if x2-r2<-sx:
		circle = plt.Circle((x2+2*sx,y2), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	if x2+r2>sx:
		circle = plt.Circle((x2-2*sx,y2), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	if y2-r2<-sy:
		circle = plt.Circle((x2,y2+2*sy), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	if y2+r2>sy:
		circle = plt.Circle((x2,y2-2*sy), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	if x2-r2<-sx and y2-r2<-sy:
		circle = plt.Circle((x2+2*sx,y2+2*sy), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	if x2-r2<-sx and y2+r2>sy:
		circle = plt.Circle((x2+2*sx,y2-2*sy), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	if x2+r2>sx and y2-r2<-sy:
		circle = plt.Circle((x2-2*sx,y2+2*sy), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	if x2+r2>sx and y2+r2>sy:
		circle = plt.Circle((x2-2*sx,y2-2*sy), radius=r2, fc='xkcd:cyan')
		plt.gca().add_patch(circle)
	plt.xlim([-sx,sx])
	plt.ylim([-sy,sy])
	plt.title('2D Elastic Collisions')
	ax=plt.gca()
	ax.set_aspect(1)
	ax.set_facecolor('xkcd:black')
	plt.subplot(223)
	ke1=0.5*m1*(v1x**2+v1y**2)
	ke2=0.5*m2*(v2x**2+v2y**2)
	ke=ke1+ke2
	kea.append(ke/ket)
	ke1a.append(ke1/ket)
	ke2a.append(ke2/ket)
	time.append(frame*dt)
	plt.plot(time,ke1a,lw=2,color='xkcd:red')
	plt.plot(time,ke2a,lw=2,color='xkcd:cyan')
	plt.plot(time,kea,lw=2,color='xkcd:emerald green')
	plt.xlim([0,frps*sec*dt])
	plt.ylim([0,1.2])
	plt.title('Normalized Kinetic Energy')
	ax=plt.gca()
	ax.set_facecolor('xkcd:black')
	plt.subplot(224)
	p1=m1*math.sqrt(v1x**2+v1y**2)
	p2=m2*math.sqrt(v2x**2+v2y**2)
	p1x=m1*v1x
	p1y=m1*v1y
	p2x=m2*v2x
	p2y=m2*v2y
	ptix=p1x+p2x
	ptiy=p1y+p2y
	pti=math.sqrt(ptix**2+ptiy**2)
	p1a.append(p1/pt)
	p2a.append(p2/pt)
	pa.append(pti/pt)
	p1xa.append(p1x/pt)
	p1ya.append(p1y/pt)
	p2xa.append(p2x/pt)
	p2ya.append(p2y/pt)
	pxa.append(ptix/pt)
	pya.append(ptiy/pt)
	list=[p1/pt,p2/pt,pti/pt,p1x/pt,p1y/pt,p2x/pt,p2y/pt,ptix/pt,ptiy/pt]
	maxl=max(list)
	minl=min(list)
	maxlist.append(maxl)
	minlist.append(minl)
	ymax=max(maxlist)
	ymin=min(minlist)
	plt.plot(time,p1a,lw=0.5,color='xkcd:red')
	plt.plot(time,p1xa,lw=0.5,color='xkcd:pink')
	plt.plot(time,p1ya,lw=0.5,color='xkcd:lavender')
	plt.plot(time,p2a,lw=0.5,color='xkcd:cyan')
	plt.plot(time,p2xa,lw=0.5,color='xkcd:light aqua')
	plt.plot(time,p2ya,lw=0.5,color='xkcd:lightish blue')
	plt.plot(time,pa,lw=0.5,color='xkcd:emerald green')
	plt.plot(time,pxa,lw=0.5,color='xkcd:teal')
	plt.plot(time,pya,lw=0.5,color='xkcd:sea green')
	plt.xlim([0,frps*sec*dt])
	plt.ylim([1.2*ymin,1.2*ymax])
	plt.title('Normalized Momentum')
	ax=plt.gca()
	ax.legend([r'$|\vec P_{1}|$',r'$P_{1x}$',r'$P_{1y}$',r'$|\vec P_{2}|$',r'$P_{2x}$',r'$P_{2y}$',r'$|\vec P|$',r'$P_{x}$','$P_{y}$'],labelcolor='w',ncol=3,fontsize='xx-small',frameon=False)
	ax.set_facecolor('xkcd:black')
		
ani=animation.FuncAnimation(fig,run,interval=1,frames=frps*sec)
writervideo = animation.FFMpegWriter(fps=frps)
ani.save('2dcol2_pbc.mp4', writer=writervideo)
plt.show()

