
import numpy as npy
import sys
import math
import wx
import datetime
import conversion

def dp(a,b):
   u=0.0
   if (not (len(a)==len(b))):
      sys.exit(1)
   for i in range(len(a)):
      u+=a[i]*b[i]
   return u
      

def der(x,y):
   d=[]
   for i in range(len(y)):
      if (i==0):
         v=(y[1]-y[0])/(x[1]-x[0])
      elif (i==(len(y)-1)):
         v=(y[len(y)-1]-y[len(y)-2])/(x[len(y)-1]-x[len(y)-2])
      else:
         v=(y[i+1]-y[i-1])/(x[i+1]-x[i-1])
      d.append(v)

   return d
   
def plateaupoint(x,y,pv,eps):
   for i in range(len(y)):
      if (abs(pv-y[i])>eps):
         return x[i] 

   return -999

def gettrunccurves(x,y, miny, maxy, starttim, endtim, ntrunc, ipass, npass):

   count=0
   ncountshow=100
   n=5+ncountshow
   dlg = wx.ProgressDialog('Progress', 'Computing truncation statistics ('+str(ipass)+'/'+str(npass)+')', n)

   # compute statistics
   thisv=[]
   for i in range(len(x)):
      tim=x[i]%1.0
      if (tim>=starttim and tim<=endtim):
         if (y[i]>=miny and y[i]<=maxy):
            thisv.append(y[i])
   count+=1
   dlg.Update(count)

   m=npy.mean(thisv)
   s=npy.std(thisv)
   thisv.sort()

   #print "m="+str(m)
   #print "s="+str(s)

   ival=0
   vblom=[]
   vprob=[]
   for val in thisv:
      ival+=1
      vblom.append((float(ival)-0.375)/(float(len(thisv)+0.25)))
      if (ival==len(thisv)):
         vprob.append(0.99999)
      else:
         vprob.append(float(ival)/float(len(thisv)))
   maxval=max(thisv)
  
   count+=1
   dlg.Update(count)

   means=[]
   stdevs=[]
   epsilons=[]
   Qbars=[]
   sigmas=[]
   iupd=0
   lastupd=0
   if (len(thisv)>0):
      for itrunc in range(ntrunc):
         iupd=int(ncountshow*float(itrunc)/float(ntrunc))
         if (iupd>lastupd):
            count+=1
            dlg.Update(count)
            lastupd=iupd
         trunc=float(itrunc)/float(ntrunc-1)
         ar=[]
         for val in thisv:
            v=(val-maxval*trunc)
            if (v<0.0):
               v=0.0
            ar.append(v)
         means.append(npy.mean(ar))
         stdevs.append(npy.std(ar))
         epsilons.append((maxval*trunc-m)/s)
         sar=sum(ar)
         Qbars.append((1.0/(float(len(thisv))*s))*sar)
         sqa=float(len(thisv))*dp(ar,ar)-sar*sar
         #print 'sqa='+str(sqa)
         sigmas.append((1.0/(float(len(thisv))*s))*math.sqrt(sqa))

   count+=1
   dlg.Update(count)
   derQbar=der(epsilons,Qbars)
   count+=1
   dlg.Update(count)
   dersigmas=der(epsilons,sigmas)
 
   dlg.Destroy()

   return Qbars, sigmas, derQbar, dersigmas, m, s, epsilons

def subsetx(x,y, minx, maxx):

   ox=[]; oy=[]
   margin=1e-6
   #margin=1e-9
   minxm=minx-margin
   maxxm=maxx+margin

   n=len(x)
   for i in range(n):
       lx=x[i]
       if (lx>=minxm):
          if (lx<=maxxm):
            ox.append(x[i])
            oy.append(y[i])

   return ox,oy

def subsetx2(x,y,minx,maxx):
    #print('minx,maxx=%f %f' % (minx,maxx))
    # find bounds:
    bounds=[0,-1]
    lo=0
    up=len(x)-1
    margin=1e-6
    
    for i in range(2):
        # lower and upper bound
        done=False
        
        if i==0:
            if minx<=x[lo]:
                bounds[0]=lo
                done=True
        elif i==1:
            if maxx>=x[up]:
                bounds[1]=up
                done=True
        while not done:
            #print('lo,up=%i %i' % (lo,up))
            ul=[up,lo]
            mi=int((lo+up)/2)
            if up-lo==0:
                #print('a%i' %i)
                bounds[i]=mi
                done=True
            elif up-lo==1:
                #print('b%i' %i)
                bounds[i]=ul[i]
                done=True
            else:
                #print('x=%i' % x[mi])
                if abs(x[mi]-minx)<margin:
                    #print('c%i' %i)
                    bounds[i]=mi
                    done=True
                elif i==0:
                    # lower bound
                    if x[mi]<minx:
                        #print('d%ia' %i)
                        lo=mi
                    else:
                        #print('e%ia' % i)
                        up=mi
                elif i==1:
                     # upper bound
                    if x[mi]>maxx:
                        #print('d%ib' %i)
                        up=mi
                    else:
                        #print('e%ib' % i)
                        lo=mi                   
        lo=bounds[0]
        up=len(x)-1
                    
    #print('bounds='+str(bounds))
    ist,ien=bounds
    return x[ist:(ien+1)], y[ist:(ien+1)]
            
def test_subsetx2():
    x=[1,2,3,4,5,6]
    y=[11,12,13,14,15,16]

    ranges=[[0,7],[1,6],[2,5],[3,4],[4,4],[2.5,4.5], [0.5,7.5]]

    print ('x=%s y=%s' % (str(x),str(y)))
    for ra in ranges:
        print('subsetx2(%f,%f)=%s' % (ra[0],ra[1],str(subsetx2(x,y,ra[0],ra[1]))))            
            
        


def subset(x,y,minx, maxx, miny, maxy, mint, maxt):
   #return subset_impl1(x,y, minx, maxx, miny, maxy, mint, maxt)
   return subset_impl2(x,y, minx, maxx, miny, maxy, mint, maxt)

def subset_impl2(x,y, minx, maxx, miny, maxy, mint, maxt):

   ox=[]; oy=[]; oi=[]
   margin=1e-6
   #margin=1e-9
   minxm=minx-margin
   maxxm=maxx+margin
   minym=miny-margin
   maxym=maxy+margin
   mintm=mint-margin 
   maxtm=maxt+margin

   n=len(x)
   a=[i for i in range(n)]   
   if mintm==0.0 and maxtm==1.0:    
      for i in a:
          lx=x[i]
          ly=y[i]
          if (lx>=minxm):
             if (lx<=maxxm):
                if (ly>=minym):
                   if (ly<=maxym):
                      ox.append(x[i])
                      oy.append(y[i])
                      oi.append(i)
   else:
      for i in a:
          lx=x[i]
          ly=y[i]
          if (lx>=minxm):
             if (lx<=maxxm):
                if (ly>=minym):
                   if (ly<=maxym):
                      t=x[i]%1.0 
                      if (t>=mintm and t<=maxtm):                      
                          ox.append(x[i])
                          oy.append(y[i]) 
                          oi.append(i)
   return ox,oy,oi

def subset_impl1(x,y, minx, maxx, miny, maxy, mint, maxt):

   ox=[]; oy=[]; oi=[]
   margin=1e-9
   for i in range(len(x)):
      if (x[i]>=minx-margin):
         if (x[i]<=maxx+margin):
            if (y[i]>=miny-margin):
               if (y[i]<=maxy+margin):
                  t=x[i]%1.0 
                  if (t>=mint-margin and t<=maxt+margin):
                     ox.append(x[i])
                     oy.append(y[i])
                     oi.append(i)

   return ox,oy,oi

def basepattern(choice, nsample, mean, sd, minx,maxx,mint, maxt, trunc):

   ox=[]; oy=[]; ot=[]
   if (choice=="normal"):
      dlg = wx.ProgressDialog('Progress', 'Generating baseline sample set', 1)

      oy=[]
      for isample in range(nsample):
         loy=npy.random.normal(mean,sd)
         if (loy>=trunc):
            oy.append(loy)
         else:
            oy.append(trunc)
      #for i in range(nsample):
          #x=mint+i*(maxt-mint)/float(nsample-1)
          #t=x
          #ox.append(x+7.0)
          #ot.append(t+7.0)
      nday=int(maxx-minx)
      nperday=int(nsample/nday)+1
      ntot=0
      for iday in range(nday):
         for jsam in range(nperday):
            ntot+=1
            if ntot>nsample:
               break
            t=mint+jsam*(maxt-mint)/float(nperday)
            x=float(int(minx)+iday)+t
            ox.append(x)
            ot.append(t)
            

      dlg.Update(1)
      dlg.Destroy()
   else:
      print ('mystat.basepattern error choice='+str(choice))
      sys.exit(1)

   return ox,oy,ot
 
def dates2weekdays(x):
   # reference day: May 2, 2011 is a Monday
   ox=[]
   for i in range(len(x)):
      d=conversion.xldate_as_datetime(x[i],0) 
      wd=d.weekday()
      ox.append(conversion.date2exceldate(datetime.datetime(2011,5,2+wd,d.hour,d.minute,d.second)))

   print( 'dates2weekdays: min,max='+str([min(ox),max(ox)]))

   return ox


      
def resample(dat,n):

   ndat=len(dat)

   outdat=[]
   a=[i for i in range(n)]
   fni=1.0/float(n)
   for i in a:
      j=i*float(ndat)*fni
      j1=int(math.floor(j)); j2=int(math.floor(j)+1)
      # linear interpolation
      if j2>=ndat:
         v=dat[j1]
      else:
         v=dat[j1]+(j-float(j1))*(dat[j2]-dat[j1])
      outdat.append(v)
  
   return outdat

def aggregate(x,y,t,blocklen):
   print ('blocklen='+str(blocklen))
   x0=min(x); x1=max(x)
   n=int((x1-x0)/blocklen+1)
   px=[]
   py=[]
   pt=[]
   for i in range(n):
      px.append([])
      py.append([])
      pt.append([])
   #print 'px='+str(px)
   #print 'py='+str(py)
   #print 'pt='+str(pt)
   for i in range(len(x)):
      j=int((x[i]-x0)/blocklen)
      #print 'i,j='+str([i,j])
      #print 'add '+str(x[i])+' to '+str(px[j])
      px[j].append(x[i])
      py[j].append(y[i])
      pt[j].append(t[i])
   #print 'pt='+str(pt)
       
   #print 'px='+str(px)
   #print 'py='+str(py)
   #print 'pt='+str(pt)
   ox=[]; oy=[]; ot=[]
   for i in range(n):
      #print 'i='+str(i)+'/'+str(n)
      if len(px[i])>0:
         ox.append(x0+i*blocklen)
         ot.append(min(pt[i])) 
         oy.append(npy.mean(py[i]))

   #print 'ox='+str(ox)
   #print 'oy='+str(oy)
   #print 'ot='+str(ot)

   #print 'input range:'
   #print 'x:'+str([min(x),max(x)])
   #print 'y:'+str([min(y),max(y)])
   #print 't:'+str([min(t),max(t)])
   #print 'n='+str(len(x))
   #print 'output range:'
   #print 'x:'+str([min(ox),max(ox)])
   #print 'y:'+str([min(oy),max(oy)])
   #print 't:'+str([min(ot),max(ot)])
   #print 'n='+str(len(ox))
   mindx=1e33
   for i in range(len(ox)-1):
      dx=ox[i+1]-ox[i]
      mindx=min(mindx,dx)
   #print 'min dx='+str(mindx)


   return ox,oy,ot
