# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 07:52:29 2021

@author: thienpe
"""

import sys
import numpy 
import conversion
import mystat
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.optimize
import scipy.stats
import os
import datetime
import operator
import itertools
import time
import pickle

tool='clusterCFPD2'
version='0.57 020210811'

# MAJOR CHANGES
# 0.57: * prepare for CFPD comparison of same period in multiple years
# 0.56: * use R2 instead of se_b for finding the best matching period
#       * check_nonzero_ts: do not allow any zero values


# disable matplotlib interactive mode
plt.ioff()

def main():
    
    single();
    #multi();

def single():
    # define datasets     

    infile='testdata/synthetic testdata clusterCFPD subset 5d.csv'
    
    perlen=1  # comparison period length
    nbar=4  # number of bars in mCFPD diagram
    maxnclust=1  # maximum number of subclusters
    suffix='_test' # output directory suffix
    
    ret, summary, istart, iend=execute(infile,perlen,nbar,maxnclust, suffix=suffix)
    
def multi():
    perlen=7 # comparison period length
    nbar=5 # number of bars in mCFPD diagram
    maxnclust=2    # maximum number of subclusters
    compare=True


    infiles=['data/column_data_2015_5s_jan-mrt.csv',
              'data/column_data_2016_5s_jan-mrt.csv',
              'data/column_data_2017_5s_jan-mrt.csv',
              'data/column_data_2018_5s_jan-mrt.csv']
        
    
    suffix='_UK'
    nowstring=datetime.datetime.now().strftime('%Y%m%d_%H%M')
    prefix='run_%s_' % nowstring
    
    outlist=[]
    periods=[]
    summarylist=[]
    for infile in infiles:
        ret, summary, istart, iend=execute(infile,perlen,nbar,maxnclust, suffix=suffix, prefix=prefix)
        outlist.append(ret)
        summarylist.append(summary)
        periods.append([istart,iend])
        
    if compare:
        # get DMA labels from the first infile
        success,xdat,ydat,DMAs=loadandparsecsv_multi(infiles[0])

        plotdata=[]
        
        print('generating comparison CFPD diagrams')
        print ('for DMAs: '+str(DMAs))

        for DMA in DMAs: 
            for ii in range(len(periods)-1):
                istart1,iend1=periods[ii]
                istart2,iend2=periods[ii+1]
                infile1=infiles[ii]
                infile2=infiles[ii+1]
                success,xdat1,ydat1,labels1=loadandparsecsv_multi(infile1)
                success,xdat2,ydat2,labels2=loadandparsecsv_multi(infile2)

                p1=conversion.exceldate2stringd(istart1)+' - '+ conversion.exceldate2stringd(iend1)
                p2=conversion.exceldate2stringd(istart2)+' - '+ conversion.exceldate2stringd(iend2)
                llab=DMA+':'+ p2 + ' vs. '+p1
                plotdata.append([llab,xdat1,ydat1[labels1.index(DMA)],istart1, iend1,
                                      xdat2,ydat2[labels2.index(DMA)],istart2, iend2])
                

        fname='output/cfpd_period_comparison_'+nowstring+'.png'
        print('outputting to file %s' % fname)
        cfpdplot3(plotdata, fname=fname, title=fname )        
    
    
    
    
    fname='output/'+prefix+'_results.txt'
    print ('writing combined results to "%s"' % fname)
    writelist(outlist,fname)
    
    fname='output/'+prefix+'_summary.txt'
    print ('writing summary of combined results to "%s"' % fname)
    writelist(summarylist,fname)        
    
def execute(infile,perlen,nbar,maxnclust, suffix='', prefix=''):
    
    nowstring=datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')
    outdir='output/'+prefix+tool+'_'+'_'+nowstring+suffix
    print('writing output to %s' % outdir)
    os.mkdir(outdir)
        
    inputs=[]
    inputs.append('tool='+str(tool))
    inputs.append('code version='+str(version))
    inputs.append('infile='+str(infile))
    inputs.append('perlen='+str(perlen))
    inputs.append('maxnclust='+str(maxnclust))
    writelist(inputs,outdir+'/inputs.txt')
    
    
    res=find_best_correspondence_multi2(infile,perlen)
    istart=res[0][0]
    iend=res[0][1]
    print('found range %i-%i' % (istart,iend))
    
    success,xdat,ydat,labels=loadandparsecsv_multi(infile)
    
    # plot time series and range box
    plot_timeseries(xdat,ydat,labels,outdir,istart,iend)

    outputs=[]
    outdat=[]
    meanflows={}
    if success: 
        for nclust in range(1,maxnclust+1):
            meanflows[nclust]={}
            print('================================================')
            print('testing nclust=%i' % nclust)
            outdat.append({})
            # generate all possible clusters
            clusters=generate_clusters2(labels, nclust)
            #print ('nclust=%i clusters: %s' % (nclust,str(clusters)))
            
            # find best fitting clustering
            res=evaluate_clusters(clusters, xdat, ydat, labels, perlen, ts_start=istart, ts_end=iend)
            #print('EC res=')
            #for item in res:
                #print(str(item))
            best=res[0][0] # clustering for best results
            mnfs=res[0][3]
            meanfs=res[0][4]
            print('best='+str(best))
            print('mnfs='+str(mnfs))
            print('meanfs='+str(meanfs))
            ii=-1
            for cluster in best:
                for lab in cluster:
                    ii+=1
                    mf=meanfs[ii]
                    #print ('nclust, lab= %s, %s' % (str(nclust), str(lab)))
                    meanflows[nclust][lab]=mf
                    
            
            # evaluate subclusters
            isubclust=-1
            letlabels=[str(chr(i)) for i in range(97,123)]
            
            outdat[-1]['best']=best
            outdat[-1]['values']={}

            for subclust in best:
                isubclust+=1
                runlabel='nclust=%i%s'% (nclust, letlabels[isubclust])
                outputs.append(runlabel)

                # location of SSR curve bend and summed MNF
                clustdat=[ runlabel, infile, istart, iend, subclust]
                res=find_range2(clustdat)
                #print('FR res='+str(res))
                [bgl_start, bgl_end, curve_x, curve_y, dataset]=res
                
                # create CFPD plot of all combinations in subcluster
                fname=outdir+'/cfpdplots_'+runlabel+'.png'
                plotdata=[]
                for llab in subclust:
                    plotdata.append([llab,xdat,ydat[labels.index(llab)],istart, iend])
                cfpdplot2(plotdata, fname=fname, title=runlabel )
                
                # create bg leakage figure and table
                #print('dataset='+str(dataset))
                [llabels,lvalues]=plot_mcfpd_range(dataset, bgl_start, bgl_end, curve_x, curve_y, outdir, runlabel,nbar=nbar)
                for i in range(len(llabels)):
                    lstr=llabels[i]+' : '+str(lvalues[i][0])
                    outputs.append(lstr)
                    print(lstr)
                    outdat[-1]['values'][llabels[i]]=lvalues[i][0]
                plot_mcfpd_range_ranks(dataset, bgl_start, bgl_end, curve_x, curve_y, outdir, runlabel)
   
        # create overview table
        letlabels=[str(chr(i)) for i in range(97,123)]       
        
        # order of labels based on clustering for largest nclust
        orderedlabels=[]
        for subclust in outdat[-1]['best']:
            orderedlabels+=subclust
      
        # values in order of labels for nclust=1
        orderedvalues=[[]]
        for label in orderedlabels:
            orderedvalues[0].append(('a', outdat[0]['values'][label]))
 
        # cluster numbers/labels
        clusternumbers=[]
        for od in outdat:
            isc=-1
            ldict={}
            for sc in od['best']:
                isc+=1
                for item in sc:
                   ldict[item]=letlabels[isc]
            clusternumbers.append(ldict)

       
        # values in order of labels for nclust>1
        for ii in range(1,maxnclust):
            orderedvalues.append([])
            for label in orderedlabels:
                v1=clusternumbers[ii][label]
                try:
                    v2=outdat[ii]['values'][label]
                except KeyError:
                    v2=None
                tup=(v1,v2)
                orderedvalues[-1].append(tup)
        
        # correction parameters for nclust>1 w.r.t. nclust=1        
           
        corrections=[[]]        
        for nclust in range(1,maxnclust):          
            corrvalues={}
            for subclust in outdat[nclust]['best']:
                llist=[]
                llist2=[]
                for label in subclust:
                    llist.append(outdat[0]['values'][label])
                    try:
                        llist2.append(outdat[nclust]['values'][label])
                    except KeyError:
                        #llist2.append(None)
                        pass
                corrval=numpy.mean(llist)-numpy.mean(llist2)
                for label in subclust:
                    corrvalues[label]=(clusternumbers[nclust-1][label],corrval)
            corrections.append(corrvalues) 
            
        correctedvalues=[[]]
        correctedvalues[0]=orderedvalues[0][:]
        for ii in range(1,maxnclust):
            correctedvalues.append([])
            for label in orderedlabels:
                cv=corrections[ii][label][1] # [0] contains the subcluster label
                v1=clusternumbers[ii][label]
                try:
                    v2=outdat[ii]['values'][label]+cv
                except KeyError:
                    v2=None
                tup=(v1,v2)
                correctedvalues[-1].append(tup)
        
        sumstr=''
        lstr='----------------------------------------------------------\n'
        lstr+=infile+' (analysis period = %i days)\n' % perlen
        lstr='----------------------------------------------------------\n'
        lstr+='\n\n%s - %s\n' % (conversion.exceldate2stringd(istart),conversion.exceldate2stringd(iend-1))
        sumstr+='\n\n%s - %s\n' % (conversion.exceldate2stringd(istart),conversion.exceldate2stringd(iend-1))
        # raw analysis output
        lstr+='\nraw analysis results \n      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\n'
        for nclust in range(maxnclust):
            lstr+='ncl=%2i' % (nclust+1)
            for tup in orderedvalues[nclust]:
                #print (str(tup)+' ('+str(type(tup))+')')
                if tup[1]==None:
                    lstr+=' | %1.1s   None' % (tup[0])
                else:
                    lstr+=' | %1.1s %6.2f' % (tup[0],tup[1])
            lstr+='\n'
 
        # corrected values
        lstr+='\nstacked analysis results - minimum background leakage \n'
        lstr+='      '
        sumstr+='\nstacked analysis results - minimum background leakage \n'
        sumstr+='      '        
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
            sumstr+=' | %8.8s' %label
        lstr+='\n'
        sumstr+='\n'
        for nclust in range(maxnclust):
            lstr+='ncl=%2i' % (nclust+1)
            sumstr+='ncl=%2i' % (nclust+1)
            for tup in correctedvalues[nclust]:
                if tup[1]==None:
                    lstr+=' | %1.1s   None' % (tup[0])
                    sumstr+=' | %1.1s   None' % (tup[0])
                else:
                    lstr+=' | %1.1s %6.2f' % (tup[0],tup[1])
                    sumstr+=' | %1.1s %6.2f' % (tup[0],tup[1])
            lstr+='\n'
            sumstr+='\n'

        # mean flows, percentages, ranking
        lstr+='\nmean flows\n'
        lstr+='      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\n'
        for nclust in range(maxnclust):
            lstr+='      '
            for label in orderedlabels:
                lstr+=' | %8.3f' % meanflows[nclust+1][label]
            lstr+='\n'
        lstr+='\nminimum background leakage (% of mean flow)\n'
        lstr+='      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\n'
        labs=[]
        perc=[]
        for nclust in range(maxnclust):
            lstr+='ncl=%2i' % (nclust+1)
            ii=-1
            perc.append([])
            labs.append([])
            for tup in correctedvalues[nclust]:
                ii+=1
                if tup[1]==None:
                    lstr+=' | %1.1s   None' % (tup[0])
                    perc[-1].append(-1.0)
                    labs[-1].append(-1.0)
                else:
                    v=100.0*tup[1]/meanflows[nclust+1][orderedlabels[ii]]
                    lstr+=' | %1.1s %4.1f %%' % (tup[0],v)
                    perc[-1].append(v)
                    labs[-1].append(tup[1])

            lstr+='\n'            
       
        lstr+='\npriotization (absolute | relative)\n'
        lstr+='      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\n'
        for nclust in range(maxnclust):
            lstr+='ncl=%2i' % (nclust+1)
          
            print ('labs='+str(labs))
            la=list(numpy.argsort(labs[nclust]))
            lp=list(numpy.argsort(perc[nclust]))
            print ('la='+str(la))
            la.reverse()
            lp.reverse()
            print ('->la='+str(la))
            la=[la.index(j) for j in range(len(la))] # from sorting index to rank
            lp=[lp.index(j) for j in range(len(lp))] # from sorting index to rank
            print ('->la='+str(la))
            la=[(v+1) for v in la]
            lp=[(v+1) for v in lp]
            print ('->la='+str(la))
            for ii in range(len(orderedlabels)):
                    lstr+=' | %2i / %2i ' % (la[ii],lp[ii])
            lstr+='\n'            
                        
            

           
        print(lstr)
        outputs.append(lstr)
        
        writelist(outputs,outdir+'/outputs.txt')
        
        return lstr, sumstr, istart, iend
    else:
        return '\nfailed\n'
            
 
def plot_mcfpd_range(datasets, l1_loc, summnf, xc, yc, outdir, runlabel, nbar=5):
    
    if not len(datasets)>1:
        print('not enough datasets for %s' % runlabel)
        return [],[]
    
    dx=(summnf-l1_loc)/float(nbar-1)
    bar_locs=[]
    for ii in range(nbar):
        bar_locs.append(l1_loc+float(ii)*dx)
        
    
    
    #print ('l1/min/l2 at %s' % str(bar_locs))


    returnvalues=[[],[]]    
    rvalues=[]
    values={}
    for v_sum_NRW in bar_locs: 
        rvalues.append(v_sum_NRW)
        success, res, ssr, dsids=analyze(datasets,v_sum_NRW, figures=False)

        for ival in range(1,len(res)+1):
            dsid=dsids[ival]
            label='b_0_%i' % ival
            val=res[label]
            
            if not dsid in values:
                values[dsid]=[val]
            else:
                values[dsid].append(val)
                    
        #print('values=%s' % str(values))
        
    # plot stacked bar charts
    lllabels=rvalues
    fig,ax=plt.subplots()
    width=0.1*summnf
    
    bottoms=[0.0]*len(bar_locs)
    
    rects=[]
    for i in range(len(values)): 
    
        lrects=ax.bar(lllabels,values[dsids[i+1]],width,label=dsids[i+1], bottom=bottoms)
        returnvalues[0].append(dsids[i+1])
        returnvalues[1].append(values[dsids[i+1]])
        for rect in lrects:
            rects.append(rect)

        for j in range(len(bar_locs)):
            bottoms[j]+=values[dsids[i+1]][j]
            
        autolabel(ax,rects)
    
    ax2 = ax.twinx()  # set up the 2nd 
    #print('ssrs='+str(ssrs))
    
    #yc=list(zip(*sorted(zip(yc,xc))))[1]
    #xc=sorted(xc)
    clist=[]
    for i in range(len(xc)):
        clist.append([xc[i],yc[i]])
    clist=sorted(clist,key=operator.itemgetter(0))
    xc=[]
    yc=[]
    for item in clist:
        xc.append(item[0])
        yc.append(item[1])

        
    
    ax2.plot(xc,yc,':k')
    #ax2.plot(xc,yc,'.k')
    
    ax.set_ylabel('background leakage (m3/h)')
    ax2.set_ylabel('log10(SSR)')
    ax.set_title('')
    ax.set_xlabel('total background leakage estimate (m3/h)')
    ax.legend()

    fname=outdir+'/clustCFPD_%s.png' % runlabel
    fig.savefig(fname, dpi=300)            
    #fig.savefig(fname, figsize=(9,6), dpi=300)            
    plt.close()
        
    olist=[]
    for i in range(len(xc)):
        olist.append([xc[i],yc[i]])
    fname=outdir+'/SSR_curve_%s.txt' % runlabel
    writelist(olist,fname)  


    return returnvalues            

def plot_mcfpd_range_ranks(datasets, l1_loc, summnf, xc, yc, outdir, runlabel, npoint=10):
    
    if not len(datasets)>1:
        print('not enough datasets for %s' % runlabel)
        return [],[]
    
    dx=(summnf-l1_loc)/float(npoint-1)
    locs=[]
    for ii in range(npoint):
        locs.append(l1_loc+float(ii)*dx)
        
    
    xs=[]
    ys=[[]]*(len(locs)) 
    #ys=[[]]*(len(datasets)-1) 
    iloc=-1
    for v_sum_NRW in locs: 
        iloc+=1
        xs.append(v_sum_NRW)
        success, res, ssr, dsids=analyze(datasets,v_sum_NRW, figures=False)

        vals=[]
        for ival in range(1,len(res)+1):
            #dsid=dsids[ival]
            label='b_0_%i' % ival
            val=res[label]
            vals.append(val)
            
        ranks=[len(vals)-(sorted(vals).index(x)) for x in vals]
        print('vals=%s' % str(vals))
        print('ranks=%s' % str(ranks))
        #ranks=scipy.stats.rankdata(vals)
        for i in range(len(ranks)):
            ys[iloc].append(ranks[i])
                    
        #print('values=%s' % str(values))
        
    # plot rank charts
    fig,ax=plt.subplots()
    
    
    for y in ys:
        ax.plot(xs,ys)
    
          
  
   
    ax.set_ylabel('Priority')
    ax.set_title('')
    ax.set_xlabel('total background leakage estimate')
    ax.legend()

    fname=outdir+'/clustCFPD_%s_ranks.png' % runlabel
    fig.savefig(fname, dpi=300)            
    #fig.savefig(fname, figsize=(9,6), dpi=300)            
    plt.close()
        


    return    
            
        
def plot_timeseries(xdat,ydat,labels,outdir, bx1,bx2):
    # create plot for the complete time series with 
    axs=[]
    nts=len(ydat)
    ysize=nts*2
    figsize=(15, ysize)
    fs2=12
    
    
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs0 = gridspec.GridSpec(nts,1, figure=fig)
    
    letlabels=[str(chr(i)) for i in range(97,123)]
    letlabels+=['a'+str(chr(i)) for i in range(97,123)]
    
    start=min(xdat)
    end=max(xdat)
    
    for its in range(nts):
        #print ('its=%i n1=%i n2=%i' % (its,len(labels),len(ydat)))
        axs.append(fig.add_subplot(gs0[its,:]))
        axs[-1].plot(xdat,ydat[its],'-k',linewidth=0.5,label=labels[its])
        axs[-1].set_xlim([start,end])
        axs[-1].set_ylabel('Q', fontsize=fs2)    
        axs[-1].tick_params(labelbottom=False)
        axs[-1].set_title(letlabels[its]+') '+labels[its], fontsize=10,loc='left')
       
        # add boxes for the selected subperiods
        maxv=max(ydat[its])
        lx=[bx1,bx1,bx2,bx2,bx1]
        ly=[0,0.95*maxv,0.95*maxv,0,0]
        axs[-1].plot(lx,ly,'-b')
    
    #fig.tight_layout(pad=0.0)
    #plt.show()
    fname=outdir+'/timeseries.png'
    fig.savefig(fname, dpi=300)  
    # fig.savefig(fname, figsize=figsize, dpi=300)  
    plt.close()


def evaluate_clusters(clusters, xdat, ydat, labels, perlen, ts_start=None, ts_end=None):

    icluster=0
    ncluster=len(clusters)
    tprev=None
    looptimes=[]
    results=[]
    
    for cluster in clusters:
        if tprev==None:
            est='unknown'
            tprev=time.time()
        else:
            now=time.time()
            looptime=now-tprev
            looptimes.append(looptime)
            mlooptime=numpy.mean(looptimes)
            tprev=now
            nrem=ncluster-icluster
            est=datetime.timedelta(seconds=mlooptime*nrem)
        print('\rprocessing cluster %i/%i estimated remaining time=%s' % (icluster, ncluster, str(est)), end="", flush=True)
        icluster+=1
        
        isubcluster=0
        cluster_se_b=[]
        mnfs=[]
        meanfs=[]
        for subcluster in cluster:
            isubcluster+=1

         
            if len(subcluster)>1:
                lydat=[]
                for lab in subcluster:
                    lydat.append(ydat[labels.index(lab)]) 
                res=performance(xdat, lydat, ts_start, perlen)
                #print('res='+str(res))
                cluster_se_b+=res[8] 
                mnfs+=list(res[9])
                # print('res[10]='+str(res[10])+' '+str(type(res[10])))
                # print('meanfs='+str(meanfs)+' '+str(type(meanfs)))
                meanfs+=list(res[10])
            elif len(subcluster)==1:
                #print('ydat='+str(ydat))
                #print('subcluster='+str(subcluster))
                lydat=ydat[labels.index(subcluster[0])]
                mnfs+=[min(lydat)]
                meanfs+=[numpy.average(lydat)]
                
        mcseb=numpy.mean(cluster_se_b)
        #print ('cluster_se_b=%s (%s)' % (str(cluster_se_b), str(type(cluster_se_b))))
        results.append([cluster,mcseb, cluster_se_b, mnfs, meanfs])
        #print ('results='+str(results))
    sresults=sorted(results,key=operator.itemgetter(1)) 
    print()
    
    return sresults
    
        
    
        
def autolabel(ax,rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        hs='%.1f' % height
#        ax.annotate('{}'.format(height),
        ax.annotate(hs,
                    xy=(rect.get_x() + rect.get_width() / 2, rect.get_y() + rect.get_height() / 2),
                    xytext=(0, 0),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

def printl(logfile,lstr):
    if logfile==None:
        #print(lstr)
        pass
    else:
        f=open(logfile,'a')
        f.write(lstr+'\n')
        f.close()        

def writelist(llist,fname):
    f=open(fname,'w')
    for item in llist:
        f.write(str(item)+'\n')
    f.close()
            
def generate_timeseries(x,y,a,b,noise_level,fname):

    outy=[]
    
    for ly in y:
       outy.append(a*ly+b+numpy.random.normal(0.0,noise_level*ly))
       

    lstr='' 
    for i in range(len(x)):
        lstr+='%f %f\n' % (x[i],outy[i]) 
    f=open(fname,'w')
    f.write(lstr)
    f.close()
    
    
def generate_timeseries2(x,y,a,b,error_model,ep1,ep2,fname):
    outy=[]
    
    for ly in y:
       outy.append(a*ly+b)
       
    if error_model=='Gaussian':
        for i in range(len(outy)):
            outy[i]+=ep2+numpy.random.normal(0.0,ep1*outy[i])
    else:
        print('error: error model %s undefined' % error_model)
        sys.exit(1)


    lstr='' 
    for i in range(len(x)):
        lstr+='%f %f\n' % (x[i],outy[i]) 
    f=open(fname,'w')
    f.write(lstr)
    f.close()   

def mean_flow(fname):
    success,x,y,t=loadandparsecsv(fname)
    if success:
        return numpy.mean(y)
    else:
        return None


    
def analyze(datasets, sum_NRW, figures=False, method='L-BFGS-B',logfile=None):
    # load data
    flowdata={}
    for dataset in datasets:
        key=dataset[0]
        #print('dataset='+str(dataset)+' ('+str(type(dataset))+')')
        if type(dataset[1])==str:
            # filename
        
            success,x,y,t=loadandparsecsv(dataset[1])
            if not success:
                print('failed to load data for key=%s' % key)
        elif type(dataset[1])==list:
            x=dataset[1][0]
            y=dataset[1][1]
        else:
            print('analyse: dataset error')
            sys.exit(1)
                
            
        flowdata[key]=[x,y]    
    
    
    
    dsids=['hypothetical base pattern'] + [a[0] for a in datasets]  #  hypothetical undisturbed pattern id=0
    #print ('dsids=%s'% str(dsids))
    

        
    n=len(dsids)     
    a_mat=numpy.zeros((n,n))
    b_mat=numpy.zeros((n,n))
    #sd_mat=numpy.zeros((n,n))

    if figures:
        fig, axs = plt.subplots(n-1,n-1)
    

    # determine CFPD fits
    for i in range(1,n):
        ilabel=dsids[i]
        for j in range(1,n):
            jlabel=dsids[j]
            #print ('CFPD %i:%s vs %i:%s' % (i,ilabel,j,jlabel))
            a,b,fiterrparams,x,y=get_cfpd_fit(flowdata[ilabel], datasets[i-1][2],datasets[i-1][3],
                                       flowdata[jlabel], datasets[j-1][2],datasets[j-1][3],
                                       ilabel,jlabel)
            if figures:
                fs1=5
                fs2=4
                if not i==j:
                        x2=[x[0],x[-1]]
                        y2=[a*x[0]+b, a*x[-1]+b]

                        axs[i-1,j-1].plot(x2,y2,'--b')                   
                        axs[i-1,j-1].plot(x,y, 'k')
                        axs[i-1,j-1].set_title('%s - %s' % (ilabel,jlabel), fontsize=fs1)

                        axs[i-1][j-1].tick_params(axis='both', which='major', labelsize=fs2)
                        axs[i-1][j-1].tick_params(axis='both', which='minor', labelsize=fs2)
                else:
                        axs[i-1][j-1].tick_params(axis='both', which='major', labelsize=0)
                        axs[i-1][j-1].tick_params(axis='both', which='minor', labelsize=0)


            a_mat[i][j]=a
            b_mat[i][j]=b
            #sd_mat[i][j]=fiterrparams
            # sd_a=fiterrparams['sd_a']
            # sd_b=fiterrparams['sd_b']
            #print ('[%i,%i] : a=%5.3f b=%5.3f sd_a=%5.3e sd_b=%5.3e' % (i,j,a,b,sd_a,sd_b))

    if figures:
        plt.setp(axs.flat, xticks=[], yticks=[])
        fig.tight_layout(pad=0.0)
        plt.show() 


           
    # set up nonlinear system of equations
    coef_idx={}
    coef_label=[]
    idx=-1
    # a coefficients
    for i in range(n):  
        for j in range(n):
           idx+=1
           label='a_%i_%i' % (i,j)
           coef_idx[label]=idx
           coef_label.append(label)
           #print ('index %i coef %s' % (idx, label))
    # b coefficients
    for i in range(n):
        for j in range(n):
           idx+=1
           label='b_%i_%i' % (i,j)
           coef_idx[label]=idx           
           coef_label.append(label)
           #print ('index %i coef %s' % (idx, label))
           
    ncoef=len(coef_idx)
    neq=n*n*n
    v_vec=numpy.zeros((ncoef))
    v_mat=numpy.zeros((ncoef,1))
    C_mat=numpy.zeros((neq,ncoef))
    D_mat=numpy.zeros((neq,ncoef))
    E_mat=numpy.zeros((neq,ncoef))
    F_mat=numpy.zeros((neq,ncoef))
    #res_vec=numpy.zeros((ncoef))
    res_mat=numpy.zeros((ncoef,1))
    
    #print ('n=%i ncoef=%i' % (n,ncoef))
    
    # free parameters
    free_params=[]
    free_params_idx=[]
    bounds=[]
    for i in range(n):
        for j in range(n):   
            if i==0 or j==0:
                if not (i==0 and j==1):  # do not include a_0_1 as free parameter
                                         # this parameter is constrained
                    label='a_%i_%i' % (i,j)
                    idx=coef_idx[label]    
                    free_params.append(label)
                    free_params_idx.append(idx)
                    # if i==0:
                    #     # lower bound for a_0_j parameters
                    #     bounds.append((0.1,None))
                    # else:
                    bounds.append((None,None))                
                    
                label='b_%i_%i' % (i,j)
                idx=coef_idx[label]    
                free_params.append(label)
                free_params_idx.append(idx)
                if i==0 and j==0:
                    # lower bound for b_0_j parameters
                    bounds.append((0.0,None))
                elif i==0:
                    # lower bound for b_0_j parameters
                    bounds.append((0.0,None))
                else:
                    bounds.append((None,None))
                    
    
                                
    nfp=len(free_params)
    
    # fill coefficients for b_ij a_jk = b_ik- b_jk
    ieq=-1
    for i in range(n):
        for j in range(n):
           # if not i==j:
                for k in range(n):
                    ieq+=1
                    # b_ij
                    label='b_%i_%i' % (i,j)
                    idx=coef_idx[label]
                    C_mat[ieq][idx]=1.0
                    # a_jk
                    label='a_%i_%i' % (j,k)
                    idx=coef_idx[label]
                    D_mat[ieq][idx]=1.0   
                    # b_ik
                    label='b_%i_%i' % (i,k)
                    idx=coef_idx[label]
                    E_mat[ieq][idx]=1.0               
                    label='b_%i_%i' % (j,k)
                    idx=coef_idx[label]
                    F_mat[ieq][idx]=1.0     
                
    # fill value vector
    for i in range(1,n):
        for j in range(1,n):                    
           label='a_%i_%i' % (i,j)                
           idx=coef_idx[label]
           v_vec[idx]=a_mat[i][j]
           label='b_%i_%i' % (i,j)                
           idx=coef_idx[label]
           v_vec[idx]=b_mat[i][j]  

    # set a_0_1 to 1, to constrain pattern 0 slope
    label='a_0_1' 
    #print('coef_idx='+str(coef_idx))                
    idx=coef_idx[label]
    v_vec[idx]=1.0
           
    # print ('v_vec:')
    # for i in range(ncoef):
    #     if coef_label[i] in free_params:
    #         flag='*'
    #     else:
    #         flag=' '
    #     print ('%i : %s%s = %f' % (i, coef_label[i], flag,v_vec[i]))

    for i in range(ncoef):
        v_mat[i][0]=v_vec[i]
    
    # plot system of equations      
    if figures:
        fig, axs = plt.subplots(1,8)
        axs[0].imshow(C_mat, cmap='binary')
        axs[0].set_title('C_mat')
        axs[1].imshow(D_mat, cmap='binary')
        axs[1].set_title('D_mat')
        axs[2].imshow(E_mat, cmap='binary')
        axs[2].set_title('E_mat')
        axs[3].imshow(F_mat, cmap='binary')
        axs[3].set_title('F_mat')
        axs[4].imshow(v_mat)
        axs[4].set_title('v_vec')
    
    # free parameters and bounds   
    fp_vec=numpy.zeros((nfp))
    i=-1
    for fp in free_params_idx:
        i+=1
        fp_vec[i]=v_vec[fp]
  
            
    # print ('free param ids = %s' % str(free_params))    
    # print ('free params = %s' % str(fp_vec))    
    
    args=(v_vec,free_params_idx, C_mat, D_mat, E_mat, F_mat, coef_idx, sum_NRW, n)

    ssr, res_vec=ssr_func(fp_vec, v_vec,free_params_idx, C_mat, D_mat, E_mat, F_mat, coef_idx, sum_NRW, n, return_vec=True )
    printl(logfile,'a priori ssr = %f' % ssr)
    
    ssr_reduction=[ssr,0.0]
    
    
    for i in range(ncoef):
        res_mat[i][0]=res_vec[i]
    prior_res_mat=numpy.copy(res_mat)

    if figures:
        axs[5].imshow(res_mat)
        axs[5].set_title('prior res_vec')


    maxfun=250000 # default 15000
    options={'maxfun': maxfun}
    
    inputs={'fp_vec':fp_vec, 'args':args, 'bounds':bounds, 'options':options}
    with open('input.pkl', 'wb') as file:
        pickle.dump(inputs, file)
    #sys.exit(0)
        
    try:
        opt_result=scipy.optimize.minimize(ssr_func,fp_vec, args, method=method, bounds=bounds, options=options)
    except: 
        return False,[],[], dsids
                
    opt_fp=opt_result['x']
    success=opt_result['success']
    
    if not success:
        message=opt_result['message']
        return False,message,[], dsids

    
    # print ('optimization result: %s' % str(success))
    # print ('a posteriori v_vec:')
    # for i in range(ncoef):
    #     if coef_label[i] in free_params:
    #         flag='*'
    #     else:
    #         flag=' '
    #     print ('%i : %s%s = %f' % (i, coef_label[i], flag,v_vec[i]))   
        
    #print ('solved b_0_j:')
    recovery_vector={}

    i=-1
    for fp in free_params_idx:
        i+=1
        cl=free_params[i]
        if cl[:3]=='b_0' and not cl=='b_0_0':
            recovery_vector[cl]=opt_fp[i]     
            
            
    if not success:
        printl(logfile,str(opt_result))
        sys.exit(1)
    
    # construct full optimize system
    
    ssr, res_vec=ssr_func(opt_fp, v_vec,free_params_idx, C_mat, D_mat, E_mat, F_mat, coef_idx, sum_NRW, n, return_vec=True )
    printl(logfile,'a posteriori ssr = %f' % ssr)
    ssr_reduction[1]=ssr

    for i in range(ncoef):
        res_mat[i][0]=res_vec[i]

    if figures:
        axs[6].imshow(res_mat)
        axs[6].set_title('post. res_vec')
    
        axs[7].imshow(prior_res_mat-res_mat)
        axs[7].set_title('delta')

        plt.show()
        
    printl(logfile,'')

    
    return True, recovery_vector, ssr_reduction, dsids
    
def ssr_func(fp_vec,v_vec, free_params_idx, C_mat, D_mat, E_mat, F_mat, coef_idx, sum_NRW, nb,return_vec=False):
    
    # set current value of free parameters in value vector
    i=-1
    for idx in free_params_idx:
        i+=1
        v_vec[idx]=fp_vec[i]
    
    # compute residuals
    res_vec=residual_func(v_vec, C_mat, D_mat, E_mat, F_mat)
    
    # compute SSR
    v=0.0
    for res in res_vec:
        v+=res*res    
    
    # residuals for the sum of b_0_? factors -> total NRW
    x=0.0
    for i in range(0,nb-1):
        label='b_0_%i'% (i+1)
        #print('balance %s' % label)
        x+=(v_vec[coef_idx[label]])
    if not (sum_NRW==None):
        v+=(x-sum_NRW)**2
    
        
    if return_vec:
        return v, res_vec
    else:
        return v
    
def residual_func(v_vec, C_mat, D_mat, E_mat, F_mat):
    # compute residuals
    p1=numpy.matmul(C_mat,numpy.transpose(v_vec))
    p2=numpy.transpose(numpy.matmul(D_mat,numpy.transpose(v_vec)))
    # note: Hadamard product p1 o p2
    res_vec=numpy.multiply(p1,p2)-numpy.matmul((E_mat-F_mat),numpy.transpose(v_vec))   

    
    return res_vec
    
            
def get_cfpd_fit(dat1, t1start, t1end,dat2, t2start,t2end,ilabel,jlabel):
    [datx1,daty1]=dat1
    [datx2,daty2]=dat2

    
    if t1start=='START':
        t1start=min(dat1[0])
    if t1end=='END':
        t1end=max(dat1[0])
        
    [sx1,sy1]=mystat.subsetx2(datx1, daty1, t1start,t1end)
    sorted_y1=sorted(sy1)

    if t2start=='START':
        t2start=min(dat2[0])
    if t2end=='END':
        t2end=max(dat2[0])
       
    [sx2,sy2]=mystat.subsetx2(datx2, daty2, t2start,t2end)

    if len(sy1)<=1 or len(sy2)<=1:
       print ('zero length dataset')
       sys.exit(1)


    mnf1=min(sy1)
    mnf2=min(sy2)
    
    meanf1=numpy.average(sy1)
    meanf2=numpy.average(sy2)
    
    #print ('mnf1=%f mnf2=%f' % (mnf1,mnf2))
  
    sorted_y2=sorted(sy2)

    # always resample to length of set 2!
    if len(sorted_y1)==len(sy2):
        newsy1=sorted_y1
    else:
        newsy1=mystat.resample(sorted_y1,len(sy2))
  
    
    x=newsy1
    y=sorted_y2
    
    
    # fig, ax = plt.subplots()
    # crv = ax.plot(x, y)
    # ax.set_xlabel(ilabel)
    # ax.set_ylabel(jlabel)
       
    # plt.show()
    
    
    
    
    [fit,residuals, rank, singular_values, rcond]=numpy.polyfit(x,y,1,full=True) # linear fit
    #print 'residuals='+str(residuals)
    #fit=numpy.polyfit(x[istart:iend],y[istart:iend],1) # linear fit
    slope=fit[0]
    intercept=fit[1]

    nused=len(x)
    # dx=(x[-1]-x[0])/float(nused-1)
    # xc=[]; yc=[]; sd=[]

    # calculation of standard deviation using:
    # http://en.wikipedia.org/wiki/Simple_linear_regression
    # with asymptotic assumption and q=1 because we just want to have the standard deviation
    xm=numpy.mean(x)
    ym=numpy.mean(y)
    sxxm=0.0
    #ssq=0.0
    syym=0.0
    for i in range(len(x)):
       sxxm+=(x[i]-xm)**2
       syym+=(y[i]-ym)**2


    # standard deviations of fit parameters:
    sdslope=math.sqrt((1.0/(float(nused)-2.0))*residuals[0]/sxxm)
 
    sdintercept=math.sqrt((1.0/(float(nused)-2.0))*residuals[0])*math.sqrt((1.0/(float(nused)-2.0))+xm*xm/sxxm) # according to Wolfram Math World


    ssyyhat=0.0
    sx2=0.0
    sx=0.0
    for i in range(nused):
        lx=x[i]
        ly=y[i]
        lyhat=intercept+slope*lx
        ssyyhat+=(ly-lyhat)**2
        
        sx2+=lx*lx
        sx+=lx
    s=math.sqrt(ssyyhat/(float(nused)-2.0))
    sdintercept_alt=s*math.sqrt(sx2/(nused*sx2-sx*sx))

    # coefficient of determination R^2
    # following http://en.wikipedia.org/wiki/Coefficient_of_determination
    R2=1.0-residuals[0]/syym
    #self.nppfitR2=1.0-residuals[0]/sxxm ! ERROR corrected 20120702 replacing sxxm with syym

    SSR=residuals[0]
    MSE=residuals[0]/float(nused)
    MSE_maxx=residuals[0]/float(nused)/max(x)
    RMSE=math.sqrt(residuals[0]/float(nused))
    RMSE_maxx=math.sqrt(residuals[0]/float(nused))/max(x)

    errpar={'sd_a' : sdslope,
            'sd_b': sdintercept,
            'sd_b_alt': sdintercept_alt,
            'R2' : R2,
            'SSR' : SSR,
            'MSE' : MSE,
            'MSE_maxx' : MSE_maxx,
            'RMSE' : RMSE,
            'RMSE_maxx' : RMSE_maxx,
            'MNF1' : mnf1,
            'MNF2' : mnf2,
            'MEANFLOW1': meanf1,
            'MEANFLOW2': meanf2}
    
    return slope,intercept,errpar,x,y    


def check_nonzero_multi_ts(xdat,ydat,istart,iend):
    res=True
    for lydat in ydat:
        res=check_nonzero_ts(xdat,lydat,istart,iend)
        if not res:
            break
    return res

def check_nonzero_ts(xdat,ydat,istart,iend):
        
    [sx1,sy1]=mystat.subsetx2(xdat, ydat, istart,iend)

    if len(sy1)<=1:
       print ('zero length dataset')
       sys.exit(1)
    
    # ssq=0.0
    # for item in sy1:
    #     ssq+=item*item
    
    # if item>0.0:
    if not 0.0 in sy1:
        return True
    else:
        return False
        
       
    
        
def loadandparsecsv_multi(fname,maxnline=2**31):
    #print('opening %s' % fname)
    f=open(fname,'r')
    lines=f.readlines()
    f.close()

    # determine seperator character based on first line
    hassemicolon=not lines[0].find(';')==-1
    hastab=not lines[0].find('\t')==-1
    hascomma=not lines[0].find(',')==-1
    hasspace=not lines[0].find(' ')==-1

    #print ('hassemicolon,hastab,hascomma,hasspace='+str([hassemicolon,hastab,hascomma,hasspace]))
 
    if hassemicolon:
       separator=';'
       cull1=','
       cull2='\t'
       cull3='\t'
    elif hascomma:
       separator=','
       cull1='\t'
       cull2=';'
       cull3='\t'
    elif hastab:
       separator='\t'
       cull1=','
       cull2=';'
       cull3=';'
    elif hasspace:
       separator=' '
       cull1='\t'
       cull2=','
       cull3=';'
    else:
       return False,[],[],[]

    # select column
    selectionstrings=lines[0].replace(cull1,'').replace(cull2,'').replace(cull3,'').replace('\r','').replace('\n','').split(separator)
 
    #print ('number of columns='+str(len(selectionstrings)))
    
    for j in range(len(selectionstrings)-1,0,-1):
        if selectionstrings[j]=='':
            selectionstrings.pop(j)

    n=len(selectionstrings)-1
    datx=[]
    daty=[]
    for i in range(n):
        daty.append([])

    labels=selectionstrings[1:]
    istart=1 # skip first 
    #print('n='+str(n))
    #print ('daty='+str(daty))

    # fill arrays
    for row_index in range(istart,min(len(lines),maxnline)):
       ldat=lines[row_index].replace(cull1,'').replace(cull2,'').replace(cull3,'').replace('\r','').replace('\n','').split(separator)
       # try:
       if not ldat[0]=='':
          #print('ldat='+str(ldat))
          x=ldat[0]
          if not x=='':
              datx.append(float(x))
          for i in range(1,n+1):
              y=ldat[i]
              if not y=='':
                  daty[i-1].append(float(y))
              else:
                  daty[i-1].append(None)
                  
       # except:
       #    print ('data error; ldat='+str(ldat))
       #    sys.exit(1)


    return True, datx, daty, labels 

def loadandparsecsv(fname,tscol=None, datcol=None,maxnline=2**31):
    #print('opening %s' % fname)
    f=open(fname,'r')
    lines=f.readlines()
    f.close()

    # determine seperator character based on first line
    hassemicolon=not lines[0].find(';')==-1
    hastab=not lines[0].find('\t')==-1
    hascomma=not lines[0].find(',')==-1
    hasspace=not lines[0].find(' ')==-1

    #print ('hassemicolon,hastab,hascomma,hasspace='+str([hassemicolon,hastab,hascomma,hasspace]))
 
    if hassemicolon:
       separator=';'
       cull1=','
       cull2='\t'
       cull3='\t'
    elif hascomma:
       separator=','
       cull1='\t'
       cull2=';'
       cull3='\t'
    elif hastab:
       separator='\t'
       cull1=','
       cull2=';'
       cull3=';'
    elif hasspace:
       separator=' '
       cull1='\t'
       cull2=','
       cull3=';'
    else:
       return False,[],[],[]

    # select column
    selectionstrings=lines[0].replace(cull1,'').replace(cull2,'').replace(cull3,'').replace('\r','').replace('\n','').split(separator)
 
    #print ('number of columns='+str(len(selectionstrings)))

    if len(selectionstrings)==2:
       # two columns -> first is timestamp and second is data\
       tscolumn=0
       datacolumn=1
       # assume first line to contain data; if not, it will be ignored by
       # the parser
       istart=0
    else:
       if not (tscol==None or datcol==None):
          tscolumn=tscol
          datacolumn=datcol
          istart=0 # if it is a label, it will be skipped anyway
       else:
          print("need to define timeseries and data columns")
          print("fname:%s" % fname)
          sys.exit(1)

    # fill arrays
    datx=[]
    daty=[]
    datt=[]
    for row_index in range(istart,min(len(lines),maxnline)):
       ldat=lines[row_index].replace(cull1,'').replace(cull2,'').replace(cull3,'').replace('\r','').replace('\n','').split(separator)
       try:
          x=ldat[tscolumn]
          y=ldat[datacolumn]
       except:
          print ('data error; ldat='+str(ldat))
          x=''; y=''
       if not (x=='' or y==''):
          try:
             try:
                lx=float(x)
             except:
                lx=conversion.string2exceldate(x)
             ly=float(y)
             lt=float(lx)%1.0
             #print 'lx='+str(lx)+' ly='+str(ly)+' lt='+str(lt)
             datx.append(lx)
             daty.append(ly)
             datt.append(lt)
          except:
             #print 'skipping row '+str(row_index)+' (not parsable)   x="'+x+'" y="'+y+'"'
             pass
       else:
          #print 'skipping row '+str(row_index)+' (empty)'
          pass

    return True, datx, daty, datt

# Print iterations progress
    # from: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()    
    

def cfpdplot2(datasets,fname=None, title=None):
    # create diagrams for time series and CFPD plots
    nts=len(datasets)
    if nts<2:
        return
    elif nts==2:
        ncfpd=1
    else:
        ncfpd=math.factorial(nts)/math.factorial(nts-2)/2
    ncr=math.ceil(math.sqrt(float(ncfpd)))
    figsize=(15, 8)
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs0 = gridspec.GridSpec(1, 2, figure=fig)
    gs01 = gs0[0].subgridspec(nts,1)
    gs02 = gs0[1].subgridspec(ncr, ncr)
    axs=[]
    ii=-1
    fsl=int(24.0/float(len(datasets))) # font size
    for dataset in datasets:
        key,x,y,tstart,tend=dataset
        ii+=1
        
        #isolate subset
        lx=[]; ly=[]
        for i in range(len(x)):
            if (tstart=='START' or x[i]>=tstart) and (tend=='END' or x[i]<=tend) :
                lx.append(x[i])
                ly.append(y[i])
        
        
        
        xmi=min(lx); ymi=min(ly)
        xm=max(lx); ym=max(ly)
        dx=xm-xmi
        dy=ym-ymi
        
        axs.append(fig.add_subplot(gs01[ii,:]))
        axs[-1].plot(lx,ly,'-k',linewidth=0.5)
        label=chr(97+ii)+') '+key 
        axs[-1].text(xmi+0.02*dx,ym, label, va="top", ha="left", fontsize=fsl)
        axs[-1].tick_params(labelbottom=False)
        axs[-1].set_xlim([xmi,xm])
        axs[-1].set_ylim([ymi-0.05*dy, ym+0.05*dy])
        
    keys=[a[0] for a in datasets]
    iplot=-1
    for ii in range(nts-1):
        k1=keys[ii]
        for jj in range(ii+1,nts):
            iplot+=1
            k2=keys[jj]
            
            x1=datasets[ii][1]
            y1=datasets[ii][2]
            ts1=datasets[ii][3]
            te1=datasets[ii][4]
    
            x2=datasets[jj][1]
            y2=datasets[jj][2]     
            ts2=datasets[jj][3]
            te2=datasets[jj][4]
            
            m1=numpy.mean(y1)
            m2=numpy.mean(y2)
            if m1<m2:
         
                a,b,fiterrparams,x,y=get_cfpd_fit([x1,y1],ts1,te1, 
                                              [x2,y2],ts2,te2,k1,k2)
                label=k2+' vs '+k1
            else:
                # swap -> always a>=1 in CFPD analysis
                a,b,fiterrparams,x,y=get_cfpd_fit([x2,y2],ts2,te2, 
                                              [x1,y1],ts1,te1,k2,k1)
                label=k1+' vs '+k2

           
            jcol=int( iplot % ncr )
            irow=int((iplot-jcol)/ncr)
            
            #print('irow,jcol,ncr=%i,%i,%i'% (irow,jcol,ncr))
            axs.append(fig.add_subplot(gs02[irow,jcol]))
            axs[-1].plot(x,y, '--k')
            x2=[0,x[-1]]
            y2=[b, a*x2[1]+b]
            axs[-1].plot(x2,y2,'-k')
            axs[-1].plot(x2,x2,':k')
            axs[-1].set_xlim([0,x2[1]])
            axs[-1].set_ylim([0,y2[1]])
            xm=max(x2)*0.03
            ym=max(y2)*0.97
            
            label=chr(97+nts+iplot)+') '+label 
            axs[-1].text(xm,ym, label, va="top", ha="left", fontsize=fsl)
            
            info='a=%4.2f\nb=%4.2f' % (a,b)
            info+='\nse for b=%4.2f' % fiterrparams['sd_b_alt']
            info+='\nR2=%4.3f' % fiterrparams['R2']
            info+='\nRMSE=%4.3f' % fiterrparams['RMSE']
            lx=min(x2)+0.6*(max(x2)-min(x2))
            ly=min(y2)+0.4*(max(y2)-min(y2))
            axs[-1].text(lx,ly, info, va="top", ha="left", fontsize=fsl)

            
            axs[-1].tick_params(labelbottom=False, labelleft=False)
            
    if not title==None:
        fig.suptitle(title)
            
    if fname==None:
        plt.show()
    if not fname==None:
        fig.savefig(fname, dpi=300)
        #fig.savefig(fname, figsize=figsize, dpi=300)
        
        plt.close()


def cfpdplot3(datasets,fname=None, title=None):
    # create a set of CFPD diagrams
    ncfpd=len(datasets)
    ncr=math.ceil(math.sqrt(float(ncfpd)))
    figsize=(15, 15)
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs0 = gridspec.GridSpec(1, 1, figure=fig)
    gs02 = gs0[0].subgridspec(ncr, ncr)
    axs=[]
    fsl=int(24.0/float(ncr)) # font size

        
    iplot=-1
    for dataset in datasets:
            iplot+=1
            
            x1=dataset[1]
            y1=dataset[2]
            ts1=dataset[3]
            te1=dataset[4]
    
            x2=dataset[5]
            y2=dataset[6]     
            ts2=dataset[7]
            te2=dataset[8]
            
            label=dataset[0]
            
            k1=conversion.exceldate2stringd(ts1)+'-'+conversion.exceldate2stringd(te1)
            k2=conversion.exceldate2stringd(ts2)+'-'+conversion.exceldate2stringd(te2)
         
            a,b,fiterrparams,x,y=get_cfpd_fit([x1,y1],ts1,te1, 
                                          [x2,y2],ts2,te2,k1,k2)


           
            jcol=int( iplot % ncr )
            irow=int((iplot-jcol)/ncr)
            
            #print('irow,jcol,ncr=%i,%i,%i'% (irow,jcol,ncr))
            axs.append(fig.add_subplot(gs02[irow,jcol]))
            axs[-1].plot(x,y, '--k')
            x2=[0,x[-1]]
            y2=[b, a*x2[1]+b]
            axs[-1].plot(x2,y2,'-k')
            axs[-1].plot(x2,x2,':k')
            axs[-1].set_xlim([0,x2[1]])
            axs[-1].set_ylim([0,y2[1]])
            xm=max(x2)*0.03
            ym=max(y2)*0.97
            
            label=chr(97+iplot)+') '+label 
            axs[-1].text(xm,ym, label, va="top", ha="left", fontsize=fsl)
            
            info='a=%4.2f\nb=%4.2f' % (a,b)
            info+='\nse for b=%4.2f' % fiterrparams['sd_b_alt']
            info+='\nR2=%4.3f' % fiterrparams['R2']
            info+='\nRMSE=%4.3f' % fiterrparams['RMSE']
            lx=min(x2)+0.6*(max(x2)-min(x2))
            ly=min(y2)+0.4*(max(y2)-min(y2))
            axs[-1].text(lx,ly, info, va="top", ha="left", fontsize=fsl)

            
            axs[-1].tick_params(labelbottom=False, labelleft=False)
            
    if not title==None:
        fig.suptitle(title)
            
    if fname==None:
        plt.show()
    if not fname==None:
        fig.savefig(fname, dpi=300)
       # fig.savefig(fname, figsize=figsize, dpi=300)
        
        plt.close()

def generate_clusters2(llist, n, sub=False):
    clusters=[[llist]] # n=1

    for csize in range(1,math.floor(float(len(llist))/2.0)+1):
       #print('csize=%i' % csize)
       lset=itertools.combinations(llist,csize)
       
       for cluster in lset:
           complement=[]
           for item in llist:
               if not item in cluster:
                   complement.append(item)
           if n>2:
               # subdivide complement
               subclusters=generate_clusters2(complement,n-1, sub=True)
               for subcluster in subclusters:
                   lclust=[list(cluster)]
                   for cl in subcluster:
                       lclust+=[cl]
                   clusters.append(lclust)
           else:
               clusters.append([list(cluster),list(complement)])
      
    # remove duplicates
    cclusters=[]
    for cluster in clusters:
        if not len(cluster)==n:
            keep=False
        else:
            #print('cluster='+str(cluster))
            permutations=list(itertools.permutations(cluster,len(cluster)))
            #print('permutations='+str(permutations))
            keep=True
            for permutation in permutations:
                #print('permutation='+str(permutation))
                #print('cclusters=')
                #for item in cclusters:
                #    print(str(item))
                #print('-lpermutation='+str(list(permutation)))
                if list(permutation) in cclusters:
                    keep=False
                    #print('ditch')
                    break
        if keep:
            cclusters.append(cluster)
    clusters=cclusters
    
    if not sub:         
        #print('clusters:')
        ncheck=0
        count=0
        for item in clusters:
            count+=1
            ln=0
            for subclust in item:
                ln+=len(subclust)
            check=(ln==len(llist))
            if check:
                ncheck+=1
        
    return clusters


 
def find_best_correspondence_multi2(infile, perlen=28):

    success,xdat,ydat,labels=loadandparsecsv_multi(infile)
    
    if success: 
        res=find_best_correspondence_multi_exec2(xdat,ydat,labels,perlen)
        return res
    else:
        print('error reading input file %s' % infile)


    
        
def performance(xdat, ydat, istart, perlen):    
    n=len(ydat)
    sum_sdb=0.0
    sum_R2=0.0
    sum_RMSE=0.0
    sum_sdb_alt=0.0
    mb=0.0
    bs=[]
    nn=0
    b_standard_errors=[]
    mnfs=[0.0]*n
    meanfs=[0.0]*n
    for ii in range(n-1):
        for jj in range(ii+1,n):
            nn+=1
            y1=ydat[ii]
            y2=ydat[jj]
            x1=xdat
            x2=xdat
    
            a,b,fiterrparams,x,y=get_cfpd_fit([x2,y2],istart, istart+perlen, 
                                              [x1,y1],istart, istart+perlen,2,1)
            
            sum_sdb+=fiterrparams['sd_b']
            sum_R2+=fiterrparams['R2']
            sum_RMSE+=fiterrparams['RMSE']
            sum_sdb_alt=fiterrparams['sd_b_alt']
            mb+=abs(b)
            bs.append(b)
            b_standard_errors.append(fiterrparams['sd_b'])
            
        mnfs[ii]=fiterrparams['MNF2'] # note order in call get_cfpd_fit
        meanfs[ii]=fiterrparams['MEANFLOW2'] # note order in call get_cfpd_fit
    mnfs[n-1]=fiterrparams['MNF1']
    meanfs[n-1]=fiterrparams['MEANFLOW1']
    
    if nn>0:
        sdb=numpy.std(bs)
        mean_sdb=sum_sdb/float(nn)
        mean_sdb_alt=sum_sdb_alt/float(nn)
        mean_R2=sum_R2/float(nn)
        mean_RMSE=sum_RMSE/float(nn)
    else:
        sdb=-1.0
        mean_sdb=-1.0
        mean_sdb_alt=-1.0
        mean_R2=-1.0
        mean_RMSE=-1.0
    return [istart,istart+perlen,mb,sdb,mean_sdb,mean_R2,mean_RMSE, mean_sdb_alt, b_standard_errors, mnfs, meanfs]
    
    
def find_best_correspondence_multi_exec2(xdat,ydat,labels,perlen=28, ts_start=None, ts_end=None):
      
    if ts_start==None:
        start=math.ceil(min(xdat))
    else:
        start=ts_start
    if ts_end==None:
        end=math.floor(max(xdat))
    else:
        end=ts_end
    
    results=[]
    for istart in range(start,end-perlen):
        
        res=check_nonzero_multi_ts(xdat,ydat,istart,istart+perlen)
        #print('looking at period %i-%i  nonzero=%s' % (istart,istart+perlen-1, str(res)))
        if res:  # only check correspondence if all timeseries in the time
                 # window have nonzero values
            res=performance(xdat, ydat, istart, perlen)
            #print('perf=%s' % str(res[4]))
            results.append(res)
    
    
    #sresults=sorted(results,key=operator.itemgetter(4))  
    sresults=sorted(results,key=operator.itemgetter(5), reverse=True) # highest R2 value first  
    
    return sresults
    
def find_range2(cluster):

    #print('cluster='+str(cluster))
    if len(cluster[4])<=1:
        return [0,0,[],[],[]]
    
    # label, datafile, start, end, datalabels
    #print('cluster='+str(cluster))
    runlabel, infile, istart, iend, sel_datalabels=cluster
    #print ('processing cluster %s' % runlabel)
    
    success,xdat,ydat,labels=loadandparsecsv_multi(infile)
    
    if not success:
        print('error reading input file %s' % infile)
        sys.exit(1)        
    
    datasets=[]
    for lab in sel_datalabels:
        # get best results for this cluster
        datasets.append([lab,[xdat,ydat[labels.index(lab)]],istart,iend]) 

     
    # data ranges
    summnf=0.0
    for llab in sel_datalabels:
        idx=labels.index(llab)
        mnf=1e10
        for i in range(len(xdat)):
            if xdat[i]>=istart and xdat[i]<=iend:
                mnf=min(mnf,ydat[idx][i])
        summnf+=mnf    
      
    #  SSR curve
    xc=[]
    yc=[]
    
    
    # get curve shape
    if True:
        for i in range(100):
            p=float(i)*0.01*summnf
            success, res, ssr, dsids=analyze(datasets,p, figures=False)
            xc.append(p)
            print('ssr='+str(ssr))
            
            yc.append(math.log10(ssr[1]))
        
 
    # bisection to find first transition - first derivative is zero - -> +
    threshold=1e-8
    minrange=0.0
    maxrange=0.5*summnf
    nbisection=20
    delta=0.01
    valsplat=0.1 # there is sometimes a slight negative slope just past the inflection point
                 # by adding this additional slope, we should detect the inflection point anyway
    for ibisection in range(nbisection):    
        #print('l1 iter %i' % ibisection)
        halfway=0.5*(minrange+maxrange)
        success, res, ssr1, dsids=analyze(datasets,(1.0-delta)*halfway, figures=False)
        #, res, ssr2, dsids=analyze(datasets,halfway, figures=False)
        success, res, ssr3, dsids=analyze(datasets,(1.0+delta)*halfway, figures=False)
        xc.append((1.0-delta)*halfway)
        #xc.append(halfway)
        xc.append((1.0+delta)*halfway)
        yc.append(math.log10(ssr1[1]))
        #yc.append(math.log10(ssr2[1]))
        yc.append(math.log10(ssr3[1]))
        #der=((ssr3[1]-ssr2[1])/(0.001*halfway)-(ssr2[1]-ssr1[1])/(0.001*halfway))/(0.002*halfway)
        der=(ssr3[1]-ssr1[1])/(2.0*delta*halfway)+ valsplat
        if abs(der)<threshold:
            # convergence
            l1_loc=halfway
            #print ('l1: convergence in iteration %i' % ibisection)
            break
        elif der<0.0:
            minrange=halfway
            l1_loc=halfway
        else:
            maxrange=halfway
            l1_loc=halfway
            
    return [l1_loc,summnf,xc,yc, datasets]        
    
  
  
main()
