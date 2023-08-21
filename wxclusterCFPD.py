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
# import matplotlib as mpl
# mpl.use('WXAgg')
from matplotlib.backends.backend_agg import FigureCanvasAgg


#from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.optimize
import scipy.stats
#import os
import datetime
import operator
import itertools
import time
import wx
#from wx.lib.plot import PlotCanvas, PolyLine, PlotGraphics
import pygame
#import locale

tool='wxclusterCFPD'
version='0.12 020210901'

# MAJOR CHANGES
# version 0.12 
#      - include CPFDtool based time series visualization
# version 0.11 
#      - clear button + functionality
#      - notebook on output window
# version 0.1 based on clusterCFPD2 version 0.57



# disable matplotlib interactive mode
plt.ioff()

class MainWindow(wx.Frame):
    """
    Main window of application
    """

    def __init__(self, *args, **kw):
        # ensure the parent's __init__ is called
        super(MainWindow, self).__init__(*args, **kw)
        
        # locale.setlocale(locale.LC_ALL, 'us_US')
        pygame.init()
        
        # data
        self.datapanels=[]
        
        

        # create a panel in the frame
        self.panel = wx.Panel(self)

        # and create a sizer to manage the layout of child widgets
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.controlsizer= wx.BoxSizer(wx.VERTICAL)
        self.datasizer= wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.controlsizer, 0,  wx.ALL, 5)
        sizer.Add(self.datasizer, 1,  wx.ALL|wx.EXPAND, 5)
        
        # fill control sizer
        self.controlsizer.Add(wx.StaticText(self.panel,-1,'Analysis'))

        self.clearbutton=wx.Button(self.panel,-1,label='Clear')
        self.clearbutton.Bind(wx.EVT_BUTTON, self.clear)
        self.controlsizer.Add(self.clearbutton,0,wx.ALL,5)
        
        plsizer=wx.BoxSizer(wx.HORIZONTAL)
        pllabel=wx.StaticText(self.panel,-1,'period:')
        plsizer.Add(pllabel,0)
        self.perlenfield=wx.TextCtrl(self.panel,-1,'7', size=(32,20))
        plsizer.Add(self.perlenfield,0,wx.ALL,5)
        self.controlsizer.Add(plsizer,0,wx.ALL,5)

        self.applybutton=wx.Button(self.panel,-1,label='Apply period')
        self.applybutton.Bind(wx.EVT_BUTTON, self.apply_period)
        self.controlsizer.Add(self.applybutton,0,wx.ALL,5)

        
        self.matchbutton=wx.Button(self.panel,-1,label='Match')
        self.matchbutton.Bind(wx.EVT_BUTTON, self.match)
        self.controlsizer.Add(self.matchbutton,0,wx.ALL,5)
        
        nbsizer=wx.BoxSizer(wx.HORIZONTAL)
        pllabel=wx.StaticText(self.panel,-1,'nbar:')
        nbsizer.Add(pllabel,0)
        self.nbarfield=wx.TextCtrl(self.panel,-1,'5', size=(32,20))
        nbsizer.Add(self.nbarfield,0,wx.ALL,5)
        self.controlsizer.Add(nbsizer,0,wx.ALL,5)        
        
        ncsizer=wx.BoxSizer(wx.HORIZONTAL)
        pllabel=wx.StaticText(self.panel,-1,'maxnclust:')
        ncsizer.Add(pllabel,0)
        self.maxnclustfield=wx.TextCtrl(self.panel,-1,'2', size=(32,20))
        ncsizer.Add(self.maxnclustfield,0,wx.ALL,5)
        self.controlsizer.Add(ncsizer,0,wx.ALL,5)         

        self.mcfpdbutton=wx.Button(self.panel,-1,label='mCFPD')
        self.mcfpdbutton.Bind(wx.EVT_BUTTON, self.mcfpd)
        self.controlsizer.Add(self.mcfpdbutton,0,wx.ALL,5)        

        self.cfpdbutton=wx.Button(self.panel,-1,label='CFPD')
        self.cfpdbutton.Bind(wx.EVT_BUTTON, self.cfpd)
        self.controlsizer.Add(self.cfpdbutton,0,wx.ALL,5)    
        
        self.datasizer.Add(wx.StaticText(self.panel,-1,'Data'))
        self.panel.SetSizer(sizer)

        # create a menu bar
        self.makeMenuBar()

        # and a status bar
        #self.CreateStatusBar()
        #self.SetStatusText("Ready")


    def makeMenuBar(self):
        """
        A menu bar is composed of menus, which are composed of menu items.
        This method builds a set of menus and binds handlers to be called
        when the menu item is selected.
        """

        # Make a file menu with Hello and Exit items
        fileMenu = wx.Menu()
        # The "\t..." syntax defines an accelerator key that also triggers
        # the same event
        loadItem = fileMenu.Append(-1, "&Load data\tCtrl-L",
                "Load data file")
        fileMenu.AppendSeparator()
        # When using a stock ID we don't need to specify the menu item's
        # label
        exitItem = fileMenu.Append(wx.ID_EXIT)

        # Now a help menu for the about item
        helpMenu = wx.Menu()
        aboutItem = helpMenu.Append(wx.ID_ABOUT)

        # Make the menu bar and add the two menus to it. The '&' defines
        # that the next letter is the "mnemonic" for the menu item. On the
        # platforms that support it those letters are underlined and can be
        # triggered from the keyboard.
        menuBar = wx.MenuBar()
        menuBar.Append(fileMenu, "&File")
        menuBar.Append(helpMenu, "&Help")

        # Give the menu bar to the frame
        self.SetMenuBar(menuBar)

        # Finally, associate a handler function with the EVT_MENU event for
        # each of the menu items. That means that when that menu item is
        # activated then the associated handler function will be called.
        self.Bind(wx.EVT_MENU, self.OnLoad, loadItem)
        self.Bind(wx.EVT_MENU, self.OnExit,  exitItem)
        self.Bind(wx.EVT_MENU, self.OnAbout, aboutItem)


    def OnExit(self, event):
        """Close the frame, terminating the application."""
        self.Close(True)


    def OnLoad(self, event):
        """Load datafile."""
        openFileDialog = wx.FileDialog(self, "Open datafile", "", "", "csv files (*.csv)|*.csv", 
        wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        openFileDialog.ShowModal()
        infile=openFileDialog.GetPath()
        openFileDialog.Destroy()
            
        success,xdat,ydat,DMAs=loadandparsecsv_multi(infile)
        
        if success:
            # first clear bitmaps from existing panels
            for datapanel in self.datapanels:
                datapanel.clear_bitmap()
            
            for i in range(len(DMAs)):
                print ('adding DMA %s' % str(DMAs[i]))
                self.datapanels.append(DataPanel(self.panel, self))
                self.datapanels[-1].setdata(DMAs[i],xdat,ydat[i])

                if i>0:
                    self.datasizer.Add(wx.StaticLine(self,-1))
                self.datasizer.Add(self.datapanels[-1], 1,  wx.ALL|wx.EXPAND, 5)
            self.datasizer.Layout()
            #self.Fit()

            for datapanel in self.datapanels:
                print('redraw datapanel %s' + str(datapanel))
                datapanel.draw()

    def OnAbout(self, event):
        """Display an About Dialog"""
        wx.MessageBox("This is "+tool+" version "+version,
                      "About "+tool,
                      wx.OK|wx.ICON_INFORMATION)

    def clear(self,evt):
        for dp in self.datapanels:
            dp.Destroy()

        self.datapanels=[]
        
    def match(self,evt):
        pd=wx.ProgressDialog('Working', 'Finding best matching period', maximum=100, parent=None,
               style=wx.PD_APP_MODAL|wx.PD_AUTO_HIDE)
        pd.Pulse()

        xdat=[]
        ydat=[]
        labels=[]
        for dp in self.datapanels:
            
            xdat.append(dp.datx)
            ydat.append(dp.daty)
            labels.append(dp.label)
            
        perlen=self.getperlen()
        
        res=find_best_correspondence_multi_exec3(xdat,ydat,labels,perlen)
        
        pd.Close()
        
        if len(res)==0:
            wx.MessageBox('No matching period found', 'Warning', wx.OK | wx.ICON_WARNING)
        
        else:
            window=res[0][0]
            
            print('window='+str(window))
            
            for dp in self.datapanels:
                dp.window=window[labels.index(dp.label)][0:2]
                dp.setstartdate(dp.window[0])
                dp.draw()
                
    def apply_period(self,evt):
        for dp in self.datapanels:
            val=dp.startdate.GetValue()
            if not val=='None':
                st=float(val)
                dp.window=[st,st+float(self.perlenfield.GetValue())]
                dp.zoom3=dp.window[:]
                dp.update_hr_dates(3)
                dp.draw()

    def mcfpd(self,evt):
        maxnclust=int(self.maxnclustfield.GetValue())
        
        xdat=[]
        ydat=[]
        labels=[]
        windows=[]
        for dp in self.datapanels:
            if dp.checked:
                xdat.append(dp.datx)
                ydat.append(dp.daty)
                labels.append(dp.label)
                windows.append(dp.window)
        
        n=0
        for i in range(1,maxnclust+1):
            clusters=generate_clusters2(labels, i)
            n+=len(clusters)    
        
            print ('n='+str(n))
    
        pd=wx.ProgressDialog('Working', 'Performing mCFPD analysis', maximum=n, parent=None,
               style=wx.PD_APP_MODAL|wx.PD_AUTO_HIDE|wx.PD_ELAPSED_TIME|wx.PD_ESTIMATED_TIME|wx.PD_REMAINING_TIME)
        

                
        nbar=int(self.nbarfield.GetValue())
        lstr, sumstr, outputs, images=execute2(xdat,ydat,labels,windows,nbar,maxnclust, progressdialog=pd)
        pd.Close()
       
        analysis_res_window=NewWindow(parent=None,id=-1)
        for item in images:
            image_string, im_dim, subcluster=item
            if not len(image_string)==0:
                analysis_res_window.add_image_tab(image_string, im_dim, str(subcluster))

        analysis_res_window.add_info_tab(lstr)
  
    def cfpd(self,evt):
        choices=[]
        mapping={}
        for dp in self.datapanels:
            if not dp.zoom1==None:
                t1=dp.label+':'+dp.daterange2.GetValue()
                mapping[t1]=[dp,dp.zoom1]
                choices.append(t1)
            if not dp.zoom2==None:
                t1=dp.label+':'+dp.daterange3.GetValue()
                mapping[t1]=[dp,dp.zoom2]
                choices.append(t1)  
                
        dialog=wx.MultiChoiceDialog(self,'Select data', 'Select data',  choices, style=wx.OK|wx.CANCEL)
        res=dialog.ShowModal()
        
        if res==wx.ID_CANCEL:
            return
        else:
            # should be OK
            sel=dialog.GetSelections()
            dialog.Destroy

            datasets=[]
            for ii in sel:
                dp,zoom=mapping[choices[ii]]
                datasets.append([choices[ii],dp.datx, dp.daty, zoom[0],zoom[1]])
             
            image_string, im_dim=self.cfpdplot4(datasets)
            
            analysis_res_window=NewWindow(parent=None,id=-1)
            if not len(image_string)==0:
                analysis_res_window.add_image_tab(image_string, im_dim, 'CFPD')
        
        
    def getperlen(self):
        val= int(self.perlenfield.GetValue())
        self.perlenfield.SetValue(str(val))
        return val
   
    
 
    
    def cfpdplot4(self,datasets):
        # create a set of CFPD diagrams
        n=len(datasets)
        ncfpd=int(n*(n-1)/2)
        ncr=math.ceil(math.sqrt(float(ncfpd)))
        figsize=(10, 10)
        fig = plt.figure(figsize=figsize, constrained_layout=True)
        gs0 = gridspec.GridSpec(1, 1, figure=fig)
        gs02 = gs0[0].subgridspec(ncr, ncr)
        axs=[]
        fsl=int(24.0/float(ncr)) # font size
    
            
        iplot=-1
        for ii in range(n):
            for jj in range(ii+1,n):
                if not ii==jj:
                    iplot+=1
                    
                    x1=datasets[ii][1]
                    y1=datasets[ii][2]
                    ts1=datasets[ii][3]
                    te1=datasets[ii][4]
            
                    x2=datasets[jj][1]
                    y2=datasets[jj][2]
                    ts2=datasets[jj][3]
                    te2=datasets[jj][4]
                
                    label=datasets[ii][0]+' - '+ datasets[jj][0]
                    
                    k1=datasets[ii][0]
                    k2=datasets[ii][1]
                 
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
                
        # get canvas, draw and get as string
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        
        image_string=canvas.tostring_rgb() 
        s, dim = canvas.print_to_buffer()           
        plt.close()
        
        return image_string, dim
            
    
class NewWindow(wx.Frame):

    def __init__(self,parent,id):
        wx.Frame.__init__(self, parent, id, 'Analysis results', size=wx.Size(1000,700))
        #wx.Frame.CenterOnScreen(self)
        
        self.panel = wx.Panel(self)
        self.notebook = wx.Notebook(self.panel)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.ALL|wx.EXPAND, 5)
        
        self.savebutton=wx.Button(self.panel,-1,"Save images and info")
        self.savebutton.Bind(wx.EVT_BUTTON, self.save)
        sizer.Add(self.savebutton,0,wx.ALL,5)
        
        
        self.panel.SetSizer(sizer)
        
        
        
        self.Layout()
        
        self.tabs=[]
        self.data={'info':'', 'images':{}}
 
        self.Show()    
        self.Raise()
        
    def add_image_tab(self, image_string, dimensions,label):
        tab=TabPanel(self.notebook)
        self.notebook.AddPage(tab, label)

        xl,yl=dimensions
        imgWx = wx.Image(xl,yl)
        self.SetClientSize((xl,yl))
        self.SetMinSize((xl,yl))
        imgWx.SetData(image_string)
        bmp = imgWx.ConvertToBitmap()
        wxbmp=wx.StaticBitmap(tab, -1, bmp, (0,0))

        tab.sizer.Add(wxbmp,1,wx.ALL|wx.EXPAND,5)
        
        self.tabs.append(tab)
        self.data['images'][label]=imgWx
 

    def add_info_tab(self, text):
        tab=TabPanel(self.notebook)
        self.notebook.AddPage(tab, 'info')

        self.text=wx.TextCtrl(tab,-1,style=wx.TE_READONLY|wx.TE_MULTILINE)
        font1 = wx.Font(10, wx.MODERN, wx.NORMAL, wx.NORMAL, False, u'Consolas')
        self.text.SetFont(font1)
        self.text.SetValue(text)
        

        tab.sizer.Add(self.text,1,wx.ALL|wx.EXPAND,5)
        tab.sizer.Layout()
  
        
        self.tabs.append(tab)
        self.data['info']=text
        
    def save(self,evt):
        saveFileDialog = wx.FileDialog(self, "Save results", "", "", "csv files (*.txt)|*.txt", 
        wx.FD_SAVE| wx.FD_OVERWRITE_PROMPT)

        saveFileDialog.ShowModal()
        outfile=saveFileDialog.GetPath()
        saveFileDialog.Destroy()
        
        if not outfile=='':
            f=open(outfile,'w')
            f.write(self.data['info'])
            f.close()


            for key,value in self.data['images'].items():
                lkey=key.replace('[','').replace(']','').replace('\'','')
                fname=outfile.replace('.txt','_')+lkey+'.png'
                value.SaveFile(fname,type=wx.BITMAP_TYPE_PNG)
                
        
 
class TabPanel(wx.Panel):
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """"""
        wx.Panel.__init__(self, parent=parent)
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)
        
class DataPanel(wx.Panel):
    def __init__(self, parent, main):
        super(DataPanel,self).__init__(parent) 
        
        dpsizer = wx.BoxSizer(wx.HORIZONTAL)
        ctrlsizer = wx.BoxSizer(wx.VERTICAL)
        self.cpsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.main=main
        self.datx=[]
        self.daty=[]
        self.label='datapanel'
        self.window=None
        
        self.zoom1=None
        self.zoom2=None
        self.zoom3=None
        self.xrange=None
        self.yrange=None
        
        self.minx=None
        self.maxx=None
        self.miny=None
        self.maxy=None
        
        self.fontsize=8
        self.linethickness=1
        
        self.click_startx=None
        
        self.labeltext=wx.StaticText(self,-1, style = wx.ALIGN_CENTER) 
        ctrlsizer.Add(self.labeltext, 0, wx.ALL , 5)
        
        # range for mCFPD comparison
        sdlabel=wx.StaticText(self,-1, label='mCFPD:', style = wx.ALIGN_LEFT) 
        self.startdate=wx.TextCtrl(self,-1,'None',style=wx.TE_PROCESS_ENTER, size=(40,20))
        sdsizer=wx.BoxSizer(wx.HORIZONTAL)
        sdsizer.Add(sdlabel,0,wx.RIGHT,5)
        sdsizer.Add(self.startdate,0,wx.RIGHT,5)
                
        self.daterange=wx.TextCtrl(self,-1,'None',style=wx.TE_READONLY|wx.TE_RICH, size=(96,14))
        f = wx.Font(6, wx.ROMAN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, True)
        self.daterange.SetStyle(0,36,wx.TextAttr("black", wx.NullColour, f))
        ctrlsizer.Add(sdsizer, 0, wx.ALL , 5 )
        ctrlsizer.Add(self.daterange, 0, wx.ALL , 5 )
        self.startdate.Bind(wx.EVT_TEXT_ENTER , self.updatewindow)
        
        # ranges for classical CFPD comparison
        sd2label=wx.StaticText(self,-1, label='CFPD set 1:', style = wx.ALIGN_LEFT)              
        self.daterange2=wx.TextCtrl(self,-1,'None',style=wx.TE_READONLY|wx.TE_RICH, size=(96,14))
        f = wx.Font(6, wx.ROMAN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, True)
        self.daterange2.SetStyle(0,36,wx.TextAttr("black", wx.NullColour, f))
        ctrlsizer.Add(sd2label, 0, wx.ALL , 5 )
        ctrlsizer.Add(self.daterange2, 0, wx.ALL , 5 )

        sd3label=wx.StaticText(self,-1, label='CFPD set 2:', style = wx.ALIGN_LEFT)              
        self.daterange3=wx.TextCtrl(self,-1,'None',style=wx.TE_READONLY|wx.TE_RICH, size=(96,14))
        f = wx.Font(6, wx.ROMAN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, True)
        self.daterange3.SetStyle(0,36,wx.TextAttr("black", wx.NullColour, f))
        ctrlsizer.Add(sd3label, 0, wx.ALL , 5 )
        ctrlsizer.Add(self.daterange3, 0, wx.ALL , 5 )        

        # and the rest
        self.includecb = wx.CheckBox(self, label = 'include') 
        self.includecb.SetValue(True)
        ctrlsizer.Add(self.includecb, 0, wx.ALL, 5 )
        dpsizer.Add(ctrlsizer, 0, wx.ALL | wx.TOP)
     
        self.curvepanel=wx.Panel(self)

        self.cpsizer.Add(self.curvepanel, 1, wx.ALL  | wx.GROW)
        dpsizer.Add(self.cpsizer, 1, wx.ALL  | wx.GROW)
        self.wx_bitmap=wx.StaticBitmap(self.curvepanel)
        self.curvepanel.Bind(wx.EVT_PAINT, lambda event: self.draw())

        wx.EVT_LEFT_DOWN(self.wx_bitmap,self.left_down)
        wx.EVT_LEFT_UP(self.wx_bitmap,self.left_up)


        self.SetSizer(dpsizer)
        self.Layout()
        print('self.curvepanel.Size=%s' % str(self.curvepanel.GetSize()))
        print('dpsizer.Size=%s' % str(dpsizer.GetSize()))
        print('self.Size=%s' % str(self.GetSize()))
        
    def clear_bitmap(self):
        self.wx_bitmap.SetBitmap(wx.NullBitmap)    

    
    def left_down(self,evt):
        print ('left down')
        x,y=evt.GetPosition()
        self.click_startx=float(self.minx)+float(x)/float(self.xl)*float(self.maxx-self.minx)
        
        
        
    def left_up(self,evt):
        print ('left up')
        x,y=evt.GetPosition()
        endx=float(self.minx)+float(x)/float(self.xl)*float(self.maxx-self.minx)
    
        st=math.floor(self.click_startx)
        en=math.ceil(endx)
        if evt.ShiftDown():
            self.zoom1=[st,en]
            self.update_hr_dates(1)
        elif evt.ControlDown():
            self.zoom2=[st,en]
            self.update_hr_dates(2)
        else:
            self.zoom3=[st,en]
            self.window=self.zoom3[:]
            self.startdate.SetValue(str(st))
            self.update_hr_dates(3)
            
        self.draw()
    
    def setdata(self,label,xdat,ydat):
        print('label=%s' % str(label))
        print('datx=%s...%s' % (str(xdat[:3]),str(xdat[-3:])))
        print('daty=%s...%s' % (str(ydat[:3]),str(ydat[-3:])))
        self.label=label
        self.datx=xdat
        self.daty=ydat

    def updatewindow(self,evt):
        start=int(self.startdate.GetValue())
        end=start+self.main.getperlen()
        self.window=[start,end]
        self.update_hr_dates(3)
        self.draw()
        
    def update_hr_dates(self, lset):
        
        if lset==1:
            st=conversion.exceldate2stringd(self.zoom1[0])
            en=conversion.exceldate2stringd(self.zoom1[1])         
            self.daterange2.SetValue(st+'-'+en)
            f = wx.Font(6, wx.ROMAN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, True)
            self.daterange2.SetStyle(0,36,wx.TextAttr("black", wx.NullColour, f))          
        elif lset==2:
            st=conversion.exceldate2stringd(self.zoom2[0])
            en=conversion.exceldate2stringd(self.zoom2[1])            
            self.daterange3.SetValue(st+'-'+en)
            f = wx.Font(6, wx.ROMAN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, True)
            self.daterange3.SetStyle(0,36,wx.TextAttr("black", wx.NullColour, f))            
        elif lset==3:
            # update human readable dates
            st=conversion.exceldate2stringd(self.window[0])
            en=conversion.exceldate2stringd(self.window[1])
            self.daterange.SetValue(st+'-'+en)
            f = wx.Font(6, wx.ROMAN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, True)
            self.daterange.SetStyle(0,36,wx.TextAttr("black", wx.NullColour, f))
        else:
            print('lset error value=%i' % lset)
            sys.exit(1)
        
        
    def setstartdate(self,val):
        self.startdate.SetValue(str(val))
        
    def checked(self):
        return self.includecb.GetValue()



    def draw(self):
        self.labeltext.SetLabel(self.label)
  
        # box
        if not self.window==None:
            self.zoom3=self.window
        
      
        self.xl, self.yl=self.cpsizer.Size
        print('time_series_panel init: curve panel size=%s, %s' % (str(self.xl), str(self.yl)))
        

     
       
        #defaultcolors=[(0,0,0),(192,0,0),(0,192,0),(0,0,192),(128,128,128),(96,0,0),(0,96,0),(0,0,96)]
    
        lminx=1e31; lmaxx=-1e31; lminy=1e31; lmaxy=-1e31
    
        #xl*=2; yl*=2 # high res (pseudo subpixels)
        
        levels=True
        suppressxmargin=False
        datelabels=False
        levellinemultiplier=0
        labelgridline=True
        checkskip=True
        

        
        xl=self.xl
        yl=self.yl
        
        if xl==0 or yl==0:
            # no real estate yet
            self.wx_bitmap.SetBitmap(wx.Bitmap(1,1))
            return
        
        xrange=self.xrange
        yrange=self.yrange
        fontsize=self.fontsize
        linethickness=self.linethickness
   
        if len(self.datx)>0 and len(self.daty)>0:
            lminx=min(self.datx)
            lmaxx=max(self.datx)
            lminy=min(self.daty)
            lmaxy=max(self.daty)
        else:
            lminx=-1
            lmaxx=1
            lminy=-1
            lmaxy=1
    
        auto=False
        if (xrange==None):
           minx=lminx
           auto=True
        else:
           minx=xrange[0]
        if (xrange==None):
           maxx=lmaxx
           auto=True
        else:
           maxx=xrange[1]
        if (yrange==None):
           miny=lminy
           auto=True
        else:
           miny=yrange[0]
        if (yrange==None):
           maxy=lmaxy
           auto=True
        else:
           maxy=yrange[1]
    
        if auto:
           dy=maxy-miny
           miny-=0.10*dy
           maxy+=0.10*dy
           if not suppressxmargin:
              dx=maxx-minx
              minx-=0.10*dx
              maxx+=0.10*dx
        
    
        if (minx==maxx):
           minx=maxx-0.1
           maxx=maxx+0.1
        if (miny==maxy):
           miny=maxy-0.1
           maxy=maxy+0.1
    
        #print 'minx,maxx,miny,maxy='+str([minx,maxx,miny,maxy])
    
        # pygame surface
        pygame_surface=pygame.Surface([xl,yl])
        pygame_surface.fill((225,225,225))
    
        self.minx=minx
        self.maxx=maxx
        self.miny=miny
        self.maxy=maxy 

        
        
        #print '----------------------------------------------'
        # day pattern
        colbg=True
        if (colbg): 
           c=1
           dx=xl/(maxx-minx)
           for ix in range(int(minx)-1, int(maxx)+2):
              #date=conversion.minimalist_xldate_as_datetime(ix+0.1,0)
              datep1=conversion.minimalist_xldate_as_datetime(ix-0.9,0) # subtract one day to make week start on Monday instead of Sunday
              w=int(datep1.strftime('%U'))%2
              # if (ix>0):
              #    ds=conversion.exceldate2string(ix+0.1)
              # else:
              #    ds='-----'
              #print 'date='+str(ix+0.1)+" "+ds+' week number='+date.strftime('%U')+' w='+str(w)+' datestr='+date.strftime('%Y %m %d %H %M %S')
              if (w==0):
                 if (c==1):
                    col=(250,250,250)
                    #col=(240,250,240)
                 else:
                    col=(240,240,240)
                    #col=(210,220,210)
              else:
                 if (c==1):
                    col=(250,252,255)
                    #col=(240,240,250)
                 else:
                    col=(240,245,250)
    
              c*=-1
              lx=(float(ix)-minx)*dx
              #print 'lx,dx,yl='+str(lx)+' '+str(dx)+' '+str(yl)
              pygame.draw.rect(pygame_surface,col,(lx,0,dx,yl),0)
              if datelabels:
                 #lfontsize=8
                 lfontsize=int(dx/3)
                 lfontsize=min(lfontsize,10)
                 d=conversion.exceldate2string(ix+0.1).split(' ')[0].split('/')[0]
                 m=conversion.exceldate2string(ix+0.1).split(' ')[0].split('/')[1]
                 ltext=d+'/'+m
                 
                 xx=lx+0.05*dx
                 yy=0.1*yl
                 col=(150,150,150)
                 #lprint 'text='+text+' lfontsize='+str(lfontsize)+' dx='+str(dx)
                 self.pygameprint(pygame_surface,ltext,xx,yy,col, fontsize=lfontsize)
    
 
        # zoom areas
        if not self.zoom1==None:
           col=(0,255,0)
           minxzoom,maxxzoom=self.zoom1
           lx=xl*float(minxzoom-minx)/(float(maxx)-float(minx))
           dx=xl*float(maxxzoom-minxzoom)/(float(maxx)-float(minx))
           pygame.draw.rect(pygame_surface,col,(lx,0,dx,yl),1)
        if not self.zoom2==None:
           col=(0,0,255)
           minxzoom,maxxzoom=self.zoom2
           lx=xl*float(minxzoom-minx)/(float(maxx)-float(minx))
           dx=xl*float(maxxzoom-minxzoom)/(float(maxx)-float(minx))
           pygame.draw.rect(pygame_surface,col,(lx,0,dx,yl),1)
        if not self.zoom3==None:
           col=(0,255,255)
           minxzoom,maxxzoom=self.zoom3
           lx=xl*float(minxzoom-minx)/(float(maxx)-float(minx))
           dx=xl*float(maxxzoom-minxzoom)/(float(maxx)-float(minx))
           pygame.draw.rect(pygame_surface,col,(lx,0,dx,yl),1)
    
        # draw standard deviations and curves
        #colar=[(255,0,0),(0,255,0),(0,0,255), (255,255,0), (0,255,255), (255,0,255), (128,0,0),(0,128,0),(0,0,128)]
        try:
           lidx=1.0/(maxx-minx)
        except:
           lidx=0
        try:
           lidy=1.0/(maxy-miny)
        except:
           lidy=1
    
    
        # draw level lines
        amaxy=max(abs(miny),abs(maxy))
        done=False
        extrafact=1
        nloop=0
        while not done:
           if (levels and amaxy>0.0):
              # orders of magnitude 
              maom=int(math.log(amaxy,10))
              #maom=int(math.log(maxy,10))
              m=int(amaxy/10**maom)+1
              #m=int(maxy/10**maom)+1
              maxlin=float(m)*10**(maom)
              nline=int(extrafact*5*2**levellinemultiplier)
              linecount=0
              for iline in range(nline+1):
                 val=(float(iline)/float(nline))*maxlin
                 y=yl*(1.0-(val-miny)/(maxy-miny))
                 if y>=miny and y<=maxy:
                    linecount+=1
                 #print 'y='+str(y)
                 col=(180,180,180)
                 pygame.draw.line(pygame_surface,col,(0,y),(xl,y),1) 

                 if labelgridline:
                    col=(0,0,0)
                    self.pygameprint(pygame_surface,str(val),30,y-int(fontsize/2),col,fontsize) 
              # also for negative y:
              for iline in range(nline):
                 val=-(float(iline+1)/float(nline))*maxlin
                 y=yl*(1.0-(val-miny)/(maxy-miny))
                 if y>=miny and y<=maxy:
                    linecount+=1
                 #print 'y='+str(y)
                 col=(180,180,180)
                 pygame.draw.line(pygame_surface,col,(0,y),(xl,y),1) 
                 if labelgridline:
                    col=(0,0,0)
                    self.pygameprint(pygame_surface,str(val),30,y-int(fontsize/2),col,fontsize) 
              if linecount>=2 or nloop>1:
                 done=True
              else:
                 #print 'refine nloop='+str(nloop)
                 extrafact*=2
                 nloop+=1
           else:
              done=True
    
         
        # curve
        if True:
           col=[255,0,0]
           lx=self.datx
           ly=self.daty
     
    
           # determine minimum point spacing
           dx=[]
           for i in range(len(lx)-1):
              ldx=lx[i+1]-lx[i]
              if (ldx>0.0):
                 dx.append(ldx)
           if (len(dx)>0):
              lmindx=min(dx)
           else:
              #print 'lmindx set to 1; no data'
              lmindx=1
           #print 'lmindx='+str(lmindx)+' n='+str(len(lx))
    
           if (len(lx)>0):
              for i in range(len(lx)-1):
                 llx=xl*(lx[i]-minx)*lidx
                 lly=yl*(1.0-(ly[i]-miny)*lidy)
                 llx2=xl*(lx[i+1]-minx)*lidx
                 lly2=yl*(1.0-(ly[i+1]-miny)*lidy)
                 #print "prepare line "+str([lx[i],ly[i],lx[i+1],ly[i+1]])
                 if (llx2>=0.0 and llx<=xl and ((checkskip and (lx[i+1]-lx[i])<2.0*lmindx) or not checkskip) and llx2>=llx):
                    try:
                       #print 'draw'
                       pygame.draw.line(pygame_surface,col,(llx,lly),(llx2,lly2),linethickness) 
                    except:
                       print ('4) draw error coor='+str([llx,lly,llx2,lly2]))
                       print ('xl,yl,i,lx[i],lx[i+1],minx,ly[i],ly[i+1],miny,lidx,lidy='+str([xl,yl,i,lx[i],lx[i+1],minx,ly[i],ly[i+1],miny,lidx,lidy]))
                       sys.exit(1)

    
        # display
        image_string = pygame.image.tostring(pygame_surface, "RGB")
        print('about to convert pygame surface of size %s' % str(pygame_surface.get_size()))
        imgWx = wx.Image(xl,yl)
        imgWx.SetData(image_string)
        #wx_bitmap = imgWx.Scale(xl/2,yl/2).ConvertToBitmap() # high res (pseudo subpixels)
        wx_bitmap = imgWx.ConvertToBitmap()
        self.wx_bitmap.SetBitmap(wx_bitmap)    
        
    
    def pygameprint(self,surf,text,xx,yy,color, fontsize=10,angle=0):
       #print('fontsize='+str(fontsize))
       font = pygame.font.SysFont("Helvetica",int(fontsize))
       #font = pygame.font.SysFont("Courier New",fontsize)
       ren = font.render(text,1,color)
       rotate = pygame.transform.rotate
       rren=rotate(ren,angle) # note: angle in degrees ROTATION DOES NOT APPEAR TO WORK
       surf.blit(rren, (xx,yy))
   
    
def execute2(xdat,ydat,labels,windows,nbar,maxnclust, progressdialog=None):
    
    nowstring=datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')
        
    inputs=[]
    inputs.append('time='+nowstring)
    inputs.append('tool='+str(tool))
    inputs.append('code version='+str(version))
    inputs.append('labels='+str(labels))
    inputs.append('windows='+str(windows))
    inputs.append('nbar='+str(nbar))
    inputs.append('maxnclust='+str(maxnclust))
    
    
  
    # plot time series and range box
    #fig_ts=plot_timeseries2(xdat,ydat,labels,outdir,istart,iend)

    outputs=[]
    outdat=[]
    meanflows={}
    images=[]
    pdcounter=-1
    if True: 
        for nclust in range(1,maxnclust+1):
            # if not progressdialog==None:
            #     progressdialog.Update(nclust-1)
            meanflows[nclust]={}
            print('================================================')
            print('testing nclust=%i' % nclust)
            outdat.append({})
            # generate all possible clusters
            clusters=generate_clusters2(labels, nclust)
            #print ('nclust=%i clusters: %s' % (nclust,str(clusters)))
            
            # find best fitting clustering
            res=evaluate_clusters2(clusters, xdat, ydat, labels, windows, progressdialog=progressdialog,
                                   pdcounter=pdcounter)
            pdcounter+=len(clusters)
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
                lxdat=[]
                lydat=[]
                lwindows=[]
                for item in subclust:
                    idx=labels.index(item)
                    lxdat.append(xdat[idx])
                    lydat.append(ydat[idx])
                    lwindows.append(windows[idx])
                
                res=find_range3(lxdat,lydat,subclust,lwindows)
                #print('FR res='+str(res))
                [bgl_start, bgl_end, curve_x, curve_y, dataset]=res
                
                # create CFPD plot of all combinations in subcluster
                # fname=outdir+'/cfpdplots_'+runlabel+'.png'
                # plotdata=[]
                # for llab in subclust:
                #     plotdata.append([llab,xdat,ydat[labels.index(llab)],istart, iend])
                # cfpdplot2(plotdata, fname=fname, title=runlabel )
                
                # create bg leakage figure and table
                #print('dataset='+str(dataset))
                [llabels,lvalues],image_string, im_dim=plot_mcfpd_range2(dataset, bgl_start, bgl_end, curve_x, curve_y, runlabel,nbar=nbar)
                images.append([image_string,im_dim,subclust])
                for i in range(len(llabels)):
                    lstr=llabels[i]+' : '+str(lvalues[i][0])
                    outputs.append(lstr)
                    print(lstr)
                    outdat[-1]['values'][llabels[i]]=lvalues[i][0]
   
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
        lstr='----------------------------------------------------------\r\n'
        #lstr+='\n\n%s - %s\n' % (conversion.exceldate2stringd(istart),conversion.exceldate2stringd(iend-1))
        #sumstr+='\n\n%s - %s\n' % (conversion.exceldate2stringd(istart),conversion.exceldate2stringd(iend-1))
        # raw analysis output
        lstr+='\r\nraw analysis results \r\n      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\r\n'
        for nclust in range(maxnclust):
            lstr+='ncl=%2i' % (nclust+1)
            for tup in orderedvalues[nclust]:
                #print (str(tup)+' ('+str(type(tup))+')')
                if tup[1]==None:
                    lstr+=' | %1.1s   None' % (tup[0])
                else:
                    lstr+=' | %1.1s %6.2f' % (tup[0],tup[1])
            lstr+='\r\n'
 
        # corrected values
        lstr+='\r\nstacked analysis results - minimum background leakage \r\n'
        lstr+='      '
        sumstr+='\r\nstacked analysis results - minimum background leakage \r\n'
        sumstr+='      '        
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
            sumstr+=' | %8.8s' %label
        lstr+='\r\n'
        sumstr+='\r\n'
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
            lstr+='\r\n'
            sumstr+='\r\n'

        # mean flows, percentages, ranking
        lstr+='\r\nmean flows\r\n'
        lstr+='      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\r\n'
        for nclust in range(maxnclust):
            lstr+='      '
            for label in orderedlabels:
                lstr+=' | %8.3f' % meanflows[nclust+1][label]
            lstr+='\r\n'
        lstr+='\r\nminimum background leakage (% of mean flow)\r\n'
        lstr+='      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\r\n'
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

            lstr+='\r\n'            
       
        lstr+='\r\npriotization (absolute | relative)\r\n'
        lstr+='      '
        for label in orderedlabels:
            lstr+=' | %8.8s' %label
        lstr+='\r\n'
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
            lstr+='\r\n'            
                        
            

           
        print(lstr)
        outputs.append(lstr)
        
        #writelist(outputs,outdir+'/outputs.txt')
        
        return lstr, sumstr, outputs, images 
    else:
        return '\r\nfailed\r\n', None, None, None
            
 
def plot_mcfpd_range2(datasets, l1_loc, summnf, xc, yc, runlabel, nbar=5):
    
    if not len(datasets)>1:
        print('not enough datasets for %s' % runlabel)
        return [[],[]],[],[]
    
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
    
    ax.set_ylabel('Background leakage')
    ax2.set_ylabel('log10(SSR)')
    ax.set_title('')
    ax.set_xlabel('total background leakage estimate')
    ax.legend()

    # get canvas, draw and get as string
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    
    image_string=canvas.tostring_rgb() 
    s, (width, height) = canvas.print_to_buffer()           
    plt.close()
        
    # olist=[]
    # for i in range(len(xc)):
    #     olist.append([xc[i],yc[i]])
    # fname=outdir+'/SSR_curve_%s.txt' % runlabel
    # writelist(olist,fname)

    return returnvalues, image_string, (width, height)            
            
        
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
    fig.savefig(fname, figsize=figsize, dpi=300)  
    plt.close()


def evaluate_clusters2(clusters, xdat, ydat, labels, windows,progressdialog=None, pdcounter=0):

    icluster=0
    ncluster=len(clusters)
    tprev=None
    looptimes=[]
    results=[]
    
    for cluster in clusters:
        if not progressdialog==None:
            pdcounter+=1
            progressdialog.Update(pdcounter)
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
                lxdat=[]
                lydat=[]
                for lab in subcluster:
                    lxdat.append(xdat[labels.index(lab)]) 
                    lydat.append(ydat[labels.index(lab)])
                res=performance2(lxdat, lydat, windows)
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


    try:
        maxfun=250000 # default 15000
        options={'maxfun': maxfun}
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
    # for i in range(ncoef):
    #     if coef_label[i] in free_params:
    #         if coef_label[i][:3]=='b_0' and not coef_label[i]=='b_0_0':
    #             #print ('%i : %s%s = %f' % (i, coef_label[i], '*',v_vec[i])) 
    #             recovery_vector[coef_label[i]]=v_vec[i]
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
    #[sx1,sy1,si1]=mystat.subset(datx1, daty1, t1start,t1end,-1e33,1e33,0.0,1.0)
    sorted_y1=sorted(sy1)

    if t2start=='START':
        t2start=min(dat2[0])
    if t2end=='END':
        t2end=max(dat2[0])
       
    [sx2,sy2]=mystat.subsetx2(datx2, daty2, t2start,t2end)
    #[sx2,sy2,si2]=mystat.subset(datx2, daty2, t2start,t2end,-1e33,1e33,0.0,1.0)

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


def check_nonzero_multi_ts2(xdat,ydat,istart,iend):
    res=True
    for i in range(len(ydat)):
        lxdat=xdat[i]
        lydat=ydat[i]
        res=check_nonzero_ts(lxdat,lydat,istart,iend)
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
        fig.savefig(fname, figsize=figsize, dpi=300)
        
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
        fig.savefig(fname, figsize=figsize, dpi=300)
        
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
            
            #if count<100 or count>(len(clusters)-100):
                #print('%i : %s check=%s ln=%i len=%i' % (count, str(item), str(check), ln, len(llist)))
        #print('length check %i/%i okay' % (ncheck, len(clusters)))
        
        
    return clusters


 
 
        
def performance2(xdat, ydat, windows):    
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
            x1=xdat[ii]
            x2=xdat[jj]
    
            istart1,iend1=windows[ii]
            istart2,iend2=windows[jj]
            
            a,b,fiterrparams,x,y=get_cfpd_fit([x2,y2],istart2, iend2, 
                                              [x1,y1],istart1, iend1,2,1)
            
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
    return [windows,windows,mb,sdb,mean_sdb,mean_R2,mean_RMSE, mean_sdb_alt, b_standard_errors, mnfs, meanfs]
    
    
def find_best_correspondence_multi_exec3(xdat,ydat,labels,perlen=28, ts_start=None, ts_end=None):
      
    if ts_start==None:
        start=-1e33
        for lxdat in xdat:
            start=max(start,math.ceil(min(lxdat)))
    else:
        start=ts_start
    if ts_end==None:
        end=1e33
        for lxdat in xdat:
            end=min(end,math.floor(max(lxdat)))        
    else:
        end=ts_end
    
    start=int(start)
    end=int(end)
    print('start-end = %i - %i' % (start,end))
    results=[]
    for istart in range(start,end-perlen):
        
        windows=[]
        for i in range(len(labels)):
            windows.append([istart,istart+perlen])
        
        res=check_nonzero_multi_ts2(xdat,ydat,istart,istart+perlen)
        #print('looking at period %i-%i  nonzero=%s' % (istart,istart+perlen-1, str(res)))
        if res:  # only check correspondence if all timeseries in the time
                 # window have nonzero values
            
            res=performance2(xdat, ydat, windows)
            #print('perf=%s' % str(res[4]))
            results.append(res)
    
    
    #sresults=sorted(results,key=operator.itemgetter(4))  
    sresults=sorted(results,key=operator.itemgetter(5), reverse=True) # highest R2 value first  
    
    return sresults
    
def find_range3(xdat,ydat,labels, windows):

        
    
    datasets=[]
    for i in range(len(labels)):
      # get best results for this cluster
        datasets.append([labels[i],[xdat[i],ydat[i]],windows[i][0],windows[i][1]]) 

     
    # data ranges
    summnf=0.0
    for idx in range(len(labels)):
        mnf=1e10
        for i in range(len(xdat[idx])):
            if xdat[idx][i]>=windows[idx][0] and xdat[idx][i]<=windows[idx][1]:
                mnf=min(mnf,ydat[idx][i])
        summnf+=mnf    
    
    print('summnf=%f' % summnf)
    
    #  SSR curve
    xc=[]
    yc=[]
    
    
    # get curve shape
    if True:
        for i in range(100):
            p=float(i)*0.01*summnf
            success, res, ssr, dsids=analyze(datasets,p, figures=False)
            xc.append(p)
            #print('ssr='+str(ssr))
            
            if not ssr[1]<=0.0:
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
        #print('ssr1='+str(ssr1))
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
    
  
  
if __name__ == '__main__':
    # When this module is run (not imported) then create the app, the
    # frame, show it, and start the event loop.
    app = wx.App()
    frm = MainWindow(None, title=tool+' '+version)
    frm.Show()
    frm.Maximize(True)
    app.MainLoop()