#!/usr/bin/env python

#Author: rosencrantz@gmx.net
#License: GPL

import wx
import wx.lib.scrolledpanel as scp
from numpy import array,mean,uint8
from PIL import Image
from scipy.ndimage import gaussian_filter1d,rotate, binary_opening, label
from scipy.stats import threshold as scthresh
import os
import pickle

class AutoWin(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, -1, "Detect Plot", size=(150,200))
        self.flag=False
        self.colour=(0,0,0)
        self.parent=wx.GetApp().TopWindow
        self.bit=wx.Panel(self,-1,size=(140,30))
        self.rb1 = wx.RadioButton(self,-1,label='line plot',style=wx.RB_GROUP)
        self.rb2 = wx.RadioButton(self,-1,label='scatter plot')
        self.pick=wx.Button(self,1,label='Pick plot colour')
        self.calc=wx.Button(self,2,label='Calculate!')
        sizer=wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.rb1,5, wx.EXPAND|wx.ALL,5)
        sizer.Add(self.rb2,5, wx.EXPAND|wx.ALL,5)
        sizer.Add(self.pick,5, wx.EXPAND|wx.ALL,5)
        sizer.Add(self.bit,5, wx.EXPAND|wx.ALL,5)
        sizer.Add(self.calc,5, wx.EXPAND|wx.ALL,5)
        self.SetSizer(sizer)
        self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
        self.Bind(wx.EVT_BUTTON, self.ColourPick,id=1)
        self.Bind(wx.EVT_BUTTON, self.Calculate,id=2)
        
    def ColourPick(self,event):
        self.flag=True
        self.CaptureMouse()    
    
    def OnLeftUp(self, event):
        if self.flag:
            s = wx.ScreenDC()
            w, h = s.Size.Get()
            b = wx.EmptyBitmap(1, 1)
            m = wx.MemoryDCFromDC(s)
            m.SelectObject(b)
            m.BlitPointSize((0,0), (1,1), s, wx.GetMousePosition())
            self.colour=m.GetPixel(0,0).Get()
            m.SelectObject(wx.NullBitmap)
            dc = wx.PaintDC(self.bit)
            dc.SetBrush(wx.Brush(self.colour))            
            dc.SetPen(wx.Pen(self.colour)) 
            dc.DrawRectangle(0,0,140,40)
            self.ReleaseMouse()
        else: pass
            
    def Calculate(self,event):
        lim=5
        im=array(Image.open(self.parent.image))
        binim=(((self.colour[0]-lim)<im[:,:,0])&((self.colour[0]+lim)>im[:,:,0]))*(((self.colour[1]-lim)<im[:,:,1])&((self.colour[1]+lim)>im[:,:,1]))*(((self.colour[2]-lim)<im[:,:,2])&((self.colour[2]+lim)>im[:,:,2]))
        #case line
        if self.rb1.GetValue():
            pt=[]
            for i in range(im.shape[1]):
                v=mean((binim[:,i]>0).nonzero())
                if v>0: pt=pt+[[i,v]]
        if self.rb2.GetValue():
            binim=binim.astype(uint8)
            blobs,nblobs=label(binim)
            pt=[]
            for i in range(1,nblobs+1):
                a=(blobs==i).nonzero()
                pt=pt+[[mean(a[1]),mean(a[0])]]
        self.parent.scp.im.points=pt
        print pt
        self.parent.sZoom()        
            
        
class HelpWin(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, -1, "Usage", size=(300,250))
        self.helptext="""Detect manually:
        -Set points by left click on image, remove previous points by right click.
        -Get higher precision by zooming in.
        -If no image is loaded, left click opens file dialog.
        -If the 'axis points' check box is activated, points will be used for axis calibration and drawn in red.
         You can set as many axis points as you like (min. 2), the program will use the minimum/maximum x/y values of the set to calculate the calibration scale.
        -Enter data values for calibration points in the comma-separated form 'xmin,ymin,xmax,ymax' into the text field.
        -Saving/loading a project restores image, data and calibration points, calibration values and zoom factor.
        -Loading a new image erases all data.
        -Saved data files are space separated ASCII (floats with two trailing digits, adjust in line 239 of source code).
        
        Detect automatically ('Detect...'):
        -select plot colour by pressing colour picker button and left clicking image in main window. RGB tolerance is +-5 /255
        -radio button 'Line plot': averages all pixels with right colour over x=const. slices. 
        -radio button 'Scatter plot': finds points by blob detections. Fails on overlapping data points.
        -The Calculate button produces a (x,y) list of data points and plots them in the main window. All previous data is erased.
        -Axis calibration data has to be entered manually."""
        self.text=wx.TextCtrl(self,-1, self.helptext, style=wx.TE_MULTILINE)
        self.text.SetEditable(False)
        sizer=wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.text,5, wx.EXPAND|wx.ALL,5)
        self.SetSizer(sizer)
        
        
    def OnResize(self,event):
        print self.GetSize()
        self.text.Wrap(self.GetSize()[0]-10)


class MyPanel(scp.ScrolledPanel):
    def __init__(self, parent):
        scp.ScrolledPanel.__init__(self,parent)
        self.SetScrollRate(1,1)
        self.EnableScrolling(True,True)
        self.parent=parent
        self.im=MyImage(self)
        imbox=wx.BoxSizer(wx.HORIZONTAL)
        imbox.Add(self.im)
        self.SetSizer(imbox)
        self.im.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))
        

class MyImage(wx.StaticBitmap):
    def __init__(self, parent):
        wx.StaticBitmap.__init__(self, parent,-1,wx.EmptyBitmap(280,280,-1), (5,5))
        self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
        self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
        self.axes=[]
        self.points=[]
        self.parent=parent
        self.pparent=parent.parent
        self.scale=1.
    
    def OnLeftUp(self,event):        
        if self.pparent.image !='':
            spt=(self.parent.GetScrollPos(wx.HORIZONTAL),self.parent.GetScrollPos(wx.VERTICAL))
            pt=event.GetPosition()
            pt=(pt[0]+spt[0], pt[1]+spt[1])
            if self.pparent.check.GetValue(): 
                self.axes=self.axes+[[pt[0]/self.scale,pt[1]/self.scale]]
            else: 
                self.points=self.points+[[pt[0]/self.scale,pt[1]/self.scale]]
            self.Redraw(1000000)    
        else:
            self.pparent.Open()
        
    def OnRightUp(self,event):
        self.Redraw(-1)
        
    def Redraw(self,num):
        if self.pparent.image !='':
            im=wx.Image(self.pparent.image)
            im.Rescale(im.GetSize()[0]*self.scale,im.GetSize()[1]*self.scale)
            bit=wx.BitmapFromImage(im)
            mydc = wx.MemoryDC(bit)
            if self.pparent.check.GetValue(): 
                self.axes=self.axes[:num]
            else:
                self.points=self.points[:num]
            
            mydc.SetPen(wx.Pen('RED'))
            for p in self.axes:     
                    mydc.DrawLine(int(p[0]*self.scale),int((p[1]-5)*self.scale), int(p[0]*self.scale),int((p[1]+5)*self.scale))
                    mydc.DrawLine(int((p[0]-5)*self.scale),int(p[1]*self.scale), int((p[0]+5)*self.scale),int(p[1]*self.scale))     
            mydc.SetPen(wx.Pen('SPRING GREEN'))
            for p in self.points:  
                mydc.DrawLine(int(p[0]*self.scale),int((p[1]-5)*self.scale), int(p[0]*self.scale),int((p[1]+5)*self.scale))
                mydc.DrawLine(int((p[0]-5)*self.scale),int(p[1]*self.scale), int((p[0]+5)*self.scale),int(p[1]*self.scale))    
            mydc.SelectObject(wx.NullBitmap)
            self.SetBitmap(bit)
       
        
class MyFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Plot digitiser", size=(450,350))
        save=wx.Button(self,1,"Save data...",size=(140,-1)) 
        fopen=wx.Button(self,2,"Image...",size=(140,-1)) 
        helpb=wx.Button(self,3,"Usage...",size=(140,-1)) 
        psave=wx.Button(self,4,"Save project...",size=(140,-1))
        pload=wx.Button(self,5,"Load project...",size=(140,-1)) 
        auto=wx.Button(self,7,"Detect...",size=(140,-1)) 
        self.check=wx.CheckBox(self,-1,label='axis points', size=(140,-1))
        self.check.SetValue(True)
        self.text=wx.TextCtrl(self,-1,'Enter axis points',size=(140,-1))
        self.cb = wx.ComboBox(self, 6, choices=['20%','50%','100%','150%','200%','300%','400%'], size=(140,-1))
        self.cb.SetValue('100%')
        self.image=''
        self.scp=MyPanel(self)
        self.cdir=os.getcwd()
        hbox=wx.BoxSizer(wx.HORIZONTAL)
        vbox=wx.BoxSizer(wx.VERTICAL)
        hbox.Add(self.scp, 1, wx.EXPAND|wx.ALL,5)
        vbox.Add(self.text, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(self.check, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(self.cb, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(save, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(fopen, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(helpb, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(psave, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(pload, 5, wx.EXPAND|wx.ALL,5)
        vbox.Add(auto, 5, wx.EXPAND|wx.ALL,5)
        hbox.Add(vbox)        
        self.SetSizer(hbox)
        self.Bind(wx.EVT_BUTTON, self.Save, id=1)
        self.Bind(wx.EVT_BUTTON, self.Bopen, id=2)
        self.Bind(wx.EVT_BUTTON, self.Help, id=3)
        self.Bind(wx.EVT_BUTTON, self.SaveProject, id=4)
        self.Bind(wx.EVT_BUTTON, self.LoadProject, id=5)
        self.Bind(wx.EVT_COMBOBOX, self.Zoom, id=6)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Bind(wx.EVT_BUTTON, self.Auto, id=7)
        
    def Save(self,event):
        flag=True
        if self.text.GetValue=='Enter axis points':
            dlg = wx.MessageDialog(self, "Please provide axis calibration values first.", "Warning", style=wx.OK)
            dlg.ShowModal()
            flag=False
        if flag:    
            dlg = wx.FileDialog(self, "Choose data file", self.cdir, style=wx.SAVE)        
            if dlg.ShowModal() == wx.ID_OK:
                self.cdir=os.path.dirname(dlg.GetPath())
                size=wx.Image(self.image).GetSize()
                ax=array(self.scp.im.axes)
                ax[:,1]=size[1]-ax[:,1]
                xmin=min(ax[:,0])
                xmax=max(ax[:,0])
                ymin=min(ax[:,1])
                ymax=max(ax[:,1])
                cal=array(self.text.GetValue().split(',')).astype(float)
                xfac=(cal[2]-cal[0])/(xmax-xmin)
                yfac=(cal[3]-cal[1])/(ymax-ymin)
                points=array(self.scp.im.points).astype(float)
                points[:,0]=(points[:,0]-xmin)*xfac+cal[0]
                points[:,1]=(size[1]-points[:,1]-ymin)*yfac+cal[1]
                with open(dlg.GetPath(),'w') as f:
                    for point in points: f.write('{0:.2f} {1:.2f}\n'.format(point[0],point[1]))

    def Bopen(self,event):
        self.Open()
        
    def Open(self):
        flag=False
        if self.image!='':
            msg = wx.MessageDialog(None, 'This will erase all axis and data points. Cancel if you want to save first.', 'Warning', wx.OK | wx.CANCEL)
            if msg.ShowModal() == wx.ID_OK:
                flag=True
        else: flag=True
        if flag:
            dlg = wx.FileDialog(self, "Choose image", self.cdir, style=wx.OPEN)        
            if dlg.ShowModal() == wx.ID_OK:
                self.image=dlg.GetPath()
                self.cdir=os.path.dirname(dlg.GetPath())
                self.cb.SetValue('100%')
                self.scp.im.scale=1.
                im=wx.Image(self.image)
                im.Rescale(im.GetSize()[0]*self.scp.im.scale,im.GetSize()[1]*self.scp.im.scale)
                ds=wx.GetDisplaySize()
                ws=(im.GetSize()[0]+120,im.GetSize()[1]+20)
                if ws[0]<ds[0] and ws[1]<ds[1]:
                    self.SetSize((im.GetSize()[0]+120,im.GetSize()[1]+20))
                else:
                    self.SetSize(ds)
                self.scp.im.SetBitmap(wx.BitmapFromImage(im))
                self.scp.im.points=[]
                self.scp.im.axes=[]
    
    def Help(self,event):
        HelpWin().Show()
    
    def Zoom(self,event): self.sZoom()
    
    def sZoom(self):
        try: 
            sc=float(self.cb.GetValue()[:-1])/100.
            self.scp.im.scale=sc
            im=wx.Image(self.image)
            ds=wx.GetDisplaySize()
            ws=(int(im.GetSize()[0]*sc+120),int(im.GetSize()[1]*sc+20))
            self.SetSize((min(ws[0],ds[0]),min(ws[1],ds[1])))
            im.Rescale(im.GetSize()[0]*sc,im.GetSize()[1]*sc)
            self.scp.im.SetBitmap(wx.BitmapFromImage(im))
            self.scp.im.Redraw(100000)
            
        except ValueError:
            self.cb.SetValue('{0:d}%'.format(int(self.scp.im.scale*100)))
    
    def SaveProject(self,event):
        self.sSaveProject()
        
    def sSaveProject(self):
        data={
            'axes':self.scp.im.axes,
            'points':self.scp.im.points,
            'scale':self.scp.im.scale,
            'image':self.image,
            'calib':self.text.GetValue()
            }
        flag=False
        dlg = wx.FileDialog(self, "Choose project file", self.cdir, style=wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            with open(dlg.GetPath(),'wb') as f:
                pickle.dump(data,f)
            print 'saved as '+dlg.GetPath()
            flag=True
        dlg.Destroy()
        return flag
        
    def LoadProject(self,event):
        dlg = wx.FileDialog(self, "Open project file", self.cdir, style=wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            with open(dlg.GetPath(),'rb') as f:
                data=pickle.load(f)
            self.image=data['image']
            self.text.SetValue(data['calib'])
            self.cb.SetValue('{0:d}%'.format(int(100*data['scale'])))
            self.scp.im.scale=data['scale']
            self.scp.im.points=data['points']
            self.scp.im.axes=data['axes']
            self.sZoom()
        dlg.Destroy()
        
    def Auto(self,event):
        AutoWin().Show()
        
    def OnClose(self,event):
        msg = wx.MessageDialog(None, 'Exit without saving?', 'Warning', wx.YES | wx.NO)
        if msg.ShowModal()==wx.ID_NO:
            if self.SaveProject():
                msg.Destroy()
                self.Destroy()
        else:
            msg.Destroy()
            self.Destroy()
            
            
app=wx.App()
frame=MyFrame()
frame.Show()
app.MainLoop()
app.Destroy()
