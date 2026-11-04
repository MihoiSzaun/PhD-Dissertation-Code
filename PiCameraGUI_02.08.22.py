########################################
# IMPORTING NECESSARY LIBRARIES AND STANDARD SETTINGS
########################################

# Time imports
from time import sleep
from time import strftime
import time

# PiCamera imports
from picamera import PiCamera
camera = PiCamera()

# GUI imports
import _tkinter
import tkinter
from tkinter import *
from tkinter import messagebox

# GPIO imports
import RPi.GPIO as GPIO

# Constant camera settings
camera.resolution = (2500, 1900)

# GPIO set up
GPIO.setwarnings(False) # Block unecessary warnings
GPIO.setmode(GPIO.BOARD) # Tell how to read pin numbers

# Setting up each pin used as output
relayPinList = [37,35,33,31,29]
GPIO.setup(relayPinList, GPIO.OUT)
motorPinList = [12,11,13,15]
GPIO.setup(motorPinList, GPIO.OUT)

# Relay pins set to output, default on wiring
GPIO.output(relayPinList,1)

# Motor pins set to no output
# GPIO.output(motorPinList,0)

# Building main GUI window
root = Tk(className = 'PiCamera GUI')

########################################
# DEFINING TERMS FOR LATER
########################################

MICROSECONDS = 1000000
FRAMERATE = 0.166666666 # 0.166666666 originally
all_lasers_settings = [] # [species,ss,iso] (NOTE: due to choices made in apply() and all individual laser functions associated with full_run(), the "species" setting for settings_ABC is obsolete.)
settings_calibration = ["",0,0]
settings_405 = ["",0,0]
settings_450 = ["",0,0]
settings_532 = ["",0,0]
settings_280 = ["",0,0]
value_calib = StringVar()
value_405 = StringVar()
value_450 = StringVar()
value_532 = StringVar()
value_280 = StringVar()
ssvalue_calib = StringVar()
ssvalue_405 = StringVar()
ssvalue_450 = StringVar()
ssvalue_532 = StringVar()
ssvalue_280 = StringVar()

###############################################################################
# hard coded procedures
p1 = [["AmbrosiaArtemisiifola",0.25,100],["AmbrosiaArtemisiifola",4.000,100],
      ["AmbrosiaArtemisiifola",4.000,100],["AmbrosiaArtemisiifola",4.000,100],["AmbrosiaArtemisiifola",4.000,100]]
p2 = ["pollen2",0.4,100]
test_run =[["test",0.1,100],["test",0.2,200],["test",0.3,320],
      ["test",0.4,400],["test",0.5,500]] # calibration, 280, 405, 450, 532
###############################################################################

# Preview definitions include self for key access
def startpreview(): # For the button
   camera.framerate = 15
   camera.start_preview()
   camera.annotate_text = 'Press = key to stop preview'
def startPreview(self): # For the keystroke, must say self
   camera.framerate = 15
   camera.start_preview()
   camera.annotate_text = 'Press = key to stop preview'
def stoppreview():
   camera.annotate_text = ''
   camera.stop_preview()
def stopPreview(self):
   camera.annotate_text = ''
   camera.stop_preview()

# GPIO ports
def in1off(): # in_off turns light off via relay flip
   GPIO.output(37,1)
def in1on(): # in_on turns light on via relay flip
   GPIO.output(37,0)
   namevar.set('405nm ')
def in2off():
   GPIO.output(35,1)
def in2on():
   GPIO.output(35,0)
   namevar.set('450nm ')
def in3off():
   GPIO.output(33,1)
def in3on():
   GPIO.output(33,0)
   namevar.set('532nm ')
def in4off():
   GPIO.output(31,1)
def in4on():
   GPIO.output(31,0)
   namevar.set('650nm ')
def in5off():
   GPIO.output(29,1)
def in5on():
   GPIO.output(29,0)
   namevar.set('280nm ')

def set_camera_settings(fr='',ss='',iso=''):
    # These if-elif clauses below are not useful for a dynamic framerate. Left around as comments in case I want to revert dynamic framerate in this function.
    # if not fr: # if no framerate value is passed, set it to the default constant
    #     camera.framerate = FRAMERATE
    # elif fr:
    #     camera.framerate = fr
    try:
        if not ss and all_lasers_settings[1]: # use same ss for all lasers if nothing for ss is passed
            ss = int((float(all_lasers_settings[1]))*1000000)
            fr = 1/(float(ss)/MICROSECONDS)
            camera.framerate = fr
            camera.shutter_speed = ss
        elif ss: # specific ss passed for a laser in full run
            fr = 1/(float(ss)/MICROSECONDS)
            camera.framerate = fr
            camera.shutter_speed = ss
        else:
            raise PressSaveError
        
        if not iso and all_lasers_settings[2]: # use same iso for all lasers if nothing for iso is passed
            iso = int(all_lasers_settings[2])
            camera.iso = iso
        elif iso:
            camera.iso = iso
        else:
            raise PressSaveError
    except PressSaveError:
        messagebox.showerror("Error","Please press save.")
        return False
    return True
    
def camera_capture(species='',ss='',iso='',laser=''):
    if not species and all_lasers_settings[0]:
        species = all_lasers_settings[0]
    if not ss and all_lasers_settings[1]:
        ss = all_lasers_settings[1]
    if not iso and all_lasers_settings[2]:
        iso = all_lasers_settings[2]
    camera.start_preview(alpha=200)
    sleep(2)
    camera.capture('/home/pi/Desktop/' + species + ' # ' + laser + strftime("%Y.%m.%d %H.%M.%S") + ' ' + str(iso) + ' ISO ' + str(ss) + ' sec.jpg')
    camera.stop_preview()
        
def take_one_image():
    if set_camera_settings(): # if save has been pressed, take image. otherwise, don't.
        camera_capture(species=species.get(), laser=namevar.get())
        print((strftime("Date: %Y/%m/%d, Time: %H:%M:%S")) + ', Shutter Speed: '
             + all_lasers_settings[1] + ', ISO: ' + all_lasers_settings[2] + ', Species: ' + (all_lasers_settings[0]))

# Definitions for filter wheel:
def clockwisecycle():
    GPIO.output(motorPinList,(1,0,0,0))
    sleep(0.01)
    GPIO.output(motorPinList,(0,0,1,0))
    sleep(0.01)
    GPIO.output(motorPinList,(0,1,0,0))
    sleep(0.01)
    GPIO.output(motorPinList,(0,0,0,1))
    sleep(0.01)

def clockwisesixth():
    for fullstep in range(2):
        clockwisecycle()

def counterclockwisecycle():
    GPIO.output(motorPinList,(0,0,0,1))
    sleep(0.01)
    GPIO.output(motorPinList,(0,1,0,0))
    sleep(0.01)
    GPIO.output(motorPinList,(0,0,1,0))
    sleep(0.01)
    GPIO.output(motorPinList,(1,0,0,0))
    sleep(0.01)

def counterclockwisesixth(): # Currently corrected by hysteresisCorrection() below.
    for fullstep in range(2):
        counterclockwisecycle()
        
def hysteresisCorrection():  # Should address hysteresis of filter wheel by making sure to go clockwise during Full_Run.
    counterclockwisesixth()
    sleep(0.01)
    clockwisesixth()
        
def calibration():
    print("Current stage: Calibration")
    if checked.get() == 1:
        if set_camera_settings():
            in2on() # 405nm
            in4on() # 650nm
            camera_capture(laser='Calib ')
    elif checked.get() == 0:
        if set_camera_settings(ss=int(float(settings_calibration[1])*1000000),iso=int(settings_calibration[2])):
            in2on() # 405nm
            in4on() # 650nm
            camera_capture(all_lasers_settings[0],settings_calibration[1],settings_calibration[2],'Calib ')
    in2off()
    in4off()

def Four05():
    print("Current stage: 405")
    if checked.get() == 1:
        if set_camera_settings():
            clockwisesixth() # Filter 2
            in1on()
            camera_capture(laser='405nm ')
    elif checked.get() == 0:
        if set_camera_settings(ss=int(float(settings_405[1])*1000000),iso=int(settings_405[2])):
            clockwisesixth() # Filter 2
            in1on()
            camera_capture(all_lasers_settings[0],settings_405[1],settings_405[2],'405nm ')
    in1off()
       
def Four50():
    print("Current stage: 450")
    if checked.get() == 1:
        if set_camera_settings():
            clockwisesixth()  # Filter 3
            in2on()
            camera_capture(laser='450nm ')
    elif checked.get() == 0:
        if set_camera_settings(ss=int(float(settings_450[1])*1000000),iso=int(settings_450[2])):
            clockwisesixth()  # Filter 3
            in2on()
            camera_capture(all_lasers_settings[0],settings_450[1],settings_450[2],'450nm ')
    in2off()
   
def Five32():
    print("Current stage: 532")
    if checked.get() == 1:
        if set_camera_settings():
            for i in range(2):
                clockwisesixth() # Filter 4
            in3on()
            camera_capture(laser='532nm ')
    elif checked.get() == 0:
        if set_camera_settings(ss=int(float(settings_532[1])*1000000),iso=int(settings_532[2])):
            for i in range(2):
                clockwisesixth() # Filter 4
            in3on()
            camera_capture(all_lasers_settings[0],settings_532[1],settings_532[2],'532nm ')
    in3off()
            
def Two80():
    print("Current stage: 280")
    if checked.get() == 1:
        if set_camera_settings():
            for i in range(4):
                counterclockwisesixth() # Filter 6
            hysteresisCorrection()
            in5on()
            camera_capture(laser='280nm ')
    elif checked.get() == 0:
        if set_camera_settings(ss=int(float(settings_280[1])*1000000),iso=int(settings_280[2])):
            for i in range(4):
              counterclockwisesixth() # Filter 6
            hysteresisCorrection()
            in5on()
            camera_capture(all_lasers_settings[0],settings_280[1],settings_280[2],'280nm ')
    in5off()
     
def check_for_errors():
    # valid entry
    try:
        temp = float(ssvar.get())
    except ValueError:
        messagebox.showerror("Error","Must enter an integer or decimal.")
        ssvar.set('')
    # max ss is 6
    try:
        if ssvar.get() and float(ssvar.get()) > 6.0:
            raise MaximumError
    except MaximumError:
        messagebox.showerror("Error","Maximum shutter speed is 6.0")
        all_lasers_settings[1] = ''
        ssvar.set('')
        
def apply():
    check_for_errors()
    # checked, apply to all lasers
    if checked.get() == 1:
        all_lasers_settings[0] = species.get()
        all_lasers_settings[1] = ssvar.get()
        all_lasers_settings[2] = isovar.get()
    # unchecked, each laser gets own values
    elif checked.get() == 0:
        all_lasers_settings[0] = species.get() # use only one name globally for all images
        if namevar.get() == "280nm ":
            settings_280[1] = ssvar.get()
            settings_280[2] = isovar.get()
        if namevar.get() == "532nm ":
            settings_532[1] = ssvar.get()
            settings_532[2] = isovar.get()
        if namevar.get() == "450nm ":
            settings_450[1] = ssvar.get()
            settings_450[2] = isovar.get()
        if namevar.get() == "405nm ":
            settings_405[1] = ssvar.get()
            settings_405[2] = isovar.get()
        if namevar.get() == "Calib ":
            settings_calibration[1] = ssvar.get()
            settings_calibration[2] = isovar.get()
            
# callback() updates visuals frame
def callback():
    if checked.get() == 1:
        value_calib.set(isovar.get())
        value_405.set(isovar.get())
        value_450.set(isovar.get())
        value_532.set(isovar.get())
        value_280.set(isovar.get())
        ssvalue_calib.set(ssvar.get())
        ssvalue_405.set(ssvar.get())
        ssvalue_450.set(ssvar.get())
        ssvalue_532.set(ssvar.get())
        ssvalue_280.set(ssvar.get())
    elif checked.get() == 0:
        if namevar.get() == "405nm ":
            settings_405[0] = species.get()
            ssvalue_405.set(ssvar.get())
            value_405.set(isovar.get())
        if namevar.get() == "450nm ":
            settings_450[0] = species.get()
            ssvalue_450.set(ssvar.get())
            value_450.set(isovar.get())
        if namevar.get() == "532nm ":
            settings_532[0] = species.get()
            ssvalue_532.set(ssvar.get())
            value_532.set(isovar.get())
        if namevar.get() == "280nm ":
            settings_280[0] = species.get()
            ssvalue_280.set(ssvar.get())
            value_280.set(isovar.get())
        if namevar.get() == "Calib ":
            settings_calibration[0] = species.get()
            ssvalue_calib.set(ssvar.get())
            value_calib.set(isovar.get())
            
def select_procedure(choice):
    if choice == "pollen1":
        all_lasers_settings[0] = p1[0] # species name is index 0
        all_lasers_settings[1] = p1[1] # ss is index 1
        all_lasers_settings[2] = p1[2] # iso is index 2, throughout
        species.set(p1[0])
        ssvar.set(p1[1])
        isovar.set(p1[2])
        namevar.set('')
            
    if choice == "pollen2":
        all_lasers_settings[0] = p2[0]
        all_lasers_settings[1]= p2[1]
        all_lasers_settings[2] = p2[2]
        species.set(p2[0])
        ssvar.set(p2[1])
        isovar.set(p2[2])
        namevar.set('')
        
    if choice == "test_run":
        # uncheck the "apply to all lasers" mark because each laser gets its own settings here
        checked.set(0)
        # set actual values
        settings_calibration[0] = test_run[0][0]
        settings_calibration[1] = test_run[0][1]
        settings_calibration[2] = test_run[0][2]
        # update visuals
        ssvalue_calib.set(test_run[0][1])
        value_calib.set(test_run[0][2])
        
        settings_280[0] = test_run[1][0]
        settings_280[1] = test_run[1][1]
        settings_280[2] = test_run[1][2]
        ssvalue_280.set(test_run[1][1])
        value_280.set(test_run[1][2])
        
        settings_405[0] = test_run[2][0]
        settings_405[1] = test_run[2][1]
        settings_405[2] = test_run[2][2]
        ssvalue_405.set(test_run[2][1])
        value_405.set(test_run[2][2])
        
        settings_450[0] = test_run[3][0]
        settings_450[1] = test_run[3][1]
        settings_450[2] = test_run[3][2]
        ssvalue_450.set(test_run[3][1])
        value_450.set(test_run[3][2])
        
        settings_532[0] = test_run[4][0]
        settings_532[1] = test_run[4][1]
        settings_532[2] = test_run[4][2]
        ssvalue_532.set(test_run[4][1])
        value_532.set(test_run[4][2])
        
        # takes the last entered laser (532nm) and displays those values as if it were a manual entry
        species.set(test_run[4][0])
        ssvar.set(test_run[4][1])
        isovar.set(test_run[4][2])
        namevar.set('')
        
    # clear choice
    if choice == "-":
        all_lasers_settings[0] = ''
        all_lasers_settings[1] = 0.0
        all_lasers_settings[2] = ''
        species.set('')
        ssvar.set(0.0)
        isovar.set('')
        namevar.set('')

# switches all lasers off
# this function needs (minimal) testing!
def all_off():
    in1off()
    in2off()
    in3off()
    in4off()
    in5off()

def full_run():
    all_off() # uncomment this if it doesn't work
    start = time.time()
    calibration()
    Five32()
    Four50()
    Four05()
    Two80()
    time_taken = time.time() - start
    print("Finished. Time: ",time_taken)
    
########################################
#BUILDING PHYSICAL ATTRIBUTES OF GUI
    #Each region of GUI has its own frame
        #The LabelFrame is named and put in location
        #Widgets then located within LabelFrame for each section
    #All locations placed using grid feature
    #Widgets call on preceeding definitions
        #Exception: filter wheel defined as own class later
    #Widgets given a name
        #Tkinter type defined
        #Command defined, calls on definition
        #Appearance defined
########################################

# Preview frame
previewframe = LabelFrame(root,text='Preview Control',bd=6,relief='solid',bg='papaya whip',labelanchor='n' )
previewframe.grid(row=0, column=1)

# Building button for start preview
startP = Button(previewframe,text='Start Preview (+)',command=startpreview,bg='yellow green',activebackground='forest green' )
startP.grid(row=0,column=0,ipadx=5,ipady=5)
# Allowing the '+' key to initiate preview
startP.bind_all('+', startPreview)

# Building button for stop preview
stopP = Button(previewframe,text='Stop Preview (=)',command=stoppreview,bg='indian red',activebackground='dark red' )
stopP.grid(row=0, column=1,ipadx=5,ipady=5)
# Allowing the '=' key to kill preview
stopP.bind_all('=', stopPreview)

# GPIO laser frame
GPIOframe = LabelFrame(root,text='Laser Control',bd=6,relief='solid',bg='papaya whip',labelanchor='n' )
GPIOframe.grid(row=1, column=0, rowspan=2)

# GPIO laser buttons
on1=Button(GPIOframe,text='405nm Laser On',bg='MediumPurple1',command=in1on,activebackground='MediumPurple3' ).grid(row=1, column=3, ipady=16,padx=1,pady=2)
off1=Button(GPIOframe,text='405nm Laser Off',bg='gray25',fg='peach puff',command=in1off,activebackground='gray12',activeforeground='peach puff' ).grid(row=1, column=4, ipady=16,padx=1,pady=2)

on2=Button(GPIOframe,text='450nm Laser On',bg='deep sky blue',command=in2on,activebackground='DeepSkyBlue3').grid(row=2, column=3, ipady=16,pady=2)
off2=Button(GPIOframe,text='450nm Laser Off',bg='gray25',fg='peach puff',command=in2off,activebackground='gray12',activeforeground='peach puff' ).grid(row=2, column=4, ipady=16,pady=2)

on3=Button(GPIOframe,text='532nm Laser On',bg='chartreuse2',command=in3on,activebackground='chartreuse3' ).grid(row=3, column=3, ipady=16,pady=2)
off3=Button(GPIOframe,text='532nm Laser Off',bg='gray25',fg='peach puff',command=in3off,activebackground='gray12',activeforeground='peach puff' ).grid(row=3, column=4, ipady=16,pady=2)

on4=Button(GPIOframe,text='650nm Laser On',bg='firebrick1',command=in4on,activebackground='firebrick2' ).grid(row=4, column=3, ipady=16,pady=2)
off4=Button(GPIOframe,text='650nm Laser Off',bg='gray25',fg='peach puff',command=in4off,activebackground='gray12',activeforeground='peach puff' ).grid(row=4, column=4, ipady=16,pady=2)

on5=Button(GPIOframe,text='280nm Laser On',bg='lavender',command=in5on,activebackground='light steel blue' ).grid(row=0, column=3, ipady=16,pady=2)
off5=Button(GPIOframe,text='280nm Laser Off',bg='gray25',fg='peach puff',command=in5off,activebackground='gray12',activeforeground='peach puff' ).grid(row=0, column=4, ipady=16,pady=2)

# Camera frame
cameraframe = LabelFrame(root,text='Camera Control',bd=6,relief='solid',pady=5,padx=5,bg='papaya whip',labelanchor='n')
cameraframe.grid(row=1,column=1,padx=15)

# Checkbox
checked = IntVar()
checked.set(1) # default is checked
checkbox = Checkbutton(cameraframe,text="Apply to all lasers (i.e. not doing Full Run)",variable=checked,onvalue=1,offvalue=0,bg='papaya whip')
checkbox.grid(row=0,column=0,pady=4)

# Entry box for species name
speciesframe=Frame(cameraframe,bg='papaya whip')
speciesframe.grid(sticky=NW)
species = StringVar(root)
Label(speciesframe, text='Species Name', bg='papaya whip').grid(row=1,column=0)
Entry(speciesframe,textvariable=species,bg='papaya whip',justify='center',highlightcolor='forest green',fg='black',width=20).grid(row=1,column=1,padx=5,pady=2)
all_lasers_settings.append(species.get())

# Building name menu and defining its variable
nameframe=Frame(cameraframe,bg='papaya whip')
nameframe.grid(sticky=W)
namevar=StringVar()
nameoptions=('-', '280nm ', '405nm ', '450nm ', '532nm ', 'Calib ')
namevar.set('')
namemenu=OptionMenu(nameframe, namevar, *nameoptions)
Label(nameframe, text='Laser', bg='papaya whip').grid(row=2,column=0)
namemenu['bg']='light grey'
namemenu.grid(row=2,column=1,pady=8)

# Entry box for shutter speed
ssvarframe=Frame(cameraframe,bg='papaya whip')
ssvarframe.grid(pady=5,sticky=W)
ssvar = StringVar(root)
Label(ssvarframe, text='Shutter Speed', bg='papaya whip').grid(row=3,column=0)
Entry(ssvarframe,textvariable=ssvar,bg='papaya whip',justify='center',highlightcolor='forest green',fg='black',width=6).grid(row=3,column=2,padx=5)
Label(ssvarframe, text='(0.015 - 6.0)', bg='papaya whip').grid(row=3,column=3)
all_lasers_settings.append(ssvar.get())

# Building ISO menu and defining its variable
isoframe=Frame(cameraframe,bg='papaya whip')
isoframe.grid(sticky=W)
isovar=StringVar(root)
isooptions=('100', '200', '320', '400', '500', '640', '800')
isovar.set('')
isomenu = OptionMenu(isoframe, isovar, *isooptions)
Label(isoframe, text='ISO', bg='papaya whip').grid(row=4,column=0)
isomenu['bg']='light grey'
isomenu.grid(row=4,column=1,pady=5)
all_lasers_settings.append(isovar.get())

# Procedure Frame
procedureframe = LabelFrame(cameraframe,text='Procedure',bd=6,relief=RAISED,bg='papaya whip',labelanchor='n')
procedureframe.grid(row=6,column=0,sticky=SW)

# Drop-down menu for saved procedures
selectframe=Frame(procedureframe,bg='papaya whip')
selectframe.grid(row=2,column=0,padx=5)
selection = StringVar(root)
selection_options=['-','pollen1','pollen2','test_run'] # list can be modified versus a tuple, which cannot
selection.set('-') # default
selectmenu = OptionMenu(selectframe, selection, *selection_options,command=select_procedure)
selectmenu['bg']='gray'
selectmenu.grid(ipadx=26,ipady=1)

# Save Variables Button
save = Button(cameraframe,text='save',bg='white',command=lambda:[apply(),callback()],activebackground='red').grid(row=4,column=1,pady=15)
 
# Image Control Frame
takeimageframe = LabelFrame(cameraframe,text='Image Control',bd=10,pady=5,padx=5,bg='papaya whip',relief=RAISED,labelanchor='n')
takeimageframe.grid(row=6,column=1)

# Capture Image Button
takeimage = Button(takeimageframe,text='1 Image',bg='white',command=take_one_image,activebackground='light blue').grid(row=0,column=0,ipadx=14,ipady=5,padx=5,pady=5)

# Full Run Button
run = Button(takeimageframe,text='Full Run',bg='white',command=full_run,activebackground='light blue').grid(row=2,column=0,ipadx=16,ipady=5,pady=5)

# Visuals frame
visuals_frame = LabelFrame(root,text='Current Values',bd=8,relief=RAISED,bg='papaya whip',labelanchor='n',padx=5)
visuals_frame.grid(row=0, column=4, padx=(100,5))
Label(visuals_frame, text='ISO', bg='papaya whip').grid(row=0,column=1)
Label(visuals_frame, text='SS', bg='papaya whip').grid(row=0,column=2)

# Entry boxes for all lasers to display current values. color coded

# Calibration
Label(visuals_frame, text='Calibration', bg='papaya whip').grid(row=1,column=0)
Entry(visuals_frame,textvariable=value_calib,bg='grey',justify='center',fg='black',width=10).grid(row=1,column=1)
Entry(visuals_frame,textvariable=ssvalue_calib,bg='grey',justify='center',fg='black',width=10).grid(row=1,column=2)

# 280
Label(visuals_frame, text='280nm', bg='papaya whip').grid(row=2,column=0)
Entry(visuals_frame,textvariable=value_280,bg='lavender',justify='center',fg='black',width=10).grid(row=2,column=1)
Entry(visuals_frame,textvariable=ssvalue_280,bg='lavender',justify='center',fg='black',width=10).grid(row=2,column=2)

# 405
Label(visuals_frame, text='405nm', bg='papaya whip').grid(row=3,column=0)
Entry(visuals_frame,textvariable=value_405,bg='MediumPurple1',justify='center',fg='black',width=10).grid(row=3,column=1)
Entry(visuals_frame,textvariable=ssvalue_405,bg='MediumPurple1',justify='center',fg='black',width=10).grid(row=3,column=2)

# 450
Label(visuals_frame, text='450nm', bg='papaya whip').grid(row=4,column=0)
Entry(visuals_frame,textvariable=value_450,bg='deep sky blue',justify='center',fg='black',width=10).grid(row=4,column=1)
Entry(visuals_frame,textvariable=ssvalue_450,bg='deep sky blue',justify='center',fg='black',width=10).grid(row=4,column=2)

# 532
Label(visuals_frame, text='532nm', bg='papaya whip').grid(row=5,column=0)
Entry(visuals_frame,textvariable=value_532,bg='chartreuse2',justify='center',fg='black',width=10).grid(row=5,column=1)
Entry(visuals_frame,textvariable=ssvalue_532,bg='chartreuse2',justify='center',fg='black',width=10).grid(row=5,column=2)

# Filter wheel frame
FWframe = LabelFrame(root,text='Filter Wheel Movement',bd=6,relief='solid',bg='papaya whip',labelanchor='n')
FWframe.grid(row=1, column=4)

class Error(Exception):
    # base class
    pass
class PressSaveError(Error):
    pass
class MaximumError(Error):
    pass
    
# Filter wheel class and buttons
class Application(Frame):
    def __init__(self,master):
        super(Application, self).__init__(master)
        self.grid()
        self.count = 0 # filter 6
        self.filter = 6
        self.create_widget()
        
    def create_widget(self):
        self.counterclockwise = Button(FWframe,text='+1/6 (counterclockwise)',bg='white',command=lambda:[counterclockwisesixth(),self.positivecount()],fg='black',activebackground='grey',activeforeground='black').grid(row=0, column=0, columnspan=2, ipadx=22, ipady=16,pady=2)
        self.clockwise = Button(FWframe,text='-1/6 (clockwise)',bg='white',command=lambda:[clockwisesixth(),self.negativecount()],fg='black',activebackground='grey25',activeforeground='black').grid(row=1, column=0, columnspan=2, ipady=16,pady=2)
        self.reset = Button(FWframe,text='Reset Wheel',bg='grey',command=self.resetcount,fg='black',activebackground='light blue',activeforeground='black').grid(row=2, column=1, ipady=1, ipadx=2)
        self.countLabel = Label(FWframe,text='Count: ' + str(self.count),bg='papaya whip').grid(row=2, column=0, ipadx=15)
        self.filterLabel = Label(FWframe,text='Current Filter: ' + str(self.filter),bg='papaya whip').grid(row=4, column=0, ipadx=15)
        
    def positivecount(self):
        if (self.count % 6 == 0): # when count is a factor of 6, reset
            self.count += 1
            self.filter = 1
        else:
            self.count += 1
            self.filter += 1
        
        self.countLabel = Label(FWframe,text='Count: ' + str(self.count),bg='papaya whip').grid(row=2, column=0, ipadx=15)
        self.filterLabel = Label(FWframe,text='Current Filter: ' + str(self.filter),bg='papaya whip').grid(row=4, column=0, ipadx=15)
        
    def negativecount(self):
        self.count -= 1
        if (self.count % 6 == 0): # when count is a factor of 6 just before next term, reset
            self.filter = 6
        else:
            self.filter -= 1
        
        self.countLabel = Label(FWframe,text='Count: ' + str(self.count),bg='papaya whip').grid(row=2, column=0, ipadx=15)
        self.filterLabel = Label(FWframe,text='Current Filter: ' + str(self.filter),bg='papaya whip').grid(row=4, column=0, ipadx=15)
        
    def resetcount(self):
        if int(self.count) > 0:
            for positivereset in range(int(self.count)):
                clockwisesixth()
            self.count = 0
            self.filter = 6 # needs to actually move to a filter which is called 6
            self.countLabel = Label(FWframe,text='Count: ' + str(self.count),bg='papaya whip').grid(row=2, column=0, ipadx=15)
            self.filterLabel = Label(FWframe,text='Current Filter: ' + str(self.filter),bg='papaya whip').grid(row=4, column=0, ipadx=15)
            
        if int(self.count) < 0:
            for negativereset in range((-1)*int(self.count)):
                counterclockwisesixth()
            self.count = 0
            self.filter = 6
            self.countLabel = Label(FWframe,text='Count: ' + str(self.count),bg='papaya whip').grid(row=2, column=0, ipadx=15)
            self.filterLabel = Label(FWframe,text='Current Filter: ' + str(self.filter),bg='papaya whip').grid(row=4, column=0, ipadx=15)
            
app = Application(root)

########################################
# CLOSING COMMANDS
########################################

# Dictate loop
root.mainloop()
app.mainloop()
GPIO.cleanup()
