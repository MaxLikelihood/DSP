SimpleTone
===========================
A keyboard interactive tone generator for generating tone at any particular frequency for specified duty cycle. 
How to Compile?
==========================
First, make sure you have all the necessary libraries. Now, run the following in the source directory `gcc -o SimpleTone SimpleTone.c -lportaudio -lm -w`
Instructions
==========================
At any time during execution, press

`q` - To quit and close audio stream

`u` - To update the frequency and associated volume

`p` - To update the duration & pause (duty cycle)

After launch, if no error occurs in initializing PortAudio, a list of audio devices, including the current default, will be displayed. For example,
  
		Device Index	Device Name
		0		HDA Intel PCH: ALC269VC Analog (hw:0,0)
		1		HDA Intel PCH: HDMI 0 (hw:0,3)
		2		sysdefault
		3		front
		4		surround40
		5		surround51
		6		surround71
		7		hdmi
		8		pulse
		9		dmix
		10		default
		Default Input Device: default
		Default Output Device: default
  
Each device is given an index, choose the appropriate output device by entering its index at the next prompt,

		Choose Output Device: 3
  
Selected device name will be displayed for verification,

		Output Device Selected: front
  
Choose the number of channels for output at the next prompt (currently only stereo is supported),

		Number of Channels: 2
  
Application will attempt to verify the default parameters for output with the device chosen,
  
		..Verifying Format:
		Device Latency: 0.011610
		Sample Rate: 44100.000000
  
If successful, you will see the following, otherwise an error will display and you will be ask for a new selection,
  
		..Format Is Supported
  
Next, input the frequency to generate, an example is,
  
		Frequency: 18000
  
Then the volume, which is out of 1.0
  
		Volume(/1.0): 1
  
The next two prompts will ask you to input the duty cycle, namely the duration & pause,
  
		Duration(ms): 1000
		Pause(ms): 1000
  
If all goes well, the application will attempt to open & start the audio stream; at any time, you can always interrupt the application by the following keys
  
`q` - To quit and close audio stream
  
`u` - To update the frequency and associated volume
  
`p` - To update the duration & pause (duty cycle)
  
That's it, enjoy!
  
Note
-----------------------------
If the device you know is not appearing in the list, try waiting for a bit and restart the application (Force quit with `Ctrl+C`). If problem still persists, try rebooting.


