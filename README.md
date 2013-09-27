DSP (Digital Signal Processing)
===============================
Inside this repository you will find varies C implemented applications for your digital signal processing needs.
How to Use?
---------------
Simply download the source file and compile it, that's it!
Do I need additional libraries?
-----------------------------
Yes, [PortAudio][1], a cross-platform API for interacting with onboard audio devices, is required. Non-standard includes will be provided in the root repository directory. Standard includes will be self explanatory in the `#include` headers.
Note
-----------------------------
Some applications in this repository require Pthread, an UNIX based multi-threading library. The standard LINUX implementation has been followed when using such library. 

[1]: http://www.portaudio.com/ "PortAudio Mainpage"