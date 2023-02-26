1. What this does
=================
As noted in the description, this code is an implementation of Algorithm A in my dissertation, which can be found at 
https://academicworks.cuny.edu/gc_etds/5187. More specifically, it is an impelementation for the building for the special linear group over Q[1/t,t]. For those unfamiliar with this topic, please keep in mind that these terms have nothing to do with architecture; they are mathematical objects. In short, this algorithm constructs a frame for an apartment that contains both of two given chambers. These chambers which are defined by coordinates and related by the change of basis matrix between their apartment frames. As it is currently written, this program uses a specific matrix in my example family in my thesis. However, this matrix and the chambers can easily be adjusted if one wishes to use this code for their own purposes. More can information on this can be found in the comments of the main file.

2. Why this was written
=======================
I originally developed a Pythonic implementation of my algorithm for the sake of testing purposes. Doing it by hand was very time consuming, and this shaved a lot of time off of the process and eventually I would use it to test large sample sets when searching for chambers / matrices satisfying certain properties. The code uploaded here is a cleaned up version of that original code, rewritten from scratch for the sake of readability and ease of use.

3. Why this is uploaded here
============================
I primarily uploaded this code here as a sample of Python code for future employers. However, if anyone wishes to use this for their own research or expand on it, they are free to do so. I only request that you do not plagiarize this work without providing credit if you use it for any publications. I would also very much appreciate it if you would reach out to me if you wish to expand on it or if you see any issues. You can contact me at michael_ferguson@alumni.brown.edu.
