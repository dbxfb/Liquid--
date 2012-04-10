Liquid--
========

Description
-----------

This is a small C++ project, all about computing a fluid simulation in 2D and
rendering it using OpenGL :)


Usage
-----

This repository contains the Eclipse CDT project files (for compilations etc.)
because this is a small project and I'm currently lazy :)

You can add some dark fluid by clicking/dragging left inside the window, and
induce movements (change velocity vectors) by dragging using the right mouse
button.

There is a small option menu which can be opened using the middle mouse button,
it can be used to:

-   Toggle the (ugly) pixel blending algorithm
-   That's all (for now)


Future
------

There is no real roadmap, but I hope I'll make it better / faster later,
starting with some cleanup and refactoring (I adapted the code from a Java
implementation of Jos Stam's method described in “Stable Fluids”).

Maybe I'll even write an OpenCL backend (or make documentation).
