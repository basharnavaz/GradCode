# Code Documentation
Hi, this is where I will put as much as I can to help you learn the tools needed to execute Quadcopter Flights uisng PX4 Firmware and a Ground Station Computer running Robot Operating System.

## Robot Operating System
The Robot Operating System (ROS) is a flexible framework for writing robot software. It is a collection of tools, libraries, and conventions that aim to simplify the task of creating complex and robust robot behavior across a wide variety of robotic platforms.

To familiarize yourself with the different concepts of ROS. Primarily you should know, what a _package_ is and other concepts like, _nodes, topics, messages_. Although there are many tutorials on YouTube, Udacity etc, I recommend following through the tutorials in the ROS Wiki.

[ROS Tutorials](http://wiki.ros.org/ROS/Tutorials)


By the end you should be able to write publisher, subscriber, create packages, create custom messaging.

## GSL: GNU Scientific Library
The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. It is free software under the GNU General Public License. GSL installation on a Linux system is pretty straight forward and I followed the steps in this video.

[GSL Installation on Linux](https://www.youtube.com/watch?v=JrvnJaj7Ogk)

Another important aspect not covered in the ROS tutorials is how to link an external library in the ROS CMakeLists. I used GNU Scientific Library in nodes and had to link it so that the compiler understood where to find the library. These two links should help you, [Link1](https://answers.ros.org/question/28491/linking-external-libraries-to-ros/),[Link2](http://www.ros.org/wiki/rosbuild/CMakeLists/Examples#Linking_against_an_external_library).

### GSL for Kernel evaluations
My project used GSL to compute kernel expressions. Most data was arranged in the form of GSL vectors. 
