# VNC remote desktop instructions and setup

### VNC server installation

Check if you have vncserver installed on the host linux machine: `which vncserver`. If you don't have it installed, install your favorite version ifit with your package manager of choice, e.g:
```bash
sudo apt install tigervnc-standalone-server
```

### Linux GUI check

Make sure you have GNOME desktop gui environment set up on the linux machine. Other desktops work great but their setup technique is different. Run the following:

```bash
$ ls /usr/bin/*session
```

You should see something like `/usr/bin/gnome-session`. If you don't have that, run:

```bash
sudo apt install ubuntu-gnome-desktop
```

### Configuring vncserver

Make a password for your vnc login by executing `vncpasswd` with superuser privileges:

```bash
demo@demoserver:~$ vncpasswd
Password:
Verify:
Would you like to enter a view-only password (y/n)? n
A view-only password is not used
```

This should create a file called `~/.vnc/passwd`.

Next, paste the following into `~/.vnc/xstartup`:

```bash
#!/bin/sh
# Start Gnome 3 Desktop 
[ -x /etc/vnc/xstartup ] && exec /etc/vnc/xstartup
[ -r $HOME/.Xresources ] && xrdb $HOME/.Xresources
vncconfig -iconic &
dbus-launch --exit-with-session gnome-session &
```

### Starting VNC server

Simply run:
```bash
vncserver
```

to start your vnc session. Note which port it is on. For example, the output may look like:


```
New 'network-DGX-Station:2 (jeff)' desktop at :2 on machine network-DGX-Station

Starting applications specified in /home/jeff/.vnc/xstartup
Log file is /home/jeff/.vnc/network-DGX-Station:2.log

Use xtigervncviewer -SecurityTypes VncAuth -passwd /home/jeff/.vnc/passwd :2 to connect to the VNC server.
```

The number `:2` here indicates the port that VNC is running on. Add 5900 to that number to get the port. For example, this example is running on port 5902. This is so different users' remote desktops don't have port collisions.

To check the current sessions:
```bash
vncserver -list
```

To kill a vnc server, find out which port your session is running on. If you're running on display :2, run the command
```bash
vncserver -kill :2
```

### Connecting to vnc session from remote machine

Regardless of operating system, you will need a way to access the port on the host machine. On mac and linux, we can hse SSH tunneling to achieve this. On windows there might be some other technique that I'm not familiar with.

mac and linux (suppose port 5902):
```
ssh -L 5902:localhost:5902 username@remotehost -N
```

and keep that running while you connect using a vnc client.

#### Mac

Mac should have a vnc client natively installed. Just run the following after your ssh tunnel is active:
```
open vnc://localhost:port
```

and enter the password you have set earlier. Now, you can use the GUI of your linux machine remotely.

#### Linux

You may need to install a vnc client. One option is `apt install tigervnc-viewer`. I haven't used it before but it should be easy to set up. You can also use the windows client since they have a version for linux too.

#### Windows

After you figure out a way to ssh tunnel on Windows, you can use [VNC Viewer](https://www.realvnc.com/en/connect/download/viewer/windows/) which should be easy to set up.

### Sources
https://www.answertopia.com/ubuntu/ubuntu-remote-desktop-access-with-vnc/
