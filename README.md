# DeformableStructure
The code for simulating a deformable structure.

The main object to be using is the Structure object. It is essentially a container of Modules. Modules then are defined by two rigid plates with Links going between them. Links are defined by their stiffness and geometry.

An example use case is:

```matlab
%first define the kinds of links to use
link = Link(10e-2, 1e9*pi*3e-3^2/4);
link_soft = Link(10e-2,1e6*pi*3e-3^2/4);

%collect the links you want
links = {link_soft,link,link_soft,link,link,link};

%define the module made of links
mod = Module(5e-2*[1,1,1],5e-2*[1,1,1],2*pi/3*(0:2),2*pi/3*(0:2),links);

%then define the structure using the modules
struct = Structure([mod,mod,mod]);
```

Then the structure may be simulated for a few different kinds of situations, for now it is limited to static and quasi-static behaviors.

```matlab
%use a cable running through all the platforms to compress the structure, will save an animation
struct.minimizeEnergyCableConstraint(struct.cableLength(0,0)-0.08,0.0,0);
```